#!/usr/bin/env bash
set -euo pipefail

# Prokaryotic WGS paired-end Illumina analysis pipeline.
# Steps:
#   1) Raw read QC (FastQC + MultiQC)
#   2) Read trimming/filtering (fastp)
#   3) Assembly (SPAdes)
#   4) Assembly QC (QUAST + BUSCO optional)
#   5) Gene prediction (Prodigal)
#
# Expected input naming convention:
#   <sample>_R1.fastq.gz and <sample>_R2.fastq.gz
# or
#   <sample>_1.fastq.gz and <sample>_2.fastq.gz

usage() {
  cat <<USAGE
Usage:
  $(basename "$0") \
    -i <input_fastq_dir> \
    -o <output_dir> \
    -t <threads> \
    [--min-len 50] \
    [--busco-lineage bacteria_odb10] \
    [--disable-busco]

Required tools:
  fastqc, multiqc, fastp, spades.py, quast.py, prodigal
Optional:
  busco (if --disable-busco is not set)

Example:
  $(basename "$0") -i ./raw_reads -o ./analysis_out -t 16 --busco-lineage bacteria_odb10
USAGE
}

INPUT_DIR=""
OUT_DIR=""
THREADS=""
MIN_LEN=50
RUN_BUSCO=true
BUSCO_LINEAGE="bacteria_odb10"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      INPUT_DIR="$2"; shift 2 ;;
    -o|--output)
      OUT_DIR="$2"; shift 2 ;;
    -t|--threads)
      THREADS="$2"; shift 2 ;;
    --min-len)
      MIN_LEN="$2"; shift 2 ;;
    --busco-lineage)
      BUSCO_LINEAGE="$2"; shift 2 ;;
    --disable-busco)
      RUN_BUSCO=false; shift ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "[ERROR] Unknown option: $1" >&2
      usage
      exit 1 ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$OUT_DIR" || -z "$THREADS" ]]; then
  echo "[ERROR] -i, -o and -t are required." >&2
  usage
  exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "[ERROR] Input directory not found: $INPUT_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"/{00_raw_qc,01_trimmed,02_trim_qc,03_assembly,04_assembly_qc,05_gene_prediction,logs}

check_tool() {
  local tool="$1"
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[ERROR] Required tool not found in PATH: $tool" >&2
    exit 1
  fi
}

for req in fastqc multiqc fastp spades.py quast.py prodigal; do
  check_tool "$req"
done

if [[ "$RUN_BUSCO" == true ]]; then
  if ! command -v busco >/dev/null 2>&1; then
    echo "[WARN] BUSCO not found; disabling BUSCO step."
    RUN_BUSCO=false
  fi
fi

# Discover paired-end samples.
# Supports *_R1.fastq.gz/*_R2.fastq.gz and *_1.fastq.gz/*_2.fastq.gz.
declare -A R1_FILES

auto_discover_samples() {
  local r1
  while IFS= read -r -d '' r1; do
    local bn sample r2
    bn=$(basename "$r1")

    if [[ "$bn" =~ ^(.+)_R1\.f(ast)?q\.gz$ ]]; then
      sample="${BASH_REMATCH[1]}"
      r2="$INPUT_DIR/${sample}_R2.fastq.gz"
    elif [[ "$bn" =~ ^(.+)_1\.f(ast)?q\.gz$ ]]; then
      sample="${BASH_REMATCH[1]}"
      r2="$INPUT_DIR/${sample}_2.fastq.gz"
    else
      continue
    fi

    if [[ -f "$r2" ]]; then
      R1_FILES["$sample"]="$r1"
    else
      echo "[WARN] Skipping $bn because matching R2 file not found."
    fi
  done < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" \) -print0)
}

auto_discover_samples

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No paired-end samples found in $INPUT_DIR" >&2
  exit 1
fi

echo "[INFO] Found ${#R1_FILES[@]} sample(s):"
printf '  - %s\n' "${!R1_FILES[@]}" | sort

for sample in "${!R1_FILES[@]}"; do
  r1="${R1_FILES[$sample]}"
  r1_bn=$(basename "$r1")

  if [[ "$r1_bn" =~ _R1\.f(ast)?q\.gz$ ]]; then
    r2="$INPUT_DIR/${sample}_R2.fastq.gz"
  else
    r2="$INPUT_DIR/${sample}_2.fastq.gz"
  fi

  echo "[INFO] ===== Processing sample: $sample ====="

  # 1) Raw read QC
  fastqc -t "$THREADS" -o "$OUT_DIR/00_raw_qc" "$r1" "$r2" \
    >"$OUT_DIR/logs/${sample}.fastqc_raw.log" 2>&1

  # 2) Trimming/filtering
  trimmed_r1="$OUT_DIR/01_trimmed/${sample}_R1.trimmed.fastq.gz"
  trimmed_r2="$OUT_DIR/01_trimmed/${sample}_R2.trimmed.fastq.gz"

  fastp \
    -i "$r1" -I "$r2" \
    -o "$trimmed_r1" -O "$trimmed_r2" \
    -w "$THREADS" \
    --length_required "$MIN_LEN" \
    --detect_adapter_for_pe \
    --html "$OUT_DIR/01_trimmed/${sample}.fastp.html" \
    --json "$OUT_DIR/01_trimmed/${sample}.fastp.json" \
    >"$OUT_DIR/logs/${sample}.fastp.log" 2>&1

  # 3) QC on trimmed reads
  fastqc -t "$THREADS" -o "$OUT_DIR/02_trim_qc" "$trimmed_r1" "$trimmed_r2" \
    >"$OUT_DIR/logs/${sample}.fastqc_trimmed.log" 2>&1

  # 4) Assembly
  sample_assembly_dir="$OUT_DIR/03_assembly/$sample"
  mkdir -p "$sample_assembly_dir"

  spades.py \
    -1 "$trimmed_r1" -2 "$trimmed_r2" \
    -o "$sample_assembly_dir" \
    --careful -t "$THREADS" \
    >"$OUT_DIR/logs/${sample}.spades.log" 2>&1

  assembly_fasta="$sample_assembly_dir/contigs.fasta"
  if [[ ! -s "$assembly_fasta" ]]; then
    echo "[ERROR] Assembly output not found for $sample: $assembly_fasta" >&2
    exit 1
  fi

  # 5) Assembly QC
  quast.py "$assembly_fasta" \
    -o "$OUT_DIR/04_assembly_qc/quast_${sample}" \
    -t "$THREADS" \
    >"$OUT_DIR/logs/${sample}.quast.log" 2>&1

  if [[ "$RUN_BUSCO" == true ]]; then
    busco \
      -i "$assembly_fasta" \
      -o "busco_${sample}" \
      -m genome \
      -l "$BUSCO_LINEAGE" \
      -c "$THREADS" \
      --out_path "$OUT_DIR/04_assembly_qc" \
      >"$OUT_DIR/logs/${sample}.busco.log" 2>&1
  fi

  # 6) Gene prediction
  prodigal \
    -i "$assembly_fasta" \
    -o "$OUT_DIR/05_gene_prediction/${sample}.genes.gbk" \
    -a "$OUT_DIR/05_gene_prediction/${sample}.proteins.faa" \
    -d "$OUT_DIR/05_gene_prediction/${sample}.cds.fna" \
    -f gbk -p single \
    >"$OUT_DIR/logs/${sample}.prodigal.log" 2>&1

  echo "[INFO] Completed: $sample"
done

# Final summaries
multiqc -o "$OUT_DIR/00_raw_qc" "$OUT_DIR/00_raw_qc" >"$OUT_DIR/logs/multiqc_raw.log" 2>&1
multiqc -o "$OUT_DIR/02_trim_qc" "$OUT_DIR/02_trim_qc" >"$OUT_DIR/logs/multiqc_trimmed.log" 2>&1
multiqc -o "$OUT_DIR/04_assembly_qc" "$OUT_DIR/04_assembly_qc" >"$OUT_DIR/logs/multiqc_assembly.log" 2>&1

echo "[INFO] Pipeline finished successfully."
echo "[INFO] Output directory: $OUT_DIR"
