# Prokaryotic WGS (Illumina Paired-End) Analysis Pipeline

This repository provides a ready-to-run Bash pipeline for **whole genome sequencing (WGS)** analysis of **prokaryotic Illumina paired-end data**, plus a **Docker image build/publish setup**.

## Pipeline steps included

1. **Raw read QC**: `FastQC` + `MultiQC`
2. **Read trimming/filtering**: `fastp`
3. **Genome assembly**: `SPAdes`
4. **Assembly/genome QC**: `QUAST` (+ optional `BUSCO`)
5. **Gene prediction**: `Prodigal`

## Files

- `illumina_prok_wgs_pipeline.sh` — analysis script
- `Dockerfile` — container image definition
- `environment.yml` — conda environment for bioinformatics tools

## Input requirements

Place paired FASTQ files in one directory using one of these naming styles:

- `<sample>_R1.fastq.gz` and `<sample>_R2.fastq.gz`
- `<sample>_1.fastq.gz` and `<sample>_2.fastq.gz`

## Local usage (without Docker)

```bash
./illumina_prok_wgs_pipeline.sh \
  -i ./raw_reads \
  -o ./analysis_out \
  -t 16 \
  --busco-lineage bacteria_odb10
```

Disable BUSCO if not installed:

```bash
./illumina_prok_wgs_pipeline.sh -i ./raw_reads -o ./analysis_out -t 16 --disable-busco
```

## Docker usage

### 1) Build image

```bash
docker build -t prok-wgs:latest .
```

### 2) Run pipeline in container

```bash
docker run --rm \
  -v "$(pwd)/raw_reads:/data/raw_reads" \
  -v "$(pwd)/analysis_out:/data/analysis_out" \
  prok-wgs:latest \
  -i /data/raw_reads \
  -o /data/analysis_out \
  -t 8 \
  --busco-lineage bacteria_odb10
```

## Publish image

### Option A: Docker Hub

```bash
# Login once
docker login

# Tag for your Docker Hub repo (replace USERNAME)
docker tag prok-wgs:latest USERNAME/prok-wgs:latest

# Push
docker push USERNAME/prok-wgs:latest
```

### Option B: GitHub Container Registry (GHCR)

```bash
# Use a GitHub PAT with write:packages scope
echo "$GITHUB_TOKEN" | docker login ghcr.io -u GITHUB_USERNAME --password-stdin

# Tag (replace OWNER and REPO)
docker tag prok-wgs:latest ghcr.io/OWNER/REPO/prok-wgs:latest

# Push
docker push ghcr.io/OWNER/REPO/prok-wgs:latest
```

## Required tools (if running outside Docker)

- fastqc
- multiqc
- fastp
- spades.py
- quast.py
- prodigal

Optional:
- busco
