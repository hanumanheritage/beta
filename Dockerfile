FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="prok-wgs-illumina-pipeline" \
      org.opencontainers.image.description="Docker image for prokaryotic Illumina paired-end WGS QC, trimming, assembly, assembly QC and gene prediction" \
      org.opencontainers.image.licenses="MIT"

USER root
WORKDIR /opt/pipeline

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create -y -n prok-wgs -f /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER illumina_prok_wgs_pipeline.sh /usr/local/bin/illumina_prok_wgs_pipeline.sh
RUN chmod +x /usr/local/bin/illumina_prok_wgs_pipeline.sh

USER $MAMBA_USER

ENV PATH=/opt/conda/envs/prok-wgs/bin:$PATH

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/usr/local/bin/illumina_prok_wgs_pipeline.sh"]
CMD ["--help"]
