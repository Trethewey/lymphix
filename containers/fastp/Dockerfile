FROM mambaorg/micromamba:1.5.8-jammy
LABEL org.opencontainers.image.source="https://github.com/trethewey/lymphix"
LABEL org.opencontainers.image.description="fastp 0.23.4 for BCR/TCR clonality pipeline"

USER root
RUN micromamba install -y -n base -c bioconda -c conda-forge \
        fastp=0.23.4 \
    && micromamba clean -ya

ENV PATH=/opt/conda/bin:$PATH
USER mambauser
