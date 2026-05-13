FROM ubuntu:22.04
LABEL org.opencontainers.image.source="https://github.com/trethewey/lymphix"
LABEL org.opencontainers.image.description="TRUST4 1.0.13 with bundled human + mouse IMGT references"

ARG TRUST4_VERSION=1.0.13
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates git wget zlib1g-dev libbz2-dev \
        liblzma-dev libcurl4-openssl-dev libncurses5-dev \
        python3 python3-pip pigz \
        samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN git clone --depth 1 --branch v${TRUST4_VERSION} https://github.com/liulab-dfci/TRUST4.git trust4 \
    && cd trust4 && make -j4

ENV PATH=/opt/trust4:$PATH

# Bundle the IMGT + hg38/GRCm38 references that ship in the repo
# (hg38_bcrtcr.fa, human_IMGT+C.fa, GRCm38_bcrtcr.fa, mouse_IMGT+C.fa)
# Already present at /opt/trust4/*.fa after the clone.

RUN pip3 install --no-cache-dir biopython==1.83 numpy==1.26.4

COPY bin/simulate_repertoire.py /opt/scripts/simulate_repertoire.py
COPY bin/airr_to_fasta.py       /opt/scripts/airr_to_fasta.py
RUN chmod +x /opt/scripts/*.py
ENV PATH=/opt/scripts:$PATH

WORKDIR /work
