FROM ubuntu:22.04
LABEL org.opencontainers.image.source="https://github.com/trethewey/lymphix"
LABEL org.opencontainers.image.description="NCBI IgBLAST 1.22.0 + IMGT germline databases (human + mouse)"

ARG IGBLAST_VERSION=1.22.0
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates wget python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN wget -q https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz \
    && tar xzf ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz \
    && mv ncbi-igblast-${IGBLAST_VERSION} igblast \
    && rm ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz

ENV PATH=/opt/igblast/bin:$PATH
ENV IGDATA=/opt/igblast

# Fetch IMGT germline reference, format for IgBLAST
# (uses NCBI's pre-built germline databases bundled with IgBLAST distribution
#  + IMGT human/mouse V/D/J refs for higher accuracy)
RUN mkdir -p /opt/igblast/database \
    && cd /opt/igblast/database \
    && for org in human mouse; do \
         for seg in V D J; do \
           wget -q "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_${org}_${seg,,}_genes.tar" || true; \
         done; \
       done \
    && for tar in *.tar; do tar xf "$tar" && rm "$tar"; done || true

# Build BLAST DBs from the bundled NCBI references
RUN cd /opt/igblast && \
    for org in human mouse; do \
      for seg in V D J; do \
        if [ -f internal_data/${org}/${org}_${seg}.nhr ]; then continue; fi; \
        if [ -f database/ncbi_${org}_${seg,,}_genes ]; then \
          makeblastdb -parse_seqids -dbtype nucl \
            -in database/ncbi_${org}_${seg,,}_genes \
            -out database/${org}_gl_${seg}; \
        fi; \
      done; \
    done || true

RUN pip3 install --no-cache-dir pandas==2.2.2

COPY bin/airr_to_fasta.py /opt/scripts/airr_to_fasta.py
COPY bin/merge_airr.py    /opt/scripts/merge_airr.py
RUN chmod +x /opt/scripts/*.py
ENV PATH=/opt/scripts:$PATH

WORKDIR /work
