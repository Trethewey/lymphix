FROM ubuntu:22.04
LABEL org.opencontainers.image.source="https://github.com/trethewey/lymphix"
LABEL org.opencontainers.image.description="NCBI IgBLAST 1.22.0 + IMGT germline databases (human + mouse)"

ARG IGBLAST_VERSION=1.22.0

# IMGT/GENE-DB reference set, without IMGT gaps, functional + ORF + in-frame
# pseudogenes. This URL is unversioned: IMGT publishes one current release, so
# the image records the checksum of whatever it fetched (see
# database/imgt_reference.provenance) rather than pretending to pin a version.
# Override with --build-arg to build against a locally mirrored copy.
ARG IMGT_REFERENCE_URL=https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP

ENV DEBIAN_FRONTEND=noninteractive

# mawk is named explicitly because build_germline_db.sh parses the IMGT
# reference with awk; relying on it being pulled in as a transitive
# dependency of the base image is how a build breaks quietly later.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates wget mawk python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN wget -q https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz \
    && tar xzf ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz \
    && mv ncbi-igblast-${IGBLAST_VERSION} igblast \
    && rm ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz

ENV PATH=/opt/igblast/bin:$PATH
ENV IGDATA=/opt/igblast
ENV IGBLASTDB=/opt/igblast/database

# Build the germline V/D/J databases. The script is strict POSIX sh and runs
# under `set -eu` with no failure suppression, so any problem with the
# reference set stops the build here rather than surfacing months later as a
# run in which every IGHV call is null. See the header of the script for why
# the reference comes from IMGT and not the NCBI germline downloads.
#
# Docker exposes ARG values to RUN as environment variables, which is how
# IMGT_REFERENCE_URL reaches the script.
COPY containers/igblast/build_germline_db.sh /opt/scripts/build_germline_db.sh
RUN chmod +x /opt/scripts/build_germline_db.sh \
    && /opt/scripts/build_germline_db.sh

RUN pip3 install --no-cache-dir pandas==2.2.2

# Both are standalone: they import pandas and nothing from bin/. If either ever
# imports lymphix_common (or any other bin/ module), that module has to be
# copied here as well or the task fails at run time with an ImportError.
COPY bin/airr_to_fasta.py /opt/scripts/airr_to_fasta.py
COPY bin/merge_airr.py    /opt/scripts/merge_airr.py
RUN chmod +x /opt/scripts/*.py
ENV PATH=/opt/scripts:$PATH

WORKDIR /work
