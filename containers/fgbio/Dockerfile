FROM eclipse-temurin:17-jre-jammy
LABEL org.opencontainers.image.source="https://github.com/trethewey/lymphix"
LABEL org.opencontainers.image.description="fgbio + bwa + samtools for UMI consensus calling (TWIST, IDT xGen UMI-UDI)"

ARG FGBIO_VERSION=2.4.0
ARG BWA_VERSION=0.7.18
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates wget unzip bzip2 \
        build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
        libncurses5-dev libssl-dev pigz \
    && rm -rf /var/lib/apt/lists/*

# --- fgbio (Java fat-jar) ----------------------------------------------------
RUN mkdir -p /opt/fgbio && cd /opt/fgbio && \
    wget -q https://github.com/fulcrumgenomics/fgbio/releases/download/${FGBIO_VERSION}/fgbio-${FGBIO_VERSION}.jar -O fgbio.jar && \
    printf '#!/usr/bin/env bash\nexec java -Xmx${FGBIO_MEM:-4g} -jar /opt/fgbio/fgbio.jar "$@"\n' > /usr/local/bin/fgbio && \
    chmod +x /usr/local/bin/fgbio

# --- bwa ---------------------------------------------------------------------
RUN cd /tmp && \
    wget -q https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 && \
    tar xjf bwa-${BWA_VERSION}.tar.bz2 && cd bwa-${BWA_VERSION} && \
    make && install -m 0755 bwa /usr/local/bin/bwa && \
    cd / && rm -rf /tmp/bwa-${BWA_VERSION}*

# --- samtools (from htslib release) -----------------------------------------
ARG HTSLIB_VERSION=1.20
RUN cd /tmp && \
    wget -q https://github.com/samtools/samtools/releases/download/${HTSLIB_VERSION}/samtools-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf samtools-${HTSLIB_VERSION}.tar.bz2 && cd samtools-${HTSLIB_VERSION} && \
    ./configure --without-curses && make && make install && \
    cd / && rm -rf /tmp/samtools-${HTSLIB_VERSION}*

WORKDIR /work
