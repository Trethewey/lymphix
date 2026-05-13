#!/usr/bin/env bash
# fetch_public_data.sh
#
# Download small, public test datasets for end-to-end validation.
# Requires: sra-toolkit (fastq-dump / fasterq-dump) or wget.
#
# Edit the accessions below for your validation scenario. Defaults are
# CONSERVATIVE PLACEHOLDERS — verify the accession describes the sample
# you want before relying on truth labels.
#
set -euo pipefail

OUTDIR="${1:-test_data_real}"
mkdir -p "$OUTDIR"

declare -A SAMPLES=(
    # Monoclonal TRB+ T-cell line — replace with verified accession for your
    # validation needs. Example T-ALL cell-line WGS available on SRA / ENA.
    ["JURKAT_TRB_clonal"]="SRRXXXXXXX"
    # Polyclonal germline donor WGS (e.g. NA12878 Platinum Genomes subset)
    ["NA12878_polyclonal"]="SRRYYYYYYY"
)

MAX_SPOTS=500000   # ~75 Mbp at 2x150 — enough for clonality QC

command -v fasterq-dump >/dev/null 2>&1 || {
    echo "ERROR: sra-toolkit (fasterq-dump) not on PATH." >&2
    echo "Install: https://github.com/ncbi/sra-tools/wiki" >&2
    exit 1
}

for name in "${!SAMPLES[@]}"; do
    srr="${SAMPLES[$name]}"
    echo "[fetch] $name -> $srr (first $MAX_SPOTS spots)"
    fasterq-dump --split-files --threads 4 \
                 --maxspotid "$MAX_SPOTS" \
                 -O "$OUTDIR" "$srr"
    mv "$OUTDIR/${srr}_1.fastq" "$OUTDIR/${name}_R1.fastq"
    mv "$OUTDIR/${srr}_2.fastq" "$OUTDIR/${name}_R2.fastq"
    pigz -f "$OUTDIR/${name}_R1.fastq" "$OUTDIR/${name}_R2.fastq"
done

# Write a samplesheet
{
    echo "sample_id,fastq_1,fastq_2"
    for name in "${!SAMPLES[@]}"; do
        echo "${name},${OUTDIR}/${name}_R1.fastq.gz,${OUTDIR}/${name}_R2.fastq.gz"
    done
} > "$OUTDIR/samplesheet_real.csv"

echo "[fetch] done. Samplesheet: $OUTDIR/samplesheet_real.csv"
