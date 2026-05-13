#!/usr/bin/env bash
# Lymphix — single-command entry point.
#
# Usage:
#     ./run.sh test                        # synthetic end-to-end on bundled mock data
#     ./run.sh --samplesheet x.csv         # local run on real data (Docker required)
#     ./run.sh dnanexus --project ID --samplesheet dx://...
#     ./run.sh --help                      # full Nextflow help

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# ---------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------
case "${1:-}" in
    test)
        # Quick smoke test — no Docker, no Nextflow required.
        # Runs the Python analysis layer on synthetic AIRR data.
        exec bash tests/test_smoke.sh
        ;;

    dnanexus)
        shift
        echo "[run.sh] Submitting to DNAnexus — see docs/DNANEXUS.md for prerequisites."
        exec dx run /applets/lymphix -i nextflow_pipeline_params="$*" --watch
        ;;

    --help|-h|help|"")
        cat <<EOF
Lymphix — BCR/TCR clonality from custom-panel NGS.

Usage:
  ./run.sh test                                     Synthetic end-to-end test
  ./run.sh --samplesheet SAMPLES.CSV [options]      Real-data run (needs Docker + Nextflow)
  ./run.sh dnanexus --samplesheet ...               DNAnexus deployment

Common options (passed through to Nextflow):
  --samplesheet PATH        CSV: sample_id,fastq_1,fastq_2,bam,umi_preset,expected_status
  --outdir PATH             Results directory (default: results/)
  --species human|mouse     Reference species (default: human)
  --umi_preset PRESET       none|twist|xgen_duplex|xgen_simplex|custom
  --total_input_reads N     For accurate background fraction in composition
  --filter_dups_in_bam      Strip flag-marked duplicates from BAM input

See README.md for the full sample-sheet schema and panel BED options.
EOF
        exit 0
        ;;
esac

# ---------------------------------------------------------------------
# Default: pass through to nextflow with Docker profile
# ---------------------------------------------------------------------
command -v nextflow >/dev/null || {
    echo "[run.sh] Nextflow not found. Install: curl -fsSL https://get.nextflow.io | bash"
    exit 1
}
command -v docker >/dev/null || {
    echo "[run.sh] WARNING: Docker not found. Use -profile singularity or install Docker Desktop."
}

exec nextflow run main.nf -profile docker "$@"
