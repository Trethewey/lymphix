#!/usr/bin/env bash
# Lymphix — no-install entry point.
#
# If Lymphix is installed (pip install .), use the `lymphix` command instead:
#     lymphix --samplesheet samples.csv --outdir results/
#
# This script does the same thing without installing anything:
#     ./lymphix.sh test                        # analysis-layer smoke test on mock AIRR
#     ./lymphix.sh --samplesheet x.csv         # local run on real data (Docker required)
#     ./lymphix.sh dnanexus --project ID --samplesheet dx://...
#     ./lymphix.sh --help                      # full Nextflow help

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# ---------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------
case "${1:-}" in
    test)
        # Exercises the Python analysis layer on synthetic AIRR fixtures.
        # Does NOT run TRUST4, IgBLAST, Docker or Nextflow.
        exec bash tests/test_smoke.sh
        ;;

    dnanexus)
        shift
        echo "[lymphix] Submitting to DNAnexus — see docs/DNANEXUS.md for prerequisites."
        exec dx run /applets/lymphix -i nextflow_pipeline_params="$*" --watch
        ;;

    --help|-h|help|"")
        cat <<EOF
Lymphix — BCR/TCR clonality from custom-panel NGS.

Installed usage (after \`pip install .\`):
  lymphix --samplesheet SAMPLES.CSV --outdir results/

No-install usage:
  ./lymphix.sh test                                 Analysis-layer smoke test (mock AIRR)
  ./lymphix.sh --samplesheet SAMPLES.CSV [options]  Real-data run (needs Docker + Nextflow)
  ./lymphix.sh dnanexus --samplesheet ...           DNAnexus deployment

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
    echo "[lymphix] Nextflow not found. Install: curl -fsSL https://get.nextflow.io | bash"
    exit 1
}
command -v docker >/dev/null || {
    echo "[lymphix] WARNING: Docker not found. Use -profile singularity or install Docker Desktop."
}

exec nextflow run main.nf -profile docker "$@"
