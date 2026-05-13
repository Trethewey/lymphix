#!/usr/bin/env bash
# DNAnexus entry point. Invoked by the platform when the applet runs.
#
# This wrapper exists for the "applet" build of the pipeline. For native
# Nextflow execution on DNAnexus, prefer:
#     dx build-nextflow-applet . --destination /applets/
# which avoids this wrapper entirely and lets the platform schedule each
# process as its own DX job.

set -euo pipefail

main() {
    echo "Starting BCR/TCR clonality pipeline on DNAnexus"

    # Download inputs declared in dxapp.json
    dx download "${samplesheet}" -o samplesheet.csv

    # Install nextflow if not present
    if ! command -v nextflow >/dev/null 2>&1; then
        curl -fsSL https://get.nextflow.io | bash
        sudo mv nextflow /usr/local/bin/
    fi

    # Run the pipeline
    nextflow run /pipeline/main.nf \
        -profile docker,dnanexus \
        --samplesheet samplesheet.csv \
        --species "${species:-human}" \
        --receptor "${receptor:-all}" \
        --igh_mutated_cutoff "${igh_mutated_cutoff:-98.0}" \
        --min_clone_count "${min_clone_count:-2}" \
        --outdir results/ \
        ${nextflow_pipeline_params:-}

    # Upload outputs
    mapfile -t result_files < <(find results/ -type f ! -name '*.fastq.gz')
    mapfile -t report_files < <(find results/ -name '*.report.html')

    results_ids=()
    for f in "${result_files[@]}"; do
        id=$(dx upload "$f" --brief)
        results_ids+=("$id")
        dx-jobutil-add-output results "$id" --class=array:file
    done

    for f in "${report_files[@]}"; do
        id=$(dx upload "$f" --brief)
        dx-jobutil-add-output reports "$id" --class=array:file
    done
}
