process REPORT {
    tag        "${sample_id}"
    label      'process_low'
    container  'ghcr.io/trethewey/lymphix/clonality:0.1.0'
    publishDir "${params.outdir}/${sample_id}", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(metrics_json), path(clonotypes_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.report.html")

    script:
    """
    generate_report.py \\
        --sample-id ${sample_id} \\
        --metrics ${metrics_json} \\
        --clonotypes ${clonotypes_tsv} \\
        --out ${sample_id}.report.html
    """
}
