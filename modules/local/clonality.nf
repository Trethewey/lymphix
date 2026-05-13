process CLONALITY {
    tag        "${sample_id}"
    label      'process_low'
    container  'ghcr.io/trethewey/lymphix/clonality:0.1.0'
    publishDir "${params.outdir}/${sample_id}/clonality", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(trust4_airr), path(igblast_airr)

    output:
    tuple val(sample_id), path("${sample_id}.metrics.json"),   emit: metrics
    tuple val(sample_id), path("${sample_id}.clonotypes.tsv"), emit: clonotypes
    tuple val(sample_id), path("${sample_id}.top_clones.tsv"), emit: top_clones
    path  "versions.yml", emit: versions

    script:
    """
    clonality_metrics.py \\
        --sample-id ${sample_id} \\
        --trust4-airr ${trust4_airr} \\
        --igblast-airr ${igblast_airr} \\
        --min-clone-count ${params.min_clone_count} \\
        --igh-mutated-cutoff ${params.igh_mutated_cutoff} \\
        --out-metrics ${sample_id}.metrics.json \\
        --out-clonotypes ${sample_id}.clonotypes.tsv \\
        --out-top ${sample_id}.top_clones.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    END_VERSIONS
    """
}
