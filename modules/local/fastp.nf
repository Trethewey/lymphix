process FASTP {
    tag        "${sample_id}"
    label      'process_medium'
    container  'ghcr.io/trethewey/lymphix/fastp:0.23.4'
    publishDir "${params.outdir}/${sample_id}/fastp", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: reads
    path  "${sample_id}.fastp.json", emit: json
    path  "${sample_id}.fastp.html", emit: html
    path  "versions.yml",            emit: versions

    script:
    """
    fastp \\
        -i ${r1} -I ${r2} \\
        -o ${sample_id}_R1.trim.fastq.gz \\
        -O ${sample_id}_R2.trim.fastq.gz \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 \\
        --length_required 50 \\
        --correction

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """
}
