/*
 * BAM_TO_FASTQ — convert input BAM to paired FASTQ
 *
 * Used when the user supplies BAM input that needs UMI consensus calling
 * (so we can't pass the BAM straight to TRUST4). Read order is preserved
 * by collate. UMI tag (RX) is propagated to read comment if present so
 * downstream tools can recover it.
 */
process BAM_TO_FASTQ {
    tag        "${sample_id}"
    label      'process_medium'
    container  'ghcr.io/trethewey/lymphix/fgbio:0.1.0'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.R1.fastq.gz"), path("${sample_id}.R2.fastq.gz"), emit: reads
    path  "versions.yml", emit: versions

    script:
    """
    samtools collate -O -@ ${task.cpus} ${bam} | \\
        samtools fastq -@ ${task.cpus} \\
            -1 ${sample_id}.R1.fastq.gz \\
            -2 ${sample_id}.R2.fastq.gz \\
            -0 /dev/null -s /dev/null \\
            -T RX,MI -n -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
