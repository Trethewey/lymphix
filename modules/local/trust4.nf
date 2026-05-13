/*
 * TRUST4 — V(D)J assembly from paired FASTQ OR BAM
 *
 * Two input modes:
 *   - FASTQ pair  : standard `-1 R1.fq -2 R2.fq`
 *   - BAM (sorted): `-b sample.bam` (TRUST4 selects V(D)J reads by alignment)
 *
 * BAM input is only used when UMI preset = 'none' (no consensus calling
 * needed). With a UMI preset, BAM is first converted to FASTQ + consensus
 * called upstream, and the consensus FASTQs land here.
 */

process TRUST4 {
    tag        "${sample_id}"
    label      'process_high'
    container  'ghcr.io/trethewey/lymphix/trust4:1.0.13'
    publishDir "${params.outdir}/${sample_id}/trust4", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), val(input_kind), path(read1_or_bam), path(read2_optional)
    val   species

    output:
    tuple val(sample_id), path("${sample_id}_airr.tsv"),  emit: airr
    tuple val(sample_id), path("${sample_id}_report.tsv"), emit: report
    path  "${sample_id}_*"
    path  "versions.yml", emit: versions

    script:
    def ref_bcrtcr = species == 'mouse' ? '/opt/trust4/GRCm38_bcrtcr.fa' : '/opt/trust4/hg38_bcrtcr.fa'
    def ref_imgt   = species == 'mouse' ? '/opt/trust4/mouse_IMGT+C.fa'  : '/opt/trust4/human_IMGT+C.fa'
    def filter_dups = params.filter_dups_in_bam ?: false
    def input_args
    def pre_cmd = ''
    if (input_kind == 'bam') {
        if (filter_dups) {
            pre_cmd = "samtools view -@ ${task.cpus} -F 1024 -b ${read1_or_bam} > _filtered.bam && samtools index _filtered.bam"
            input_args = "-b _filtered.bam"
        } else {
            input_args = "-b ${read1_or_bam}"
        }
    } else {
        input_args = "-1 ${read1_or_bam} -2 ${read2_optional}"
    }
    """
    ${pre_cmd}
    run-trust4 \\
        -f ${ref_bcrtcr} \\
        --ref ${ref_imgt} \\
        ${input_args} \\
        -o ${sample_id} \\
        -t ${task.cpus}

    # Make sure the AIRR table exists even if no clonotypes were detected
    if [[ ! -s ${sample_id}_airr.tsv ]]; then
        echo -e "sequence_id\\tsequence\\tv_call\\td_call\\tj_call\\tjunction\\tjunction_aa\\tlocus\\tproductive\\tconsensus_count\\tduplicate_count" > ${sample_id}_airr.tsv
    fi
    ls -la ${sample_id}*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | head -1 | sed 's/.*v//' || echo "1.0.13")
        input_kind: ${input_kind}
    END_VERSIONS
    """
}
