/*
 * UMI_CONSENSUS — fgbio-based UMI grouping + consensus calling
 *
 * Supports three input scenarios:
 *   1. paired FASTQ + inline UMI read structure (TWIST, xGen)
 *   2. BAM with raw reads + UMI in RX tag (most common from upstream demux)
 *   3. BAM with reads but no UMI tags (preset=none — skip this module entirely)
 *
 * Outputs paired FASTQs with one consensus read per molecule.
 *
 * Read structures per preset (fgbio convention; M=molecular UMI, S=spacer, T=template):
 *   twist        5M2S+T 5M2S+T   (TWIST UMI Adapter System, 5 bp UMI + 2 bp spacer)
 *   xgen_duplex  9M+T   9M+T     (IDT xGen UMI-UDI, 9 bp UMI on each strand → duplex)
 *   xgen_simplex 9M+T   9M+T     (xGen UMI without duplex calling)
 *   custom       <user-supplied via --umi_read_structure>
 *
 * The duplex consensus caller (CallDuplexConsensusReads) requires xGen-style
 * dual-strand UMIs. Other presets use CallMolecularConsensusReads.
 */

process UMI_CONSENSUS {
    tag        "${sample_id}"
    label      'process_high'
    container  'ghcr.io/trethewey/lymphix/fgbio:0.1.0'
    publishDir "${params.outdir}/${sample_id}/umi_consensus", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(r1), path(r2), val(umi_preset)
    path  bwa_index_dir          // directory containing the bwa index files

    output:
    tuple val(sample_id), path("${sample_id}.consensus_R1.fastq.gz"),
                          path("${sample_id}.consensus_R2.fastq.gz"), emit: reads
    path  "${sample_id}.umi_metrics.tsv",                              emit: metrics
    path  "versions.yml",                                              emit: versions

    script:
    def read_struct = (umi_preset == 'twist')        ? '5M2S+T 5M2S+T' :
                      (umi_preset == 'xgen_duplex')  ? '9M+T 9M+T'     :
                      (umi_preset == 'xgen_simplex') ? '9M+T 9M+T'     :
                      params.umi_read_structure ?: '8M+T 8M+T'
    def is_duplex   = (umi_preset == 'xgen_duplex')
    def consensus_cmd = is_duplex ? 'CallDuplexConsensusReads' : 'CallMolecularConsensusReads'
    def min_reads   = is_duplex ? '1 1 0' : '1'
    def bwa_prefix  = "${bwa_index_dir}/genome"
    """
    # 1) FASTQ -> unmapped BAM with UMI extraction
    fgbio --tmp-dir=. FastqToBam \\
        --input ${r1} ${r2} \\
        --read-structures ${read_struct} \\
        --sample ${sample_id} --library ${sample_id} \\
        --output ${sample_id}.unmapped.bam

    # 2) Align (BWA-MEM) preserving UMI tag
    samtools fastq ${sample_id}.unmapped.bam | \\
        bwa mem -t ${task.cpus} -p ${bwa_prefix} - | \\
        fgbio --tmp-dir=. ZipperBams \\
            --unmapped ${sample_id}.unmapped.bam \\
            --ref ${bwa_prefix}.fa \\
            --output ${sample_id}.mapped.bam

    # 3) Group reads by UMI (adjacency, allows 1 edit). For duplex use 'paired'.
    fgbio --tmp-dir=. GroupReadsByUmi \\
        --input ${sample_id}.mapped.bam \\
        --strategy ${is_duplex ? 'paired' : 'adjacency'} \\
        --edits 1 \\
        --output ${sample_id}.grouped.bam \\
        --family-size-histogram ${sample_id}.umi_metrics.tsv

    # 4) Call consensus (duplex or molecular)
    fgbio --tmp-dir=. ${consensus_cmd} \\
        --input ${sample_id}.grouped.bam \\
        --output ${sample_id}.consensus.bam \\
        --min-reads ${min_reads}

    # 5) Filter consensus (keep families with sufficient support)
    fgbio --tmp-dir=. FilterConsensusReads \\
        --input ${sample_id}.consensus.bam \\
        --output ${sample_id}.consensus.filtered.bam \\
        --ref ${bwa_prefix}.fa \\
        --min-reads ${min_reads} \\
        --max-read-error-rate 0.05 \\
        --min-base-quality 30

    # 6) Back to FASTQ for TRUST4
    samtools collate -O ${sample_id}.consensus.filtered.bam | \\
        samtools fastq -1 ${sample_id}.consensus_R1.fastq.gz \\
                       -2 ${sample_id}.consensus_R2.fastq.gz \\
                       -0 /dev/null -s /dev/null -n -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio:    \$(fgbio --version 2>&1 | head -1 | sed 's/^.*Version: //')
        bwa:      \$(bwa 2>&1 | grep "Version" | sed 's/Version: //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
