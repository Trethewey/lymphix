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
 *
 * REFERENCE HANDLING
 * ------------------
 * UMI grouping needs alignment positions, so this is the one module that
 * requires a BWA-indexed reference. It used to assume the index basename was
 * 'genome' and never checked that --bwa_index had been given at all. Both
 * failure modes were poor: with no --bwa_index the reference path resolved to
 * nothing and bwa reported only "fail to locate the index files", naming no
 * pipeline parameter; with an index present under any other basename the run
 * failed just as opaquely, several expensive steps in.
 *
 * The index is now resolved and validated before any work is done. The
 * basename is discovered from the staged directory rather than assumed, so
 * both `bwa index reference.fa` and `bwa index -p <prefix> reference.fa`
 * layouts work.
 */

process UMI_CONSENSUS {
    tag        "${sample_id}"
    label      'process_high'
    container  'ghcr.io/trethewey/lymphix/fgbio:0.1.0'
    publishDir "${params.outdir}/${sample_id}/umi_consensus", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(r1), path(r2), val(umi_preset)
    path  bwa_index              // directory holding the reference FASTA + its bwa index

    output:
    tuple val(sample_id), path("${sample_id}.consensus_R1.fastq.gz"),
                          path("${sample_id}.consensus_R2.fastq.gz"), emit: reads
    path  "${sample_id}.umi_metrics.tsv",                              emit: metrics
    path  "versions.yml",                                              emit: versions

    script:
    // Fail here, not inside bwa. This runs when the task is created, so the
    // message names the missing parameter before anything is staged or run.
    if (!params.bwa_index) {
        error """
        Sample '${sample_id}' uses umi_preset='${umi_preset}', which needs a BWA-indexed
        reference for UMI grouping and consensus filtering, but --bwa_index was not supplied.

          Pass  --bwa_index <dir>  where <dir> contains the reference FASTA together with
          its bwa index (.amb .ann .bwt .pac .sa), its .fai and its sequence dictionary.
          It must be a directory: Nextflow stages only what it is given, so pointing at
          the FASTA alone would leave the index files behind.

          If this sample has no UMIs, set umi_preset=none for it instead.
        """.stripIndent()
    }

    // Optional escape hatch for a directory that holds more than one index.
    def pinned_prefix = params.containsKey('bwa_index_prefix') ? (params.bwa_index_prefix ?: '') : ''

    def read_struct = (umi_preset == 'twist')        ? '5M2S+T 5M2S+T' :
                      (umi_preset == 'xgen_duplex')  ? '9M+T 9M+T'     :
                      (umi_preset == 'xgen_simplex') ? '9M+T 9M+T'     :
                      params.umi_read_structure ?: '8M+T 8M+T'
    def is_duplex   = (umi_preset == 'xgen_duplex')
    def consensus_cmd = is_duplex ? 'CallDuplexConsensusReads' : 'CallMolecularConsensusReads'
    def min_reads   = is_duplex ? '1 1 0' : '1'
    """
    set -o pipefail

    # -----------------------------------------------------------------
    # 0) Resolve and validate the BWA index
    # -----------------------------------------------------------------
    if [ ! -d "${bwa_index}" ]; then
        echo "ERROR: --bwa_index must be a directory containing the reference FASTA and" >&2
        echo "       its bwa index files. Got: '${bwa_index}'" >&2
        exit 1
    fi

    PINNED_PREFIX='${pinned_prefix}'
    if [ -n "\${PINNED_PREFIX}" ]; then
        BWA_PREFIX="${bwa_index}/\${PINNED_PREFIX}"
    else
        # bwa writes <prefix>.bwt, so the .bwt file identifies the index
        # whichever way it was built.
        BWT_FILES=\$(find -L ${bwa_index} -maxdepth 2 -name '*.bwt' | sort)
        N_BWT=\$(printf '%s\\n' "\${BWT_FILES}" | grep -c . || true)

        if [ "\${N_BWT}" -eq 0 ]; then
            echo "ERROR: no BWA index (no *.bwt file) found under --bwa_index '${bwa_index}'." >&2
            echo "       Build one with: bwa index <reference.fa>" >&2
            exit 1
        fi
        if [ "\${N_BWT}" -gt 1 ]; then
            echo "ERROR: --bwa_index '${bwa_index}' holds \${N_BWT} BWA indexes:" >&2
            printf '         %s\\n' \${BWT_FILES} >&2
            echo "       Choose one with --bwa_index_prefix <basename>." >&2
            exit 1
        fi
        BWA_PREFIX="\${BWT_FILES%.bwt}"
    fi

    for ext in amb ann bwt pac sa; do
        if [ ! -s "\${BWA_PREFIX}.\${ext}" ]; then
            echo "ERROR: incomplete BWA index — \${BWA_PREFIX}.\${ext} is missing or empty." >&2
            echo "       Re-run: bwa index <reference.fa>" >&2
            exit 1
        fi
    done

    # fgbio needs the reference FASTA itself, not just the bwa index. With
    # 'bwa index reference.fa' the prefix *is* the FASTA; with 'bwa index -p'
    # it is a bare basename, so try the usual extensions too.
    REF_FASTA=""
    for candidate in "\${BWA_PREFIX}" "\${BWA_PREFIX}.fa" "\${BWA_PREFIX}.fasta" "\${BWA_PREFIX}.fna"; do
        if [ -s "\${candidate}" ]; then REF_FASTA="\${candidate}"; break; fi
    done
    if [ -z "\${REF_FASTA}" ]; then
        echo "ERROR: found the BWA index '\${BWA_PREFIX}' but not its reference FASTA." >&2
        echo "       Expected \${BWA_PREFIX}, or that name with .fa/.fasta/.fna." >&2
        exit 1
    fi

    # ZipperBams and FilterConsensusReads both read the sequence dictionary.
    # Checking now saves discovering it after the alignment and grouping.
    if [ ! -s "\${REF_FASTA}.fai" ]; then
        echo "ERROR: \${REF_FASTA}.fai is missing. Create it with: samtools faidx \${REF_FASTA}" >&2
        exit 1
    fi
    if [ ! -s "\${REF_FASTA%.*}.dict" ] && [ ! -s "\${REF_FASTA}.dict" ]; then
        echo "ERROR: no sequence dictionary beside \${REF_FASTA}." >&2
        echo "       Create it with: samtools dict -o \${REF_FASTA%.*}.dict \${REF_FASTA}" >&2
        exit 1
    fi

    echo "[umi_consensus] ${sample_id}: bwa index \${BWA_PREFIX}, reference \${REF_FASTA}"

    # -----------------------------------------------------------------
    # 1) FASTQ -> unmapped BAM with UMI extraction
    # -----------------------------------------------------------------
    fgbio --tmp-dir=. FastqToBam \\
        --input ${r1} ${r2} \\
        --read-structures ${read_struct} \\
        --sample ${sample_id} --library ${sample_id} \\
        --output ${sample_id}.unmapped.bam

    # 2) Align (BWA-MEM) preserving UMI tag
    samtools fastq ${sample_id}.unmapped.bam | \\
        bwa mem -t ${task.cpus} -p "\${BWA_PREFIX}" - | \\
        fgbio --tmp-dir=. ZipperBams \\
            --unmapped ${sample_id}.unmapped.bam \\
            --ref "\${REF_FASTA}" \\
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
        --ref "\${REF_FASTA}" \\
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
