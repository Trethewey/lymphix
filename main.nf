#!/usr/bin/env nextflow
/*
 * =====================================================================
 * BCR/TCR Clonality Pipeline
 * TRUST4 (V(D)J assembly) + IgBLAST (re-annotation) + clonality metrics
 *
 * Input modes:
 *   1. FASTQ + umi_preset=none           → fastp → TRUST4 (FASTQ)
 *   2. FASTQ + umi_preset=twist/xgen_*   → fastp → UMI_CONSENSUS → TRUST4
 *   3. BAM   + umi_preset=none           → TRUST4 (BAM, -b)
 *
 * Sample sheet columns:
 *   sample_id, fastq_1, fastq_2, bam, umi_preset, expected_status
 *
 * umi_preset:       none | twist | xgen_duplex | xgen_simplex | custom
 * expected_status:  clonal | polyclonal | negative   (optional, drives QC)
 *
 * BAM + UMI is rejected with a clear error — convert BAM → FASTQ first
 * (preserving the inline UMI bases) and rerun as case 2.
 * =====================================================================
 */
nextflow.enable.dsl = 2

include { FASTP                  } from './modules/local/fastp.nf'
include { UMI_CONSENSUS          } from './modules/local/umi_consensus.nf'
include { TRUST4                 } from './modules/local/trust4.nf'
include { IGBLAST                } from './modules/local/igblast.nf'
include { CLONALITY              } from './modules/local/clonality.nf'
include { REPORT                 } from './modules/local/report.nf'
include { SIMULATE_TEST_DATA     } from './modules/local/simulate_test_data.nf'


def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf -profile docker --samplesheet samples.csv --outdir results/

    Sample sheet (CSV):
      sample_id,fastq_1,fastq_2,bam,umi_preset,expected_status

      umi_preset: none | twist | xgen_duplex | xgen_simplex | custom
      expected_status: clonal | polyclonal | negative (optional QC assertion)

    Optional flags:
      --species              human | mouse                      (default: human)
      --receptor             all | bcr | tcr                    (default: all)
      --umi_preset           Default preset if a row leaves it blank (default: none)
      --umi_read_structure   fgbio read-structure (for umi_preset=custom)
      --bwa_index            Path to BWA index dir (required if any sample uses a UMI preset)
      --igh_mutated_cutoff   V-identity %% for IGHV-unmutated     (default: 98.0)
      --min_clone_count      Min reads per clonotype             (default: 2)
      --collapse_clonotypes  Collapse TRUST4 assembly-variant rows into clones
                             before any metric is computed       (default: false)
      --collapse_key         locus_junction_nt | locus_junction_nt_hamming1
                                                                 (default: locus_junction_nt)
      --collapse_minor_fraction
                             Abundance gate for the hamming1 key (default: 0.02)

    Profiles:
      -profile docker | singularity | test | dnanexus
    """.stripIndent()
}


def normalise_row(row) {
    def preset = (row.umi_preset && row.umi_preset.trim()) ? row.umi_preset.trim()
                                                           : params.umi_preset ?: 'none'
    def has_fastq = row.fastq_1 && row.fastq_1.trim()
    def has_bam   = row.bam && row.bam.trim()
    if (!has_fastq && !has_bam) {
        exit 1, "ERROR: sample ${row.sample_id} has neither fastq_1 nor bam"
    }
    if (has_bam && preset != 'none') {
        exit 1, "ERROR: sample ${row.sample_id} is BAM + umi_preset='${preset}'. " +
                "Convert BAM → FASTQ preserving inline UMI bases first, then rerun."
    }
    return [
        sample_id       : row.sample_id,
        fastq_1         : has_fastq ? row.fastq_1.trim() : null,
        fastq_2         : has_fastq ? (row.fastq_2 ?: '').trim() : null,
        bam             : has_bam   ? row.bam.trim() : null,
        umi_preset      : preset,
        expected_status : (row.expected_status ?: '').trim() ?: null,
        kind            : has_bam ? 'bam'
                                  : (preset == 'none' ? 'fastq_no_umi' : 'fastq_umi'),
    ]
}


workflow {
    if (params.help) { helpMessage(); exit 0 }

    // ---------- Build the input channel ---------------------------------
    if (params.containsKey('run_simulate_test_data') && params.run_simulate_test_data) {
        SIMULATE_TEST_DATA()
        ch_rows = SIMULATE_TEST_DATA.out.samplesheet.splitCsv(header: true)
            .map { row -> normalise_row(row + [bam: null, umi_preset: 'none', expected_status: null]) }
    } else {
        if (!params.samplesheet) {
            exit 1, "ERROR: --samplesheet is required (or use -profile test)"
        }
        ch_rows = Channel.fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row -> normalise_row(row) }
    }

    // Capture expected_status for QC assertions
    ch_expected = ch_rows.map { r -> tuple(r.sample_id, r.expected_status) }

    // ---------- Route 1: FASTQ no UMI -----------------------------------
    ch_fastq_no_umi = ch_rows
        .filter { it.kind == 'fastq_no_umi' }
        .map { tuple(it.sample_id, file(it.fastq_1), file(it.fastq_2)) }

    // ---------- Route 2: FASTQ + UMI ------------------------------------
    ch_fastq_umi = ch_rows
        .filter { it.kind == 'fastq_umi' }
        .map { tuple(it.sample_id, file(it.fastq_1), file(it.fastq_2), it.umi_preset) }

    // ---------- Route 3: BAM no UMI -------------------------------------
    ch_bam = ch_rows
        .filter { it.kind == 'bam' }
        .map { tuple(it.sample_id, 'bam', file(it.bam), file('NO_FILE_R2')) }

    // ---------- Trim FASTQ inputs (routes 1 & 2) ------------------------
    ch_fastp_input = ch_fastq_no_umi.mix(ch_fastq_umi.map { sid, r1, r2, _p -> tuple(sid, r1, r2) })
    FASTP(ch_fastp_input)
    ch_fastp_out = FASTP.out.reads      // tuple(sample_id, r1_trim, r2_trim)

    // Split fastp output into UMI and no-UMI lanes by re-joining the original row info
    ch_umi_lanes = ch_rows
        .filter { it.kind == 'fastq_umi' }
        .map { tuple(it.sample_id, it.umi_preset) }
    ch_fastp_for_umi = ch_fastp_out
        .join(ch_umi_lanes)
        .map { sid, r1, r2, preset -> tuple(sid, r1, r2, preset) }
    ch_fastp_no_umi  = ch_fastp_out
        .join(ch_rows.filter { it.kind == 'fastq_no_umi' }.map { tuple(it.sample_id, true) })
        .map { sid, r1, r2, _flag -> tuple(sid, r1, r2) }

    // ---------- UMI consensus calling (route 2) -------------------------
    bwa_dir = params.bwa_index ? file(params.bwa_index)
                               : file("${projectDir}/assets/empty_bwa_index", checkIfExists: false)
    UMI_CONSENSUS(ch_fastp_for_umi, bwa_dir)
    ch_consensus_fastq = UMI_CONSENSUS.out.reads

    // ---------- TRUST4 (routes 1+2 as FASTQ, route 3 as BAM) ------------
    ch_trust4_fastq = ch_fastp_no_umi.mix(ch_consensus_fastq)
        .map { sid, r1, r2 -> tuple(sid, 'fastq', r1, r2) }
    TRUST4(ch_trust4_fastq.mix(ch_bam), params.species)

    // ---------- Re-annotation + metrics + report -----------------------
    IGBLAST(TRUST4.out.airr, params.species)
    CLONALITY(TRUST4.out.airr.join(IGBLAST.out.airr))
    REPORT(CLONALITY.out.metrics.join(CLONALITY.out.clonotypes))

    // ---------- Negative-control assertion ------------------------------
    CLONALITY.out.metrics
        .join(ch_expected)
        .map { sid, metrics_json, expected ->
            def ok = true
            def msg = ''
            if (expected == 'negative') {
                def m = new groovy.json.JsonSlurper().parseText(metrics_json.text)
                def n = (m.aggregate?.n_clonotypes ?: 0) as int
                if (n > (params.negative_max_clonotypes ?: 0)) {
                    ok = false
                    msg = "expected negative, got ${n} clonotypes"
                }
            }
            log.info "[QC] ${sid}: expected=${expected ?: 'n/a'} ${ok ? 'PASS' : 'FAIL — ' + msg}"
            return tuple(sid, expected, ok, msg)
        }
        .collectFile(name: 'qc_assertions.tsv',
                     storeDir: "${params.outdir}/pipeline_info",
                     newLine: true)
            { it -> "${it[0]}\t${it[1] ?: ''}\t${it[2] ? 'PASS' : 'FAIL'}\t${it[3]}" }

    // ---------- Versions -------------------------------------------------
    Channel.empty()
        .mix(FASTP.out.versions)
        .mix(TRUST4.out.versions)
        .mix(IGBLAST.out.versions)
        .mix(CLONALITY.out.versions)
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}

workflow.onComplete {
    log.info "Pipeline complete — ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Outputs: ${params.outdir}"
}
