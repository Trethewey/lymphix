process SIMULATE_TEST_DATA {
    label      'process_low'
    container  'ghcr.io/trethewey/lymphix/trust4:1.0.13'
    publishDir "${params.outdir}/test_data", mode: params.publish_dir_mode

    output:
    path "samplesheet_test.csv",     emit: samplesheet
    path "*.fastq.gz"

    script:
    """
    # Generate two synthetic samples using TRUST4's bundled IMGT reference:
    #   - clonal:      one dominant IGH + one dominant TRB clone at >70%
    #   - polyclonal:  thousands of distinct clones, no dominance
    simulate_repertoire.py \\
        --imgt-ref /opt/trust4/human_IMGT+C.fa \\
        --bcrtcr-ref /opt/trust4/hg38_bcrtcr.fa \\
        --mode clonal \\
        --n-reads ${params.test_clonal_reads} \\
        --seed ${params.test_seed} \\
        --out-prefix CLONAL_TEST

    simulate_repertoire.py \\
        --imgt-ref /opt/trust4/human_IMGT+C.fa \\
        --bcrtcr-ref /opt/trust4/hg38_bcrtcr.fa \\
        --mode polyclonal \\
        --n-reads ${params.test_polyclonal_reads} \\
        --seed ${params.test_seed} \\
        --out-prefix POLYCLONAL_TEST

    cat > samplesheet_test.csv <<EOF
sample_id,fastq_1,fastq_2
CLONAL_TEST,\$PWD/CLONAL_TEST_R1.fastq.gz,\$PWD/CLONAL_TEST_R2.fastq.gz
POLYCLONAL_TEST,\$PWD/POLYCLONAL_TEST_R1.fastq.gz,\$PWD/POLYCLONAL_TEST_R2.fastq.gz
EOF
    """
}
