process IGBLAST {
    tag        "${sample_id}"
    label      'process_medium'
    container  'ghcr.io/trethewey/lymphix/igblast:1.22.0'
    publishDir "${params.outdir}/${sample_id}/igblast", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(trust4_airr)
    val   species

    output:
    tuple val(sample_id), path("${sample_id}.igblast.airr.tsv"), emit: airr
    path  "versions.yml", emit: versions

    script:
    def organism = species == 'mouse' ? 'mouse' : 'human'
    """
    # 1) Extract full V(D)J nucleotide sequences from TRUST4 AIRR table -> FASTA
    python3 /opt/scripts/airr_to_fasta.py \\
        --airr ${trust4_airr} \\
        --out ${sample_id}.input.fa

    # 2) Run IgBLAST in AIRR mode (handles both BCR and TCR via auto loci)
    export IGDATA=/opt/igblast
    export IGBLASTDB=/opt/igblast/database

    # Auto-route by locus: if any IGH/IGK/IGL present -> ig; if TRA/TRB/TRG/TRD -> tcr
    # Simplest approach: run both and concatenate, then dedup.
    for ig_seqtype in Ig TCR; do
        igblastn \\
            -germline_db_V \${IGBLASTDB}/${organism}_gl_V \\
            -germline_db_D \${IGBLASTDB}/${organism}_gl_D \\
            -germline_db_J \${IGBLASTDB}/${organism}_gl_J \\
            -auxiliary_data /opt/igblast/optional_file/${organism}_gl.aux \\
            -domain_system imgt \\
            -ig_seqtype \${ig_seqtype} \\
            -organism ${organism} \\
            -outfmt 19 \\
            -num_threads ${task.cpus} \\
            -query ${sample_id}.input.fa \\
            -out ${sample_id}.\${ig_seqtype}.airr.tsv || true
    done

    # 3) Merge Ig + TCR AIRR tables (header from first, drop duplicate header rows)
    python3 /opt/scripts/merge_airr.py \\
        --inputs ${sample_id}.Ig.airr.tsv ${sample_id}.TCR.airr.tsv \\
        --out ${sample_id}.igblast.airr.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblast: \$(igblastn -version 2>&1 | head -1 | sed 's/.*: //')
    END_VERSIONS
    """
}
