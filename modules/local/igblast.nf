/*
 * IGBLAST — re-annotate the TRUST4 contigs with NCBI IgBLAST
 *
 * TRUST4 assembles the V(D)J contigs and makes its own gene calls, but its
 * V-gene assignment is coarser than IgBLAST's. The %V-identity that drives
 * the IGHV mutational status call therefore comes from here, not TRUST4.
 *
 * A capture panel returns immunoglobulin and T-cell receptor rearrangements
 * in the same library, so igblastn is run twice over the same germline
 * databases — once with -ig_seqtype Ig and once with TCR — and the two AIRR
 * tables are merged.
 *
 * FAILURE SEMANTICS
 * -----------------
 * There are two very different outcomes that both end in "few or no rows",
 * and this process must never confuse them:
 *
 *   IgBLAST ran and found nothing — the sample has no assembled
 *   rearrangements. Legitimate for a negative control or a poor library.
 *   The process succeeds and emits a header-only AIRR file, which
 *   clonality_metrics.py reads as a genuine zero.
 *
 *   IgBLAST did not run — a missing germline database, a corrupt image, an
 *   out-of-memory kill. This used to be swallowed by `|| true`, so the
 *   process exited 0, the merge produced an empty file and every downstream
 *   IGHV call silently became null. It now fails the process.
 *
 * The discriminator is the number of query sequences handed to igblastn: if
 * there are none there is nothing to run and nothing can be concluded from
 * igblastn's silence, so the empty result is produced directly. If there is
 * at least one query then igblastn must run to completion and write output,
 * and any other outcome aborts the task.
 */
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
    export IGDATA=/opt/igblast
    export IGBLASTDB=/opt/igblast/database

    # The image captures the -outfmt 19 column header at build time, because
    # igblastn cannot be asked for a header when there is nothing to query.
    # Writing an empty file instead is not equivalent: read_airr() rejects a
    # zero-byte AIRR file, and rightly so — it cannot tell a truncated
    # download from an honest absence of rearrangements.
    AIRR_HEADER=/opt/igblast/airr_header.tsv
    if [ ! -s "\${AIRR_HEADER}" ]; then
        echo "ERROR: \${AIRR_HEADER} is missing from the IgBLAST image." >&2
        echo "       Rebuild it with: make build-igblast" >&2
        exit 1
    fi

    # 1) Extract the V(D)J nucleotide sequences from the TRUST4 AIRR table
    python3 /opt/scripts/airr_to_fasta.py \\
        --airr ${trust4_airr} \\
        --out ${sample_id}.input.fa

    # grep -c exits 1 on no match, which would abort the task under `set -e`;
    # the count it prints is still correct, hence the guard.
    n_query=\$(grep -c '^>' ${sample_id}.input.fa || true)
    echo "[igblast] ${sample_id}: \${n_query} contig(s) to annotate"

    if [ "\${n_query}" -eq 0 ]; then
        # IgBLAST has nothing to do. Emit the schema so downstream reads a
        # real, empty repertoire rather than a missing file.
        echo "[igblast] ${sample_id}: no assembled rearrangements — writing a header-only AIRR table"
        cp "\${AIRR_HEADER}" ${sample_id}.igblast.airr.tsv
    else
        # 2) Run IgBLAST in AIRR mode over both receptor classes.
        #    No failure suppression here: a non-zero exit must fail the task.
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
                -out ${sample_id}.\${ig_seqtype}.airr.tsv

            # -outfmt 19 always writes the header, so an empty file means
            # igblastn reported success without producing a table.
            if [ ! -s ${sample_id}.\${ig_seqtype}.airr.tsv ]; then
                echo "ERROR: igblastn exited 0 for \${ig_seqtype} but wrote no output." >&2
                exit 1
            fi
        done

        # 3) Merge the Ig and TCR tables (header from the first, duplicate
        #    sequence_ids resolved in favour of the row that has a V call)
        python3 /opt/scripts/merge_airr.py \\
            --inputs ${sample_id}.Ig.airr.tsv ${sample_id}.TCR.airr.tsv \\
            --out ${sample_id}.igblast.airr.tsv

        # merge_airr.py writes a zero-byte file when it finds nothing to
        # merge. That cannot happen now that both inputs are checked above,
        # but the output contract is "always a readable AIRR table", so it is
        # enforced here rather than assumed.
        if [ ! -s ${sample_id}.igblast.airr.tsv ]; then
            cp "\${AIRR_HEADER}" ${sample_id}.igblast.airr.tsv
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igblast: \$(igblastn -version 2>&1 | head -1 | sed 's/.*: //')
    END_VERSIONS
    """
}
