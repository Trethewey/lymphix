#!/bin/sh
# =====================================================================
# build_germline_db.sh — build the IgBLAST germline V/D/J databases
# =====================================================================
#
# WHY this is a script and not a few inline RUN lines:
#
#   The version this replaces used the bash-only lowercase expansion
#   ${seg,,}. Docker runs RUN through /bin/sh, which on Ubuntu is dash,
#   and dash treats that expansion as a parse error — so the layer
#   aborted and the image could not be built at all. Every command in
#   those layers was also wrapped in `|| true`, so had the parse error
#   been fixed the image would simply have shipped with no databases:
#   IgBLAST would have failed on every sample and every downstream IGHV
#   call would have come back null with nothing to explain why.
#
#   Keeping the logic in a real script lets it stay strict POSIX sh, run
#   under `set -eu`, and — most importantly — assert its own output so a
#   broken reference set fails the build instead of the analysis.
#
# WHY IMGT rather than the NCBI germline downloads:
#
#   The old loop fetched
#     .../igblast/release/database/ncbi_<organism>_<segment>_genes.tar
#   Those files do not exist. That directory publishes only
#   mouse_gl_VDJ.tar, rhesus_monkey_VJ.tar and the human *constant*
#   region set ncbi_human_c_genes.tar — there is no human germline
#   V/D/J download from NCBI. The IgBLAST setup documentation directs
#   users to IMGT for human germline sequences, so both species are
#   taken from the one IMGT reference set. Using a single source also
#   keeps human and mouse annotation mutually comparable, which matters
#   when the same %V-identity cutoff is applied to both.
#
#   IMGT/GENE-DB is free for academic use and carries citation
#   requirements; see https://www.imgt.org/ for the current terms.
#
# WHY the databases mix IG and TR loci:
#
#   modules/local/igblast.nf runs igblastn twice over the same three
#   databases, once with -ig_seqtype Ig and once with TCR, because a
#   capture panel returns both. So each <organism>_gl_<segment> database
#   holds the immunoglobulin *and* T-cell receptor genes for that
#   segment; the seqtype flag, not the database, selects the locus.
#
# Environment:
#   IMGT_REFERENCE_URL  (required) IMGT reference FASTA to build from
#   IGBLAST_HOME        (optional) IgBLAST install root, default /opt/igblast
# =====================================================================
set -eu

IGBLAST_HOME="${IGBLAST_HOME:-/opt/igblast}"
DB_DIR="${IGBLAST_HOME}/database"
REF_FASTA="${DB_DIR}/imgt_reference.fasta"

if [ -z "${IMGT_REFERENCE_URL:-}" ]; then
    echo "ERROR: IMGT_REFERENCE_URL is not set." >&2
    exit 1
fi

mkdir -p "${DB_DIR}"

# ---------------------------------------------------------------------
# 1. Fetch the reference set
# ---------------------------------------------------------------------
echo "==> Downloading the IMGT reference set"
wget --tries=3 --timeout=60 --quiet -O "${REF_FASTA}" "${IMGT_REFERENCE_URL}"

if [ ! -s "${REF_FASTA}" ]; then
    echo "ERROR: the download from ${IMGT_REFERENCE_URL} produced no data." >&2
    exit 1
fi

# Record what was actually pulled. The URL is unversioned, so the only
# honest statement about which release an image contains is the checksum
# taken at build time — do not replace this with a hard-coded digest.
{
    echo "url    ${IMGT_REFERENCE_URL}"
    echo "date   $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "sha256 $(sha256sum "${REF_FASTA}" | cut -d ' ' -f 1)"
    echo "bytes  $(wc -c < "${REF_FASTA}")"
} > "${DB_DIR}/imgt_reference.provenance"

cat "${DB_DIR}/imgt_reference.provenance"

# ---------------------------------------------------------------------
# 2. Split into per-organism, per-segment FASTA
# ---------------------------------------------------------------------
# IMGT headers are pipe-delimited:
#   >accession|allele|species|functionality|label|coordinates|...
# so field 2 is the allele, field 3 the species and field 5 the label.
echo "==> Splitting the reference into per-organism, per-segment FASTA"
awk -v db_dir="${DB_DIR}" '
    BEGIN { FS = "|"; keep = 0 }

    /^>/ {
        keep = 0
        allele = $2; species = $3; label = $5

        # IMGT appends the strain or breed to the species field
        # ("Mus musculus_C57BL/6"), so match the binomial prefix.
        if      (species ~ /^Homo sapiens/) { org = "human" }
        else if (species ~ /^Mus musculus/) { org = "mouse" }
        else                                 { next }

        if      (label == "V-REGION") { seg = "V" }
        else if (label == "D-REGION") { seg = "D" }
        else if (label == "J-REGION") { seg = "J" }
        else                            { next }

        # Immunoglobulin and T-cell receptor loci only. The reference set
        # also carries MHC and other IMGT entries, which would pollute the
        # germline search space.
        if (allele !~ /^(IG[HKL]|TR[ABGD])/) { next }

        # makeblastdb -parse_seqids rejects duplicate identifiers, and IMGT
        # lists an allele once per source accession, so keep the first.
        tag = org "/" seg "/" allele
        if (tag in seen) { next }
        seen[tag] = 1

        current = db_dir "/" org "_gl_" seg ".fasta"
        print ">" allele > current
        keep = 1
        next
    }

    keep { print toupper($0) > current }
' "${REF_FASTA}"

# ---------------------------------------------------------------------
# 3. Build the BLAST databases
# ---------------------------------------------------------------------
echo "==> Building the BLAST databases"
for org in human mouse; do
    for seg in V D J; do
        fasta="${DB_DIR}/${org}_gl_${seg}.fasta"
        db="${DB_DIR}/${org}_gl_${seg}"

        if [ ! -s "${fasta}" ]; then
            echo "ERROR: no ${org} ${seg}-REGION sequences were extracted." >&2
            echo "       The IMGT header layout has most likely changed." >&2
            echo "       Do not ship this image — IgBLAST would annotate nothing." >&2
            exit 1
        fi

        makeblastdb -parse_seqids -dbtype nucl \
            -in "${fasta}" \
            -out "${db}" \
            -title "${org}_gl_${seg}" > /dev/null

        # nhr/nin/nsq are the three core BLAST volume files and are written
        # by every BLAST version; the extra -parse_seqids indices vary.
        for ext in nhr nin nsq; do
            if [ ! -s "${db}.${ext}" ]; then
                echo "ERROR: makeblastdb did not write ${db}.${ext}." >&2
                exit 1
            fi
        done

        n_alleles=$(grep -c '^>' "${fasta}" || true)
        printf '    %-14s %6s alleles\n' "${org}_gl_${seg}" "${n_alleles}"
    done
done

# ---------------------------------------------------------------------
# 4. Check the files IgBLAST needs alongside the databases
# ---------------------------------------------------------------------
# The .aux files carry the J-gene coding frame and CDR3 end position, and
# internal_data carries the domain annotation. Both ship in the IgBLAST
# tarball, but a missing one degrades silently into unannotated output
# rather than an error, so assert them here.
echo "==> Checking the IgBLAST auxiliary data"
for org in human mouse; do
    aux="${IGBLAST_HOME}/optional_file/${org}_gl.aux"
    if [ ! -s "${aux}" ]; then
        echo "ERROR: ${aux} is missing from the IgBLAST distribution." >&2
        exit 1
    fi
    if [ ! -d "${IGBLAST_HOME}/internal_data/${org}" ]; then
        echo "ERROR: ${IGBLAST_HOME}/internal_data/${org} is missing." >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------
# 5. Prove the databases are usable, and capture the AIRR header
# ---------------------------------------------------------------------
# Running igblastn once per organism is the only build-time evidence that
# the databases, the auxiliary data and the domain system actually work
# together. The query is the first record of the organism's own V set, so
# no sequence has to be invented or vendored.
#
# The header line of that run is kept at ${IGBLAST_HOME}/airr_header.tsv.
# modules/local/igblast.nf needs it to emit a header-only AIRR file for a
# sample that genuinely has no rearrangements: IgBLAST cannot be asked for
# a header when there is nothing to query, and an empty file is not the
# same thing as an empty result.
echo "==> Test-running igblastn against the new databases"
export IGDATA="${IGBLAST_HOME}"

for org in human mouse; do
    probe_fa="${DB_DIR}/.probe_${org}.fasta"
    probe_tsv="${DB_DIR}/.probe_${org}.tsv"

    head -n 2 "${DB_DIR}/${org}_gl_V.fasta" > "${probe_fa}"

    igblastn \
        -germline_db_V "${DB_DIR}/${org}_gl_V" \
        -germline_db_D "${DB_DIR}/${org}_gl_D" \
        -germline_db_J "${DB_DIR}/${org}_gl_J" \
        -auxiliary_data "${IGBLAST_HOME}/optional_file/${org}_gl.aux" \
        -domain_system imgt \
        -ig_seqtype Ig \
        -organism "${org}" \
        -outfmt 19 \
        -query "${probe_fa}" \
        -out "${probe_tsv}"

    if [ ! -s "${probe_tsv}" ]; then
        echo "ERROR: igblastn exited 0 for ${org} but wrote no AIRR output." >&2
        exit 1
    fi

    if [ "${org}" = "human" ]; then
        head -n 1 "${probe_tsv}" > "${IGBLAST_HOME}/airr_header.tsv"
    fi

    rm -f "${probe_fa}" "${probe_tsv}"
done

if [ ! -s "${IGBLAST_HOME}/airr_header.tsv" ]; then
    echo "ERROR: failed to capture the AIRR column header." >&2
    exit 1
fi

echo "==> Germline databases built and verified in ${DB_DIR}"
