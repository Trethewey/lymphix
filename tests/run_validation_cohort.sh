#!/usr/bin/env bash
#
# Run the published validation cohort end-to-end:
#   1. Download paired FASTQ for each cell line from ENA.
#   2. Run TRUST4 + clonality_metrics + generate_report per sample.
#   3. Grade the cohort against tests/validation_expected.json.
#   4. Build the cohort overview HTML.
#
# Requirements: TRUST4 (run-trust4 in PATH or TRUST4_DIR set), Python 3.10+,
# pandas, numpy, plotly, jinja2, wget. ~25 GB free disk for downloads.
#
# Usage:
#   tests/run_validation_cohort.sh [DATA_DIR]
# DATA_DIR defaults to ./validation_data and is reused on re-runs (idempotent).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
DATA_DIR="${1:-${LYMPHIX_VALIDATION_DATA:-$ROOT/validation_data}}"
TRUST4_DIR="${TRUST4_DIR:-$HOME/repertoire_tools/TRUST4}"
THREADS="${THREADS:-8}"

# ---- Progress helpers ------------------------------------------------------
# Colour only if stdout is a TTY.
if [[ -t 1 ]]; then
    C_RESET=$'\033[0m'; C_DIM=$'\033[2m'
    C_BLUE=$'\033[34m'; C_GREEN=$'\033[32m'; C_YELLOW=$'\033[33m'; C_RED=$'\033[31m'
    C_BOLD=$'\033[1m'
else
    C_RESET=""; C_DIM=""; C_BLUE=""; C_GREEN=""; C_YELLOW=""; C_RED=""; C_BOLD=""
fi
STEP_I=0; STEP_N=0
banner() {
    local sample="$1" idx="$2" total="$3"
    printf "\n${C_BOLD}${C_BLUE}==> [%d/%d] %s${C_RESET}\n" "$idx" "$total" "$sample"
}
step() {
    STEP_I=$((STEP_I + 1))
    printf "  ${C_BLUE}[%d/%d]${C_RESET} %s\n" "$STEP_I" "$STEP_N" "$1"
}
ok()  { printf "       ${C_GREEN}done${C_RESET} %s\n" "$1"; }
warn(){ printf "       ${C_YELLOW}warn${C_RESET} %s\n" "$1"; }
err() { printf "       ${C_RED}fail${C_RESET} %s\n" "$1"; }

# Run a command and show an elapsed-time spinner. stdout/stderr go to logfile.
spin() {
    local label="$1" logfile="$2"; shift 2
    "$@" >"$logfile" 2>&1 &
    local pid=$! spin='|/-\\' i=0 t0=$SECONDS
    if [[ -t 1 ]]; then
        while kill -0 "$pid" 2>/dev/null; do
            printf "\r       ${C_DIM}%s${C_RESET} %s  %ds " \
                "${spin:i++%${#spin}:1}" "$label" "$((SECONDS - t0))"
            sleep 0.2
        done
        printf "\r%80s\r" ""
    else
        wait "$pid" 2>/dev/null || true
    fi
    wait "$pid"
}

require() {
    command -v "$1" >/dev/null || { err "required tool not in PATH: $1"; exit 1; }
}
require run-trust4
require python3
require wget
require zcat
if [[ ! -f "$TRUST4_DIR/hg38_bcrtcr.fa" || ! -f "$TRUST4_DIR/human_IMGT+C.fa" ]]; then
    err "TRUST4 references not found in $TRUST4_DIR (override with TRUST4_DIR=...)"
    exit 1
fi

mkdir -p "$DATA_DIR"
LOG="$DATA_DIR/_cohort.log"
echo "[$(date -u +%FT%TZ)] Cohort validation start" > "$LOG"
printf "${C_BOLD}Lymphix validation cohort${C_RESET}\n"
printf "  data dir : %s\n" "$DATA_DIR"
printf "  TRUST4   : %s\n" "$TRUST4_DIR"
printf "  threads  : %s\n" "$THREADS"

declare -A SRA=(
    ["JURKAT"]="ERR3931301"       ["MOLT-4"]="SRR25601594"
    ["KARPAS-299"]="DRR505206"    ["NAMALWA"]="SRR387396"
    ["DAUDI"]="SRR17188123"       ["RAJI"]="SRR387394"
    ["OCI-LY1"]="SRR17084931"     ["U-266"]="SRR25601616"
    ["MM.1S"]="SRR13272091"       ["PBMC_HEALTHY"]="SRR26965686"
)
SAMPLE_ORDER=(MM.1S JURKAT MOLT-4 U-266 OCI-LY1 RAJI KARPAS-299 DAUDI NAMALWA PBMC_HEALTHY)
TOTAL=$((${#SAMPLE_ORDER[@]} + 1))   # +1 for POLYCLONAL_SIM

# ENA path: vol1/fastq/<prefix6>/<acc> for 9-char accessions, otherwise
# vol1/fastq/<prefix6>/0XX/<acc> where XX = last (length-6) digits.
ena_url() {
    python3 - "$1" "$2" <<'PY'
import sys
acc, mate = sys.argv[1], sys.argv[2]
prefix = acc[:6]; numeric = acc[3:]
extra = max(0, len(numeric) - 6)
if extra == 0:
    path = f"ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{acc}"
else:
    path = f"ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{acc[-extra:].zfill(3)}/{acc}"
print(f"https://{path}/{acc}_{mate}.fastq.gz")
PY
}

# ---- Real cell-line samples ------------------------------------------------
idx=0
for sample in "${SAMPLE_ORDER[@]}"; do
    idx=$((idx + 1))
    banner "$sample" "$idx" "$TOTAL"
    STEP_I=0; STEP_N=4

    acc="${SRA[$sample]}"
    sdata="$DATA_DIR/$sample"
    sout="$DATA_DIR/${sample}_results"
    mkdir -p "$sdata" "$sout"
    cd "$sout"

    if [[ -s "$sample.metrics.json" && -s "$sample.report.html" ]]; then
        ok "cached - skipping"
        continue
    fi

    R1="$sdata/${acc}_1.fastq.gz"; R2="$sdata/${acc}_2.fastq.gz"
    U1=$(ena_url "$acc" 1);        U2=$(ena_url "$acc" 2)

    step "download FASTQ pairs (ENA: $acc)"
    for pair in "$R1|$U1" "$R2|$U2"; do
        f="${pair%%|*}"; u="${pair##*|}"
        if [[ -s "$f" ]]; then ok "have $(basename "$f")"; continue; fi
        wget --show-progress -q -O "$f" "$u"
        ok "$(basename "$f") $(du -h "$f" | cut -f1)"
    done

    step "count input reads"
    n_r1=$(zcat "$R1" | wc -l); n_r1=$((n_r1 / 4)); total=$((n_r1 * 2))
    ok "$n_r1 pairs ($total reads)"

    step "TRUST4 assembly"
    spin "running TRUST4 on $THREADS threads" "trust4.log" \
        run-trust4 -f "$TRUST4_DIR/hg38_bcrtcr.fa" --ref "$TRUST4_DIR/human_IMGT+C.fa" \
        -1 "$R1" -2 "$R2" -o "$sample" -t "$THREADS"
    ok "TRUST4 finished ($(wc -l < ${sample}_airr.tsv) AIRR rows)"

    step "clonality metrics + HTML report"
    python3 "$ROOT/bin/clonality_metrics.py" \
        --sample-id "$sample" \
        --trust4-airr "${sample}_airr.tsv" --igblast-airr "${sample}_airr.tsv" \
        --min-clone-count 2 --total-input-reads "$total" \
        --composition-denominator vdj --read-length 150 \
        --out-metrics    "$sample.metrics.json" \
        --out-clonotypes "$sample.clonotypes.tsv" \
        --out-top        "$sample.top_clones.tsv" >/dev/null
    python3 "$ROOT/bin/generate_report.py" \
        --sample-id "$sample" \
        --metrics    "$sample.metrics.json" \
        --clonotypes "$sample.clonotypes.tsv" \
        --out        "$sample.report.html" >/dev/null
    ok "report: $sample.report.html"
    echo "[$(date -u +%FT%TZ)] $sample done" >> "$LOG"
done

# ---- Synthetic polyclonal positive control ---------------------------------
banner "POLYCLONAL_SIM (synthetic)" "$TOTAL" "$TOTAL"
STEP_I=0; STEP_N=3
sample=POLYCLONAL_SIM
sdata="$DATA_DIR/$sample"; sout="$DATA_DIR/${sample}_results"
mkdir -p "$sdata" "$sout"; cd "$sout"
if [[ -s "$sample.metrics.json" ]]; then
    ok "cached - skipping"
else
    step "simulate polyclonal repertoire (400k reads, seed 7)"
    python3 "$ROOT/bin/simulate_repertoire.py" \
        --imgt-ref "$TRUST4_DIR/human_IMGT+C.fa" \
        --mode polyclonal --n-reads 400000 --seed 7 \
        --out-prefix "$sdata/$sample" >/dev/null
    ok "wrote $sdata/${sample}_R{1,2}.fastq.gz"

    step "TRUST4 assembly"
    spin "running TRUST4 on $THREADS threads" "trust4.log" \
        run-trust4 -f "$TRUST4_DIR/hg38_bcrtcr.fa" --ref "$TRUST4_DIR/human_IMGT+C.fa" \
        -1 "$sdata/${sample}_R1.fastq.gz" -2 "$sdata/${sample}_R2.fastq.gz" \
        -o "$sample" -t "$THREADS"
    ok "TRUST4 finished ($(wc -l < ${sample}_airr.tsv) AIRR rows)"

    step "clonality metrics + HTML report"
    python3 "$ROOT/bin/clonality_metrics.py" \
        --sample-id "$sample" \
        --trust4-airr "${sample}_airr.tsv" --igblast-airr "${sample}_airr.tsv" \
        --min-clone-count 2 --total-input-reads 400000 \
        --composition-denominator vdj --read-length 150 \
        --out-metrics    "$sample.metrics.json" \
        --out-clonotypes "$sample.clonotypes.tsv" \
        --out-top        "$sample.top_clones.tsv" >/dev/null
    python3 "$ROOT/bin/generate_report.py" \
        --sample-id "$sample" \
        --metrics    "$sample.metrics.json" \
        --clonotypes "$sample.clonotypes.tsv" \
        --out        "$sample.report.html" >/dev/null
    ok "report: $sample.report.html"
fi

# ---- Grade + cohort overview ----------------------------------------------
printf "\n${C_BOLD}${C_BLUE}==> Cohort grading + overview${C_RESET}\n"
python3 "$ROOT/bin/grade_validation.py" \
    --results-root "$DATA_DIR" \
    --expected     "$ROOT/tests/validation_expected.json" \
    --out-json     "$DATA_DIR/_validation_grading.json" \
    --out-tsv      "$DATA_DIR/_validation_grading.tsv" | tee -a "$LOG"

python3 "$ROOT/bin/cohort_report.py" \
    --results-root "$DATA_DIR" \
    --expected     "$ROOT/tests/validation_expected.json" \
    --grading      "$DATA_DIR/_validation_grading.json" \
    --out          "$DATA_DIR/cohort_overview.html"

printf "\n${C_BOLD}${C_GREEN}Cohort overview:${C_RESET} %s\n" "$DATA_DIR/cohort_overview.html"
