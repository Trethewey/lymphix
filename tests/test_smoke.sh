#!/usr/bin/env bash
# End-to-end smoke test for Lymphix analysis layer.
# Generates mock AIRR data → clonality_metrics.py → generate_report.py.
# Does NOT require TRUST4, IgBLAST, Docker, or Nextflow.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
OUT="$ROOT/results_smoke_test"

if command -v py >/dev/null 2>&1;       then PY="py -3"
elif command -v python3 >/dev/null 2>&1; then PY="python3"
else PY="python"
fi

cleanup_and_fail() {
    echo "[FAIL] smoke test failed"
    rm -rf "$OUT"
    exit 1
}
trap cleanup_and_fail ERR

echo "[1/4] Generate mock AIRR fixtures"
rm -rf "$OUT"; mkdir -p "$OUT"
$PY "$ROOT/tests/make_mock_airr.py" --outdir "$OUT" --seed 42 >/dev/null

for sample in CLONAL_TEST POLYCLONAL_TEST; do
    sdir="$OUT/$sample"
    echo "[2/4] Compute clonality metrics for $sample"
    $PY "$ROOT/bin/clonality_metrics.py" \
        --sample-id "$sample" \
        --trust4-airr  "$sdir/$sample.trust4.airr.tsv" \
        --igblast-airr "$sdir/$sample.igblast.airr.tsv" \
        --min-clone-count 2 \
        --total-input-reads 25000 \
        --clonal-dominance-threshold 0.05 \
        --out-metrics    "$sdir/$sample.metrics.json" \
        --out-clonotypes "$sdir/$sample.clonotypes.tsv" \
        --out-top        "$sdir/$sample.top_clones.tsv"

    echo "[3/4] Generate HTML report for $sample"
    $PY "$ROOT/bin/generate_report.py" \
        --sample-id "$sample" \
        --metrics    "$sdir/$sample.metrics.json" \
        --clonotypes "$sdir/$sample.clonotypes.tsv" \
        --out        "$sdir/$sample.report.html"
done

echo "[4/4] Verify expected outputs"
errors=0
for sample in CLONAL_TEST POLYCLONAL_TEST; do
    for f in metrics.json clonotypes.tsv top_clones.tsv report.html; do
        path="$OUT/$sample/$sample.$f"
        if [[ ! -s "$path" ]]; then
            echo "  MISSING: $path"
            errors=$((errors + 1))
        fi
    done
done

# Verify the clonal sample is correctly classified.
# cd into OUT first so relative paths work (sidesteps Git Bash → Windows Python
# path translation: /d/... isn't a valid path for native Windows Python).
pushd "$OUT" >/dev/null
clonal_idx=$($PY -c "import json; print(json.load(open('CLONAL_TEST/CLONAL_TEST.metrics.json'))['aggregate']['clonality_index'])")
poly_idx=$($PY -c "import json; print(json.load(open('POLYCLONAL_TEST/POLYCLONAL_TEST.metrics.json'))['aggregate']['clonality_index'])")
popd >/dev/null

awk -v c="$clonal_idx" -v p="$poly_idx" 'BEGIN {
    if (c <= 0.3) { print "  FAIL: CLONAL aggregate clonality (" c ") is not > 0.3"; exit 1 }
    if (p >= 0.1) { print "  FAIL: POLYCLONAL aggregate clonality (" p ") is not < 0.1"; exit 1 }
    print "  OK: clonal=" c "   polyclonal=" p
}' || errors=$((errors + 1))

trap - ERR
if [[ $errors -gt 0 ]]; then
    echo "[FAIL] $errors error(s) detected"
    exit 1
fi

echo
echo "[PASS] smoke test passed"
echo "Outputs in: $OUT"
