#!/usr/bin/env python3
"""
cohort_compare.py — cross-sample IGH clonotype comparison for within-patient clonal-relationship's
transformation cohorts.

Groups samples by patient code, lists each sample's dominant IGH clone, and
calls clonal relationship between same-patient samples:

    same          identical CDR3 nt
    related       CDR3 aa identical, CDR3 nt edit distance ≤ 3   (likely SHM descendant)
    different     CDR3 differs by more than 3 nt or different V/J
    no_dominant   sample has no IGH clone above threshold (germline-like)

Output: an HTML cohort summary + a TSV with the same data.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from collections import defaultdict

import pandas as pd

DOMINANT_FRACTION_THRESHOLD = 0.05   # ≥5% of IGH locus reads to count as a "dominant" clone
EDIT_DISTANCE_RELATED       = 3       # ≤3 nt edits = treat as SHM-related


# ---------------------------------------------------------------------------
# Edit distance (Levenshtein, simple DP)
# ---------------------------------------------------------------------------
def edit_distance(a, b) -> int:
    # Defensive: handle None / NaN
    if a is None or (isinstance(a, float) and a != a): a = ""
    if b is None or (isinstance(b, float) and b != b): b = ""
    a, b = str(a), str(b)
    if a == b: return 0
    if not a or not b: return max(len(a), len(b))
    if abs(len(a) - len(b)) > 10: return 999    # short-circuit unrelated lengths
    m, n = len(a), len(b)
    prev = list(range(n + 1))
    for i in range(1, m + 1):
        curr = [i] + [0] * n
        for j in range(1, n + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            curr[j] = min(prev[j] + 1, curr[j - 1] + 1, prev[j - 1] + cost)
        prev = curr
    return prev[n]


# ---------------------------------------------------------------------------
# Load each sample's metrics + dominant IGH clone
# ---------------------------------------------------------------------------
def dominant_igh_clone(clonotypes_path: Path) -> dict | None:
    if not clonotypes_path.exists():
        return None
    df = pd.read_csv(clonotypes_path, sep="\t")
    if df.empty or "locus" not in df.columns:
        return None
    igh = df[df["locus"] == "IGH"].copy()
    if igh.empty:
        return None
    igh["read_count"]     = pd.to_numeric(igh["read_count"], errors="coerce").fillna(0).astype(int)
    igh["locus_fraction"] = igh["read_count"] / igh["read_count"].sum()
    top = igh.sort_values("read_count", ascending=False).iloc[0]
    if top["locus_fraction"] < DOMINANT_FRACTION_THRESHOLD:
        return None
    return dict(
        v_call       = top.get("v_call", ""),
        d_call       = top.get("d_call", "") or "",
        j_call       = top.get("j_call", ""),
        cdr3_aa      = top.get("junction_aa", "") or "",
        cdr3_nt      = top.get("junction", "") or "",
        reads        = int(top["read_count"]),
        fraction     = float(top["locus_fraction"]),
        v_identity   = top.get("igblast_v_identity", None),
    )


def _gene(call) -> str:
    """Defensive gene-name extraction; tolerates NaN / float / None."""
    if call is None or (isinstance(call, float) and call != call):
        return ""
    return str(call).split("*")[0].split(",")[0]


def classify_pair(c1: dict | None, c2: dict | None) -> tuple[str, str]:
    """Return (call, note) — see docstring."""
    if c1 is None and c2 is None:
        return "no_dominant_both", ""
    if c1 is None or c2 is None:
        return "no_dominant", "one sample has no dominant IGH clone"
    v1, v2 = _gene(c1.get("v_call")), _gene(c2.get("v_call"))
    j1, j2 = _gene(c1.get("j_call")), _gene(c2.get("j_call"))
    if not v1 or not v2:
        return "no_dominant", "missing V call"

    if c1["cdr3_nt"] == c2["cdr3_nt"] and v1 == v2 and j1 == j2:
        return "same", f"identical CDR3 nt, same V/J"
    d_nt = edit_distance(c1["cdr3_nt"], c2["cdr3_nt"])
    d_aa = edit_distance(c1["cdr3_aa"], c2["cdr3_aa"])
    if v1 == v2 and j1 == j2 and d_aa == 0 and d_nt <= EDIT_DISTANCE_RELATED:
        return "related", f"same V/J, identical CDR3 aa, {d_nt} nt diffs (SHM)"
    if v1 == v2 and j1 == j2 and d_nt <= EDIT_DISTANCE_RELATED:
        return "related", f"same V/J, {d_nt} nt / {d_aa} aa diffs"
    return "different", f"V {v1} vs {v2}, J {j1} vs {j2}, CDR3 nt edit={d_nt}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True, type=Path,
                    help="Directory containing one subfolder per sample with {sample}.clonotypes.tsv and {sample}.metrics.json")
    ap.add_argument("--samplesheet", required=True, type=Path,
                    help="CSV with sample_id,patient columns")
    ap.add_argument("--out-html",    required=True, type=Path)
    ap.add_argument("--out-tsv",     required=True, type=Path)
    args = ap.parse_args()

    # Load sample → patient mapping
    sheet = pd.read_csv(args.samplesheet)
    if "patient" not in sheet.columns:
        # Infer: first alpha-only token of sample_id
        sheet["patient"] = sheet["sample_id"].apply(
            lambda s: (s.split("_")[0] if s.split("_")[0].isalpha()
                       else s.split("_")[1] if len(s.split("_")) > 1 else s))

    # Load every sample's dominant IGH clone + aggregate clonality
    samples = []
    for _, r in sheet.iterrows():
        sid = r["sample_id"]; pat = r["patient"]
        sdir = args.results_dir / sid
        clone = dominant_igh_clone(sdir / f"{sid}.clonotypes.tsv")
        metrics_path = sdir / f"{sid}.metrics.json"
        clonality = None; verdict = None
        if metrics_path.exists():
            m = json.loads(metrics_path.read_text())
            agg = m.get("aggregate") or {}
            clonality = agg.get("clonality_index")
            igh_m = (m.get("per_locus") or {}).get("IGH") or {}
            igh_clonality = igh_m.get("clonality_index")
            ighv_status = (m.get("ighv_status") or {}).get("dominant_status")
            samples.append(dict(
                sample_id=sid, patient=pat,
                aggregate_clonality=clonality,
                igh_clonality=igh_clonality,
                ighv_status=ighv_status,
                dominant=clone,
            ))
        else:
            samples.append(dict(sample_id=sid, patient=pat,
                                aggregate_clonality=None, igh_clonality=None,
                                ighv_status=None, dominant=None))

    # Group by patient and pairwise-classify
    by_patient = defaultdict(list)
    for s in samples:
        by_patient[s["patient"]].append(s)

    # Build long-form rows for the HTML and TSV
    pair_rows = []
    for pat, group in sorted(by_patient.items()):
        if len(group) < 2:
            pair_rows.append(dict(patient=pat,
                                  sample_a=group[0]["sample_id"], sample_b="—",
                                  call="singleton", note="no within-patient comparison"))
            continue
        # Sort by sample_id for stable output
        group = sorted(group, key=lambda x: x["sample_id"])
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                a, b = group[i], group[j]
                call, note = classify_pair(a["dominant"], b["dominant"])
                pair_rows.append(dict(
                    patient=pat, sample_a=a["sample_id"], sample_b=b["sample_id"],
                    call=call, note=note,
                    cdr3_a=(a["dominant"] or {}).get("cdr3_aa", "—"),
                    cdr3_b=(b["dominant"] or {}).get("cdr3_aa", "—"),
                    v_a=(a["dominant"] or {}).get("v_call", "—"),
                    v_b=(b["dominant"] or {}).get("v_call", "—"),
                    j_a=(a["dominant"] or {}).get("j_call", "—"),
                    j_b=(b["dominant"] or {}).get("j_call", "—"),
                ))

    # TSV out
    pd.DataFrame(pair_rows).to_csv(args.out_tsv, sep="\t", index=False)

    # HTML out
    CALL_COLOR = {
        "same":              "#c0392b",
        "related":           "#e67e22",
        "different":         "#27ae60",
        "no_dominant":       "#95a5a6",
        "no_dominant_both":  "#bdc3c7",
        "singleton":         "#7f8c8d",
    }
    html_parts = ["""<!DOCTYPE html><html><head><meta charset='utf-8'>
<title>Lymphix cohort — within-patient clonal relationship</title>
<style>
body { font-family: -apple-system, Segoe UI, Roboto, sans-serif;
       max-width:1300px; margin:0 auto; padding:24px; color:#222; background:#fafafa; }
h1 { font-size:22px; }
h2 { font-size:16px; border-left:4px solid #4a7; padding-left:8px; margin-top:24px; }
table { border-collapse:collapse; width:100%; font-size:12px; }
th,td { padding:6px 10px; border-bottom:1px solid #ddd; text-align:left; }
th { background:#f4f4f4; }
.tag { font-size:10px; padding:2px 8px; border-radius:99px; color:#fff; text-transform:uppercase; }
code { font-size:11px; }
.cohort-note { color:#666; font-size:12px; margin-bottom:18px; }
</style></head><body>"""]
    html_parts.append("<h1>Lymphix cohort — within-patient clonal relationship</h1>")
    html_parts.append(f"<div class='cohort-note'>"
                       f"{len(samples)} samples across {len(by_patient)} patients. "
                       f"A pair of same-patient samples is called <b>same</b> "
                       f"if CDR3 nt is identical, <b>related</b> if same V/J + CDR3 nt edit ≤ {EDIT_DISTANCE_RELATED}, "
                       f"otherwise <b>different</b>.</div>")

    # Per-sample table
    html_parts.append("<h2>Per-sample dominant IGH clone</h2>")
    html_parts.append("<table><tr><th>Sample</th><th>Patient</th><th>Clonality (agg)</th>"
                       "<th>IGH clonality</th><th>IGHV status</th><th>V</th><th>J</th>"
                       "<th>CDR3 aa</th><th>Top clone %</th></tr>")
    for s in sorted(samples, key=lambda x: (x["patient"], x["sample_id"])):
        d = s["dominant"] or {}
        html_parts.append(
            f"<tr><td><b>{s['sample_id']}</b></td>"
            f"<td>{s['patient']}</td>"
            f"<td>{(s['aggregate_clonality'] or 0):.3f}</td>"
            f"<td>{(s['igh_clonality'] or 0):.3f}</td>"
            f"<td>{s['ighv_status'] or '—'}</td>"
            f"<td>{d.get('v_call','—')}</td>"
            f"<td>{d.get('j_call','—')}</td>"
            f"<td><code>{d.get('cdr3_aa','—')}</code></td>"
            f"<td>{(d.get('fraction',0))*100:.1f}%</td></tr>")
    html_parts.append("</table>")

    # Pairwise comparison
    html_parts.append("<h2>Within-patient pairwise clonal relationship</h2>")
    html_parts.append("<table><tr><th>Patient</th><th>Sample A</th><th>Sample B</th>"
                       "<th>Call</th><th>CDR3 A</th><th>CDR3 B</th><th>V A → V B</th><th>Note</th></tr>")
    for r in pair_rows:
        color = CALL_COLOR.get(r["call"], "#95a5a6")
        html_parts.append(
            f"<tr><td><b>{r['patient']}</b></td>"
            f"<td>{r['sample_a']}</td>"
            f"<td>{r['sample_b']}</td>"
            f"<td><span class='tag' style='background:{color}'>{r['call']}</span></td>"
            f"<td><code>{r.get('cdr3_a','')}</code></td>"
            f"<td><code>{r.get('cdr3_b','')}</code></td>"
            f"<td>{r.get('v_a','—')} → {r.get('v_b','—')}</td>"
            f"<td>{r['note']}</td></tr>")
    html_parts.append("</table></body></html>")

    args.out_html.write_text("\n".join(html_parts), encoding="utf-8")
    print(f"[cohort_compare] wrote {args.out_html} and {args.out_tsv}")


if __name__ == "__main__":
    main()
