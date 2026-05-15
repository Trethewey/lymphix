#!/usr/bin/env python3
"""
cohort_report.py

Render a self-contained HTML cohort overview from per-sample
<sid>.metrics.json and <sid>.clonotypes.tsv outputs:

  * Verdict table
  * Lineage composition (stacked bar)
  * Per-locus clonality-index heatmap
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
try:
    from plotly.io import get_plotlyjs
except ImportError:
    from plotly.offline import get_plotlyjs


COMPOSITION_POOLS = [
    ("clonal_IGH",            "Clonal IGH",        "#1f4e79"),
    ("clonal_IGK_kappa",      "Clonal IGK (kappa)", "#2e75b6"),
    ("clonal_IGL_lambda",     "Clonal IGL (lambda)","#5b9bd5"),
    ("polyclonal_B",          "Polyclonal B",       "#bdd7ee"),
    ("clonal_TRB",            "Clonal TRB",         "#c00000"),
    ("clonal_TRG_gamma_delta","Clonal TRG/TRD",     "#e6804f"),
    ("polyclonal_T",          "Polyclonal T",       "#f8cbad"),
    ("background",            "Background",         "#bfbfbf"),
]
LOCI_ORDER = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]


def load_sample(root: Path, sample_id: str) -> tuple[dict, pd.DataFrame] | None:
    results = root / f"{sample_id}_results"
    metrics_path = results / f"{sample_id}.metrics.json"
    clones_path  = results / f"{sample_id}.clonotypes.tsv"
    if not metrics_path.exists():
        return None
    metrics = json.loads(metrics_path.read_text())
    clones = pd.read_csv(clones_path, sep="\t") if clones_path.exists() else pd.DataFrame()
    return metrics, clones


def derive_verdict(metrics: dict) -> tuple[str, list[str]]:
    """Mirror the report/grading logic for verdict + clonal loci."""
    comp = metrics.get("composition") or {}
    agg  = metrics.get("aggregate") or {}
    vdj_reads = (comp or {}).get("vdj_assigned_reads", 0) or 0
    n_clones  = agg.get("n_clonotypes", 0) or 0
    if vdj_reads == 0 or n_clones == 0:
        return "no_signal", []
    clonal_loci = []
    for L, m in (metrics.get("per_locus") or {}).items():
        if not m:
            continue
        ci    = m.get("clonality_index") or 0
        top   = m.get("top_clone_fraction") or 0
        n     = m.get("n_clonotypes") or 0
        reads = m.get("n_reads") or 0
        if (ci >= 0.30 and top >= 0.20) or (n == 1 and reads >= 20):
            clonal_loci.append(L)
    if clonal_loci:
        return "clonal", clonal_loci
    if n_clones < 5:
        return "indeterminate", []
    return "no_clonal", []


def top_clone_row(clones: pd.DataFrame) -> dict:
    if clones.empty or "read_count" not in clones.columns:
        return {}
    c = clones.sort_values("read_count", ascending=False).iloc[0]
    return {
        "locus":   c.get("locus"),
        "v":       c.get("v_call"),
        "j":       c.get("j_call"),
        "cdr3":    c.get("junction_aa"),
        "reads":   int(c.get("read_count") or 0),
    }


def build_table(samples: dict[str, dict]) -> go.Figure:
    rows = []
    for sid, s in samples.items():
        m = s["metrics"]
        v, loci = derive_verdict(m)
        tc = top_clone_row(s["clones"])
        ighv = (m.get("ighv_status") or {})
        rows.append({
            "Sample":      sid,
            "Biology":     s["expected"].get("biology", ""),
            "Verdict":     v,
            "Clonal loci": ", ".join(loci) or "—",
            "Top V":       tc.get("v") or "—",
            "Top J":       tc.get("j") or "—",
            "Top CDR3":    tc.get("cdr3") or "—",
            "Top reads":   tc.get("reads") or 0,
            "Top %":       round(((m.get("aggregate") or {}).get("top_clone_fraction") or 0) * 100, 1),
            "IGHV":        ighv.get("dominant_status") or "—",
            "Pass":        "PASS" if s.get("pass") else "FAIL",
        })
    df = pd.DataFrame(rows)
    colors = [["#d9f2d9" if p == "PASS" else "#fde0e0" for p in df["Pass"]]] * len(df.columns)
    fig = go.Figure(go.Table(
        header=dict(values=list(df.columns),
                    fill_color="#1f4e79", font=dict(color="white", size=12),
                    align="left", height=32),
        cells=dict(values=[df[c] for c in df.columns],
                   fill_color=colors,
                   align="left", height=28, font=dict(size=11))))
    fig.update_layout(margin=dict(l=10, r=10, t=10, b=10), height=420)
    return fig


def build_composition(samples: dict[str, dict]) -> go.Figure:
    sample_ids = list(samples.keys())
    fig = go.Figure()
    for key, label, color in COMPOSITION_POOLS:
        ys = [(samples[sid]["metrics"].get("composition") or {})
              .get("fractions", {}).get(key, 0) * 100 for sid in sample_ids]
        fig.add_bar(name=label, x=sample_ids, y=ys, marker_color=color,
                    hovertemplate=label + ": %{y:.1f}%<extra>%{x}</extra>")
    fig.update_layout(
        barmode="stack",
        title="Lineage composition (% of total input reads)",
        yaxis=dict(title="% of reads", range=[0, 100]),
        xaxis=dict(title=""),
        legend=dict(orientation="h", y=-0.25),
        margin=dict(l=60, r=20, t=60, b=80),
        height=420)
    return fig


def build_clonality_heatmap(samples: dict[str, dict]) -> go.Figure:
    sample_ids = list(samples.keys())
    z, hover = [], []
    for L in LOCI_ORDER:
        row, hrow = [], []
        for sid in sample_ids:
            pl = (samples[sid]["metrics"].get("per_locus") or {}).get(L) or {}
            ci = pl.get("clonality_index")
            n  = pl.get("n_clonotypes") or 0
            reads = pl.get("n_reads") or 0
            row.append(ci if ci is not None else None)
            hrow.append(f"{sid} · {L}<br>CI={ci}<br>n={n}, reads={reads}")
        z.append(row)
        hover.append(hrow)
    fig = go.Figure(go.Heatmap(
        z=z, x=sample_ids, y=LOCI_ORDER,
        colorscale="RdYlBu_r", zmin=0, zmax=1,
        colorbar=dict(title="Clonality<br>index"),
        text=hover, hovertemplate="%{text}<extra></extra>"))
    fig.update_layout(
        title="Clonality index by sample × locus (1 = monoclonal, 0 = polyclonal)",
        margin=dict(l=60, r=20, t=60, b=60),
        height=360)
    return fig


def render(samples: dict[str, dict], out: Path) -> None:
    n_pass = sum(1 for s in samples.values() if s.get("pass"))
    n_total = len(samples)
    table_html = build_table(samples).to_html(include_plotlyjs=False, full_html=False, div_id="tbl")
    comp_html  = build_composition(samples).to_html(include_plotlyjs=False, full_html=False, div_id="comp")
    heat_html  = build_clonality_heatmap(samples).to_html(include_plotlyjs=False, full_html=False, div_id="heat")
    plotly_js = get_plotlyjs()

    html = f"""<!doctype html><meta charset="utf-8">
<title>Lymphix validation cohort overview</title>
<style>
  body {{ font-family: -apple-system, Segoe UI, sans-serif; margin: 24px; color: #222; }}
  h1   {{ margin: 0 0 4px 0; font-size: 22px; }}
  h2   {{ margin: 28px 0 8px 0; font-size: 16px; color: #1f4e79; }}
  .meta {{ color: #666; font-size: 13px; margin-bottom: 18px; }}
  .pass {{ color: #2e7d32; font-weight: 600; }}
  .fail {{ color: #c62828; font-weight: 600; }}
  .card {{ border: 1px solid #e1e4e8; border-radius: 6px; padding: 12px; margin-bottom: 18px; background: #fff; }}
</style>
<script>{plotly_js}</script>

<h1>Lymphix — validation cohort overview</h1>
<div class="meta">
  {n_total} samples ·
  <span class="{'pass' if n_pass == n_total else 'fail'}">{n_pass}/{n_total} pass</span>
</div>

<div class="card"><h2>Verdict table</h2>{table_html}</div>
<div class="card"><h2>Lineage composition</h2>{comp_html}</div>
<div class="card"><h2>Per-locus clonality index</h2>{heat_html}</div>
"""
    out.write_text(html, encoding="utf-8")
    print(f"[cohort_report] wrote {out}  ({n_pass}/{n_total} pass)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", required=True, type=Path,
                    help="Folder containing one <sample>_results subfolder per sample")
    ap.add_argument("--expected", required=True, type=Path,
                    help="validation_expected.json")
    ap.add_argument("--grading", type=Path,
                    help="Optional _validation_grading.json from grade_validation.py.")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    spec = json.loads(args.expected.read_text())
    grading_pass: dict[str, bool] = {}
    if args.grading and args.grading.exists():
        for r in json.loads(args.grading.read_text()).get("per_sample", []):
            grading_pass[r["sample"]] = bool(r.get("pass"))

    samples: dict[str, dict] = {}
    for sid, sspec in spec["samples"].items():
        loaded = load_sample(args.results_root, sid)
        if loaded is None:
            continue
        metrics, clones = loaded
        samples[sid] = {
            "metrics":  metrics,
            "clones":   clones,
            "expected": sspec,
            "pass":     grading_pass.get(sid, True),
        }
    if not samples:
        raise SystemExit("No sample results found.")
    render(samples, args.out)


if __name__ == "__main__":
    main()
