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
import argparse, json, datetime
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

# Shared definitions live in bin/lymphix_common.py. This file previously kept
# its own composition palette and its own copy of the clonality rule; the
# palette disagreed with the per-sample report's (so a sample changed colour
# between the two documents) and the rule was a third transcription of the
# same five literals.
from lymphix_common import (            # noqa: E402
    COMP_COLORS,
    COMP_LABELS,
    COMP_ORDER,
    LOCI,
    inline_plotly_js,
    load_logo_svg,
    verdict_category,
)

COMPOSITION_POOLS = [(k, COMP_LABELS[k], COMP_COLORS[k]) for k in COMP_ORDER]
LOCI_ORDER = LOCI


def load_sample(root: Path, sample_id: str) -> tuple[dict, pd.DataFrame] | None:
    """Load one sample's metrics + clonotypes from the cohort root.

    Tries the canonical pipeline layout first:
        <root>/<sample_id>/<sample_id>.metrics.json
    Falls back to the legacy `<sample_id>_results/` suffix for older runs.
    """
    candidates = [
        root / sample_id / f"{sample_id}.metrics.json",
        root / f"{sample_id}_results" / f"{sample_id}.metrics.json",  # legacy
    ]
    metrics_path = next((p for p in candidates if p.exists()), None)
    if metrics_path is None:
        return None
    clones_path = metrics_path.parent / f"{sample_id}.clonotypes.tsv"
    metrics = json.loads(metrics_path.read_text())
    clones = pd.read_csv(clones_path, sep="\t") if clones_path.exists() else pd.DataFrame()
    return metrics, clones


def derive_verdict(metrics: dict) -> tuple[str, list[str]]:
    """Verdict + clonal loci for one sample, from the shared rule."""
    comp = metrics.get("composition") or {}
    agg  = metrics.get("aggregate") or {}
    return verdict_category(metrics.get("per_locus") or {},
                            comp.get("vdj_assigned_reads", 0) or 0,
                            agg.get("n_clonotypes", 0) or 0)


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
            "IGHV":        ighv.get("dominant_clone_status") or "—",
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


def render(samples: dict[str, dict], out: Path, plotly_mode: str = "inline") -> None:
    n_pass = sum(1 for s in samples.values() if s.get("pass"))
    n_total = len(samples)
    table_html = build_table(samples).to_html(include_plotlyjs=False, full_html=False, div_id="tbl")
    comp_html  = build_composition(samples).to_html(include_plotlyjs=False, full_html=False, div_id="comp")
    heat_html  = build_clonality_heatmap(samples).to_html(include_plotlyjs=False, full_html=False, div_id="heat")
    if plotly_mode == "cdn":
        plotly_block = '<script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>'
    else:
        plotly_block = inline_plotly_js()
    logo_svg  = load_logo_svg()
    stamp     = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    status_cls = "pass" if n_pass == n_total else "fail"
    badge_label = "PASS" if n_pass == n_total else "FAIL"

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Lymphix - validation cohort overview</title>
<style>
  :root {{
    --ink: #1f2937;
    --ink-soft: #475569;
    --brand: #1f4e79;
    --brand-soft: #2e75b6;
    --pass: #166534;
    --fail: #991b1b;
    --rule: #e2e8f0;
    --card: #ffffff;
    --bg: #f7fafc;
  }}
  * {{ box-sizing: border-box; }}
  body {{ font-family: -apple-system, "Segoe UI", Roboto, sans-serif;
         margin: 0; padding: 32px 40px 48px; color: var(--ink);
         background: var(--bg); }}
  header {{ display: flex; align-items: center; gap: 20px;
            padding-bottom: 20px; border-bottom: 1px solid var(--rule);
            margin-bottom: 24px; }}
  header .logo svg {{ width: 72px; height: auto; display: block; }}
  header h1 {{ margin: 0; font-size: 24px; color: var(--brand); }}
  header .subtitle {{ margin-top: 4px; color: var(--ink-soft); font-size: 14px; }}
  .badge {{ display: inline-block; padding: 4px 10px; border-radius: 4px;
           font-weight: 600; font-size: 12px; letter-spacing: 0.4px; margin-left: 8px; }}
  .badge.pass {{ background: #dcfce7; color: var(--pass); }}
  .badge.fail {{ background: #fee2e2; color: var(--fail); }}
  h2 {{ margin: 0 0 12px 0; font-size: 15px; color: var(--brand);
        text-transform: uppercase; letter-spacing: 0.6px; }}
  .card {{ background: var(--card); border: 1px solid var(--rule);
          border-radius: 8px; padding: 18px 20px;
          margin-bottom: 18px; box-shadow: 0 1px 2px rgba(15,23,42,0.04); }}
  footer {{ margin-top: 24px; color: var(--ink-soft); font-size: 12px;
            text-align: right; }}
</style>
{plotly_block}
</head><body>

<header>
  <div class="logo">{logo_svg}</div>
  <div>
    <h1>Validation cohort overview
      <span class="badge {status_cls}">{badge_label}</span>
    </h1>
    <div class="subtitle">{n_total} samples - {n_pass}/{n_total} pass - generated {stamp}</div>
  </div>
</header>

<div class="card"><h2>Verdict table</h2>{table_html}</div>
<div class="card"><h2>Lineage composition</h2>{comp_html}</div>
<div class="card"><h2>Per-locus clonality index</h2>{heat_html}</div>

<footer>Lymphix - github.com/Trethewey/lymphix</footer>

</body></html>
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
    ap.add_argument("--plotly", choices=["inline", "cdn"], default="inline",
                    help="Embed Plotly.js (~4.8 MB, offline-safe) or load from "
                         "CDN (~50 KB file, needs internet). Default: inline.")
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
    render(samples, args.out, plotly_mode=args.plotly)


if __name__ == "__main__":
    main()
