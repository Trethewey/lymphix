#!/usr/bin/env python3
"""
cohort_summary.py — cohort-level Lymphix overview HTML.

Aggregates every per-sample metrics.json under a results directory and emits a
single self-contained HTML matching the per-sample report styling (gradient
banner, verdict palette, BCR/TCR colour scheme, inline Plotly).

Usage:
    python bin/cohort_summary.py \
        --results-dir results/2026-05-15_cappseq \
        --cohort-id cappseq \
        --out results/2026-05-15_cappseq/_cohort_summary.html

Discovers metrics by glob `<results-dir>/*/<sample>.metrics.json`. If the
results dir contains a `data_runinfo.json`, its fields are surfaced in the
header.

Sections:
  - Cohort KPIs (n samples, total clones, verdict counts, IGHV split)
  - Verdict + locus matrix table (one row per sample)
  - Stacked composition bar chart (one bar per sample, eight pools)
  - Per-locus clonotype-count heatmap (samples x loci)
  - Aggregate clonality index strip plot
  - IGHV mutation status distribution
  - Top recurring CDR3 (cross-sample) — flags possible shared clones / artefacts
  - Run provenance (from data_runinfo.json)
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from jinja2 import Template


# ---------------------------------------------------------------------------
# Constants — shared with generate_report.py via lymphix_common
# ---------------------------------------------------------------------------
# These were previously re-declared here with the comment "mirror
# generate_report.py". They did not mirror it: the labels and the verdict rule
# had both drifted, so the same sample could carry one colour and one call in
# its own report and a different pair in the cohort view.
from lymphix_common import (            # noqa: E402
    BCR_LOCI,
    COMP_COLORS,
    COMP_LABELS,
    COMP_ORDER,
    LOCI,
    VERDICT_COLORS,
    VERDICT_LABELS,
    inline_plotly_js,
    lineage_verdict,
    load_logo_svg,
    per_locus_from_flat_row,
    safe_float as _safe_float,
    verdict_category,
)


# ---------------------------------------------------------------------------
# Metrics ingestion
# ---------------------------------------------------------------------------
def discover_metrics(results_dir: Path) -> list[Path]:
    """Find all <sample>.metrics.json under <results_dir>/<sample>/ (one level)."""
    paths = sorted(results_dir.glob("*/*.metrics.json"))
    return paths


def load_labels(labels_csv: Path | None) -> dict[str, dict]:
    """Map sample_id (or any prefix-matching token) -> label metadata.

    The CSV must have header row 'cmdl_id' (or 'sample_id') plus optional
    columns: sample_label, sample_type (germline|tumour|cell_line|unknown),
    patient_or_line, germline_cmdl_id, visit. Any other columns are ignored
    but kept on the row for templating.
    """
    if labels_csv is None:
        return {}
    import csv
    out = {}
    with open(labels_csv, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            key = r.get("cmdl_id") or r.get("sample_id")
            if key:
                out[key] = r
    return out


def load_vidjil(vidjil_csv: Path | None) -> dict[str, str]:
    """Map cmdl_short -> a one-line Vidjil V(D)J summary (top real-clone row)."""
    if vidjil_csv is None or not vidjil_csv.exists():
        return {}
    import csv
    seen = {}
    with vidjil_csv.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            cmdl = r.get("cmdl_short", "").strip()
            if not cmdl or cmdl in seen:
                continue
            v = r.get("v") or "—"
            d = r.get("d") or "—"
            j = r.get("j") or "—"
            reads = r.get("reads") or "?"
            seen[cmdl] = f"{v} / {d} / {j}  ({reads} reads)"
    return seen


def lookup_label(labels: dict[str, dict], sample_id: str) -> dict | None:
    """Match sample_id against the labels dict, allowing prefix match
    (so 'CMDL20000596_S137_L004' matches 'CMDL20000596'). Returns None
    if no match."""
    if not labels:
        return None
    if sample_id in labels:
        return labels[sample_id]
    # Try splitting on common separators
    for sep in ("_S", "_L", "."):
        if sep in sample_id:
            prefix = sample_id.split(sep, 1)[0]
            if prefix in labels:
                return labels[prefix]
    return None


def load_runinfo(results_dir: Path) -> dict | None:
    p = results_dir / "data_runinfo.json"
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


def load_clonotypes(results_dir: Path, sample_id: str) -> pd.DataFrame | None:
    """Try to load the per-sample clonotype TSV (next to metrics.json)."""
    tsv = results_dir / sample_id / f"{sample_id}.clonotypes.tsv"
    if not tsv.exists():
        return None
    try:
        df = pd.read_csv(tsv, sep="\t", dtype=str)
        return df
    except Exception:
        return None


def summarise_sample(metrics_path: Path) -> dict:
    """Pull cohort-relevant fields out of a single sample metrics.json."""
    m = json.loads(metrics_path.read_text(encoding="utf-8"))
    agg = m.get("aggregate", {}) or {}
    pl = m.get("per_locus", {}) or {}
    comp = m.get("composition", {}) or {}
    ighv = m.get("ighv_status", {}) or {}

    row = {
        "sample_id":            m.get("sample_id", metrics_path.parent.name),
        "wgs_mode":             bool(m.get("wgs_mode")),
        "min_clone_count":      m.get("min_clone_count"),
        # n_clonotypes and every diversity column below mean different things
        # depending on these two, so they travel with the numbers rather than
        # living only in the per-sample metrics.json. Absent from metrics
        # written before collapsing existed, which is the same as "off".
        "collapse_clonotypes":  bool(m.get("collapse_clonotypes")),
        "collapse_key":         m.get("collapse_key"),
        "n_clonotypes":         int(agg.get("n_clonotypes") or 0),
        "vdj_reads":            int(agg.get("n_reads") or 0),
        "top_clone_fraction":   _safe_float(agg.get("top_clone_fraction")),
        "clonality_index":      _safe_float(agg.get("clonality_index")),
        "shannon_H":            _safe_float(agg.get("shannon_H")),
        "D50":                  agg.get("D50"),
        "ighv_reads_total":     int(ighv.get("reads_total") or 0) if ighv else 0,
        "ighv_dominant_status": ighv.get("dominant_clone_status") if ighv else None,
        "ighv_fraction_unmutated": _safe_float(ighv.get("repertoire_unmutated_read_fraction")) if ighv else None,
        "kappa_lambda_ratio":   _safe_float(comp.get("kappa_lambda_ratio")),
        "kappa_lambda_call":    comp.get("kappa_lambda_call"),
        "composition_denominator": comp.get("denominator_mode"),
        "metrics_path":         str(metrics_path),
    }
    for L in LOCI:
        Lm = pl.get(L, {}) or {}
        row[f"{L}_n_clones"] = int(Lm.get("n_clonotypes") or 0)
        row[f"{L}_n_reads"]   = int(Lm.get("n_reads") or 0)
        row[f"{L}_clonality"] = _safe_float(Lm.get("clonality_index"))
        row[f"{L}_top_fraction"] = _safe_float(Lm.get("top_clone_fraction"))

    fractions = comp.get("fractions") or {}
    for k in COMP_ORDER:
        row[f"comp_{k}"] = _safe_float(fractions.get(k))

    # Derived verdict-like label (mirrors generate_report.py logic at a high level)
    row["verdict"] = derive_verdict(row)
    return row


def derive_verdict(r: dict) -> str:
    """Cohort-table verdict for one flattened row.

    Delegates to the shared rule so the cohort table and the per-sample report
    cannot disagree about the same metrics.json — which is precisely what the
    hand-written copy that used to live here did. It tested the V(D)J yield
    *before* the clonality and returned 'indeterminate' for any sample under
    LOW_VDJ_YIELD_ABSOLUTE reads, so a low-input sample with an unambiguous
    dominant clone appeared as clonal in its own report and as indeterminate
    in the cohort. The shared rule assesses clonality first, so a positive
    survives low input — but a *negative* on low input stays indeterminate,
    because this view attaches no low-yield caveat of its own and
    "No clonal expansion" would otherwise read as a confident exclusion the
    read depth cannot support.

    The clonal_B / clonal_T split is kept because the cohort figures colour by
    it; it comes from lineage_verdict() applied to the shared clonal loci.
    """
    per_locus = per_locus_from_flat_row(r)
    category, loci = verdict_category(per_locus,
                                      r.get("vdj_reads") or 0,
                                      r.get("n_clonotypes") or 0)
    return lineage_verdict(category, loci)


# ---------------------------------------------------------------------------
# Recurring CDR3 detection — cross-sample clonotype overlap
# ---------------------------------------------------------------------------
def find_recurring_clones(results_dir: Path, rows: list[dict],
                          locus_filter: tuple = ("IGH",),
                          min_samples: int = 2,
                          top_n: int = 20) -> list[dict]:
    """For each (locus, CDR3aa) combination, count how many samples it appears in.
    Useful for spotting cross-sample artefacts (sterile V-J products) or genuine
    shared clones in same-patient cohorts."""
    seen = defaultdict(list)
    for r in rows:
        df = load_clonotypes(results_dir, r["sample_id"])
        if df is None or df.empty:
            continue
        if "locus" not in df.columns or "junction_aa" not in df.columns:
            continue
        sub = df[df["locus"].isin(locus_filter)].copy()
        if sub.empty:
            continue
        # de-duplicate within a sample
        for _, row in sub.iterrows():
            key = (row["locus"], row.get("junction_aa") or "", row.get("v_call") or "",
                   row.get("j_call") or "")
            if not key[1]:
                continue
            try:
                rc = int(row.get("read_count") or 0)
            except (TypeError, ValueError):
                rc = 0
            seen[key].append((r["sample_id"], rc))

    recurring = []
    for key, hits in seen.items():
        sample_set = set(s for s, _ in hits)
        if len(sample_set) >= min_samples:
            recurring.append({
                "locus":    key[0],
                "cdr3_aa":  key[1],
                "v":        key[2],
                "j":        key[3],
                "n_samples": len(sample_set),
                "total_reads": sum(rc for _, rc in hits),
                "samples":  ", ".join(sorted(sample_set)),
            })
    recurring.sort(key=lambda x: (-x["n_samples"], -x["total_reads"]))
    return recurring[:top_n]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
PLOTLY_CONFIG = {"displayModeBar": False, "responsive": True}

def _to_html(fig: go.Figure, height: int = 380) -> str:
    fig.update_layout(
        margin=dict(l=40, r=20, t=30, b=40),
        height=height,
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="-apple-system, Segoe UI, Roboto, sans-serif", size=11),
    )
    return fig.to_html(include_plotlyjs=False, full_html=False, config=PLOTLY_CONFIG)


def fig_composition_bar(df: pd.DataFrame) -> str:
    """Stacked bar of composition fractions, one bar per sample."""
    if df.empty:
        return ""
    sids = df["sample_id"].tolist()
    fig = go.Figure()
    for pool in COMP_ORDER:
        col = f"comp_{pool}"
        if col not in df.columns:
            continue
        vals = (df[col].fillna(0) * 100).tolist()
        fig.add_trace(go.Bar(
            name=COMP_LABELS[pool], x=sids, y=vals,
            marker_color=COMP_COLORS[pool],
            hovertemplate="<b>%{x}</b><br>" + COMP_LABELS[pool] + ": %{y:.1f}%<extra></extra>",
        ))
    fig.update_layout(
        barmode="stack",
        xaxis=dict(tickangle=-45, automargin=True),
        yaxis=dict(title="% of denominator", range=[0, 100]),
        legend=dict(orientation="h", y=-0.25, x=0),
        title="Lineage composition per sample",
    )
    return _to_html(fig, height=max(360, 60 + 20 * 8))


def fig_locus_heatmap(df: pd.DataFrame) -> str:
    """Heatmap: samples x loci, value = n_clonotypes per locus."""
    if df.empty:
        return ""
    sids = df["sample_id"].tolist()
    z = np.array([[df.loc[df["sample_id"] == s, f"{L}_n_clones"].iloc[0] for L in LOCI]
                  for s in sids])
    fig = go.Figure(go.Heatmap(
        z=z, x=LOCI, y=sids,
        colorscale="Blues",
        colorbar=dict(title="n_clones"),
        hovertemplate="sample=%{y}<br>locus=%{x}<br>n_clones=%{z}<extra></extra>",
    ))
    fig.update_layout(
        title="Clonotype count per locus, per sample",
        xaxis=dict(side="top"),
        yaxis=dict(autorange="reversed"),
    )
    return _to_html(fig, height=max(300, 40 + 16 * len(sids)))


def fig_clonality_strip(df: pd.DataFrame) -> str:
    """Per-sample aggregate clonality index strip plot, coloured by verdict."""
    if df.empty:
        return ""
    fig = go.Figure()
    for verdict, group in df.groupby("verdict"):
        fig.add_trace(go.Scatter(
            x=group["sample_id"], y=group["clonality_index"],
            mode="markers",
            marker=dict(color=VERDICT_COLORS.get(verdict, "#888"), size=10,
                        line=dict(width=0.5, color="#333")),
            name=VERDICT_LABELS.get(verdict, verdict),
            hovertemplate="<b>%{x}</b><br>clonality=%{y:.3f}<extra></extra>",
        ))
    fig.update_layout(
        title="Aggregate clonality index per sample (coloured by verdict)",
        xaxis=dict(tickangle=-45, automargin=True),
        yaxis=dict(title="Clonality index (1 − H/log N)", range=[0, 1]),
        showlegend=True,
    )
    return _to_html(fig, height=380)


def fig_ighv_distribution(df: pd.DataFrame) -> str:
    """Stacked bar: each sample's IGHV mutated/unmutated/unknown read split."""
    if df.empty:
        return ""
    rows = df[df["ighv_reads_total"] > 0].copy()
    if rows.empty:
        return ""
    sids = rows["sample_id"].tolist()
    # Compute per-sample mutated/unmutated counts from ighv_dominant + fraction
    # (we don't have the raw counts on the cohort row — recompute from metrics)
    mutated, unmutated = [], []
    for sid, p in zip(sids, rows["metrics_path"].tolist()):
        m = json.loads(Path(p).read_text(encoding="utf-8"))
        ighv = m.get("ighv_status") or {}
        mutated.append(int(ighv.get("reads_mutated") or 0))
        unmutated.append(int(ighv.get("reads_unmutated") or 0))
    fig = go.Figure()
    fig.add_trace(go.Bar(name="Mutated (<98% V-id)", x=sids, y=mutated,
                         marker_color="#27ae60"))
    fig.add_trace(go.Bar(name="Unmutated (≥98%)", x=sids, y=unmutated,
                         marker_color="#c0392b"))
    fig.update_layout(
        barmode="stack",
        xaxis=dict(tickangle=-45, automargin=True),
        yaxis=dict(title="IGHV-mapped reads"),
        legend=dict(orientation="h", y=-0.25, x=0),
        title="IGHV mutation status per sample (mutated vs unmutated reads)",
    )
    return _to_html(fig, height=380)


def fig_verdict_donut(df: pd.DataFrame) -> str:
    counts = df["verdict"].value_counts()
    fig = go.Figure(go.Pie(
        labels=[VERDICT_LABELS.get(v, v) for v in counts.index],
        values=counts.values,
        hole=0.55,
        marker=dict(colors=[VERDICT_COLORS.get(v, "#888") for v in counts.index]),
        textinfo="label+value",
    ))
    fig.update_layout(title="Verdict distribution", showlegend=False)
    return _to_html(fig, height=340)


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------
TEMPLATE = """
<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>Lymphix cohort summary — {{ cohort_id }}</title>
{{ plotly_js|safe }}
<style>
:root { --bcr:#1B4F72; --tcr:#922B21; --pos:#c0392b; --neg:#27ae60; --warn:#e67e22; --muted:#666; }
body { font-family: -apple-system, Segoe UI, Roboto, sans-serif; max-width: 1300px;
       margin: 0 auto; color: #222; background: #fafafa; }

header { background: linear-gradient(135deg, #1B4F72 0%, #2C3E50 60%, #5D2018 100%);
         padding: 20px 32px 22px;
         display: grid; grid-template-columns: auto 1fr; column-gap: 28px;
         align-items: center;
         border-bottom: 4px solid #E67E22; color: #ECF0F1; }
header .logo { width: 72px; height: auto; flex-shrink: 0; }
header .titles { min-width: 0; }
header .brand { font-size: 11px; letter-spacing: 2px; text-transform: uppercase;
                color: #BDC3C7; }
header h1 { margin: 4px 0 6px; font-size: 24px; color: #FFFFFF; font-weight: 600; }
header .meta { color: #BDC3C7; font-size: 12px; line-height: 1.6; }
header .logo text[fill="#2C3E50"] { fill: #FFFFFF; }
header .logo text[fill="#7F8C8D"] { fill: #BDC3C7; }
header .logo line[stroke="#2C3E50"] { stroke: #FFFFFF; }

section { padding: 16px 28px 24px; border-top: 1px solid #ddd; }
section h2 { font-size: 18px; border-left: 4px solid #4a7; padding-left: 10px; margin: 0 0 8px; }
section p.intro { color: var(--muted); font-size: 12px; margin: 0 0 8px; }
.card { background: white; border: 1px solid #e0e0e0; border-radius: 6px;
        padding: 8px 14px; margin-bottom: 10px; }
.kpi-grid { display: grid; grid-template-columns: repeat(6, 1fr); gap: 8px; margin: 8px 0; }
.kpi { background: #fff; border: 1px solid #e0e0e0; border-radius: 4px;
       padding: 8px 10px; }
.kpi .l { font-size: 10px; color: var(--muted); text-transform: uppercase;
          letter-spacing: 0.5px; }
.kpi .v { font-size: 18px; font-weight: 600; }

.row { display: grid; gap: 14px; }
.row.two { grid-template-columns: 1fr 1fr; }

table.cohort { width: 100%; border-collapse: collapse; font-size: 11px; background: white; }
table.cohort th, table.cohort td { padding: 5px 8px; border-bottom: 1px solid #eee;
                                    text-align: left; vertical-align: top; }
table.cohort th { background: #f4f4f4; font-weight: 600; position: sticky; top: 0; }
table.cohort td.num { text-align: right; font-variant-numeric: tabular-nums; }
table.cohort td.verdict-pill { font-weight: 600; }
table.cohort tr:hover { background: #fafbfc; }
.pill { display: inline-block; padding: 2px 8px; border-radius: 10px;
        font-size: 10px; color: white; }

.runinfo { background: #f4f6f7; border-left: 4px solid #95a5a6;
           padding: 10px 14px; font-size: 12px; line-height: 1.6;
           border-radius: 0 4px 4px 0; }
.runinfo dt { font-weight: 600; color: var(--muted); }
.runinfo dd { margin: 0 0 6px 12px; font-family: SF Mono, Consolas, monospace; }

footer { background: linear-gradient(135deg, #1B4F72 0%, #2C3E50 60%, #5D2018 100%);
         color: #ECF0F1; padding: 16px 32px 18px; border-top: 4px solid #E67E22;
         margin-top: 24px; font-size: 11px; text-align: center; line-height: 1.6; }
footer .footer-meta { color: #BDC3C7; }
footer .footer-meta .name { color: #ECF0F1; }
</style>
</head><body>

<header>
  {{ logo_svg|safe }}
  <div class="titles">
    <div class="brand">Cohort summary</div>
    <h1>{{ cohort_id }}</h1>
    <div class="meta">
      Generated {{ generated_on }} &nbsp;|&nbsp;
      Samples: {{ n_samples }} &nbsp;|&nbsp;
      Total V(D)J reads: {{ '{:,}'.format(total_vdj_reads) }} &nbsp;|&nbsp;
      Total clonotypes: {{ '{:,}'.format(total_clones) }}
    </div>
  </div>
</header>

<section id="kpis">
  <h2>Cohort at a glance</h2>
  <div class="kpi-grid">
    {% for label, val in kpis %}
      <div class="kpi"><div class="l">{{ label }}</div><div class="v">{{ val }}</div></div>
    {% endfor %}
  </div>
  <div class="row two">
    <div class="card">{{ fig_verdict_donut|safe }}</div>
    <div class="card">{{ fig_clonality_strip|safe }}</div>
  </div>
</section>

<section id="cohort_table">
  <h2>Per-sample table</h2>
  <p class="intro">Sortable in your browser: click a column header. Verdict is derived from the
     same per-locus rules as generate_report.py.</p>
  <div class="card" style="overflow-x: auto; max-height: 600px; overflow-y: auto">
    <table class="cohort">
      <thead><tr>
        <th>Sample</th><th>Label</th><th>DL</th><th>Verdict</th><th class="num">V(D)J reads</th><th class="num">Clones</th>
        <th class="num">Clonality</th><th class="num">Top %</th>
        <th class="num">IGH</th><th class="num">IGK</th><th class="num">IGL</th>
        <th class="num">TRB</th><th>IGHV</th><th>κ:λ</th><th>Vidjil top V(D)J</th>
      </tr></thead>
      <tbody>
      {% for r in rows %}
      <tr>
        <td><b>{{ r.sample_id }}</b></td>
        <td>{% if r.sample_label %}<b>{{ r.sample_label }}</b>{% if r.modality %} <span style="color:#999; font-size:10px">({{ r.modality.replace('_',' ') }})</span>{% elif r.sample_type_lab %} <span style="color:#999; font-size:10px">({{ r.sample_type_lab }})</span>{% endif %}{% endif %}</td>
        <td>{% if r.dl_code %}<span class="pill" style="background:#34495e">{{ r.dl_code }}</span>{% endif %}</td>
        <td class="verdict-pill"><span class="pill" style="background:{{ verdict_colors[r.verdict] }}">{{ verdict_labels[r.verdict] }}</span></td>
        <td class="num">{{ '{:,}'.format(r.vdj_reads) }}</td>
        <td class="num">{{ r.n_clonotypes }}</td>
        <td class="num">{{ '%.3f'|format(r.clonality_index) if r.clonality_index is not none else '—' }}</td>
        <td class="num">{{ '%.1f%%'|format(100 * r.top_clone_fraction) if r.top_clone_fraction is not none else '—' }}</td>
        <td class="num">{{ r.IGH_n_clones }}</td>
        <td class="num">{{ r.IGK_n_clones }}</td>
        <td class="num">{{ r.IGL_n_clones }}</td>
        <td class="num">{{ r.TRB_n_clones }}</td>
        <td>{{ r.ighv_dominant_status or '—' }}</td>
        <td>{{ r.kappa_lambda_call or '—' }}</td>
        <td style="font-family: SF Mono, Consolas, monospace; font-size: 10px">{% if r.vidjil_top_clone %}<b style="color:#27ae60">{{ r.vidjil_top_clone }}</b>{% endif %}</td>
      </tr>
      {% endfor %}
      </tbody>
    </table>
  </div>
</section>

<section id="composition">
  <h2>Lineage composition across samples</h2>
  <p class="intro">Each bar = one sample. Pools sum to 100 % of the configured denominator
     (per sample: total reads or V(D)J reads).</p>
  <div class="card">{{ fig_composition_bar|safe }}</div>
</section>

<section id="loci">
  <h2>Clonotype-count heatmap (samples × loci)</h2>
  <p class="intro">Darker = more clonotypes called at that locus. Empty cells highlight panels
     that miss specific loci.</p>
  <div class="card">{{ fig_locus_heatmap|safe }}</div>
</section>

<section id="ighv">
  <h2>IGHV mutation status</h2>
  <p class="intro">Per-sample IGHV V-identity classification. Mutated = &lt;98 % V-identity
     (favourable CLL prognosis); unmutated = ≥98 % (poor prognosis).</p>
  <div class="card">{{ fig_ighv_distribution|safe }}</div>
</section>

{% if recurring_clones %}
<section id="recurring">
  <h2>Recurring IGH clonotypes across samples</h2>
  <p class="intro">CDR3 amino-acid sequences appearing in ≥2 samples. In same-patient cohorts
     these are shared clones (transformation pairs). Across unrelated samples these are
     usually artefacts — germline-rearrangement / sterile V-J products that escaped the filter,
     or low-yield single-read clonotypes.</p>
  <div class="card" style="overflow-x: auto">
    <table class="cohort">
      <thead><tr>
        <th>Locus</th><th>CDR3 (aa)</th><th>V</th><th>J</th>
        <th class="num">Samples</th><th class="num">Reads</th><th>Sample list</th>
      </tr></thead>
      <tbody>
      {% for c in recurring_clones %}
      <tr>
        <td><b>{{ c.locus }}</b></td>
        <td><code style="font-size:10px">{{ c.cdr3_aa }}</code></td>
        <td>{{ c.v }}</td><td>{{ c.j }}</td>
        <td class="num">{{ c.n_samples }}</td>
        <td class="num">{{ '{:,}'.format(c.total_reads) }}</td>
        <td style="font-size:10px; color:#666">{{ c.samples }}</td>
      </tr>
      {% endfor %}
      </tbody>
    </table>
  </div>
</section>
{% endif %}

{% if runinfo %}
<section id="provenance">
  <h2>Run provenance</h2>
  <p class="intro">From <code>data_runinfo.json</code> in the cohort directory.</p>
  <div class="runinfo">
    <dl>
      <dt>Pipeline</dt><dd>{{ runinfo.get('pipeline', {}).get('name', 'lymphix') }} {{ runinfo.get('pipeline', {}).get('version', '') }}</dd>
      <dt>Run date</dt><dd>{{ runinfo.get('run_date', '—') }}</dd>
      {% if runinfo.get('input') %}
      <dt>Input type</dt><dd>{{ runinfo.input.get('type', '—') }}</dd>
      {% if runinfo.input.get('samplesheet') %}
      <dt>Samplesheet</dt><dd>{{ runinfo.input.samplesheet }}</dd>
      {% endif %}
      {% endif %}
      {% if runinfo.get('filters_applied') %}
      <dt>Germline-rearrangement filter</dt>
      <dd>{{ 'ON' if runinfo.filters_applied.get('germline_rearrangement_filter') else 'OFF' }}</dd>
      {% endif %}
      {% if runinfo.get('caveats') %}
      <dt>Caveats</dt>
      <dd><ul style="margin:0; padding-left:18px">{% for c in runinfo.caveats %}<li>{{ c }}</li>{% endfor %}</ul></dd>
      {% endif %}
    </dl>
  </div>
</section>
{% endif %}

<footer>
  <div class="footer-meta">
    Lymphix cohort summary &middot; Generated {{ generated_on }} &middot;
    <span class="name">Dr C.S. Trethewey</span>
  </div>
  <div class="footer-meta">
    {{ n_samples }} sample(s) &middot; {{ '{:,}'.format(total_vdj_reads) }} V(D)J reads &middot;
    {{ '{:,}'.format(total_clones) }} clonotypes
  </div>
</footer>

</body></html>
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def build_cohort_summary(results_dir: Path, cohort_id: str, out_path: Path,
                          inline_plotly: bool = True,
                          labels_csv: Path | None = None,
                          vidjil_csv: Path | None = None) -> dict:
    metrics_paths = discover_metrics(results_dir)
    if not metrics_paths:
        raise SystemExit(f"No metrics.json files found under {results_dir}/*/")

    rows = [summarise_sample(p) for p in metrics_paths]

    # Attach optional labels (sample label, patient/line, tissue type, visit)
    labels = load_labels(labels_csv)
    vidjil = load_vidjil(vidjil_csv) if vidjil_csv else {}
    for r in rows:
        lab = lookup_label(labels, r["sample_id"])
        r["sample_label"]    = (lab or {}).get("sample_label") or ""
        r["sample_type_lab"] = (lab or {}).get("sample_type") or ""
        r["modality"]        = (lab or {}).get("modality") or ""
        r["patient_or_line"] = (lab or {}).get("patient_or_line") or ""
        r["visit"]           = (lab or {}).get("visit") or ""
        r["dl_code"]         = (lab or {}).get("dl_code") or ""
        # Vidjil lookup uses the CMDL short ID (split on _S to drop the lane suffix)
        cmdl_short = r["sample_id"].split("_S")[0] if "_S" in r["sample_id"] else r["sample_id"].split("_")[0]
        r["vidjil_top_clone"] = vidjil.get(cmdl_short, "")

    df = pd.DataFrame(rows)

    # Sort rows: clonal_B first, then clonal_T, no_clonal, indeterminate, no_signal
    verdict_rank = {"clonal_B": 0, "clonal_T": 1, "no_clonal": 2,
                    "indeterminate": 3, "no_signal": 4}
    df["_rank"] = df["verdict"].map(verdict_rank).fillna(99)
    df = df.sort_values(["_rank", "sample_id"]).drop(columns=["_rank"]).reset_index(drop=True)

    # KPIs
    total_clones = int(df["n_clonotypes"].sum())
    total_vdj    = int(df["vdj_reads"].sum())
    verdict_counts = df["verdict"].value_counts().to_dict()
    n_mutated = int((df["ighv_dominant_status"] == "mutated").sum())
    n_unmutated = int((df["ighv_dominant_status"] == "unmutated").sum())
    n_wgs = int(df["wgs_mode"].sum())

    kpis = [
        ("Samples",          f"{len(df)}"),
        ("Clonal B",         f"{verdict_counts.get('clonal_B', 0)}"),
        ("Clonal T",         f"{verdict_counts.get('clonal_T', 0)}"),
        ("No clonal",        f"{verdict_counts.get('no_clonal', 0)}"),
        ("Indeterminate",    f"{verdict_counts.get('indeterminate', 0)}"),
        ("No V(D)J signal",  f"{verdict_counts.get('no_signal', 0)}"),
        ("IGHV mutated",     f"{n_mutated}"),
        ("IGHV unmutated",   f"{n_unmutated}"),
        ("Total clones",     f"{total_clones:,}"),
        ("Total V(D)J reads", f"{total_vdj:,}"),
        ("Median V(D)J reads", f"{int(df['vdj_reads'].median()):,}"),
        ("--wgs samples",     f"{n_wgs}"),
    ]

    recurring = find_recurring_clones(results_dir, rows, locus_filter=tuple(BCR_LOCI))
    runinfo = load_runinfo(results_dir)

    html = Template(TEMPLATE).render(
        cohort_id     = cohort_id,
        generated_on  = _dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_samples     = len(df),
        total_vdj_reads = total_vdj,
        total_clones  = total_clones,
        kpis          = kpis,
        rows          = df.to_dict(orient="records"),
        verdict_colors = VERDICT_COLORS,
        verdict_labels = VERDICT_LABELS,
        fig_verdict_donut    = fig_verdict_donut(df),
        fig_clonality_strip  = fig_clonality_strip(df),
        fig_composition_bar  = fig_composition_bar(df),
        fig_locus_heatmap    = fig_locus_heatmap(df),
        fig_ighv_distribution = fig_ighv_distribution(df),
        recurring_clones = recurring,
        runinfo          = runinfo,
        plotly_js        = inline_plotly_js() if inline_plotly else
                           '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>',
        logo_svg         = load_logo_svg(),
    )

    out_path.write_text(html, encoding="utf-8")
    return {"n_samples": len(df), "out_path": str(out_path),
            "verdict_counts": verdict_counts,
            "recurring_clones": len(recurring)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True, type=Path,
                    help="Directory containing per-sample subdirs with metrics.json.")
    ap.add_argument("--cohort-id", required=True,
                    help="Cohort label (shown in the report header).")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output HTML path.")
    ap.add_argument("--cdn-plotly", action="store_true",
                    help="Use plotly.js from CDN (small file) instead of inline. "
                         "Default is inline so the report is offline-friendly.")
    ap.add_argument("--labels-csv", type=Path, default=None,
                    help="CSV mapping sample_id (or CMDL ID prefix) -> sample_label / "
                         "patient_or_line / sample_type / visit. When supplied, the "
                         "cohort table shows readable labels (e.g. 'EB5941_V02') instead "
                         "of just the BAM-derived sample_id.")
    ap.add_argument("--vidjil-csv", type=Path, default=None,
                    help="CSV with per-sample top Vidjil V(D)J calls (columns: "
                         "cmdl_short, v, d, j, reads). Adds a 'Vidjil top V(D)J' column "
                         "to the cohort table so independent IGHV calls can be compared "
                         "side by side with lymphix's verdicts.")
    args = ap.parse_args(argv)

    summary = build_cohort_summary(args.results_dir, args.cohort_id, args.out,
                                    inline_plotly=not args.cdn_plotly,
                                    labels_csv=args.labels_csv,
                                    vidjil_csv=args.vidjil_csv)
    print(f"[cohort_summary] {summary['n_samples']} samples -> {summary['out_path']}")
    print(f"  verdicts: {summary['verdict_counts']}")
    print(f"  recurring clones (≥2 samples): {summary['recurring_clones']}")


if __name__ == "__main__":
    main()
