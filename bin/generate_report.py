#!/usr/bin/env python3
"""
generate_report.py — main per-sample BCR/TCR clonality report.

Locked-in figure set (curated from the brochure):

  1. Composition thermometer bar          — headline
  2. Composition donut                    — alternate headline view
  4. Composition Sankey                   — Total reads → Lineage → Pool
  5. κ:λ ratio gauge                      — light-chain restriction flag
  6. Top clones — BCR | TCR side-by-side
 10. CDR3 length — faceted histograms (BCR row / TCR row)
 18. Per-locus clonality bar (BCR | TCR separator)
 22. IGHV mutation donut
 23. V-identity histogram

Computes a single-line clinical verdict string from the metrics block.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from jinja2 import Template


# ---------------------------------------------------------------------------
# Lineage definitions and palettes (mirror build_brochure.py)
# ---------------------------------------------------------------------------
LOCI     = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
BCR_LOCI = ["IGH", "IGK", "IGL"]
TCR_LOCI = ["TRA", "TRB", "TRG", "TRD"]

BCR_COLORS = {"IGH": "#1B4F72", "IGK": "#2E86C1", "IGL": "#85C1E9"}
TCR_COLORS = {"TRB": "#922B21", "TRA": "#C0392B", "TRG": "#E67E22", "TRD": "#F5B041"}
LOCUS_COLORS = {**BCR_COLORS, **TCR_COLORS}
LOCUS_ORDER  = BCR_LOCI + TCR_LOCI
LINEAGE_OF   = {**{l: "BCR" for l in BCR_LOCI}, **{l: "TCR" for l in TCR_LOCI}}

COMP_LABELS = {
    "clonal_IGH":              "Clonal B-cell (IGH)",
    "clonal_IGK_kappa":        "Clonal B-cell (κ-restricted, IGK)",
    "clonal_IGL_lambda":       "Clonal B-cell (λ-restricted, IGL)",
    "polyclonal_B":            "Polyclonal B-cell",
    "clonal_TRB":              "Clonal T-cell (αβ, TRB)",
    "clonal_TRG_gamma_delta":  "Clonal T-cell (γδ, TRG/TRD)",
    "polyclonal_T":            "Polyclonal T-cell",
    "background":              "Background / germline",
}
COMP_COLORS = {
    "clonal_IGH":              "#1B4F72",
    "clonal_IGK_kappa":        "#2E86C1",
    "clonal_IGL_lambda":       "#85C1E9",
    "polyclonal_B":            "#D6EAF8",
    "clonal_TRB":              "#922B21",
    "clonal_TRG_gamma_delta":  "#E67E22",
    "polyclonal_T":            "#FAD7A0",
    "background":              "#7F8C8D",
}
COMP_ORDER = list(COMP_LABELS.keys())

# Clonal thresholds used to derive the headline verdict
LOCUS_CLONAL_INDEX_THRESHOLD = 0.30
TOP_CLONE_FRACTION_THRESHOLD = 0.20
LOW_VDJ_YIELD_FRACTION       = 0.005


# ---------------------------------------------------------------------------
# Verdict logic
# ---------------------------------------------------------------------------
def compute_verdict(metrics: dict) -> dict:
    """Return {'headline': str, 'severity': 'positive'|'negative'|'uncertain',
              'details': list[str], 'warnings': list[str]}."""
    per_locus = metrics.get("per_locus") or {}
    comp      = metrics.get("composition") or {}
    ighv      = metrics.get("ighv_status") or {}
    fractions = comp.get("fractions", {})

    clonal_loci = []
    for L in LOCI:
        m = per_locus.get(L) or {}
        if (m.get("clonality_index") or 0) >= LOCUS_CLONAL_INDEX_THRESHOLD and \
           (m.get("top_clone_fraction") or 0) >= TOP_CLONE_FRACTION_THRESHOLD:
            clonal_loci.append(L)
    bcr_clonal = [L for L in clonal_loci if L in BCR_LOCI]
    tcr_clonal = [L for L in clonal_loci if L in TCR_LOCI]

    details, warnings = [], []
    severity = "negative"
    if bcr_clonal and tcr_clonal:
        headline = f"Bi-clonal — B-cell ({'+'.join(bcr_clonal)}) and T-cell ({'+'.join(tcr_clonal)})"
        severity = "positive"
    elif bcr_clonal:
        headline = f"Clonal B-cell — {'+'.join(bcr_clonal)} dominant"
        severity = "positive"
    elif tcr_clonal:
        headline = f"Clonal T-cell — {'+'.join(tcr_clonal)} dominant"
        severity = "positive"
    else:
        headline = "Polyclonal repertoire — no dominant clone detected"

    # Light-chain restriction supporting flag
    kl_call = comp.get("kappa_lambda_call")
    klr     = comp.get("kappa_lambda_ratio")
    if kl_call and kl_call not in ("balanced", "no_lambda_reads"):
        details.append(f"Light-chain restriction: <b>{kl_call.replace('_',' ')}</b> "
                       f"(κ:λ = {klr:.2f}; normal range 0.5–2.5)")
        if severity == "negative":
            severity = "uncertain"

    # IGHV mutation status (if IGH present)
    if ighv and ighv.get("reads_total", 0) > 0:
        details.append(f"IGHV mutation status (CLL prognostic): "
                       f"<b>{ighv['dominant_status']}</b> "
                       f"({100*ighv['fraction_unmutated']:.1f}% unmutated)")

    # Warnings
    if comp:
        vdj_frac = comp.get("vdj_assigned_reads", 0) / max(1, comp.get("total_input_reads", 1))
        if vdj_frac < LOW_VDJ_YIELD_FRACTION:
            warnings.append(f"Low V(D)J yield: only {100*vdj_frac:.2f}% of input reads "
                            "assembled into clonotypes. Interpret cautiously — possible "
                            "capture failure or non-immune sample.")
        if not comp.get("total_input_reads_known"):
            warnings.append("Total input read count was not supplied — background fraction is "
                            "treated as 0. Pass --total-input-reads to make composition accurate.")

    return dict(headline=headline, severity=severity,
                details=details, warnings=warnings,
                clonal_loci=clonal_loci)


# ---------------------------------------------------------------------------
# Figure builders (all return Plotly HTML fragments)
# ---------------------------------------------------------------------------
def _to_html(fig: go.Figure, height: int = 400) -> str:
    fig.update_layout(margin=dict(l=40, r=20, t=46, b=40), height=height)
    return fig.to_html(include_plotlyjs=False, full_html=False,
                        config={"displaylogo": False, "responsive": True})


def fig_composition_bar(comp: dict, sample_id: str) -> str:
    fig = go.Figure()
    for k in COMP_ORDER:
        v = comp["fractions"][k] * 100
        if v <= 0: continue
        fig.add_trace(go.Bar(
            y=[sample_id], x=[v], name=COMP_LABELS[k], orientation="h",
            marker=dict(color=COMP_COLORS[k], line=dict(color="white", width=1)),
            text=f"{v:.1f}%" if v >= 3 else "",
            textposition="inside", insidetextanchor="middle",
            hovertemplate=f"<b>{COMP_LABELS[k]}</b><br>%{{x:.2f}}%<br>%{{customdata:,}} reads<extra></extra>",
            customdata=[comp["reads"][k]],
        ))
    fig.update_layout(barmode="stack", xaxis=dict(range=[0, 100], title="% of total reads"),
                      yaxis=dict(title=""),
                      legend=dict(orientation="h", y=-0.5))
    return _to_html(fig, height=210)


def fig_composition_donut(comp: dict, sample_id: str) -> str:
    labels = [COMP_LABELS[k] for k in COMP_ORDER]
    vals   = [comp["fractions"][k] * 100 for k in COMP_ORDER]
    reads  = [comp["reads"][k] for k in COMP_ORDER]
    colors = [COMP_COLORS[k] for k in COMP_ORDER]
    fig = go.Figure(go.Pie(
        labels=labels, values=vals, marker=dict(colors=colors),
        hole=0.55, sort=False, direction="clockwise",
        textinfo="percent", textposition="inside",
        hovertemplate="<b>%{label}</b><br>%{percent}<br>%{customdata:,} reads<extra></extra>",
        customdata=reads,
    ))
    fig.update_layout(annotations=[dict(text=sample_id, showarrow=False, font=dict(size=13))])
    return _to_html(fig, height=380)


def fig_composition_sankey(comp: dict) -> str:
    reads = comp["reads"]
    lineage_groups = {
        "B-cell":     ["clonal_IGH", "clonal_IGK_kappa", "clonal_IGL_lambda", "polyclonal_B"],
        "T-cell":     ["clonal_TRB", "clonal_TRG_gamma_delta", "polyclonal_T"],
        "Background": ["background"],
    }
    nodes = (["Total input"]
             + list(lineage_groups.keys())
             + [COMP_LABELS[k] for k in COMP_ORDER])
    node_colors = (["#cccccc", "#1B4F72", "#922B21", "#7F8C8D"]
                   + [COMP_COLORS[k] for k in COMP_ORDER])
    sources, targets, values = [], [], []
    for li, (lin, keys) in enumerate(lineage_groups.items()):
        v = sum(reads[k] for k in keys)
        if v > 0:
            sources.append(0); targets.append(1 + li); values.append(v)
    pool_idx = {k: 4 + i for i, k in enumerate(COMP_ORDER)}
    for li, (lin, keys) in enumerate(lineage_groups.items()):
        for k in keys:
            if reads[k] > 0:
                sources.append(1 + li); targets.append(pool_idx[k]); values.append(reads[k])
    fig = go.Figure(go.Sankey(
        node=dict(label=nodes, color=node_colors, pad=14, thickness=14),
        link=dict(source=sources, target=targets, value=values),
    ))
    return _to_html(fig, height=420)


def fig_kappa_lambda_gauge(comp: dict) -> str | None:
    klr  = comp.get("kappa_lambda_ratio")
    call = comp.get("kappa_lambda_call")
    if klr is None:
        return None
    color = "#27ae60" if call == "balanced" else "#c0392b"
    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta", value=klr,
        number=dict(valueformat=".2f"),
        delta=dict(reference=1.5, valueformat=".2f"),
        gauge=dict(
            axis=dict(range=[0, 5], tickvals=[0, 0.5, 1.5, 2.5, 5],
                       ticktext=["0", "0.5", "1.5", "2.5", "≥5"]),
            bar=dict(color=color, thickness=0.3),
            steps=[dict(range=[0,   0.5], color="#fadbd8"),
                   dict(range=[0.5, 2.5], color="#d4efdf"),
                   dict(range=[2.5, 5.0], color="#fadbd8")],
            threshold=dict(line=dict(color="black", width=3), value=klr),
        ),
        title=dict(text=f"κ:λ ratio<br><b>{call.replace('_',' ')}</b>", font=dict(size=14)),
    ))
    return _to_html(fig, height=320)


def fig_top_clones_split(df: pd.DataFrame) -> str:
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("BCR (IGH/IGK/IGL)", "TCR (TRA/TRB/TRG/TRD)"),
                        horizontal_spacing=0.12)
    for col, loci, palette in [(1, BCR_LOCI, BCR_COLORS), (2, TCR_LOCI, TCR_COLORS)]:
        sub = df[df["locus"].isin(loci)].sort_values("read_count", ascending=False).head(15).copy()
        if sub.empty:
            fig.add_annotation(text="No clonotypes", xref=f"x{col}", yref=f"y{col}",
                               x=0.5, y=0.5, showarrow=False)
            continue
        sub["rank"] = range(1, len(sub) + 1)
        for locus in loci:
            ss = sub[sub["locus"] == locus]
            if ss.empty: continue
            fig.add_trace(go.Bar(
                x=ss["rank"], y=ss["read_count"], name=locus,
                marker_color=palette[locus],
                customdata=np.stack([ss["v_gene"], ss["j_gene"], ss["junction_aa"].fillna("")], axis=-1),
                hovertemplate=("Rank %{x}<br>V: %{customdata[0]}<br>"
                                "J: %{customdata[1]}<br>CDR3: %{customdata[2]}<br>"
                                "Reads: %{y:,}<extra>" + locus + "</extra>"),
            ), row=1, col=col)
        fig.update_yaxes(title_text="Reads", row=1, col=col)
        fig.update_xaxes(title_text="Rank", row=1, col=col)
    fig.update_layout(barmode="stack")
    return _to_html(fig, height=400)


def fig_cdr3_faceted(df: pd.DataFrame) -> str:
    fig = make_subplots(rows=2, cols=4,
                        subplot_titles=BCR_LOCI + [""] + TCR_LOCI,
                        vertical_spacing=0.18, horizontal_spacing=0.06)
    for i, locus in enumerate(BCR_LOCI):
        sub = df[(df["locus"] == locus) & (df["cdr3_len"] > 0)]
        if sub.empty: continue
        fig.add_trace(go.Histogram(x=sub["cdr3_len"], marker_color=BCR_COLORS[locus],
                                    nbinsx=25), row=1, col=i + 1)
    for i, locus in enumerate(TCR_LOCI):
        sub = df[(df["locus"] == locus) & (df["cdr3_len"] > 0)]
        if sub.empty: continue
        fig.add_trace(go.Histogram(x=sub["cdr3_len"], marker_color=TCR_COLORS[locus],
                                    nbinsx=25), row=2, col=i + 1)
    fig.update_layout(showlegend=False)
    fig.update_xaxes(title="CDR3 length (aa)", row=2)
    return _to_html(fig, height=460)


def fig_locus_clonality(per_locus: dict) -> str:
    df = pd.DataFrame([dict(locus=L, lineage=LINEAGE_OF[L],
                            **(per_locus.get(L) or {})) for L in LOCUS_ORDER]).fillna(0)
    fig = px.bar(df, x="locus", y="clonality_index", color="locus",
                  color_discrete_map=LOCUS_COLORS,
                  category_orders={"locus": LOCUS_ORDER}, text_auto=".3f")
    fig.add_vline(x=2.5, line_dash="dot", line_color="black",
                   annotation_text="BCR | TCR", annotation_position="top right")
    fig.add_hline(y=LOCUS_CLONAL_INDEX_THRESHOLD, line_dash="dash", line_color="red",
                   annotation_text=f"Clonal threshold ({LOCUS_CLONAL_INDEX_THRESHOLD})",
                   annotation_position="top left")
    fig.update_layout(yaxis=dict(range=[0, 1], title="Clonality index"),
                       xaxis_title="Locus", showlegend=False)
    return _to_html(fig, height=380)


def _unused_fig_ighv_donut(ighv: dict) -> str:
    labels = ["Unmutated (≥98%)", "Mutated (<98%)", "Unknown"]
    values = [ighv["reads_unmutated"], ighv["reads_mutated"], ighv["reads_unknown"]]
    fig = go.Figure(go.Pie(
        labels=labels, values=values, hole=0.55,
        marker=dict(colors=["#c0392b", "#27ae60", "#95a5a6"]),
        textinfo="label+percent",
    ))
    return _to_html(fig, height=340)


def fig_v_identity_hist(df_igh: pd.DataFrame) -> str | None:
    if "igblast_v_identity" not in df_igh.columns:
        return None
    vids = pd.to_numeric(df_igh["igblast_v_identity"], errors="coerce").dropna()
    if vids.empty:
        return None
    if vids.max() <= 1.0:
        vids = vids * 100
    fig = px.histogram(vids, nbins=30, labels=dict(value="V-identity (%)"))
    fig.add_vline(x=98.0, line_dash="dash", line_color="red",
                   annotation_text="98% cut-off")
    fig.update_layout(showlegend=False, xaxis_title="V-identity (%)", yaxis_title="Clonotypes")
    return _to_html(fig, height=340)


# ---------------------------------------------------------------------------
# Top-clones tables (per lineage)
# ---------------------------------------------------------------------------
def top_clones_table(df: pd.DataFrame, lineage: str, n: int = 8) -> list[dict]:
    loci = BCR_LOCI if lineage == "BCR" else TCR_LOCI
    sub = df[df["locus"].isin(loci)].sort_values("read_count", ascending=False).head(n)
    return [dict(
        locus = r["locus"],
        v     = r.get("v_call", ""),
        d     = r.get("d_call", "") or "",
        j     = r.get("j_call", ""),
        cdr3  = r.get("junction_aa", "") or "",
        reads = int(r["read_count"]),
        pct   = f"{100 * r['locus_fraction']:.1f}%",
    ) for _, r in sub.iterrows()]


# ---------------------------------------------------------------------------
# HTML template (Jinja2)
# ---------------------------------------------------------------------------
TEMPLATE = """
<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>BCR/TCR Clonality — {{ sample_id }}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
:root { --bcr:#1B4F72; --tcr:#922B21; --pos:#c0392b; --neg:#27ae60; --warn:#e67e22; --muted:#666; }
body { font-family: -apple-system, Segoe UI, Roboto, sans-serif; max-width: 1200px;
       margin: 0 auto; color: #222; background: #fafafa; }
header { padding: 20px 28px 8px; border-bottom: 1px solid #ddd; }
header .brand { font-size: 11px; letter-spacing: 1.5px; text-transform: uppercase;
                color: var(--muted); }
header h1 { margin: 4px 0 2px; font-size: 24px; }
header .meta { color: var(--muted); font-size: 12px; }
.verdict { margin: 16px 28px; padding: 14px 18px; border-radius: 6px;
           border-left: 5px solid var(--neg); background: #f4faf6; }
.verdict.positive  { border-color: var(--pos); background: #fdf3f1; }
.verdict.uncertain { border-color: var(--warn); background: #fdf6ec; }
.verdict h2 { margin: 0; font-size: 18px; }
.verdict ul { margin: 8px 0 0 18px; padding: 0; font-size: 13px; color: #333; }
.verdict .warn { color: var(--warn); font-size: 12px; margin-top: 6px; }
section { padding: 16px 28px 24px; border-top: 1px solid #ddd; }
section h2 { font-size: 18px; border-left: 4px solid #4a7; padding-left: 10px; margin: 0 0 8px; }
section p.intro { color: var(--muted); font-size: 12px; margin: 0 0 8px; }
.row { display: grid; gap: 14px; }
.row.two   { grid-template-columns: 1fr 1fr; }
.row.equal { grid-template-columns: 1fr 1fr; }
.row.three { grid-template-columns: 1fr 1fr 1fr; }
.card { background: white; border: 1px solid #e0e0e0; border-radius: 6px;
        padding: 8px 14px; }
.card h3 { font-size: 14px; margin: 4px 0 6px; }
.bcr-only { border-top: 3px solid var(--bcr); }
.tcr-only { border-top: 3px solid var(--tcr); }
table.clones { width: 100%; border-collapse: collapse; font-size: 11px; }
table.clones th, table.clones td { padding: 4px 8px; border-bottom: 1px solid #eee;
                                    text-align: left; }
table.clones th { background: #f4f4f4; font-weight: 600; }
table.clones code { font-size: 10px; }
.kpi-grid { display: grid; grid-template-columns: repeat(4, 1fr); gap: 8px; margin: 8px 0; }
.kpi { background: #fff; border: 1px solid #e0e0e0; border-radius: 4px;
       padding: 8px 10px; }
.kpi .l { font-size: 10px; color: var(--muted); text-transform: uppercase;
          letter-spacing: 0.5px; }
.kpi .v { font-size: 18px; font-weight: 600; }
footer { color: var(--muted); font-size: 10px; padding: 16px 28px 30px; }
</style>
</head><body>

<header>
  <div class="brand">Lymphix &nbsp;|&nbsp; BCR / TCR Clonality Report</div>
  <h1>{{ sample_id }}</h1>
  <div class="meta">
    Generated {{ generated_on }} &nbsp;|&nbsp;
    Total input reads: {{ '{:,}'.format(comp.total_input_reads) if comp else 'n/a' }} &nbsp;|&nbsp;
    V(D)J assembled: {{ '{:,}'.format(comp.vdj_assigned_reads) if comp else 'n/a' }} &nbsp;|&nbsp;
    Min clone count: {{ min_clone_count }}
  </div>
</header>

<div class="verdict {{ verdict.severity }}">
  <h2>{{ verdict.headline }}</h2>
  {% if verdict.details %}<ul>{% for d in verdict.details %}<li>{{ d|safe }}</li>{% endfor %}</ul>{% endif %}
  {% for w in verdict.warnings %}<div class="warn">{{ w }}</div>{% endfor %}
</div>

<section id="summary">
  <h2>Aggregate summary</h2>
  <p class="intro">Sample-level repertoire metrics (all seven loci pooled).</p>
  <div class="kpi-grid" style="grid-template-columns: repeat(4, 1fr)">
    {% for label, val in agg_kpis %}
      <div class="kpi"><div class="l">{{ label }}</div><div class="v">{{ val }}</div></div>
    {% endfor %}
  </div>
</section>

<section id="composition">
  <h2>1. Lineage composition</h2>
  <p class="intro">B-cell vs T-cell vs background, clonal vs polyclonal sub-pools.
     Eight pools sum to 100% of total input reads.</p>
  <div class="card">{{ fig_comp_bar|safe }}</div>
  <div class="row equal" style="margin-top:10px">
    <div class="card">{{ fig_comp_donut|safe }}</div>
    <div class="card">{{ fig_comp_sankey|safe }}</div>
  </div>
</section>

<section id="topclones">
  <h2>2. Top clones — BCR | TCR</h2>
  <p class="intro">Top-ranked clonotypes per lineage, with V/D/J/CDR3 and read fraction.</p>
  <div class="card">{{ fig_top_clones|safe }}</div>
  <div class="row" style="grid-template-columns: 1fr 1fr; margin-top:10px">
    <div class="card bcr-only">
      <h3>Top BCR clonotypes</h3>
      <table class="clones">
        <tr><th>Locus</th><th>V</th><th>D</th><th>J</th><th>CDR3 (aa)</th><th>Reads</th><th>% locus</th></tr>
        {% for c in bcr_clones %}<tr>
          <td>{{ c.locus }}</td><td>{{ c.v }}</td><td>{{ c.d }}</td><td>{{ c.j }}</td>
          <td><code>{{ c.cdr3 }}</code></td><td>{{ '{:,}'.format(c.reads) }}</td><td>{{ c.pct }}</td>
        </tr>{% endfor %}
        {% if not bcr_clones %}<tr><td colspan="7" style="color:#999">No BCR clonotypes detected.</td></tr>{% endif %}
      </table>
    </div>
    <div class="card tcr-only">
      <h3>Top TCR clonotypes</h3>
      <table class="clones">
        <tr><th>Locus</th><th>V</th><th>D</th><th>J</th><th>CDR3 (aa)</th><th>Reads</th><th>% locus</th></tr>
        {% for c in tcr_clones %}<tr>
          <td>{{ c.locus }}</td><td>{{ c.v }}</td><td>{{ c.d }}</td><td>{{ c.j }}</td>
          <td><code>{{ c.cdr3 }}</code></td><td>{{ '{:,}'.format(c.reads) }}</td><td>{{ c.pct }}</td>
        </tr>{% endfor %}
        {% if not tcr_clones %}<tr><td colspan="7" style="color:#999">No TCR clonotypes detected.</td></tr>{% endif %}
      </table>
    </div>
  </div>
</section>

<section id="cdr3">
  <h2>3. CDR3 length spectratype</h2>
  <p class="intro">Classical clonality fingerprint. Polyclonal repertoires show
     Gaussian-shaped curves; clonal expansions produce sharp single-length spikes.
     BCR loci in the top row; TCR in the bottom row.</p>
  <div class="card">{{ fig_cdr3|safe }}</div>
</section>

<section id="locus_table">
  <h2>4. Per-locus summary table</h2>
  <p class="intro">Per-locus clonality metrics. BCR loci first, then TCR.</p>
  <div class="card">
    <table class="clones">
      <tr>
        <th>Locus</th><th>Lineage</th><th>Clonotypes</th><th>Reads</th>
        <th>Top clone %</th><th>Clonality</th><th>Shannon H</th>
        <th>Simpson D</th><th>Gini</th><th>D50</th>
      </tr>
      {% for L in loci_order %}
      {% set m = per_locus.get(L) or {} %}
      <tr>
        <td style="border-left:4px solid {{ locus_colors[L] }}; padding-left:8px">
          <b>{{ L }}</b></td>
        <td>{{ lineage_of[L] }}</td>
        <td>{{ '{:,}'.format(m.get('n_clonotypes') or 0) }}</td>
        <td>{{ '{:,}'.format(m.get('n_reads') or 0) }}</td>
        <td>{{ '%.1f%%'|format(100*(m.get('top_clone_fraction') or 0)) if m.get('top_clone_fraction') is not none else '—' }}</td>
        <td>{{ '%.3f'|format(m.get('clonality_index') or 0) if m.get('clonality_index') is not none else '—' }}</td>
        <td>{{ '%.3f'|format(m.get('shannon_H') or 0) if m.get('shannon_H') is not none else '—' }}</td>
        <td>{{ '%.3f'|format(m.get('simpson_D') or 0) if m.get('simpson_D') is not none else '—' }}</td>
        <td>{{ '%.3f'|format(m.get('gini') or 0) if m.get('gini') is not none else '—' }}</td>
        <td>{{ m.get('D50') or '—' }}</td>
      </tr>
      {% endfor %}
    </table>
  </div>
</section>

{% if ighv and ighv.reads_total > 0 %}
<section id="ighv_table">
  <h2>5. IGHV mutation status</h2>
  <p class="intro">CLL prognostic indicator. ≥{{ ighv.cutoff_percent_v_identity }}% V-identity
     = unmutated (poor prognosis); below = mutated (favourable).</p>
  <div class="card">
    <table class="clones">
      <tr><th>Reads total</th><th>Unmutated</th><th>Mutated</th><th>Unknown</th>
          <th>% unmutated</th><th>Dominant call</th></tr>
      <tr>
        <td>{{ '{:,}'.format(ighv.reads_total) }}</td>
        <td>{{ '{:,}'.format(ighv.reads_unmutated) }}</td>
        <td>{{ '{:,}'.format(ighv.reads_mutated) }}</td>
        <td>{{ '{:,}'.format(ighv.reads_unknown) }}</td>
        <td>{{ '%.1f%%'|format(100*(ighv.fraction_unmutated or 0)) }}</td>
        <td style="font-weight:bold; color:{{ '#c0392b' if ighv.dominant_status=='unmutated' else '#27ae60' }}">
          {{ ighv.dominant_status }}</td>
      </tr>
    </table>
  </div>
</section>
{% endif %}

<footer>
  Pipeline: Lymphix v0.1.0 (TRUST4 + IgBLAST) ·
  Clonal threshold for verdict: per-locus clonality index ≥ {{ '%.2f'|format(locus_threshold) }},
  top-clone fraction ≥ {{ '%.0f'|format(top_clone_threshold * 100) }}% ·
  Composition clonal threshold: per-clone ≥ {{ '%.0f'|format(composition_threshold * 100) }}% of locus.
</footer>

</body></html>
"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    import datetime
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id",  required=True)
    ap.add_argument("--metrics",    required=True, type=Path)
    ap.add_argument("--clonotypes", required=True, type=Path)
    ap.add_argument("--out",        required=True, type=Path)
    args = ap.parse_args(argv)

    metrics = json.loads(args.metrics.read_text())
    df = pd.read_csv(args.clonotypes, sep="\t")

    if not df.empty:
        df["read_count"]     = pd.to_numeric(df["read_count"], errors="coerce").fillna(0).astype(int)
        df["cdr3_len"]       = df["junction_aa"].fillna("").str.len()
        df["v_gene"]         = df["v_call"].fillna("").str.split("*").str[0].str.split(",").str[0]
        df["j_gene"]         = df["j_call"].fillna("").str.split("*").str[0].str.split(",").str[0]
        df["locus_fraction"] = df.groupby("locus")["read_count"].transform(
            lambda s: s / s.sum() if s.sum() else 0)

    comp      = metrics.get("composition")
    per_locus = metrics.get("per_locus") or {}
    ighv      = metrics.get("ighv_status")

    verdict = compute_verdict(metrics)

    fig_comp_bar    = fig_composition_bar(comp, args.sample_id) if comp else ""
    fig_comp_donut  = fig_composition_donut(comp, args.sample_id) if comp else ""
    fig_comp_sankey = fig_composition_sankey(comp) if comp else ""
    fig_kl_gauge    = fig_kappa_lambda_gauge(comp) if comp else None
    fig_top_clones  = fig_top_clones_split(df) if not df.empty else ""
    fig_cdr3        = fig_cdr3_faceted(df) if not df.empty else ""
    fig_locus       = fig_locus_clonality(per_locus)

    # Build aggregate KPI cards from metrics.aggregate
    agg = metrics.get("aggregate") or {}
    def _fmt(v, kind="num"):
        if v is None: return "—"
        if kind == "pct":   return f"{v*100:.1f}%"
        if kind == "float": return f"{v:.3f}"
        return f"{int(v):,}"
    agg_kpis = [
        ("Clonotypes",         _fmt(agg.get("n_clonotypes"))),
        ("Reads",              _fmt(agg.get("n_reads"))),
        ("Top clone %",        _fmt(agg.get("top_clone_fraction"), "pct")),
        ("Clonality index",    _fmt(agg.get("clonality_index"), "float")),
        ("Shannon H",          _fmt(agg.get("shannon_H"), "float")),
        ("Simpson D",          _fmt(agg.get("simpson_D"), "float")),
        ("Gini",               _fmt(agg.get("gini"), "float")),
        ("D50",                _fmt(agg.get("D50"))),
    ]

    html = Template(TEMPLATE).render(
        sample_id       = args.sample_id,
        generated_on    = datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        min_clone_count = metrics.get("min_clone_count", 2),
        comp            = comp,
        verdict         = verdict,
        per_locus       = per_locus,
        ighv            = ighv,
        loci_order      = LOCUS_ORDER,
        locus_colors    = LOCUS_COLORS,
        lineage_of      = LINEAGE_OF,
        locus_threshold = LOCUS_CLONAL_INDEX_THRESHOLD,
        top_clone_threshold = TOP_CLONE_FRACTION_THRESHOLD,
        composition_threshold = (comp or {}).get("clonal_dominance_threshold", 0.05),
        bcr_clones      = top_clones_table(df, "BCR") if not df.empty else [],
        tcr_clones      = top_clones_table(df, "TCR") if not df.empty else [],
        fig_comp_bar    = fig_comp_bar,
        fig_comp_donut  = fig_comp_donut,
        fig_comp_sankey = fig_comp_sankey,
        fig_kl_gauge    = fig_kl_gauge,
        fig_top_clones  = fig_top_clones,
        fig_cdr3        = fig_cdr3,
        agg_kpis        = agg_kpis,
    )
    args.out.write_text(html, encoding="utf-8")
    print(f"[report] wrote {args.out}")


if __name__ == "__main__":
    main()
