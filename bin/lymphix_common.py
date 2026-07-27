#!/usr/bin/env python3
"""
lymphix_common.py — the single definition of everything the bin/ scripts share.

WHY THIS FILE EXISTS
--------------------
The locus list, the clonality thresholds, the colour palette, the logo loader
and — most importantly — the rule that decides whether a locus is clonal were
each copy-pasted into five or six scripts. The copies drifted. A sample could
be reported as clonal by generate_report.py and indeterminate by
cohort_summary.py at the same time, from the same metrics.json, because one
copy tested the V(D)J yield before the clonality and the other tested it
after. That class of defect is invisible in review: every individual file is
internally consistent and reads correctly.

So the rule lives here once, and every consumer imports it. If the rule is
wrong it is now wrong in one place and fixable in one place.

IMPORT MECHANICS
----------------
Nextflow places <projectDir>/bin on the PATH of every task, and the scripts in
bin/ import each other by bare module name. This file must therefore sit in
bin/ alongside them, and the container images must copy the whole bin/
directory rather than an enumerated list of scripts — otherwise this module is
simply absent inside the image and the failure appears only at run time, in a
task, on a real sample.
"""
from __future__ import annotations

import math
from pathlib import Path


# ---------------------------------------------------------------------------
# Loci and lineage groupings
# ---------------------------------------------------------------------------
# LOCI is the canonical iteration order. Every consumer walks it rather than
# walking whatever keys happen to be present in a metrics.json, so the order of
# the "clonal loci" column is stable across samples and across reports, and a
# locus missing from an older metrics file is treated as absent rather than
# silently skipped.
LOCI     = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
BCR_LOCI = ["IGH", "IGK", "IGL"]
TCR_LOCI = ["TRA", "TRB", "TRG", "TRD"]
LOCUS_ORDER = BCR_LOCI + TCR_LOCI
LINEAGE_OF  = {**{L: "BCR" for L in BCR_LOCI}, **{L: "TCR" for L in TCR_LOCI}}


# ---------------------------------------------------------------------------
# Clonality thresholds
# ---------------------------------------------------------------------------
# These were bare literals scattered through the verdict code (0.30, 0.20, 20,
# 200, 5). Naming them is not cosmetic: the literals had already diverged in
# how they were *applied*, and an unnamed 200 in one file gives a reviewer no
# way to tell it is meant to be the same 200 as in another.
LOCUS_CLONAL_INDEX_THRESHOLD = 0.30   # normalised clonality index of the locus
TOP_CLONE_FRACTION_THRESHOLD = 0.20   # dominant clone's share of locus reads
SINGLE_CLONE_READS_MIN       = 20     # reads needed to trust an n=1 locus

LOW_VDJ_YIELD_FRACTION = 0.005   # V(D)J reads as a fraction of total input
LOW_VDJ_YIELD_ABSOLUTE = 200     # absolute V(D)J read count
NO_VDJ_SIGNAL_ABSOLUTE = 0       # at or below this, there is nothing to assess

# Below this many clonotypes in total, a non-clonal result is uninformative
# rather than genuinely polyclonal, and is reported as indeterminate.
INDETERMINATE_MAX_CLONOTYPES = 5

# cohort_compare.py: the share of locus reads a clone needs before it is worth
# comparing between two samples of the same patient. Deliberately lower than
# TOP_CLONE_FRACTION_THRESHOLD — the question there is "is there a clone worth
# tracking", not "is this sample clonal".
DOMINANT_FRACTION_THRESHOLD = 0.05
EDIT_DISTANCE_RELATED       = 3      # CDR3 nt edits still consistent with SHM


# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
# One palette for the per-sample report, the cohort summary and the cohort
# overview. Three palettes previously existed for the same eight composition
# pools, so the same sample changed colour depending on which report you opened.
BCR_COLORS = {"IGH": "#1B4F72", "IGK": "#2E86C1", "IGL": "#85C1E9"}
TCR_COLORS = {"TRB": "#922B21", "TRA": "#C0392B", "TRG": "#E67E22", "TRD": "#F5B041"}
LOCUS_COLORS = {**BCR_COLORS, **TCR_COLORS}

COMP_LABELS = {
    "clonal_IGH":              "Clonal IGH",
    "clonal_IGK_kappa":        "Clonal IGK (κ)",
    "clonal_IGL_lambda":       "Clonal IGL (λ)",
    "polyclonal_B":            "Polyclonal B",
    "clonal_TRB":              "Clonal TRB",
    "clonal_TRG_gamma_delta":  "Clonal TRG/TRD",
    "polyclonal_T":            "Polyclonal T",
    "background":              "Background",
}
# The per-sample report spells the pools out in full; the cohort views need
# them to fit in a legend. Same keys, same colours, longer wording.
COMP_LABELS_LONG = {
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

VERDICT_COLORS = {
    "clonal":        "#1B4F72",
    "clonal_B":      "#1B4F72",
    "clonal_T":      "#922B21",
    "no_clonal":     "#27ae60",
    "indeterminate": "#e67e22",
    "no_signal":     "#7F8C8D",
}
VERDICT_LABELS = {
    "clonal":        "Clonal",
    "clonal_B":      "Clonal B",
    "clonal_T":      "Clonal T",
    "no_clonal":     "No clonal expansion",
    # Fires on too few clonotypes to judge diversity OR on yield too low to
    # exclude a clone, so the label must not name only one of the two.
    "indeterminate": "Indeterminate",
    "no_signal":     "No V(D)J signal",
}


# ---------------------------------------------------------------------------
# Small defensive helpers
# ---------------------------------------------------------------------------
# pandas hands us NaN where JSON would give null, so every consumer had grown
# its own `isinstance(v, float) and v != v` incantation. One copy each of these
# is enough.
def safe_float(value) -> float | None:
    """Coerce to float, mapping None, NaN and unparseable values to None.

    Returning None rather than 0.0 matters: a missing clonality index is not a
    clonality index of zero, and conflating them is how an unassessable locus
    starts reporting as confidently polyclonal.
    """
    if value is None:
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def safe_str(value) -> str:
    """Coerce to str, mapping None and NaN to '' so they never reach a report
    as the literal text 'nan'."""
    if value is None or (isinstance(value, float) and value != value):
        return ""
    return str(value)


def gene_name(call) -> str:
    """Reduce an AIRR V/D/J call to a bare gene name.

    Calls arrive as comma-separated ambiguity lists with allele suffixes
    ('IGHV3-23*01,IGHV3-23*04'). Comparing those verbatim between two samples
    of the same patient reports clonal relationship as 'different' whenever
    IgBLAST happened to break a tie differently, so the allele and the
    ambiguity list are both dropped.
    """
    text = safe_str(call)
    return text.split("*")[0].split(",")[0] if text else ""


# ---------------------------------------------------------------------------
# The clonality rule — one definition, every consumer
# ---------------------------------------------------------------------------
def is_locus_clonal(clonality_index=None,
                    top_clone_fraction=None,
                    n_clonotypes: int = 0,
                    n_reads: int = 0,
                    *,
                    require_dominance: bool = True,
                    clonality_min: float | None = None) -> bool:
    """Decide whether one locus carries a clonal expansion.

    Two routes to a clonal call:

      (a) the multi-clone route — the locus is skewed (clonality index at or
          above LOCUS_CLONAL_INDEX_THRESHOLD) *and* one clone dominates it
          (at or above TOP_CLONE_FRACTION_THRESHOLD). Both are needed: a
          skewed locus with no single dominant clone is an oligoclonal or
          reactive picture, not a diagnosable clone.

      (b) the single-clone route — exactly one clonotype at the locus with at
          least SINGLE_CLONE_READS_MIN reads behind it. The clonality index is
          undefined for n=1 (there is no distribution to measure), so route
          (a) can never fire on a genuinely monoclonal sample. Without this
          fallback a cell line with one IGH clone and hundreds of reads is
          reported as 'no clonal expansion', which is the exact opposite of
          the truth.

    `require_dominance=False` drops the top-clone test from route (a). It
    exists for one caller: clonality_metrics.compute_composition(), which is
    answering a different question — "should this locus's reads go in a clonal
    pool" — and has already applied its own per-clone dominance gate
    (`clonal_threshold`, 5% by default) to every read it is about to bin.
    Applying the 20% test there as well would double-count dominance and move
    reads out of the clonal pools of samples that the verdict still calls
    clonal. No caller other than compute_composition() should pass it.

    `clonality_min` overrides LOCUS_CLONAL_INDEX_THRESHOLD. It exists only so
    compute_composition() can keep exposing its own tunable; nothing on the
    command line currently sets it.
    """
    ci  = safe_float(clonality_index)
    top = safe_float(top_clone_fraction)
    n     = int(n_clonotypes or 0)
    reads = int(n_reads or 0)
    ci_min = (LOCUS_CLONAL_INDEX_THRESHOLD if clonality_min is None
              else float(clonality_min))

    multi_clone = ci is not None and ci >= ci_min
    if require_dominance:
        multi_clone = multi_clone and top is not None and top >= TOP_CLONE_FRACTION_THRESHOLD

    single_clone = n == 1 and reads >= SINGLE_CLONE_READS_MIN
    return bool(multi_clone or single_clone)


def clonal_loci(per_locus: dict) -> list[str]:
    """Return the clonal loci from a metrics.json `per_locus` block, in
    canonical LOCI order."""
    per_locus = per_locus or {}
    out = []
    for locus in LOCI:
        m = per_locus.get(locus) or {}
        if is_locus_clonal(m.get("clonality_index"),
                           m.get("top_clone_fraction"),
                           m.get("n_clonotypes"),
                           m.get("n_reads")):
            out.append(locus)
    return out


def verdict_category(per_locus: dict,
                     vdj_reads: int,
                     n_clonotypes: int) -> tuple[str, list[str]]:
    """Return (category, clonal_loci) for one sample.

    Categories:
        no_signal      — no V(D)J reads at all, or no clonotypes survived the
                         minimum-count filter. Clonality is not assessable.
        clonal         — at least one locus meets is_locus_clonal().
        indeterminate  — no clonal locus and too few clonotypes to call the
                         repertoire diverse.
        no_clonal      — a diverse repertoire with no dominance.

    ORDER MATTERS, AND THE TREATMENT OF LOW YIELD IS DELIBERATELY ASYMMETRIC.

    A positive finding stands on low input: a sample with 150 V(D)J reads all
    belonging to one clone is 'clonal', not 'indeterminate'. Short-circuiting
    to indeterminate on yield alone would discard a real positive on exactly
    the low-input samples that are hardest to repeat, and bury it behind a
    category that reads as "we found nothing".

    A negative finding does not. Below LOW_VDJ_YIELD_ABSOLUTE the assay has
    not sampled the repertoire deeply enough to exclude a clone, so the
    absence of one is 'indeterminate', never 'no_clonal'. Reporting
    "No clonal expansion" on 150 reads states a confident negative the data
    cannot support, and the cohort table renders it as a green pill with no
    caveat attached — a false-negative risk on precisely the samples that
    most need repeating.

    So 'no_clonal' means "adequately sampled and genuinely diverse", and is
    the only category that carries that guarantee.
    """
    vdj = int(vdj_reads or 0)
    n   = int(n_clonotypes or 0)

    if vdj <= NO_VDJ_SIGNAL_ABSOLUTE or n == 0:
        return "no_signal", []

    loci = clonal_loci(per_locus)
    if loci:
        return "clonal", loci

    if n < INDETERMINATE_MAX_CLONOTYPES or vdj < LOW_VDJ_YIELD_ABSOLUTE:
        return "indeterminate", []

    return "no_clonal", []


def lineage_verdict(category: str, loci: list[str]) -> str:
    """Collapse a verdict into the five labels the cohort views colour by,
    splitting `clonal` into clonal_B / clonal_T.

    A sample clonal at both a BCR and a TCR locus is labelled clonal_B. That
    is a presentation compromise, not a biological claim: the cohort strip
    plot has one colour per sample. The per-sample report says 'bi-clonal' and
    names both loci, and `loci` is carried alongside this label everywhere it
    is used so the detail is never lost.
    """
    if category != "clonal":
        return category
    if any(L in BCR_LOCI for L in loci):
        return "clonal_B"
    return "clonal_T"


def per_locus_from_flat_row(row: dict) -> dict:
    """Rebuild a metrics.json-shaped `per_locus` block from a flattened cohort
    table row (keys like 'IGH_clonality', 'IGH_top_fraction').

    The cohort table is flat because it becomes a DataFrame. Rather than give
    the clonality rule a second entry point that understands flat rows — which
    is how the two rules diverged in the first place — the row is reshaped
    here and fed to the single rule.
    """
    return {
        L: {
            "clonality_index":    row.get(f"{L}_clonality"),
            "top_clone_fraction": row.get(f"{L}_top_fraction"),
            "n_clonotypes":       row.get(f"{L}_n_clones"),
            "n_reads":            row.get(f"{L}_n_reads"),
        }
        for L in LOCI
    }


# ---------------------------------------------------------------------------
# Shared report assets
# ---------------------------------------------------------------------------
def load_logo_svg() -> str:
    """Return inline SVG markup carrying class='logo', ready to embed.

    Inlined rather than linked because the reports are single files that get
    emailed and archived; an <img src> would render as a broken icon the
    moment the file leaves the results directory. Prefers the mark-only logo
    (the report headers already carry the wordmark as text), falls back to the
    older logo.svg, and finally to a text logo so a missing assets directory
    degrades the branding rather than the report.
    """
    asset_dirs = [
        Path(__file__).resolve().parents[1] / "assets",
        Path(__file__).resolve().parent / "assets",
    ]
    for directory in asset_dirs:
        for name in ("lymphix-mark.svg", "logo.svg"):
            path = directory / name
            if path.exists():
                svg = path.read_text(encoding="utf-8")
                return svg.replace("<svg ", '<svg class="logo" ', 1)
    return '<div class="logo" style="font-weight:700; font-size:22px">Lymphix</div>'


def inline_plotly_js() -> str:
    """Return a <script> block containing the whole of plotly.js.

    Same reason as the logo: the reports must work with no network. The import
    path moved between plotly versions, hence the fallback. There is
    deliberately no CDN fallback — a report that silently needs the internet to
    draw its figures is worse than one that fails while it is being written.
    """
    try:
        from plotly.io import get_plotlyjs
    except ImportError:
        from plotly.offline import get_plotlyjs
    return f'<script type="text/javascript">{get_plotlyjs()}</script>'
