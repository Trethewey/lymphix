#!/usr/bin/env python3
"""
clonality_metrics.py

Merge TRUST4 AIRR + IgBLAST AIRR, compute per-locus and aggregate clonality
metrics, and emit:
    - metrics.json    (per-locus + aggregate)
    - clonotypes.tsv  (merged clonotype table)
    - top_clones.tsv  (top N per locus with V/D/J/CDR3/SHM%)

Clonality conventions follow Adaptive Biotechnologies / immunoSEQ:
    H            = Shannon entropy over clonotype read fractions
    clonality    = 1 - H / log(N_clonotypes)        (0 = polyclonal, 1 = monoclonal)
    D50          = number of clonotypes accounting for 50% of reads
    Simpson D    = sum(p_i^2)
    Gini         = standard inequality coefficient on read counts

IGHV mutation status (CLL convention):
    unmutated if v_identity >= 98.0 %, mutated otherwise.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


import re

LOCI     = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
BCR_LOCI = ["IGH", "IGK", "IGL"]
TCR_LOCI = ["TRA", "TRB", "TRG", "TRD"]

# Clinical κ:λ ratio reference range (light-chain restriction flag)
KAPPA_LAMBDA_NORMAL_LOW  = 0.5
KAPPA_LAMBDA_NORMAL_HIGH = 2.5

# Approximate J-segment contribution to CDR3 (nt). Used by the germline-
# rearrangement filter to estimate how much of the junction is "novel" vs
# explained by germline. Values are conservative — slight overestimates
# would only reduce sensitivity for very short real CDR3s.
J_STEM_NT = {
    "IGHJ1": 30, "IGHJ2": 27, "IGHJ3": 30, "IGHJ4": 21, "IGHJ5": 21, "IGHJ6": 33,
    "IGKJ1": 27, "IGKJ2": 27, "IGKJ3": 27, "IGKJ4": 27, "IGKJ5": 27,
    "IGLJ1": 27, "IGLJ2": 27, "IGLJ3": 27, "IGLJ7": 27,
    "TRAJ":  21, "TRBJ":  30, "TRGJ":  21, "TRDJ":  24,
}


def _cigar_total_match(cigar: str | float | None) -> int:
    """Sum all 'M' (match) operations in a CIGAR string."""
    if not cigar or (isinstance(cigar, float) and math.isnan(cigar)):
        return 0
    return sum(int(m.group(1)) for m in re.finditer(r"(\d+)M", str(cigar)))


def _j_stem_estimate(j_call: str | None) -> int:
    """Return approximate J-segment nt contribution to CDR3."""
    if not isinstance(j_call, str) or not j_call:
        return 21
    name = j_call.split("*")[0].split(",")[0]
    if name in J_STEM_NT:
        return J_STEM_NT[name]
    # Fallback to family prefix
    for prefix, length in J_STEM_NT.items():
        if name.startswith(prefix):
            return length
    return 21


def is_germline_rearrangement(row: dict, *,
                              min_v_match: int   = 100,
                              min_v_identity: float = 85.0,
                              min_junction_diversity: int = 3) -> tuple[bool, str]:
    """
    Detect sterile V-J fusions (germline-rearrangement artefacts).

    A clone is flagged only if all three signatures are present together:
      * V alignment shorter than min_v_match nt
      * V identity at or above min_v_identity %  (near-germline)
      * Junction diversity below min_junction_diversity nt

    The conjunction is required because somatic hypermutation in real B-cell
    lymphomas drops v_identity below 85% on its own.
    """
    junction = row.get("junction") or ""
    j_stem   = _j_stem_estimate(row.get("j_call"))
    diversity = len(junction) - 3 - j_stem
    v_match = _cigar_total_match(row.get("v_cigar"))
    try:
        v_id = float(row.get("v_identity") or 0)
    except (ValueError, TypeError):
        v_id = 0.0

    short_v       = v_match < min_v_match
    sterile_junc  = diversity < min_junction_diversity
    germline_v_id = v_id >= min_v_identity

    if short_v and sterile_junc and germline_v_id:
        return True, (f"V={v_match}nt id={v_id:.1f}% junc_div={diversity}nt "
                      f"(short V + germline-identity + sterile junction)")
    return False, "ok"


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------
def read_airr(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    return df


def infer_locus(v_call: str) -> str | None:
    if not isinstance(v_call, str) or not v_call:
        return None
    for locus in LOCI:
        if v_call.startswith(locus):
            return locus
    return None


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def shannon(counts: np.ndarray) -> float:
    counts = counts[counts > 0]
    if counts.size == 0:
        return 0.0
    p = counts / counts.sum()
    return float(-np.sum(p * np.log(p)))


def simpson_d(counts: np.ndarray) -> float:
    counts = counts[counts > 0]
    if counts.size == 0:
        return 0.0
    p = counts / counts.sum()
    return float(np.sum(p ** 2))


def gini(counts: np.ndarray) -> float:
    counts = np.sort(counts[counts > 0].astype(float))
    n = counts.size
    if n == 0:
        return 0.0
    cum = np.cumsum(counts)
    return float((n + 1 - 2 * np.sum(cum) / cum[-1]) / n)


def d50(counts: np.ndarray) -> int:
    counts = np.sort(counts[counts > 0])[::-1]
    if counts.size == 0:
        return 0
    cumsum = np.cumsum(counts)
    threshold = 0.5 * counts.sum()
    return int(np.searchsorted(cumsum, threshold) + 1)


def clonality_index(counts: np.ndarray) -> float:
    """Adaptive's normalised clonality: 1 - H/log(N). NaN if N<2."""
    counts = counts[counts > 0]
    n = counts.size
    if n < 2:
        return float("nan")
    h = shannon(counts)
    return float(1.0 - h / math.log(n))


def summarise(counts: np.ndarray) -> dict:
    counts = counts[counts > 0]
    if counts.size == 0:
        return {
            "n_clonotypes": 0, "n_reads": 0,
            "top_clone_fraction": None, "shannon_H": None,
            "shannon_H_normalised": None, "clonality_index": None,
            "simpson_D": None, "gini": None, "D50": None,
        }
    total = int(counts.sum())
    h = shannon(counts)
    return {
        "n_clonotypes":          int(counts.size),
        "n_reads":               total,
        "top_clone_fraction":    float(counts.max() / total),
        "shannon_H":             h,
        "shannon_H_normalised":  h / math.log(counts.size) if counts.size > 1 else None,
        "clonality_index":       clonality_index(counts),
        "simpson_D":             simpson_d(counts),
        "gini":                  gini(counts),
        "D50":                   d50(counts),
    }


# ---------------------------------------------------------------------------
# Merge TRUST4 + IgBLAST
# ---------------------------------------------------------------------------
def build_clonotype_table(trust4: pd.DataFrame, igblast: pd.DataFrame) -> pd.DataFrame:
    if trust4.empty:
        return pd.DataFrame()

    # TRUST4 columns we rely on
    keep_cols = [
        "sequence_id", "v_call", "d_call", "j_call", "c_call",
        "junction", "junction_aa", "productive",
        "consensus_count", "duplicate_count",
        "v_cigar", "v_identity",
    ]
    for col in keep_cols:
        if col not in trust4.columns:
            trust4[col] = None
    df = trust4[keep_cols].copy()

    df["consensus_count"] = pd.to_numeric(df["consensus_count"], errors="coerce").fillna(0).astype(int)
    df["duplicate_count"] = pd.to_numeric(df["duplicate_count"], errors="coerce").fillna(0).astype(int)
    df["read_count"]      = df[["consensus_count", "duplicate_count"]].max(axis=1)
    df["locus"]           = df["v_call"].map(infer_locus)

    # Merge IgBLAST v_identity & productive flags by sequence_id
    if not igblast.empty and "sequence_id" in igblast.columns:
        ig_cols = ["sequence_id"]
        for c in ["v_identity", "productive", "complete_vdj", "stop_codon", "vj_in_frame"]:
            if c in igblast.columns:
                ig_cols.append(c)
        ig = igblast[ig_cols].drop_duplicates("sequence_id")
        ig = ig.rename(columns={
            "v_identity":  "igblast_v_identity",
            "productive":  "igblast_productive",
            "complete_vdj": "igblast_complete_vdj",
            "stop_codon":  "igblast_stop_codon",
            "vj_in_frame": "igblast_vj_in_frame",
        })
        df = df.merge(ig, on="sequence_id", how="left")
        if "igblast_v_identity" in df.columns:
            df["igblast_v_identity"] = pd.to_numeric(df["igblast_v_identity"], errors="coerce")

    # Drop rows with no locus or no junction (keep out-of-frame: they are real
    # clones for clonality purposes even when not productive).
    df = df[df["locus"].notna()]
    has_junction = (df["junction"].fillna("").str.len() > 0) | \
                   (df["junction_aa"].fillna("").str.len() > 0)
    df = df[has_junction]

    # Numeric coercion for fields used by the germline-rearrangement filter
    df["v_identity"] = pd.to_numeric(df.get("v_identity"), errors="coerce")
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Read-length awareness
# ---------------------------------------------------------------------------
# A CDR3 is "spanned by a single read" if read_length covers the junction
# plus a V- and J-anchor on either side (15 nt each is a common threshold).
SINGLE_READ_SPAN_ANCHOR_NT = 15


def annotate_single_read_spanning(df: pd.DataFrame, read_length: int) -> pd.DataFrame:
    """Add `cdr3_spanned_by_single_read` and `assembly_inferred` columns.

    `cdr3_spanned_by_single_read = True` when read_length is large enough to
    cover the entire junction plus V- and J-anchors → the CDR3 was directly
    observed on at least one read (high-confidence by construction).
    `assembly_inferred = True` is the inverse — the CDR3 was reconstructed by
    TRUST4's overlap-based assembly across multiple reads and should be flagged
    for further confidence assessment downstream (see bin/infer_long_cdr3.py).
    """
    if df.empty:
        return df
    df = df.copy()
    junc_len = df["junction"].fillna("").str.len()
    span_needed = junc_len + 2 * SINGLE_READ_SPAN_ANCHOR_NT
    df["cdr3_spanned_by_single_read"] = span_needed <= read_length
    df["assembly_inferred"]           = ~df["cdr3_spanned_by_single_read"]
    return df


def adaptive_min_v_match(read_length: int) -> int:
    """V CIGAR match threshold for the germline-rearrangement filter scales
    with read length so the filter doesn't over-reject short-read data.

      75 bp reads  → ~49 nt min V match
     150 bp reads  → ~98 nt    (close to the original 100 default)
     250 bp reads  → ~163 nt
    1000 bp reads  → 200 nt    (capped)
    """
    return int(max(20, min(200, round(read_length * 0.65))))


_CIGAR_M_RE = re.compile(r"(\d+)M")


def _cigar_match_vec(s: pd.Series) -> np.ndarray:
    """Vectorised sum of all 'M' operations in each CIGAR string."""
    return s.fillna("").astype(str).map(
        lambda c: sum(int(n) for n in _CIGAR_M_RE.findall(c)) if c else 0
    ).to_numpy(dtype=np.int32)


def _j_stem_vec(s: pd.Series) -> np.ndarray:
    """Vectorised J-stem length lookup. Strips allele and gene-list suffixes,
    then falls back to a family-prefix scan, defaulting to 21 nt."""
    names = s.fillna("").astype(str).str.split("*", n=1).str[0].str.split(",", n=1).str[0]
    out = names.map(J_STEM_NT)
    if out.isna().any():
        prefixes = sorted(J_STEM_NT.keys(), key=len, reverse=True)
        def fam(name: str) -> int:
            for p in prefixes:
                if name.startswith(p):
                    return J_STEM_NT[p]
            return 21
        out = out.where(out.notna(), names.map(fam))
    return out.to_numpy(dtype=np.int32)


def apply_germline_rearrangement_filter(df: pd.DataFrame, *,
                                         min_v_match: int = 100,
                                         min_v_identity: float = 85.0,
                                         min_junction_diversity: int = 3) -> tuple[pd.DataFrame, dict]:
    """Drop sterile V-J artefacts. Returns (filtered_df, summary_stats)."""
    thresholds = {
        "min_v_match":            min_v_match,
        "min_v_identity":         min_v_identity,
        "min_junction_diversity": min_junction_diversity,
    }
    if df.empty:
        return df, {"n_input": 0, "n_dropped": 0, "n_kept": 0,
                    "by_reason": {}, "thresholds": thresholds}

    v_match  = _cigar_match_vec(df["v_cigar"])
    j_stem   = _j_stem_vec(df["j_call"])
    junc_len = df["junction"].fillna("").astype(str).str.len().to_numpy(dtype=np.int32)
    diversity = junc_len - 3 - j_stem
    v_id = pd.to_numeric(df["v_identity"], errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)

    short_v       = v_match < min_v_match
    sterile_junc  = diversity < min_junction_diversity
    germline_v_id = v_id >= min_v_identity
    is_artefact   = short_v & sterile_junc & germline_v_id

    n_dropped = int(is_artefact.sum())
    by_reason = ({f"short V + germline-identity + sterile junction "
                  f"(<{min_v_match}nt / >={min_v_identity}% / <{min_junction_diversity}nt)": n_dropped}
                 if n_dropped else {})
    kept = df.loc[~is_artefact].reset_index(drop=True)
    return kept, {
        "n_input":    len(df),
        "n_dropped":  n_dropped,
        "n_kept":     int(len(kept)),
        "by_reason":  by_reason,
        "thresholds": thresholds,
    }


def assign_ighv_mutation_status(df: pd.DataFrame, cutoff: float) -> pd.DataFrame:
    """For IGH clonotypes, classify unmutated (>=cutoff) vs mutated."""
    if "igblast_v_identity" not in df.columns:
        df["ighv_status"] = None
        return df
    status = []
    for _, row in df.iterrows():
        if row["locus"] != "IGH":
            status.append(None); continue
        ident = row["igblast_v_identity"]
        if pd.isna(ident):
            status.append("unknown")
        else:
            # IgBLAST v_identity is a fraction (0-100) or 0-1 depending on version; normalise
            v = ident * 100 if ident <= 1.0 else ident
            status.append("unmutated" if v >= cutoff else "mutated")
    df["ighv_status"] = status
    return df


# ---------------------------------------------------------------------------
# Lineage composition (BCR vs TCR, clonal vs polyclonal, vs background)
# ---------------------------------------------------------------------------
def compute_composition(df: pd.DataFrame,
                        total_input_reads: int | None,
                        clonal_threshold: float,
                        denominator: str = "total",
                        locus_clonality_min: float = 0.30) -> dict:
    """
    Partition total input reads into eight mutually-exclusive pools:
        clonal_IGH, clonal_IGK_kappa, clonal_IGL_lambda, polyclonal_B,
        clonal_TRB, clonal_TRG_gamma_delta, polyclonal_T, background.

    `clonal_threshold` is the per-clone fraction of locus reads above which
    a clonotype counts as clonal (default 5%). TRD reads are bundled with
    TRG (γδ T-cell rearrangement); TRA reads contribute to polyclonal T
    only, as TRA is rarely diagnostic of clonality on its own.

    Background = total_input_reads − sum(clonotype reads). If
    total_input_reads is None it falls back to the sum of clonotype reads
    (background = 0); the metrics emit a flag so the report can warn.
    """
    pools = {
        "clonal_IGH": 0, "clonal_IGK_kappa": 0, "clonal_IGL_lambda": 0,
        "polyclonal_B": 0,
        "clonal_TRB": 0, "clonal_TRG_gamma_delta": 0,
        "polyclonal_T": 0,
        "background": 0,
    }

    if df.empty:
        vdj_reads = 0
    else:
        vdj_reads = int(df["read_count"].sum())
        # Per-clone fraction of its locus (recomputed here so this function is self-contained)
        df = df.copy()
        df["_locus_total"]    = df.groupby("locus")["read_count"].transform("sum")
        df["_locus_fraction"] = df["read_count"] / df["_locus_total"]

        # Locus is clonal if either:
        #   (a) clonality_index >= locus_clonality_min, or
        #   (b) n_clonotypes == 1 AND n_reads >= SINGLE_CLONE_READS_MIN.
        # (b) covers the monoclonal case where clonality_index is undefined.
        SINGLE_CLONE_READS_MIN = 20
        locus_clonal_call = {}
        for locus in LOCI:
            cnts = df.loc[df["locus"] == locus, "read_count"].to_numpy()
            ci = clonality_index(cnts)
            multi_clone_clonal = (ci is not None
                                  and not (isinstance(ci, float) and ci != ci)
                                  and ci >= locus_clonality_min)
            single_clone_clonal = (cnts.size == 1 and int(cnts.sum()) >= SINGLE_CLONE_READS_MIN)
            locus_clonal_call[locus] = multi_clone_clonal or single_clone_clonal

        for _, r in df.iterrows():
            locus  = r["locus"]
            reads  = int(r["read_count"])
            is_clonal = (r["_locus_fraction"] >= clonal_threshold
                          and locus_clonal_call.get(locus, False))
            if locus == "IGH":
                pools["clonal_IGH" if is_clonal else "polyclonal_B"] += reads
            elif locus == "IGK":
                pools["clonal_IGK_kappa" if is_clonal else "polyclonal_B"] += reads
            elif locus == "IGL":
                pools["clonal_IGL_lambda" if is_clonal else "polyclonal_B"] += reads
            elif locus == "TRB":
                pools["clonal_TRB" if is_clonal else "polyclonal_T"] += reads
            elif locus in ("TRG", "TRD"):
                pools["clonal_TRG_gamma_delta" if is_clonal else "polyclonal_T"] += reads
            elif locus == "TRA":
                pools["polyclonal_T"] += reads  # TRA rarely diagnostic on its own

    denom_known = total_input_reads is not None and total_input_reads > 0

    # Two denominator modes:
    #   "total" — fractions are % of total input reads. Useful for repertoire
    #             panels where V(D)J reads dominate; on cancer-gene panels
    #             (e.g. CAPP-seq) it makes the background pool ~99% and the
    #             clonal pools disappear visually.
    #   "vdj"   — fractions are % of V(D)J-assigned reads only. Drops the
    #             background pool. Use when IG is a small panel target.
    if denominator == "vdj":
        denom = vdj_reads
        pools["background"] = 0
    else:
        denom = int(total_input_reads) if denom_known else vdj_reads
        pools["background"] = max(0, denom - vdj_reads)

    fractions = {k: (v / denom if denom else 0.0) for k, v in pools.items()}

    # κ:λ ratio — clinical light-chain restriction indicator
    if df.empty:
        igk_reads = igl_reads = 0
    else:
        igk_reads = int(df.loc[df["locus"] == "IGK", "read_count"].sum())
        igl_reads = int(df.loc[df["locus"] == "IGL", "read_count"].sum())
    kappa_lambda_ratio = (igk_reads / igl_reads) if igl_reads > 0 else None
    if kappa_lambda_ratio is None:
        kappa_lambda_call = "no_lambda_reads"
    elif kappa_lambda_ratio > KAPPA_LAMBDA_NORMAL_HIGH:
        kappa_lambda_call = "kappa_restricted"
    elif kappa_lambda_ratio < KAPPA_LAMBDA_NORMAL_LOW:
        kappa_lambda_call = "lambda_restricted"
    else:
        kappa_lambda_call = "balanced"

    return {
        "total_input_reads_known":   denom_known,
        "total_input_reads":         int(total_input_reads) if denom_known else vdj_reads,
        "vdj_assigned_reads":        vdj_reads,
        "denominator_mode":          denominator,
        "denominator_value":         denom,
        "clonal_dominance_threshold": clonal_threshold,
        "reads":                      pools,
        "fractions":                  fractions,
        "kappa_reads":                igk_reads,
        "lambda_reads":               igl_reads,
        "kappa_lambda_ratio":         kappa_lambda_ratio,
        "kappa_lambda_call":          kappa_lambda_call,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id",        required=True)
    ap.add_argument("--trust4-airr",      required=True, type=Path)
    ap.add_argument("--igblast-airr",     required=True, type=Path)
    ap.add_argument("--min-clone-count",  type=int,   default=2)
    ap.add_argument("--igh-mutated-cutoff", type=float, default=98.0)
    ap.add_argument("--top-n",            type=int,   default=50)
    ap.add_argument("--clonal-dominance-threshold", type=float, default=0.05,
                    help="Per-clone fraction of locus reads required to count as 'clonal' (default 0.05 = 5%%).")
    ap.add_argument("--total-input-reads", type=int, default=None,
                    help="Total reads in the input FASTQ/BAM (denominator for background fraction). "
                         "If omitted, background defaults to 0.")
    ap.add_argument("--composition-denominator", choices=["total", "vdj"], default="total",
                    help="Composition fraction denominator: 'total' (default — %% of total input "
                         "reads, suitable for repertoire panels) or 'vdj' (%% of V(D)J-assigned "
                         "reads only, suitable for cancer-gene panels where IG is a small target).")
    ap.add_argument("--filter-germline-rearrangements",
                    dest="filter_germline_rearrangements",
                    action="store_true", default=True,
                    help="Drop clones that look like germline-rearrangement / sterile V-J "
                         "transcripts. Detected by short V alignment, low V identity, or "
                         "junction = V Cys + J stem with no novel nt. Default ON.")
    ap.add_argument("--no-filter-germline-rearrangements",
                    dest="filter_germline_rearrangements",
                    action="store_false",
                    help="Disable the germline-rearrangement filter (keep raw TRUST4 output).")
    ap.add_argument("--read-length", type=int, default=150,
                    help="Sequencing read length (nt). Used to (a) compute the "
                         "cdr3_spanned_by_single_read flag and (b) auto-scale the "
                         "germline-filter V-match threshold. Default 150.")
    ap.add_argument("--germline-min-v-match",      type=int,   default=None,
                    help="V CIGAR match length below which a clone is flagged as "
                         "germline-rearrangement. If omitted, scales automatically "
                         "with --read-length (e.g. 100 nt for 150 bp reads).")
    ap.add_argument("--germline-min-v-identity",   type=float, default=85.0,
                    help="V identity (%%) at or above which a clone counts as "
                         "near-germline for the artefact filter. Default 85.")
    ap.add_argument("--germline-min-junction-diversity", type=int, default=3,
                    help="Min nt of junction not explained by V Cys + J stem (default 3).")
    ap.add_argument("--out-metrics",      required=True, type=Path)
    ap.add_argument("--out-clonotypes",   required=True, type=Path)
    ap.add_argument("--out-top",          required=True, type=Path)
    args = ap.parse_args(argv)

    trust4  = read_airr(args.trust4_airr)
    igblast = read_airr(args.igblast_airr)

    df = build_clonotype_table(trust4, igblast)

    # Apply min-count filter (guard for the empty case — TRUST4 may have
    # produced no clonotypes at all, e.g. on a 3'-biased library)
    if not df.empty and "read_count" in df.columns:
        df = df[df["read_count"] >= args.min_clone_count].copy()
    else:
        # Make sure df has the columns downstream code expects
        for col in ["locus", "read_count", "v_call", "d_call", "j_call",
                    "junction", "junction_aa", "v_cigar", "v_identity"]:
            if col not in df.columns:
                df[col] = pd.Series(dtype=object)

    # Resolve adaptive V-match threshold
    min_v_match = args.germline_min_v_match
    if min_v_match is None:
        min_v_match = adaptive_min_v_match(args.read_length)

    # Germline-rearrangement filter
    germline_filter_stats = None
    if args.filter_germline_rearrangements:
        df, germline_filter_stats = apply_germline_rearrangement_filter(
            df,
            min_v_match=min_v_match,
            min_v_identity=args.germline_min_v_identity,
            min_junction_diversity=args.germline_min_junction_diversity)
        germline_filter_stats["thresholds"]["min_v_match_adaptive"] = (args.germline_min_v_match is None)

    # Read-length-aware annotation: which CDR3s were directly observed on a
    # single read, vs reconstructed by assembly
    df = annotate_single_read_spanning(df, args.read_length)

    df = assign_ighv_mutation_status(df, args.igh_mutated_cutoff)

    # ---- per-locus metrics ------------------------------------------------
    per_locus = {}
    for locus in LOCI:
        sub = df[df["locus"] == locus]
        per_locus[locus] = summarise(sub["read_count"].to_numpy())

    # ---- aggregate (all loci pooled) -------------------------------------
    aggregate = summarise(df["read_count"].to_numpy())

    # ---- IGHV mutation summary -------------------------------------------
    igh = df[df["locus"] == "IGH"]
    if not igh.empty:
        total = igh["read_count"].sum()
        unmut = igh[igh["ighv_status"] == "unmutated"]["read_count"].sum()
        mut   = igh[igh["ighv_status"] == "mutated"]["read_count"].sum()
        unk   = igh[igh["ighv_status"] == "unknown"]["read_count"].sum()
        ighv_summary = {
            "cutoff_percent_v_identity": args.igh_mutated_cutoff,
            "reads_total":      int(total),
            "reads_unmutated":  int(unmut),
            "reads_mutated":    int(mut),
            "reads_unknown":    int(unk),
            "fraction_unmutated": float(unmut / total) if total else None,
            "dominant_status":  ("unmutated" if unmut > mut
                                 else "mutated" if mut > unmut
                                 else "indeterminate"),
        }
    else:
        ighv_summary = None

    composition = compute_composition(df, args.total_input_reads,
                                       args.clonal_dominance_threshold,
                                       denominator=args.composition_denominator)

    # Long-CDR3 / assembly-inferred summary
    if df.empty:
        cdr3_inference = {"n_clonotypes": 0, "n_single_read_spanned": 0, "n_assembly_inferred": 0}
    else:
        cdr3_inference = {
            "read_length":              args.read_length,
            "anchor_required_nt":       SINGLE_READ_SPAN_ANCHOR_NT,
            "n_clonotypes":             int(len(df)),
            "n_single_read_spanned":    int(df["cdr3_spanned_by_single_read"].sum()),
            "n_assembly_inferred":      int(df["assembly_inferred"].sum()),
            "dominant_clone_spanned":   bool(df.sort_values("read_count", ascending=False)
                                              .iloc[0]["cdr3_spanned_by_single_read"]) if len(df) else None,
        }

    metrics = {
        "sample_id":      args.sample_id,
        "min_clone_count": args.min_clone_count,
        "read_length":    args.read_length,
        "germline_rearrangement_filter": germline_filter_stats,
        "cdr3_inference": cdr3_inference,
        "per_locus":      per_locus,
        "aggregate":      aggregate,
        "ighv_status":    ighv_summary,
        "composition":    composition,
    }

    args.out_metrics.write_text(json.dumps(metrics, indent=2, default=str))

    # ---- clonotype tables -------------------------------------------------
    df_sorted = df.sort_values(["locus", "read_count"], ascending=[True, False])
    df_sorted.to_csv(args.out_clonotypes, sep="\t", index=False)

    top_rows = []
    for locus in LOCI:
        sub = df_sorted[df_sorted["locus"] == locus].head(args.top_n)
        top_rows.append(sub)
    pd.concat(top_rows, ignore_index=True).to_csv(args.out_top, sep="\t", index=False)

    print(f"[clonality_metrics] sample={args.sample_id} "
          f"clonotypes={len(df)} "
          f"aggregate_clonality={aggregate.get('clonality_index')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
