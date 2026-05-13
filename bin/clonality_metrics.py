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


LOCI     = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
BCR_LOCI = ["IGH", "IGK", "IGL"]
TCR_LOCI = ["TRA", "TRB", "TRG", "TRD"]

# Clinical κ:λ ratio reference range (light-chain restriction flag)
KAPPA_LAMBDA_NORMAL_LOW  = 0.5
KAPPA_LAMBDA_NORMAL_HIGH = 2.5


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
    return df.reset_index(drop=True)


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
                        clonal_threshold: float) -> dict:
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
        for _, r in df.iterrows():
            locus  = r["locus"]
            reads  = int(r["read_count"])
            is_clonal = (r["_locus_fraction"] >= clonal_threshold)
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
        "total_input_reads":         denom,
        "vdj_assigned_reads":        vdj_reads,
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
    ap.add_argument("--out-metrics",      required=True, type=Path)
    ap.add_argument("--out-clonotypes",   required=True, type=Path)
    ap.add_argument("--out-top",          required=True, type=Path)
    args = ap.parse_args(argv)

    trust4  = read_airr(args.trust4_airr)
    igblast = read_airr(args.igblast_airr)

    df = build_clonotype_table(trust4, igblast)

    # Apply min-count filter
    df = df[df["read_count"] >= args.min_clone_count].copy()

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
                                       args.clonal_dominance_threshold)

    metrics = {
        "sample_id":      args.sample_id,
        "min_clone_count": args.min_clone_count,
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
