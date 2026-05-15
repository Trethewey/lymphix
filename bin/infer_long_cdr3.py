#!/usr/bin/env python3
"""
infer_long_cdr3.py — confidence inference for assembly-reconstructed CDR3s.

For clonotypes where the CDR3 could not be spanned by a single sequencing read
(flagged `assembly_inferred=True` by clonality_metrics.py), this script
computes a 0–1 confidence score from six components and assigns a clinical
band (HIGH / MODERATE / LOW / VERY_LOW).

Components and weights:
    coverage           25%  reads supporting the assembled junction
    paired_bridging    20%  any R1/R2 pair where R1 anchors in V and R2 in J
    consensus_quality  20%  % agreement between contributing reads at the junction
    assembly_unique    15%  contig assembly produced one solution, not multiple
    junction_novelty   10%  N-region length (proxy for biological recombination)
    cohort_unique      10%  CDR3 nt is not shared with unrelated samples in cohort

Reads-realignment (coverage, paired_bridging, consensus_quality, assembly_unique)
requires BWA-aligning the source FASTQ to each clone's assembled contig. This
scaffold ships a "lightweight" mode that uses TRUST4's already-reported read
counts and identity scores instead — useful for cohorts where running BWA
per-clone is impractical. The full BWA-backed mode is gated by --mode full
and will be wired into the Nextflow pipeline behind --enable-cdr3-inference.

Outputs:
    <sample>.cdr3_confidence.tsv   per-clone confidence components + band
    <sample>.cdr3_confidence.json  per-sample summary stats
"""
from __future__ import annotations
import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Component weights (must sum to 1.0)
# ---------------------------------------------------------------------------
WEIGHTS = {
    "coverage":          0.25,
    "paired_bridging":   0.20,
    "consensus_quality": 0.20,
    "assembly_unique":   0.15,
    "junction_novelty":  0.10,
    "cohort_unique":     0.10,
}
assert abs(sum(WEIGHTS.values()) - 1.0) < 1e-9

BANDS = [
    (0.85, "HIGH"),
    (0.65, "MODERATE"),
    (0.40, "LOW"),
    (0.00, "VERY_LOW"),
]


def confidence_band(score: float) -> str:
    for thresh, label in BANDS:
        if score >= thresh:
            return label
    return "VERY_LOW"


# ---------------------------------------------------------------------------
# Lightweight component scorers — derive each metric from existing data
# (no per-clone BWA realignment). The full BWA-backed implementations live
# in the `score_*_full` functions further down (currently placeholders).
# ---------------------------------------------------------------------------
def score_coverage_lite(row: dict, k: int = 20) -> float:
    """Coverage proxy: log10(read_count) saturating at log10(k+1) = ~3 for k=1000."""
    n = max(0, int(row.get("read_count") or 0))
    return min(1.0, math.log10(n + 1) / math.log10(1000 + 1))


def score_paired_bridging_lite(row: dict, read_length: int) -> float:
    """Bridging proxy: 1 if the junction is short enough that a 250 bp PE
    fragment (the typical capture-prep insert) could plausibly span V-anchor
    → J-anchor; 0 otherwise. Real implementation uses paired-end alignment
    against the assembled contig."""
    junc_len = len(str(row.get("junction") or ""))
    fragment_span = 2 * read_length - 30  # PE overlap allowance
    return 1.0 if (junc_len + 30 <= fragment_span) else 0.0


def score_consensus_quality_lite(row: dict) -> float:
    """Use TRUST4's reported V identity as a proxy: high V identity means
    contributing reads agree well at every base."""
    try:
        v_id = float(row.get("v_identity") or 0.0)
    except (TypeError, ValueError):
        v_id = 0.0
    # 100% identity → 1.0, 80% → 0.0
    return max(0.0, min(1.0, (v_id - 80.0) / 20.0))


def score_assembly_unique_lite(row: dict) -> float:
    """Lightweight stand-in: penalise clones with very low junction diversity
    (likely assembled from sparse/ambiguous evidence)."""
    junc_len = len(str(row.get("junction") or ""))
    # Junctions much shorter than typical (~45-75 nt) suggest ambiguous assembly
    if junc_len < 30: return 0.3
    if junc_len < 45: return 0.7
    return 1.0


def score_junction_novelty(row: dict, j_stem_estimator) -> float:
    """N-region length: junction length minus the V Cys codon (3 nt) minus
    the estimated J-stem contribution. Real recombinations have ≥3 nt
    N-region; minimal-N junctions are biologically suspect."""
    junc = str(row.get("junction") or "")
    j_stem = j_stem_estimator(row.get("j_call"))
    n_region = max(0, len(junc) - 3 - j_stem)
    # 0 nt → 0.0,  3 nt → 0.5,  ≥10 nt → 1.0
    if n_region <= 0: return 0.0
    if n_region >= 10: return 1.0
    return n_region / 10.0


def score_cohort_unique(row: dict, junction_counts: Counter,
                        n_samples_in_cohort: int) -> float:
    """1.0 if the CDR3 nt is unique to one sample in the cohort; falls toward
    0 as the same nt sequence appears in more samples (suggesting an
    artefact). Vacuously 1.0 for single-sample runs."""
    if n_samples_in_cohort <= 1:
        return 1.0
    j = str(row.get("junction") or "")
    if not j:
        return 0.0
    n_with = junction_counts.get(j, 1)
    # 1 sample → 1.0, all samples → 0.0
    return max(0.0, 1.0 - (n_with - 1) / max(1, n_samples_in_cohort - 1))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
J_STEM_NT = {
    "IGHJ1": 30, "IGHJ2": 27, "IGHJ3": 30, "IGHJ4": 21, "IGHJ5": 21, "IGHJ6": 33,
    "IGKJ1": 27, "IGKJ2": 27, "IGKJ3": 27, "IGKJ4": 27, "IGKJ5": 27,
    "IGLJ1": 27, "IGLJ2": 27, "IGLJ3": 27, "IGLJ7": 27,
    "TRAJ":  21, "TRBJ":  30, "TRGJ":  21, "TRDJ":  24,
}


def j_stem_estimator(j_call) -> int:
    if not isinstance(j_call, str) or not j_call:
        return 21
    name = j_call.split("*")[0].split(",")[0]
    if name in J_STEM_NT:
        return J_STEM_NT[name]
    for prefix, length in J_STEM_NT.items():
        if name.startswith(prefix):
            return length
    return 21


def load_cohort_junction_counts(cohort_dir: Path | None) -> tuple[Counter, int]:
    """Scan every sample's clonotypes.tsv in cohort_dir to count how many
    samples carry each CDR3 nt sequence (cross-sample uniqueness signal)."""
    if cohort_dir is None or not cohort_dir.exists():
        return Counter(), 0
    counts = Counter()
    samples_seen = set()
    for p in cohort_dir.glob("*/*.clonotypes.tsv"):
        sid = p.parent.name
        try:
            df = pd.read_csv(p, sep="\t", dtype=str)
        except Exception:
            continue
        if df.empty or "junction" not in df.columns:
            continue
        for j in df["junction"].dropna().unique():
            counts[j] += 1
        samples_seen.add(sid)
    return counts, len(samples_seen)


# ---------------------------------------------------------------------------
# Full mode (BWA-realignment) — placeholder
# ---------------------------------------------------------------------------
def score_full_mode(*args, **kwargs):  # pragma: no cover
    """Full BWA-backed scoring. Not yet implemented in this scaffold.
    Requires:
      * TRUST4 *_assembled_reads.fa  (contigs per clone)
      * Source FASTQ
      * BWA on PATH
    For each clonotype:
      1. Align raw reads to the contig with bwa mem
      2. Count reads where the alignment spans the entire junction
      3. Check for paired-end pairs where R1 and R2 anchor in V/J
      4. Compute per-base consensus from aligned reads
      5. Detect alternative contig solutions (multi-path graph)
    Returns the same fields as the lightweight scorer but with empirical
    values. Designed to plug in as a drop-in replacement.
    """
    raise NotImplementedError("Full BWA-backed scoring is not yet wired in. "
                              "Use --mode lite for now (default).")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--clonotypes",    required=True, type=Path,
                    help="Path to <sample>.clonotypes.tsv produced by clonality_metrics.py.")
    ap.add_argument("--sample-id",     required=True)
    ap.add_argument("--read-length",   type=int, default=150)
    ap.add_argument("--cohort-dir",    type=Path, default=None,
                    help="Optional: directory of sibling sample folders to pull "
                         "cross-sample CDR3 uniqueness signal from.")
    ap.add_argument("--mode", choices=["lite", "full"], default="lite",
                    help="lite (default): score from existing TRUST4 fields. "
                         "full: BWA-realignment based scoring (NOT YET IMPLEMENTED).")
    ap.add_argument("--out-tsv",       required=True, type=Path)
    ap.add_argument("--out-json",      required=True, type=Path)
    args = ap.parse_args(argv)

    if args.mode == "full":
        score_full_mode()

    df = pd.read_csv(args.clonotypes, sep="\t", dtype=str)
    if df.empty:
        args.out_tsv.write_text("sample_id\tn_clonotypes\tnotes\n")
        args.out_json.write_text(json.dumps({"sample_id": args.sample_id,
                                              "n_clonotypes": 0,
                                              "mode": args.mode}, indent=2))
        print(f"[infer_long_cdr3] no clonotypes for {args.sample_id}")
        return 0

    df["read_count"] = pd.to_numeric(df.get("read_count"), errors="coerce").fillna(0).astype(int)
    df["v_identity"] = pd.to_numeric(df.get("v_identity"), errors="coerce")
    df["assembly_inferred"] = df.get("assembly_inferred", "False").astype(str).str.lower() == "true"

    # Cohort-wide junction counts (for cross-sample uniqueness)
    cohort_junctions, n_cohort_samples = load_cohort_junction_counts(args.cohort_dir)

    out_rows = []
    for _, row in df.iterrows():
        r = row.to_dict()
        components = {
            "coverage":          score_coverage_lite(r),
            "paired_bridging":   score_paired_bridging_lite(r, args.read_length),
            "consensus_quality": score_consensus_quality_lite(r),
            "assembly_unique":   score_assembly_unique_lite(r),
            "junction_novelty":  score_junction_novelty(r, j_stem_estimator),
            "cohort_unique":     score_cohort_unique(r, cohort_junctions, n_cohort_samples),
        }
        score = sum(WEIGHTS[k] * v for k, v in components.items())
        out_rows.append({
            "sample_id":         args.sample_id,
            "clone_v":           r.get("v_call"),
            "clone_j":           r.get("j_call"),
            "cdr3_aa":           r.get("junction_aa"),
            "read_count":        r["read_count"],
            "assembly_inferred": r["assembly_inferred"],
            "confidence_score":  round(score, 3),
            "confidence_band":   confidence_band(score),
            **{f"score_{k}": round(v, 3) for k, v in components.items()},
        })

    out_df = pd.DataFrame(out_rows).sort_values("read_count", ascending=False)
    out_df.to_csv(args.out_tsv, sep="\t", index=False)

    # Per-sample summary
    band_counts = Counter(out_df["confidence_band"])
    summary = {
        "sample_id":   args.sample_id,
        "mode":        args.mode,
        "read_length": args.read_length,
        "n_clonotypes": int(len(out_df)),
        "n_assembly_inferred": int(out_df["assembly_inferred"].sum()),
        "band_counts": dict(band_counts),
        "dominant_clone": (out_df.iloc[0].to_dict() if not out_df.empty else None),
        "weights": WEIGHTS,
        "cohort_dir": str(args.cohort_dir) if args.cohort_dir else None,
        "n_cohort_samples": n_cohort_samples,
    }
    args.out_json.write_text(json.dumps(summary, indent=2, default=str))

    print(f"[infer_long_cdr3] {args.sample_id} mode={args.mode}  "
          f"clones={len(out_df)}  bands={dict(band_counts)}")
    return 0


if __name__ == "__main__":
    main()
