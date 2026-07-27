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

from lymphix_common import (
    LOCI,
    LOCUS_CLONAL_INDEX_THRESHOLD,
    is_locus_clonal,
)

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
    """Read an AIRR TSV.

    A missing or zero-byte file is an upstream failure, not an empty
    repertoire. The two are indistinguishable once the frame is empty, and
    conflating them produced signed "no V(D)J signal detected" reports from
    runs where the input never arrived. A header-only file is the legitimate
    way to say "this sample genuinely yielded nothing".
    """
    if not path.exists():
        raise FileNotFoundError(
            f"AIRR input not found: {path}. If the sample genuinely yielded no "
            f"rearrangements, pass a header-only file rather than omitting it."
        )
    if path.stat().st_size == 0:
        raise ValueError(
            f"AIRR input is zero bytes: {path}. A completed run writes at least "
            f"a header line, so this indicates the upstream step failed. Pass a "
            f"header-only file if the empty repertoire is genuine."
        )
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
    for further confidence assessment downstream (see the unfinished scaffold
    at docs/scaffolds/infer_long_cdr3.py).
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


# ---------------------------------------------------------------------------
# Clonotype collapsing
# ---------------------------------------------------------------------------
# WHY THIS EXISTS
# ---------------
# TRUST4 emits one AIRR row per (assembly, CDR3 variant), not one row per clone.
# A single rearrangement therefore arrives as a dominant row plus a tail of
# near-identical rows carrying a handful of reads each. Counting those rows as
# clonotypes inflates N in every diversity metric — Shannon, Simpson, Gini, D50
# and the clonality index all take N or the read distribution directly — so a
# monoclonal sample reads as more diverse than it is, and the "dominant clone"
# may be a fragment rather than the clone.
#
# Two things were established on the real cohorts before this was written, and
# they shape the keys offered here:
#
#   1. Keying on V + J + junction nt merges nothing at all. Across 40 samples
#      with clonotypes, the number of distinct (V, J, junction nt) triples
#      equalled the row count in every single sample. The rows that need
#      merging differ in the junction itself, so any key containing the full
#      junction nt is a no-op against them.
#   2. Keying on the junction amino acid sequence is actively unsafe. Rows with
#      an untranslatable junction (out-of-frame, stop codon) carry a blank
#      junction_aa and all key to the empty string; 22 of those 40 samples had
#      more than one such row, one had nine. Merging them lumps unrelated
#      rearrangements together under "unknown".
#
# So collapsing is done in two stages, and the caller chooses how far to go.
COLLAPSE_KEY_EXACT    = "locus_junction_nt"
COLLAPSE_KEY_HAMMING1 = "locus_junction_nt_hamming1"
COLLAPSE_KEYS = (COLLAPSE_KEY_EXACT, COLLAPSE_KEY_HAMMING1)

# Abundance gate for the stage-2 near-neighbour merge. A minor row is absorbed
# only if it holds no more than this fraction of the anchor row's reads.
# Deliberately a parameter rather than a constant: nothing in the sequence
# distinguishes a sequencing error from a genuine somatic-hypermutation
# variant of the same clone — only relative abundance does, and where the line
# sits is a judgement about the assay, not a fact about the data. The observed
# error tail sits at 0.1-0.3% of its parent, so 2% clears it with room to spare
# while leaving a subclone at 5% of the dominant clone standing.
DEFAULT_COLLAPSE_MINOR_FRACTION = 0.02

# How a collapsed clone's read count is derived. Recorded in metrics.json
# because the choice is not neutral — see _aggregate_read_count().
COLLAPSE_READ_AGGREGATION = "sum"


def _assembly_of(sequence_id) -> str:
    """Return the TRUST4 assembly (contig) a row came from.

    TRUST4 names AIRR rows `<assembly>_<variant index>`, e.g. `assemble5_0`,
    `assemble5_1`. Everything after the final underscore is the variant index,
    so stripping it recovers the contig the variants were called against.
    Anything not shaped that way is treated as its own assembly, which is the
    conservative reading for non-TRUST4 input.
    """
    text = "" if sequence_id is None else str(sequence_id)
    if isinstance(sequence_id, float) and math.isnan(sequence_id):
        text = ""
    head, sep, tail = text.rpartition("_")
    return head if (sep and tail.isdigit() and head) else text


def _aggregate_read_count(sub: pd.DataFrame) -> int:
    """Read support for a collapsed clone: sum within an assembly, and across.

    Summing within one assembly is what TRUST4 intends. The per-variant
    abundances in `_cdr3.out` are fractional read weights that partition the
    reads spanning that CDR3 — a read ambiguous between two variants is split
    between them, not counted twice — so the sum reconstructs the assembly's
    read support.

    Summing ACROSS assemblies was initially rejected as unsafe, on the grounds
    that two contigs built from one rearrangement might share reads and adding
    them would invent support. That is now measured rather than assumed, on
    CMDL20001026_S127_L004, the sample the collapse affects most:

      * The two contigs really are one rearrangement. assemble18_0
        (IGKV3D-15*01, 315 reads) and assemble32_0 (IGKV3-15*01, 170 reads)
        carry an identical 33 nt junction, an identical junction_aa and the
        same IGKJ1*01; their sequences share a 334 nt block and are 80%
        identical overall. IGKV3-15 and IGKV3D-15 are the two copies of a
        duplicated gene, so this is one clone assembled twice against two
        near-identical references.
      * Reads are only lightly reused. Of the 158,660 read records in
        `_assembled_reads.fa`, 152,631 are distinct — 3.8% appear against more
        than one contig.

    So taking the maximum discards 170 of 485 reads, understating this clone by
    35%, to avoid an overcount bounded at roughly 4%. That trade is the wrong
    way round, and the understatement is not harmless: the discarded reads also
    disappear from aggregate.n_reads, the per-locus totals and the composition
    pools, so the sample loses reads it genuinely has.

    The residual overcount is real but small, and it is recorded: the collapse
    audit block in metrics.json reports reads in and out, so any inflation is
    visible rather than silent. TRUST4 writes no read-to-contig map, so an
    exact figure is not recoverable from its outputs.
    """
    per_assembly = sub.groupby("_assembly", sort=False)["read_count"].sum()
    return int(per_assembly.sum())


def _best_identity(sub: pd.DataFrame, column: str):
    """Pick one V identity for a collapsed clone from its member rows.

    `sub` must already be sorted with the dominant row first.

    Preference order: the dominant row's own value, then the value from the
    highest-read member that carries one. Never a mean, and never a
    read-weighted mean: averaging a 94% variant with a 100% variant lands on
    97%, a number belonging to neither row and sitting on the mutated side of
    the 98% IGHV cutoff. That number goes straight into a prognostic call.

    Borrowing a minor row's identity when the dominant row has none is sound
    because the members of a collapsed clone differ in the junction, not in the
    V region — they are the same rearrangement seen through assembly noise.
    """
    if column not in sub.columns:
        return np.nan
    values = pd.to_numeric(sub[column], errors="coerce")
    present = values.dropna()
    return float(present.iloc[0]) if not present.empty else np.nan


def _cell_text(value) -> str:
    """Normalise a table cell to a string, mapping missing values to ''.

    pandas leaves an absent junction or gene call as NaN, which is truthy and
    stringifies to the literal 'nan'. Left unnormalised that is a
    three-character "sequence" that can be compared, bucketed and merged like
    any other.
    """
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value)
    return "" if text == "nan" else text


def _hamming_is_one(a: str, b: str) -> bool:
    """True if a and b are the same length and differ at exactly one position."""
    if len(a) != len(b):
        return False
    diff = 0
    for x, y in zip(a, b):
        if x != y:
            diff += 1
            if diff > 1:
                return False
    return diff == 1


def _exact_groups(df: pd.DataFrame) -> list[list[int]]:
    """Stage 1 — group rows sharing a locus and an identical junction nt.

    The V call is deliberately excluded. The only exact duplicates in the real
    data are one rearrangement assembled twice against duplicate V references
    (IGKV3-15 vs IGKV3D-15, IGKV2-28 vs IGKV2D-29, the IGLV5-37/45/48/52
    family): same junction nt, same junction aa, same J, different V paralogue.
    Keying on V keeps them apart, which is exactly the split we are here to
    close — in one sample they are the top two "clones".

    Rows with a blank junction nt are never grouped. They would otherwise all
    key to the empty string and merge unrelated rearrangements together, the
    same trap that rules out keying on junction_aa.
    """
    groups: dict[tuple, list[int]] = {}
    singletons: list[list[int]] = []
    for pos, (locus, junction) in enumerate(zip(df["locus"], df["junction"])):
        junc = _cell_text(junction)
        if not junc:
            singletons.append([pos])
            continue
        groups.setdefault((locus, junc), []).append(pos)
    return list(groups.values()) + singletons


def _absorb_near_neighbours(df: pd.DataFrame,
                            groups: list[list[int]],
                            reads: list[int],
                            minor_fraction_max: float) -> list[list[int]]:
    """Stage 2 — absorb single-substitution minor variants into their parent.

    Candidates are compared only within the same locus, J gene and junction
    length; a length change is an indel or a different rearrangement, not the
    substitution noise this addresses.

    Two rules keep this from over-merging:

      * the abundance gate — a group is absorbed only if it holds at most
        `minor_fraction_max` of the anchor's reads. Somatic hypermutation
        produces genuine lineage members one nucleotide from the parent, and
        nothing in the sequence tells them apart from an error; only the
        relative abundance does.
      * every absorbed group must be within one substitution OF THE ANCHOR,
        not of any group already absorbed. Single-linkage chaining walks
        A-B-C and merges sequences two or more substitutions apart; anchored
        comparison cannot.

    Anchors are taken in descending read order, so the largest group in a
    neighbourhood always acts as the parent and an absorbed group can never
    become an anchor itself.
    """
    order = sorted(range(len(groups)), key=lambda i: (-reads[i], i))
    anchor_of: dict[int, int] = {}
    absorbed: set[int] = set()
    anchors: set[int] = set()

    # Bucket by (locus, J gene, junction length) so the pairwise scan stays
    # local rather than quadratic over the whole table.
    buckets: dict[tuple, list[int]] = {}
    facts: dict[int, tuple] = {}
    for gi, members in enumerate(groups):
        top = members[0]                      # members are read-sorted already
        junc = _cell_text(df["junction"].iloc[top])
        j_gene = _cell_text(df["j_call"].iloc[top]).split("*")[0].split(",")[0]
        locus = df["locus"].iloc[top]
        facts[gi] = (junc, (locus, j_gene, len(junc)))
        buckets.setdefault((locus, j_gene, len(junc)), []).append(gi)

    for gi in order:
        if gi in absorbed:
            continue
        junc, bucket_key = facts[gi]
        if not junc:
            continue
        anchors.add(gi)
        for other in buckets.get(bucket_key, ()):
            # A group that has already acted as an anchor is never absorbed:
            # otherwise its own absorbed members would be orphaned and their
            # reads would vanish from the table.
            if other == gi or other in absorbed or other in anchors:
                continue
            if reads[other] > reads[gi]:
                continue
            if reads[other] > minor_fraction_max * reads[gi]:
                continue
            if _hamming_is_one(junc, facts[other][0]):
                anchor_of[other] = gi
                absorbed.add(other)

    merged: list[list[int]] = []
    for gi, members in enumerate(groups):
        if gi in absorbed:
            continue
        combined = list(members)
        for other, anchor in anchor_of.items():
            if anchor == gi:
                combined.extend(groups[other])
        merged.append(combined)
    return merged


def collapse_clonotype_rows(df: pd.DataFrame, *,
                            key: str = COLLAPSE_KEY_EXACT,
                            minor_fraction_max: float = DEFAULT_COLLAPSE_MINOR_FRACTION
                            ) -> tuple[pd.DataFrame, dict]:
    """Collapse TRUST4 assembly-variant rows into clones.

    Returns (collapsed_df, stats). Every column of the input survives, taken
    from the dominant (highest-read) member of each clone, except:
      * `read_count`, replaced by the aggregate described in
        _aggregate_read_count();
      * `v_identity` / `igblast_v_identity`, filled from the best available
        member as described in _best_identity();
      * `n_collapsed_rows`, added — how many AIRR rows the clone came from.

    The dominant member also supplies the V, D, J, junction and junction_aa,
    on the same reasoning: it is the best-supported assembly of the clone, so
    where the variants disagree its call is the one with the most reads behind
    it. Ties are broken on sequence_id so the output is deterministic.
    """
    if key not in COLLAPSE_KEYS:
        raise ValueError(f"unknown collapse key: {key!r} (expected one of {COLLAPSE_KEYS})")

    if df.empty:
        return df, collapse_stats(df, applied=True, key=key,
                                  minor_fraction_max=minor_fraction_max)

    work = df.reset_index(drop=True).copy()
    work["_assembly"] = work["sequence_id"].map(_assembly_of) if "sequence_id" in work.columns \
        else pd.Series([""] * len(work), index=work.index)

    # Sort dominant-first, deterministically, then work on positions. Every
    # downstream step relies on members[0] being the dominant row.
    work["_seq_id_str"] = (work["sequence_id"].astype(str)
                           if "sequence_id" in work.columns else "")
    work = work.sort_values(["read_count", "_seq_id_str"], ascending=[False, True],
                            kind="mergesort").reset_index(drop=True)

    groups = _exact_groups(work)
    groups = [sorted(g) for g in groups]      # positions are already read-sorted
    group_reads = [_aggregate_read_count(work.iloc[g]) for g in groups]

    if key == COLLAPSE_KEY_HAMMING1:
        groups = _absorb_near_neighbours(work, groups, group_reads, minor_fraction_max)
        groups = [sorted(g) for g in groups]

    rows = []
    for members in groups:
        sub = work.iloc[members]
        row = sub.iloc[0].to_dict()           # dominant member supplies V/D/J/junction
        row["read_count"] = _aggregate_read_count(sub)
        for column in ("v_identity", "igblast_v_identity"):
            if column in work.columns:
                row[column] = _best_identity(sub, column)
        row["n_collapsed_rows"] = len(members)
        rows.append(row)

    out = pd.DataFrame(rows, columns=list(work.columns) + ["n_collapsed_rows"])
    out = out.drop(columns=["_assembly", "_seq_id_str"])
    out = out.sort_values(["read_count"], ascending=False,
                          kind="mergesort").reset_index(drop=True)

    stats = collapse_stats(df, applied=True, key=key,
                           minor_fraction_max=minor_fraction_max,
                           collapsed=out)
    return out, stats


def collapse_stats(df: pd.DataFrame, *,
                   applied: bool,
                   key: str | None = None,
                   minor_fraction_max: float | None = None,
                   collapsed: pd.DataFrame | None = None) -> dict:
    """Describe a collapse (or the absence of one) for metrics.json.

    Emitted whether or not collapsing ran, so that a metrics.json is never
    ambiguous about which convention produced its clonotype counts.
    """
    if collapsed is None:
        collapsed = df
    rows_in = 0 if df is None or df.empty else int(len(df))
    clones_out = 0 if collapsed is None or collapsed.empty else int(len(collapsed))

    # Reads in vs out. Collapsing must not quietly change a sample's read
    # total: a shortfall means support was discarded, a surplus means reads
    # shared between contigs were counted twice. Neither is visible from the
    # clonotype counts alone, so both are recorded here and expected to be
    # equal under the "sum" rule.
    reads_in = 0 if df is None or df.empty else int(df["read_count"].sum())
    reads_out = (0 if collapsed is None or collapsed.empty
                 else int(collapsed["read_count"].sum()))

    per_locus = {}
    for locus in LOCI:
        n_in = 0 if df is None or df.empty else int((df["locus"] == locus).sum())
        n_out = 0 if collapsed is None or collapsed.empty else int((collapsed["locus"] == locus).sum())
        per_locus[locus] = {"rows_in": n_in, "clones_out": n_out,
                            "rows_merged": n_in - n_out}
    return {
        "applied":            bool(applied),
        "key":                key if applied else None,
        "read_aggregation":   COLLAPSE_READ_AGGREGATION if applied else None,
        "minor_fraction_max": (float(minor_fraction_max)
                               if applied and key == COLLAPSE_KEY_HAMMING1 else None),
        "n_rows_in":          rows_in,
        "n_clones_out":       clones_out,
        "n_rows_merged":      rows_in - clones_out,
        "n_reads_in":         reads_in,
        "n_reads_out":        reads_out,
        "n_reads_delta":      reads_out - reads_in,
        "per_locus":          per_locus,
    }


def _json_safe(obj):
    """Recursively replace NaN/Infinity with None so json.dumps stays strict."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, float) and not math.isfinite(obj):
        return None
    if isinstance(obj, (np.floating, np.integer)):
        return _json_safe(obj.item())
    return obj


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
                        locus_clonality_min: float = LOCUS_CLONAL_INDEX_THRESHOLD) -> dict:
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

        # Which loci carry a clonal expansion. The rule is the shared one in
        # lymphix_common, with require_dominance=False: the per-clone gate on
        # `clonal_threshold` below already tests dominance for every read being
        # binned, so applying the repertoire-level top-clone test here as well
        # would demand dominance twice and empty the clonal pools of samples
        # the verdict still calls clonal. See is_locus_clonal() for the full
        # reasoning — this is the only caller allowed to pass that flag.
        locus_clonal_call = {}
        for locus in LOCI:
            cnts = df.loc[df["locus"] == locus, "read_count"].to_numpy()
            locus_clonal_call[locus] = is_locus_clonal(
                clonality_index=clonality_index(cnts),
                n_clonotypes=cnts.size,
                n_reads=int(cnts.sum()),
                require_dominance=False,
                clonality_min=locus_clonality_min,
            )

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
    # κ:λ call. Three sources of an undefined ratio used to all land in
    # "no_lambda_reads" (misleading: implied a λ deficit when in fact most
    # of the time it meant "panel doesn't cover light chains at all").
    # Discriminate them now so cohort-level QC can distinguish "uninterpretable"
    # from "kappa-only with no lambda counter-reads".
    if igk_reads == 0 and igl_reads == 0:
        kappa_lambda_ratio = None
        kappa_lambda_call  = "no_light_chain_reads"      # panel or signal limit
    elif igl_reads == 0:
        kappa_lambda_ratio = None
        kappa_lambda_call  = "kappa_only"                # strong κ skew OR IGL-uncaptured panel
    elif igk_reads == 0:
        kappa_lambda_ratio = 0.0
        kappa_lambda_call  = "lambda_restricted"         # extreme λ dominance
    else:
        kappa_lambda_ratio = igk_reads / igl_reads
        if   kappa_lambda_ratio > KAPPA_LAMBDA_NORMAL_HIGH: kappa_lambda_call = "kappa_restricted"
        elif kappa_lambda_ratio < KAPPA_LAMBDA_NORMAL_LOW:  kappa_lambda_call = "lambda_restricted"
        else:                                              kappa_lambda_call = "balanced"

    return {
        "total_input_reads_known":   denom_known,
        # None, not vdj_reads: substituting the V(D)J count for the library
        # size makes the V(D)J fraction identically 1.0, which silently
        # disables the capture-underperformance warning downstream.
        "total_input_reads":         int(total_input_reads) if denom_known else None,
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
_SENTINEL = object()  # marker for "user did not set this argument on CLI"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id",        required=True)
    ap.add_argument("--trust4-airr",      required=True, type=Path)
    ap.add_argument("--igblast-airr",     required=True, type=Path)
    ap.add_argument("-c", "--clones", "--min-clone-count",
                    dest="min_clone_count", type=int, default=2,
                    help="Drop clonotypes with fewer than N supporting reads (default 2). "
                         "At WGS depth, 1 surfaces sub-threshold clonotypes but raises the "
                         "noise floor; 3+ is stricter than the default CAPP-seq tuning. "
                         "Set what fits your data — the pipeline does not decide for you.")
    ap.add_argument("--igh-mutated-cutoff", type=float, default=98.0)
    ap.add_argument("--top-n",            type=int,   default=50)
    ap.add_argument("--clonal-dominance-threshold", type=float, default=0.05,
                    help="Per-clone fraction of locus reads required to count as 'clonal' (default 0.05 = 5%%).")
    ap.add_argument("--total-input-reads", type=int, default=None,
                    help="Total reads in the input FASTQ/BAM (denominator for background fraction). "
                         "If omitted, background defaults to 0.")
    ap.add_argument("--composition-denominator", choices=["total", "vdj"], default=_SENTINEL,
                    help="Composition fraction denominator: 'total' (default — %% of total input "
                         "reads, suitable for repertoire panels) or 'vdj' (%% of V(D)J-assigned "
                         "reads only, suitable for cancer-gene panels where IG is a small target; "
                         "default under --wgs).")
    ap.add_argument("--wgs", action="store_true",
                    help="WGS preset for whole-genome data (~30-40x per-position; "
                         "yields ~150-300 V(D)J reads vs CAPP-seq's thousands). "
                         "Sets germline-min-v-match=60 (WGS reads less reliably span the V "
                         "region) and composition-denominator=vdj. Does NOT change "
                         "--clones / min-clone-count — set that explicitly with -c.")
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
    ap.add_argument("--collapse-clonotypes", action="store_true", default=False,
                    help="Collapse TRUST4 assembly-variant rows into clones before any "
                         "metric is computed. DEFAULT OFF: it changes n_clonotypes and "
                         "therefore every diversity metric, so it must be switched on "
                         "deliberately and the choice is recorded in metrics.json. "
                         "Measure the effect on your own data before adopting it.")
    ap.add_argument("--collapse-key", choices=list(COLLAPSE_KEYS),
                    default=COLLAPSE_KEY_EXACT,
                    help="Which rows count as the same clone. "
                         "'locus_junction_nt' (default) merges only rows with an "
                         "identical junction nt at the same locus — exact, no "
                         "judgement, catches one rearrangement assembled twice against "
                         "paralogous V references. "
                         "'locus_junction_nt_hamming1' additionally absorbs "
                         "single-substitution minor variants into their parent clone "
                         "subject to --collapse-minor-fraction; that is the merge that "
                         "addresses the assembly-noise tail, and it is a judgement "
                         "about error versus subclone. Ignored without "
                         "--collapse-clonotypes.")
    ap.add_argument("--collapse-minor-fraction", type=float,
                    default=DEFAULT_COLLAPSE_MINOR_FRACTION,
                    help="Abundance gate for the near-neighbour merge: a single-"
                         "substitution variant is absorbed only if it holds at most "
                         "this fraction of its parent's reads (default 0.02 = 2%%). "
                         "Nothing in the sequence separates a sequencing error from a "
                         "real hypermutated subclone — only abundance does — so raising "
                         "this absorbs genuine subclones and lowering it leaves noise "
                         "behind. Only used with --collapse-key "
                         "locus_junction_nt_hamming1.")
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

    # Resolve --wgs preset defaults. Sentinel == "not set on CLI".
    wgs = bool(args.wgs)
    if args.composition_denominator is _SENTINEL:
        args.composition_denominator = "vdj" if wgs else "total"
    if wgs and args.germline_min_v_match is None:
        args.germline_min_v_match = 60   # WGS reads less reliably span the V region

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

    # Clonotype collapsing. Deliberately placed here — after the row-level
    # filters, before everything else — so that the min-count and germline
    # filters still see raw TRUST4 rows (which is what their thresholds were
    # tuned against) while every metric below sees the same collapsed set. No
    # metric is computed on the uncollapsed table.
    if args.collapse_clonotypes:
        df, collapse_summary = collapse_clonotype_rows(
            df,
            key=args.collapse_key,
            minor_fraction_max=args.collapse_minor_fraction)
    else:
        collapse_summary = collapse_stats(df, applied=False)

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
    # IGHV mutation status is a property of the tumour clone, not of the
    # repertoire: the CLL convention grades the dominant IGH rearrangement
    # against the germline-identity cutoff. A read-weighted majority over every
    # IGH clonotype answers a different question, and a polyclonal naive-B tail
    # (which is unmutated by definition) can outvote a genuinely mutated
    # tumour clone and invert the prognostic call.
    igh = df[df["locus"] == "IGH"]
    if not igh.empty:
        total = int(igh["read_count"].sum())
        unmut = int(igh[igh["ighv_status"] == "unmutated"]["read_count"].sum())
        mut   = int(igh[igh["ighv_status"] == "mutated"]["read_count"].sum())
        unk   = int(igh[igh["ighv_status"] == "unknown"]["read_count"].sum())
        assessed = unmut + mut

        dominant = igh.sort_values("read_count", ascending=False).iloc[0]
        dom_status = dominant["ighv_status"] or "unknown"
        dom_identity = dominant.get("igblast_v_identity")
        dom_identity = None if pd.isna(dom_identity) else float(dom_identity)
        dom_reads = int(dominant["read_count"])

        ighv_summary = {
            "cutoff_percent_v_identity": args.igh_mutated_cutoff,
            "reads_total":      total,
            "reads_unmutated":  unmut,
            "reads_mutated":    mut,
            "reads_unknown":    unk,
            "reads_assessed":   assessed,

            # The clinical call. "unknown" means the dominant clone carries no
            # IgBLAST identity and its status was not established — it must not
            # be reported as either mutated or unmutated.
            "dominant_clone_status":      dom_status,
            "dominant_clone_v_identity":  dom_identity,
            "dominant_clone_reads":       dom_reads,
            "dominant_clone_fraction":    float(dom_reads / total) if total else None,

            # Descriptive repertoire-wide tally over assessed reads only.
            # Not the clinical call; unknown reads are excluded from the
            # denominator so the fraction cannot be diluted by unassessed ones.
            "repertoire_unmutated_read_fraction":
                float(unmut / assessed) if assessed else None,
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
        "wgs_mode":       wgs,
        "min_clone_count": args.min_clone_count,
        # Flat mirrors of clonotype_collapse.applied / .key. n_clonotypes and
        # every diversity metric mean different things depending on these two
        # values, so they sit alongside the other run-defining settings where a
        # reader (and the flattened cohort table) cannot miss them.
        "collapse_clonotypes": bool(collapse_summary["applied"]),
        "collapse_key":        collapse_summary["key"],
        "clonotype_collapse":  collapse_summary,
        "read_length":    args.read_length,
        "germline_rearrangement_filter": germline_filter_stats,
        "cdr3_inference": cdr3_inference,
        "per_locus":      per_locus,
        "aggregate":      aggregate,
        "ighv_status":    ighv_summary,
        "composition":    composition,
    }

    # json.dumps emits bare NaN/Infinity tokens by default, which no strict
    # JSON parser accepts — and a monoclonal sample legitimately has an
    # undefined clonality_index. Map them to null so metrics.json is always
    # valid JSON, and so consumers see "not defined" rather than 0.0.
    metrics = _json_safe(metrics)
    args.out_metrics.write_text(
        json.dumps(metrics, indent=2, default=str, allow_nan=False)
    )

    # ---- clonotype tables -------------------------------------------------
    df_sorted = df.sort_values(["locus", "read_count"], ascending=[True, False])
    df_sorted.to_csv(args.out_clonotypes, sep="\t", index=False)

    top_rows = []
    for locus in LOCI:
        sub = df_sorted[df_sorted["locus"] == locus].head(args.top_n)
        top_rows.append(sub)
    pd.concat(top_rows, ignore_index=True).to_csv(args.out_top, sep="\t", index=False)

    collapse_note = ""
    if collapse_summary["applied"]:
        collapse_note = (f"collapsed={collapse_summary['n_rows_in']}"
                         f"->{collapse_summary['n_clones_out']} "
                         f"key={collapse_summary['key']} ")
    print(f"[clonality_metrics] sample={args.sample_id} "
          f"clonotypes={len(df)} "
          f"{collapse_note}"
          f"aggregate_clonality={aggregate.get('clonality_index')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
