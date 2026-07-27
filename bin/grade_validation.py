#!/usr/bin/env python3
"""
grade_validation.py — compare a Lymphix result against expected ground truth.

Reads:
  * <sample>.metrics.json
  * <sample>.clonotypes.tsv
  * validation_expected.json (provides per-sample biology + expected calls)

Emits:
  * Per-sample TSV row + verdict-vs-expected pass/fail
  * Aggregate JSON summary

Designed to be run after each pipeline invocation in the validation loop.
"""
from __future__ import annotations
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

# The grader must apply exactly the rule the reports apply, or a green
# validation run says nothing about what a clinician will actually see. It
# previously carried two separate transcriptions of that rule — one to derive
# the category and a second, a few lines later, to list the clonal loci.
from lymphix_common import verdict_category   # noqa: E402


def grade_sample(sample_id: str, metrics_path: Path, clonotypes_path: Path,
                  expected: dict) -> dict:
    if not metrics_path.exists():
        return {"sample":                sample_id,
                "pass":                  False,
                "observed_category":     "(not_yet_processed)",
                "observed_clonal_loci":  [],
                "expected_category":     expected.get("expected_verdict_category"),
                "expected_clonal_loci_any_of": expected.get("expected_clonal_loci_any_of") or [],
                "reasons_pass":          [],
                "reasons_fail":          [f"metrics.json missing at {metrics_path}"],
                "n_clonotypes":          None,
                "top_clone_pct":         None,
                "vdj_reads":             None}
    metrics = json.loads(metrics_path.read_text())
    df = pd.DataFrame() if not clonotypes_path.exists() else \
         pd.read_csv(clonotypes_path, sep="\t")

    agg = metrics.get("aggregate") or {}
    comp = metrics.get("composition") or {}
    ighv = metrics.get("ighv_status") or {}

    # Observed verdict, from the same rule the reports use.
    vdj_reads = (comp or {}).get("vdj_assigned_reads", 0) or 0
    n_clones  = agg.get("n_clonotypes", 0) or 0
    obs_category, obs_clonal_loci = verdict_category(
        metrics.get("per_locus") or {}, vdj_reads, n_clones)

    reasons_pass, reasons_fail = [], []

    # ---- Check category ---------------------------------------------------
    exp_cat_single = expected.get("expected_verdict_category")
    exp_cat_any    = expected.get("expected_verdict_category_any_of") or []
    if exp_cat_single:
        if obs_category == exp_cat_single:
            reasons_pass.append(f"verdict category {obs_category} matches expected {exp_cat_single}")
        else:
            reasons_fail.append(f"verdict category {obs_category} but expected {exp_cat_single}")
    elif exp_cat_any:
        if obs_category in exp_cat_any:
            reasons_pass.append(f"verdict category {obs_category} matches one of expected {exp_cat_any}")
        else:
            reasons_fail.append(f"verdict category {obs_category} but expected any of {exp_cat_any}")
    exp_cat = exp_cat_single or (exp_cat_any[0] if exp_cat_any else None)

    # ---- Check clonal-loci subset matches expected ------------------------
    exp_any_of = expected.get("expected_clonal_loci_any_of") or []
    if exp_any_of:
        hit = [L for L in exp_any_of if L in obs_clonal_loci]
        if hit:
            reasons_pass.append(f"clonal loci {obs_clonal_loci} include expected {hit}")
        else:
            reasons_fail.append(f"clonal loci {obs_clonal_loci} do not include any of expected {exp_any_of}")

    # ---- Check top clones (V/J/CDR3 patterns) -----------------------------
    exp_clones = expected.get("expected_top_clones") or []
    for spec in exp_clones:
        locus_pat = spec.get("locus")
        v_pat = spec.get("v_call_pattern")
        j_pat = spec.get("j_call_pattern")
        cdr3_contains = spec.get("cdr3_aa_contains")
        if df.empty:
            reasons_fail.append(f"expected clone for {locus_pat} but no clonotypes in output")
            continue
        sub = df[df["locus"] == locus_pat] if locus_pat else df
        match = sub
        if v_pat:
            match = match[match["v_call"].fillna("").str.contains(v_pat, na=False, regex=False)]
        if j_pat:
            match = match[match["j_call"].fillna("").str.contains(j_pat, na=False, regex=False)]
        if cdr3_contains:
            match = match[match["junction_aa"].fillna("").str.contains(cdr3_contains, na=False, regex=False)]
        if not match.empty:
            top = match.iloc[0]
            reasons_pass.append(
                f"found expected clone @ {locus_pat}: V={top.get('v_call')} J={top.get('j_call')} "
                f"CDR3={top.get('junction_aa')}")
        else:
            reasons_fail.append(
                f"expected clone @ {locus_pat} V~{v_pat} J~{j_pat} CDR3~{cdr3_contains} not found")

    # ---- IGHV reportable when expected -----------------------------------
    exp_ighv = expected.get("expected_ighv_reportable")
    if exp_ighv is True:
        if ighv and ighv.get("reads_total", 0) > 0 and "IGH" in obs_clonal_loci:
            reasons_pass.append(f"IGHV reportable ({ighv.get('dominant_clone_status')})")
        else:
            reasons_fail.append("IGHV expected to be reportable but was not")
    elif exp_ighv is False:
        if ighv and ighv.get("reads_total", 0) > 0 and "IGH" in obs_clonal_loci:
            reasons_pass.append("IGHV reported (not expected, but not a failure)")

    # ---- Optional explicit IGHV mutation status check --------------------
    exp_ighv_status = expected.get("expected_ighv_status")
    if exp_ighv_status:
        actual = (ighv or {}).get("dominant_clone_status")
        if actual == exp_ighv_status:
            reasons_pass.append(f"IGHV mutation status {actual} matches expected {exp_ighv_status}")
        else:
            reasons_fail.append(f"IGHV mutation status {actual} but expected {exp_ighv_status}")

    # An empty checklist used to come back pass=True: if no expectation key
    # matched this sample, nothing was actually verified, and reporting that as
    # a pass green-lights a null result.
    if not reasons_pass and not reasons_fail:
        reasons_fail.append(
            "No expectation matched this sample — nothing was checked, so this "
            "cannot count as a pass."
        )

    overall_pass = not reasons_fail

    return {
        "sample":              sample_id,
        "pass":                overall_pass,
        "observed_category":   obs_category,
        "observed_clonal_loci": obs_clonal_loci,
        "expected_category":   exp_cat,
        "expected_clonal_loci_any_of": exp_any_of,
        "reasons_pass":        reasons_pass,
        "reasons_fail":        reasons_fail,
        "n_clonotypes":        n_clones,
        "top_clone_pct":       round((agg.get("top_clone_fraction") or 0) * 100, 2),
        "vdj_reads":           vdj_reads,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", required=True, type=Path,
                    help="Root containing one folder per sample with <sid>.metrics.json + .clonotypes.tsv")
    ap.add_argument("--expected",     required=True, type=Path,
                    help="validation_expected.json")
    ap.add_argument("--out-json",     required=True, type=Path)
    ap.add_argument("--out-tsv",      required=True, type=Path)
    args = ap.parse_args()

    spec = json.loads(args.expected.read_text())
    rows = []
    for sample_id, sample_spec in spec["samples"].items():
        # Find this sample's results directory
        candidate = args.results_root / f"{sample_id}_results"
        if not candidate.exists():
            candidate = args.results_root / sample_id
        m_path = candidate / f"{sample_id}.metrics.json"
        c_path = candidate / f"{sample_id}.clonotypes.tsv"
        rows.append(grade_sample(sample_id, m_path, c_path, sample_spec))

    pd.DataFrame(rows).to_csv(args.out_tsv, sep="\t", index=False)
    args.out_json.write_text(json.dumps({
        "n_samples":    len(rows),
        "n_pass":       sum(1 for r in rows if r["pass"]),
        "n_fail":       sum(1 for r in rows if not r["pass"]),
        "per_sample":   rows,
    }, indent=2))

    print(f"=== Validation grading ===")
    print(f"  {sum(1 for r in rows if r['pass'])} / {len(rows)} samples pass")
    for r in rows:
        flag = "PASS" if r["pass"] else "FAIL"
        print(f"  [{flag}] {r['sample']:14s} observed={r['observed_category']:14s} loci={r['observed_clonal_loci']}")
        for reason in r["reasons_fail"]:
            print(f"           FAIL: {reason}")
    return 0 if all(r["pass"] for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
