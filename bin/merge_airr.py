#!/usr/bin/env python3
"""Concatenate multiple AIRR TSVs, harmonising columns."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True, type=Path)
    ap.add_argument("--out",    required=True, type=Path)
    args = ap.parse_args(argv)

    frames = []
    for p in args.inputs:
        if p.exists() and p.stat().st_size > 0:
            try:
                df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
                frames.append(df)
            except pd.errors.EmptyDataError:
                continue

    if not frames:
        args.out.write_text("")
        return 0

    merged = pd.concat(frames, ignore_index=True, sort=False)
    # Prefer rows where v_call is set when duplicate sequence_id
    if "sequence_id" in merged.columns:
        merged["__has_v"] = merged.get("v_call", "").fillna("").astype(bool).astype(int)
        merged = (merged.sort_values("__has_v", ascending=False)
                        .drop_duplicates("sequence_id")
                        .drop(columns=["__has_v"]))
    merged.to_csv(args.out, sep="\t", index=False)
    print(f"[merge_airr] wrote {len(merged)} rows -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
