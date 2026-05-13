#!/usr/bin/env python3
"""
airr_to_fasta.py

Extract per-clonotype nucleotide sequences from a TRUST4 AIRR table and
emit a FASTA file for IgBLAST input. Uses `sequence` if available, else
falls back to `sequence_alignment` / `junction` flanks.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--airr", required=True, type=Path)
    ap.add_argument("--out",  required=True, type=Path)
    ap.add_argument("--seq-col", default="sequence",
                    help="Column with the V(D)J nucleotide sequence")
    args = ap.parse_args(argv)

    if not args.airr.exists() or args.airr.stat().st_size == 0:
        args.out.write_text("")  # empty FASTA
        return 0

    df = pd.read_csv(args.airr, sep="\t", dtype=str, low_memory=False)
    seq_col = args.seq_col if args.seq_col in df.columns else None
    if seq_col is None:
        for cand in ("sequence", "sequence_alignment", "junction"):
            if cand in df.columns:
                seq_col = cand; break
    if seq_col is None:
        print("[airr_to_fasta] No sequence column found", file=sys.stderr)
        args.out.write_text("")
        return 0

    n = 0
    with open(args.out, "w") as fh:
        for _, row in df.iterrows():
            seq = row.get(seq_col)
            sid = row.get("sequence_id")
            if not isinstance(seq, str) or not seq or not sid:
                continue
            fh.write(f">{sid}\n{seq.replace('-', '').replace('.', '')}\n")
            n += 1
    print(f"[airr_to_fasta] wrote {n} sequences to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
