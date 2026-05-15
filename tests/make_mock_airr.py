"""
Generate AIRR-formatted TSV fixtures for testing the analysis layer
without needing TRUST4 / IgBLAST binaries.

Produces two samples per call:
    CLONAL_TEST     — one dominant IGH clone + one dominant TRB clone
    POLYCLONAL_TEST — ~600 distinct clones across all seven loci
"""
from __future__ import annotations
import argparse
import random
from pathlib import Path

import pandas as pd

LOCI = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]

V_GENES = {
    "IGH": ["IGHV1-2*02", "IGHV1-69*01", "IGHV3-23*01", "IGHV3-30*03",
            "IGHV4-34*01", "IGHV4-39*01", "IGHV5-51*01"],
    "IGK": ["IGKV1-5*01", "IGKV1-39*01", "IGKV3-20*01", "IGKV3-15*01"],
    "IGL": ["IGLV2-14*01", "IGLV1-44*01", "IGLV3-21*01", "IGLV1-40*01"],
    "TRA": ["TRAV1-1*01", "TRAV12-1*01", "TRAV21*01", "TRAV29/DV5*01"],
    "TRB": ["TRBV20-1*01", "TRBV5-1*01", "TRBV6-5*01", "TRBV28*01", "TRBV9*01"],
    "TRG": ["TRGV9*01", "TRGV2*01", "TRGV4*01"],
    "TRD": ["TRDV1*01", "TRDV2*01", "TRDV3*01"],
}
D_GENES = {
    "IGH": ["IGHD3-10*01", "IGHD2-2*01", "IGHD6-19*01"],
    "TRB": ["TRBD1*01", "TRBD2*01"],
    "TRD": ["TRDD1*01", "TRDD2*01", "TRDD3*01"],
}
J_GENES = {
    "IGH": ["IGHJ4*02", "IGHJ6*02", "IGHJ3*02"],
    "IGK": ["IGKJ1*01", "IGKJ2*01", "IGKJ4*01"],
    "IGL": ["IGLJ2*01", "IGLJ3*02"],
    "TRA": ["TRAJ23*01", "TRAJ42*01", "TRAJ49*01"],
    "TRB": ["TRBJ1-1*01", "TRBJ2-1*01", "TRBJ2-7*01"],
    "TRG": ["TRGJ1*01", "TRGJ2*01", "TRGJP1*01"],
    "TRD": ["TRDJ1*01", "TRDJ2*01"],
}


def rand_cdr3_aa(rng, length=None):
    if length is None:
        length = rng.randint(10, 22)
    return "C" + "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=length - 2)) + rng.choice("WF")


def make_clonotype(rng, locus, seq_id, read_count, v_identity_pct=None):
    aa  = rand_cdr3_aa(rng)
    nt  = "".join(rng.choices("ACGT", k=len(aa) * 3))
    v   = rng.choice(V_GENES[locus])
    j   = rng.choice(J_GENES[locus])
    d   = rng.choice(D_GENES.get(locus, [""])) if locus in D_GENES else ""
    v_stub = "".join(rng.choices("ACGT", k=240))
    j_stub = "".join(rng.choices("ACGT", k=40))
    full_seq = v_stub + nt + j_stub
    if v_identity_pct is None:
        v_identity_pct = round(rng.uniform(85.0, 100.0), 2)
    # Synthesize a realistic-looking V CIGAR so the germline-rearrangement
    # filter recognises this as a real clone (≥100 nt V alignment).
    v_cigar = f"0S250M{len(j_stub) + len(nt)}S2N"
    j_cigar = f"{240 + len(nt)}S50M0S"
    return dict(
        sequence_id=seq_id, sequence=full_seq, rev_comp="F",
        productive="T", vj_in_frame="T", stop_codon="F", complete_vdj="T",
        locus=locus, v_call=v, d_call=d, j_call=j, c_call="",
        junction=nt, junction_aa=aa,
        v_cigar=v_cigar, d_cigar="", j_cigar=j_cigar,
        consensus_count=read_count, duplicate_count=read_count,
        v_identity=v_identity_pct,
    )


def generate(mode, rng):
    seq_id = 0
    if mode == "clonal":
        # Dominant IGH clone (UNMUTATED)
        seq_id += 1
        yield make_clonotype(rng, "IGH", f"seq_{seq_id}", 12500, v_identity_pct=99.2)
        # Dominant TRB clone
        seq_id += 1
        yield make_clonotype(rng, "TRB", f"seq_{seq_id}", 6200)
        # Some mid-frequency clones
        for _ in range(8):
            seq_id += 1
            yield make_clonotype(rng, rng.choice(LOCI), f"seq_{seq_id}",
                                 rng.randint(50, 300))
        # Background
        for _ in range(400):
            seq_id += 1
            yield make_clonotype(rng, rng.choice(LOCI), f"seq_{seq_id}",
                                 rng.randint(2, 20))
    elif mode == "polyclonal":
        for _ in range(600):
            seq_id += 1
            yield make_clonotype(rng, rng.choice(LOCI), f"seq_{seq_id}",
                                 rng.randint(10, 80))


def write_trust4(records, out):
    cols = ["sequence_id", "sequence", "rev_comp", "productive", "vj_in_frame",
            "stop_codon", "complete_vdj", "locus", "v_call", "d_call", "j_call",
            "c_call", "junction", "junction_aa", "v_cigar", "d_cigar", "j_cigar",
            "v_identity", "consensus_count", "duplicate_count"]
    pd.DataFrame(records)[cols].to_csv(out, sep="\t", index=False)


def write_igblast(records, out):
    cols = ["sequence_id", "productive", "vj_in_frame", "stop_codon",
            "complete_vdj", "locus", "v_call", "d_call", "j_call",
            "junction", "junction_aa", "v_identity"]
    pd.DataFrame(records)[cols].to_csv(out, sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--seed",   type=int, default=42)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    for sample, mode in [("CLONAL_TEST", "clonal"),
                         ("POLYCLONAL_TEST", "polyclonal")]:
        records = list(generate(mode, rng))
        d = args.outdir / sample
        d.mkdir(exist_ok=True)
        write_trust4(records,  d / f"{sample}.trust4.airr.tsv")
        write_igblast(records, d / f"{sample}.igblast.airr.tsv")
        print(f"[mock] {sample}: {len(records)} clonotypes -> {d}/")


if __name__ == "__main__":
    main()
