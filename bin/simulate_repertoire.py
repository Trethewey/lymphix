#!/usr/bin/env python3
"""
simulate_repertoire.py

Generate synthetic paired-end 2x150bp FASTQ data exercising the BCR/TCR
clonality pipeline end-to-end. Uses TRUST4's bundled IMGT+C reference so
the V/D/J segments are real germline sequences that TRUST4 will recognise.

Two modes:
    --mode clonal       One dominant IGH + one dominant TRB clone at ~70%,
                        rest distributed across ~50 background clones.
    --mode polyclonal   ~2000 distinct clones uniformly across all loci.

The output is two gzipped FASTQs (R1, R2) at the requested total read count.
"""
from __future__ import annotations
import argparse
import gzip
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

from lymphix_common import LOCI      # noqa: E402  (see bin/lymphix_common.py)

SEG_TYPES  = ("V", "D", "J", "C")
QUAL_STR   = chr(33 + 35)  # Q35 ASCII
COMP       = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------------------
# IMGT reference parsing
# ---------------------------------------------------------------------------
GENE_RE = re.compile(r"^([A-Z]{2,4})([VDJC])")   # e.g. IGHV1-2*02 -> ("IGH","V")


def parse_imgt_fasta(path: Path) -> dict:
    """Return {locus: {segment_type: [(name, seq), ...]}}."""
    out: dict = {L: {S: [] for S in SEG_TYPES} for L in LOCI}
    name, locus, seg, seq = None, None, None, []

    def flush():
        if name and locus in out and seg in out[locus]:
            s = "".join(seq).upper().replace(".", "").replace("-", "")
            if len(s) >= 20:
                out[locus][seg].append((name, s))

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
                first  = header.split("|")[0]
                name   = first
                m      = GENE_RE.match(first)
                if m:
                    locus, seg = m.group(1), m.group(2)
                else:
                    locus, seg = None, None
                seq = []
            else:
                seq.append(line)
        flush()
    return out


# ---------------------------------------------------------------------------
# Synthetic V(D)J construct
# ---------------------------------------------------------------------------
def random_ndn(rng: random.Random, min_len: int = 4, max_len: int = 20) -> str:
    n = rng.randint(min_len, max_len)
    return "".join(rng.choices("ACGT", k=n))


def build_clone(ref: dict, locus: str, rng: random.Random) -> tuple[str, dict]:
    """Return a single recombined V[+D]+NDN+J[+C] sequence and metadata."""
    v_pool = ref[locus]["V"]
    j_pool = ref[locus]["J"]
    if not v_pool or not j_pool:
        return "", {}
    v_name, v_seq = rng.choice(v_pool)
    j_name, j_seq = rng.choice(j_pool)

    # Trim some from V 3' and J 5' (exonuclease activity)
    v_trim = rng.randint(0, 6)
    j_trim = rng.randint(0, 6)
    v_used = v_seq[: len(v_seq) - v_trim] if v_trim else v_seq
    j_used = j_seq[j_trim:] if j_trim else j_seq

    # IGH/TRB/TRD/TRG use D segments
    d_used = ""
    d_name = None
    if locus in ("IGH", "TRB", "TRD") and ref[locus]["D"]:
        d_name, d_seq = rng.choice(ref[locus]["D"])
        # nibble D ends + flanking NDN
        d_used = d_seq[rng.randint(0,3): len(d_seq) - rng.randint(0,3)]

    n1 = random_ndn(rng, 2, 12)
    n2 = random_ndn(rng, 2, 12) if d_used else ""

    # Add C region if available (gives reads a 3' anchor like real capture)
    c_used = ""
    if ref[locus]["C"]:
        _, c_seq = rng.choice(ref[locus]["C"])
        c_used = c_seq[:80]      # first 80 bp of constant

    seq = v_used + n1 + d_used + n2 + j_used + c_used
    meta = dict(locus=locus, v=v_name, d=d_name, j=j_name,
                n1=n1, n2=n2, length=len(seq))
    return seq, meta


# ---------------------------------------------------------------------------
# Paired-end read sampler
# ---------------------------------------------------------------------------
def sample_paired_reads(template: str, n_reads: int, rng: random.Random,
                        read_len: int = 150,
                        frag_mean: int = 280, frag_sd: int = 40):
    """Yield (r1_seq, r2_seq) tuples sampled from a template."""
    L = len(template)
    if L < read_len + 20:
        # Pad short templates with random bases either side
        pad = "".join(rng.choices("ACGT", k=read_len))
        template = pad + template + pad
        L = len(template)

    for _ in range(n_reads):
        frag = max(read_len + 10, int(rng.gauss(frag_mean, frag_sd)))
        frag = min(frag, L)
        start = rng.randint(0, L - frag)
        fragment = template[start:start + frag]
        r1 = fragment[:read_len]
        r2 = revcomp(fragment[-read_len:])
        yield r1, r2


# ---------------------------------------------------------------------------
# FASTQ writing
# ---------------------------------------------------------------------------
def write_fastq_pair(out_prefix: str, reads, sample_id: str):
    r1_path = f"{out_prefix}_R1.fastq.gz"
    r2_path = f"{out_prefix}_R2.fastq.gz"
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for i, (r1, r2) in enumerate(reads, 1):
            qid = f"@{sample_id}_{i}"
            q1  = QUAL_STR * len(r1)
            q2  = QUAL_STR * len(r2)
            f1.write(f"{qid}/1\n{r1}\n+\n{q1}\n")
            f2.write(f"{qid}/2\n{r2}\n+\n{q2}\n")
    return r1_path, r2_path


# ---------------------------------------------------------------------------
# Repertoire generators
# ---------------------------------------------------------------------------
def generate_clonal(ref: dict, n_reads: int, rng: random.Random):
    """One dominant IGH clone (~50%) + one dominant TRB clone (~25%) + background."""
    plans: list[tuple[str, int]] = []   # list of (template_seq, n_reads)
    log: list[dict] = []

    # Dominant IGH
    igh_seq, igh_meta = build_clone(ref, "IGH", rng)
    if igh_seq:
        n_dom_igh = int(n_reads * 0.50)
        plans.append((igh_seq, n_dom_igh)); igh_meta["abundance"] = n_dom_igh
        log.append({"role": "dominant_IGH", **igh_meta})

    # Dominant TRB
    trb_seq, trb_meta = build_clone(ref, "TRB", rng)
    if trb_seq:
        n_dom_trb = int(n_reads * 0.25)
        plans.append((trb_seq, n_dom_trb)); trb_meta["abundance"] = n_dom_trb
        log.append({"role": "dominant_TRB", **trb_meta})

    # Background: 50 random clones across all loci sharing remainder
    remainder = n_reads - sum(c for _, c in plans)
    bg_n = 50
    per_bg = max(1, remainder // bg_n)
    for _ in range(bg_n):
        locus = rng.choice(LOCI)
        s, m = build_clone(ref, locus, rng)
        if not s: continue
        plans.append((s, per_bg))
        m["abundance"] = per_bg
        log.append({"role": "background", **m})

    return plans, log


def generate_polyclonal(ref: dict, n_reads: int, rng: random.Random):
    plans = []; log = []
    n_clones = max(200, n_reads // 25)
    per = max(1, n_reads // n_clones)
    for _ in range(n_clones):
        locus = rng.choice(LOCI)
        s, m = build_clone(ref, locus, rng)
        if not s: continue
        plans.append((s, per))
        m["abundance"] = per
        log.append({"role": "polyclonal", **m})
    return plans, log


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--imgt-ref",   required=True, type=Path,
                    help="TRUST4 bundled human_IMGT+C.fa (or mouse)")
    ap.add_argument("--bcrtcr-ref", required=False, type=Path,
                    help="(unused, kept for nf interface)")
    ap.add_argument("--mode",       required=True, choices=["clonal", "polyclonal"])
    ap.add_argument("--n-reads",    type=int, default=20000)
    ap.add_argument("--seed",       type=int, default=42)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args(argv)

    rng = random.Random(args.seed)
    ref = parse_imgt_fasta(args.imgt_ref)
    seg_counts = {L: {S: len(ref[L][S]) for S in SEG_TYPES} for L in LOCI}
    print(f"[simulate] parsed reference: {seg_counts}", file=sys.stderr)

    if args.mode == "clonal":
        plans, log = generate_clonal(ref, args.n_reads, rng)
    else:
        plans, log = generate_polyclonal(ref, args.n_reads, rng)

    def all_reads():
        for tmpl, n in plans:
            for r1, r2 in sample_paired_reads(tmpl, n, rng):
                yield r1, r2

    r1, r2 = write_fastq_pair(args.out_prefix, all_reads(), args.out_prefix)
    print(f"[simulate] wrote {r1}, {r2}")

    # Dump truth set
    import json
    Path(f"{args.out_prefix}.truth.json").write_text(json.dumps(log, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
