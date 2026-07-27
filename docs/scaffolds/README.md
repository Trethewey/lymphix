# Scaffolds — unfinished, not part of the pipeline

Nothing in this directory is wired into `main.nf`, the `lymphix` command or the
test suite. It is kept for the design work it records, not because it runs.
Do not cite output from anything here.

## `infer_long_cdr3.py`

Confidence scoring for CDR3s that no single read spans, so an assembly-inferred
junction can be graded rather than silently trusted. `clonality_metrics.py`
already flags those clonotypes with `assembly_inferred`; this script was the
sketch of what to do with the flag.

Two modes were planned. The lightweight one, which scores from TRUST4's own
read counts and identities, is written but has never been validated against
anything. The full mode — BWA-realigning the source reads to each assembled
contig to measure junction-spanning coverage, paired bridging and consensus
agreement — raises `NotImplementedError`.

Moved here from `bin/` in July 2026. `bin/clonality_metrics.py` still carries a
docstring cross-reference to the old `bin/infer_long_cdr3.py` path.
