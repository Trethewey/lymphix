<p align="center"><img src="assets/lymphix-stacked.svg" alt="Lymphix" width="340"/></p>

<p align="center">
BCR / TCR clonality and V(D)J rearrangement calling from 2×150 bp Illumina
DNA capture NGS. Plain, TWIST UMI, and IDT xGen UMI-UDI libraries.
Runs as a Python analysis layer over TRUST4 output, or as a containerised
Nextflow workflow.
</p>

<p align="center">
<a href="LICENSE"><img src="https://img.shields.io/badge/licence-AGPL--3.0--or--later-blue.svg" alt="Licence: AGPL-3.0-or-later"/></a>
</p>

---

## Two layers, two levels of support

**The analysis layer** — clonality metrics, lineage composition, κ:λ, IGHV
status, per-sample and cohort HTML — is the supported path. Every real cohort
processed with Lymphix, the validation cohort below included, was run this way:
TRUST4 invoked natively, then the `lymphix` command over its AIRR output.

**The Nextflow workflow** (`main.nf`) wraps that same analysis layer together
with fastp, fgbio UMI consensus, TRUST4 and IgBLAST in containers, for local
Docker, HPC Singularity and DNAnexus. It has not been run end to end on real
data, and it has known gaps — for instance the `CLONALITY` process never passes
`--total-input-reads`, so composition falls back to a V(D)J-read denominator and
the background pool is always zero. Read `main.nf` before you rely on it.

## Quickstart — analysis layer

```bash
pip install .                # installs the `lymphix` command; Python ≥ 3.9

lymphix test                 # smoke test on mock AIRR: no TRUST4, IgBLAST, Docker or Nextflow

# 1. Assemble V(D)J natively with TRUST4 (installed separately)
run-trust4 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa \
    -1 S001_R1.fastq.gz -2 S001_R2.fastq.gz -o S001 -t 8

# 2. Clonality metrics + clonotype tables
#    TOTAL_READS = reads in the input FASTQ/BAM; omit for a V(D)J-only denominator
lymphix metrics --sample-id S001 \
    --trust4-airr S001_airr.tsv --igblast-airr S001_airr.tsv \
    --total-input-reads "$TOTAL_READS" --composition-denominator vdj \
    --out-metrics    S001.metrics.json \
    --out-clonotypes S001.clonotypes.tsv \
    --out-top        S001.top_clones.tsv

# 3. Per-sample HTML report
lymphix report --sample-id S001 \
    --metrics S001.metrics.json --clonotypes S001.clonotypes.tsv \
    --out S001.report.html

# 4. Cohort overview across a results directory
lymphix cohort --results-dir results/2026-07-26_mycohort --cohort-id mycohort \
    --out results/2026-07-26_mycohort/_cohort_summary.html
```

`--igblast-airr` is required. If you have not run IgBLAST, pass the TRUST4 AIRR
table to both flags: V-identity, productivity and in-frame flags then come from
TRUST4's own annotation. That is what the validation cohort does.

`lymphix --help` lists the rest — `compare`, `grade`, `simulate`, `merge-airr`,
`airr-to-fasta`.

## Quickstart — Nextflow workflow (not yet validated)

```bash
lymphix --samplesheet samples.csv --outdir results/        # nextflow run -profile docker
lymphix dnanexus --samplesheet dx://project:/samples.csv   # DNAnexus
```

Without installing, `./lymphix.sh` takes the same arguments.

## Install

- Lymphix — `pip install .` from a clone, or `pipx install .` to keep it
  isolated. Provides the `lymphix` command. Python ≥ 3.9. The Python
  dependencies (pandas, numpy, plotly, jinja2) come with it; add
  `pip install .[test]` for pytest.
- TRUST4 — needed for real data; Lymphix does not install it.
  [github.com/liulab-dfci/TRUST4](https://github.com/liulab-dfci/TRUST4)
- IgBLAST — optional, as above.
- Nextflow ≥ 23.10 (`curl -fsSL https://get.nextflow.io | bash`) plus Docker
  Desktop *or* Singularity — only for the containerised workflow.

DNAnexus: [`docs/DNANEXUS.md`](docs/DNANEXUS.md). The applet runs the Nextflow
workflow, so it carries the same caveat.

## Sample sheet

Only the Nextflow workflow reads a sample sheet; the analysis-layer commands
take file paths directly. A CSV with one row per sample. Pick the template
matching your data, copy it, replace the file names with your own.

| Your data | Template |
|---|---|
| Paired FASTQ, plain (no UMI) | [FASTQ template](assets/samplesheet_fastq.csv) |
| Paired FASTQ + UMIs (TWIST / xGen) | [UMI template](assets/samplesheet_umi.csv) |
| Pre-aligned BAMs | [BAM template](assets/samplesheet_bam.csv) |
| Validation run with controls | [Controls template](assets/samplesheet_controls.csv) |

```csv
sample_id,fastq_1,fastq_2
PT001,PT001_R1.fastq.gz,PT001_R2.fastq.gz
PT002,PT002_R1.fastq.gz,PT002_R2.fastq.gz
```

#### Columns

| Column | When | Values |
|---|---|---|
| `sample_id` | always | short label, no spaces |
| `fastq_1`,`fastq_2` | FASTQ input | paths (`.gz` OK) |
| `bam` | BAM input | path (sorted + indexed) |
| `umi_preset` | UMI library | `twist`, `xgen_duplex`, `xgen_simplex`, `custom` |
| `expected_status` | controls only | `clonal`, `polyclonal`, `negative` |

Each row uses **FASTQ or BAM, not both**. UMI presets require FASTQ.

## What it does

Steps 1–4 belong to the Nextflow workflow; on the supported route you run
TRUST4 yourself, and IgBLAST only if you want its annotation. Steps 5–9 are the
analysis layer, and are what `lymphix metrics` and `lymphix report` do. Step 10
is Nextflow-only.

1. **Trim & QC** — fastp.
2. **UMI consensus** (when `umi_preset != none`) — fgbio: `FastqToBam` →
   `bwa mem` → `GroupReadsByUmi` → consensus calling (duplex for xGen) →
   `FilterConsensusReads`.
3. **V(D)J assembly** — TRUST4. Recovers clonotypes for IGH, IGK, IGL,
   TRA, TRB, TRG, TRD.
4. **Re-annotation** — IgBLAST. Adds V-identity %, productivity, in-frame
   flags.
5. **Clonality metrics** — per locus and aggregate: Shannon H, normalised
   H, Simpson D, Gini, D50, top-clone fraction, **clonality index**
   (1 − H/log N).
6. **Composition call** — eight mutually-exclusive read pools:
   - Clonal B-cell (IGH-, κ-, λ-restricted) — three pools
   - Polyclonal B-cell
   - Clonal T-cell (αβ TRB, γδ TRG/TRD) — two pools
   - Polyclonal T-cell (TRA reads always land here; TRA alone is not
     diagnostic of clonality)
   - Background / germline

   They sum to 100% of total input reads only when `--total-input-reads` is
   supplied and `--composition-denominator total` is in force. Otherwise the
   denominator is V(D)J-assigned reads and background is zero.
7. **κ:λ ratio** — light-chain restriction flag (clinical range 0.5–2.5).
8. **IGHV mutation status** — CLL prognostic call on the *dominant* IGH clone
   (≥ 98% V-identity = unmutated), plus a descriptive repertoire-wide tally.
   Reported as `unknown` when the dominant clone carries no V-identity.
9. **HTML report** per sample — Plotly + Jinja2, self-contained.
10. **Cohort QC** (Nextflow only) — per-sample pass/fail vs `expected_status`
    written to `pipeline_info/qc_assertions.tsv`. Only `expected_status=negative`
    is actually asserted on.

## UMI presets

UMI consensus calling lives in the Nextflow workflow, so it inherits that
caveat.

| Preset | Read structure (fgbio) | Consensus | Typical use |
|---|---|---|---|
| `none`         | — | none | Plain Illumina, or pre-deduplicated BAMs |
| `twist`        | `5M2S+T 5M2S+T` | molecular | TWIST UMI Adapter System |
| `xgen_duplex`  | `9M+T 9M+T`     | duplex | IDT xGen UMI-UDI (MRD-grade) |
| `xgen_simplex` | `9M+T 9M+T`     | molecular | xGen UMI without duplex calling |
| `custom`       | `--umi_read_structure ...` | molecular | Other chemistries |

UMI samples require `--bwa_index <dir>` (a BWA-indexed reference).

## WGS samples

Lymphix was tuned on CAPP-seq / capture data (500–2000× per-position depth at
IG loci, 5,000–20,000 V(D)J reads per sample). WGS at 30–40× gives ~30 reads
per IG-locus position and yields only **~150–300 V(D)J reads per sample**,
around or below the 200-read low-yield threshold at which the report raises a
reduced-confidence warning.

For WGS inputs, pass `--wgs` to relax the germline-filter V-match floor
(98 → 60 nt at the default 150 bp read length; the floor otherwise scales as
0.65 × read length) and switch the composition denominator to V(D)J reads only.
Decide your supporting-read sensitivity explicitly with `-c/--clones`
(default 2; use 1 to recover sub-threshold clonotypes at the cost of a higher
noise floor; 3+ to be stricter):

```bash
lymphix metrics --wgs -c 1 \
    --sample-id S001 --trust4-airr S001_airr.tsv --igblast-airr S001_airr.tsv \
    --out-metrics S001.metrics.json \
    --out-clonotypes S001.clonotypes.tsv --out-top S001.top_clones.tsv
```

WGS samples are still depth-limited. The dominant IGH clone (if present) and
IGHV mutation status are robust; per-locus clonotype counts in the long tail
are noise. Translocation-disrupted IGH (common in DLBCL) cannot be assembled
by TRUST4 regardless of threshold — that's an out-of-scope limitation.

## Clonotype collapsing (opt-in, default OFF)

TRUST4 writes one AIRR row per **(assembly, CDR3 variant)**, not one row per
clone. A single rearrangement therefore arrives as a dominant row plus a tail
of near-identical rows carrying a handful of reads each — on the validation
cohorts, 132 of 195 minor rows were exactly one substitution from the dominant
row of their group, and 159 of them carried five reads or fewer. Counted as
clonotypes, that tail inflates N in Shannon, Simpson, Gini, D50 and the
clonality index simultaneously, and can hand "dominant clone" to a fragment.

Collapsing is available and **off by default**. It changes `n_clonotypes` and
therefore every diversity metric for every sample, so it must be an explicit
choice, and every `metrics.json` records which convention produced its numbers
(`collapse_clonotypes`, `collapse_key`, and a `clonotype_collapse` block with
rows in / clones out per locus). Measure the effect on your own data before
adopting it.

```bash
lymphix metrics --collapse-clonotypes \
    --collapse-key locus_junction_nt_hamming1 --collapse-minor-fraction 0.02 \
    --sample-id S001 --trust4-airr S001_airr.tsv --igblast-airr S001_airr.tsv \
    --out-metrics S001.metrics.json \
    --out-clonotypes S001.clonotypes.tsv --out-top S001.top_clones.tsv
```

Two keys:

| `--collapse-key` | Merges | Judgement involved |
|---|---|---|
| `locus_junction_nt` (default) | Rows with an identical junction nt at the same locus | None — exact |
| `locus_junction_nt_hamming1`  | The above, plus single-substitution variants below `--collapse-minor-fraction` of their parent | Yes — where error ends and subclone begins |

The V call is deliberately **not** part of either key. The only exact
duplicates in real data are one rearrangement assembled twice against
paralogous V references (IGKV3-15 / IGKV3D-15, IGKV2-28 / IGKV2D-29, the
IGLV5-37/45/48/52 family); keying on V keeps them apart, which is the split
that most needs closing. Nor is the junction amino acid sequence a key: rows
with an untranslatable junction all carry a blank `junction_aa` and would pool
unrelated rearrangements under "unknown".

Two things to know before you switch it on:

* **Read counts are summed within an assembly, and maxed across assemblies.**
  TRUST4's per-variant abundances partition the reads spanning a CDR3, so
  summing them within one contig reconstructs its support. No TRUST4 output
  maps read identifiers to contigs, so whether two contigs share reads is
  unverifiable — taking the larger never invents depth, which is the failure
  that matters in a report.
* **The direction of travel is not uniform.** `clonality_index = 1 − H/log N`
  has N in the denominator, so removing the noise tail raises clonality on a
  single-clone sample but lowers it on a sample with two or more real clones.
  Collapsing also pushes low-yield samples below `INDETERMINATE_MAX_CLONOTYPES`
  (5), so some genuinely polyclonal samples will start reading as
  indeterminate.

Under Nextflow: `--collapse_clonotypes true --collapse_key ... --collapse_minor_fraction ...`.

## Outputs

The Nextflow workflow lays results out like this. On the analysis-layer route
you choose the paths yourself; name run directories `YYYY-MM-DD_description`.

```
results/
├── <sample>/
│   ├── fastp/                 trimmed reads + JSON + HTML
│   ├── umi_consensus/         (UMI runs only) consensus FASTQ + family-size histogram
│   ├── trust4/                TRUST4 raw outputs incl. AIRR
│   ├── igblast/               AIRR table with V-identity %, productivity
│   ├── clonality/
│   │   ├── <sample>.metrics.json
│   │   ├── <sample>.clonotypes.tsv
│   │   └── <sample>.top_clones.tsv
│   └── <sample>.report.html
└── pipeline_info/
    ├── execution_report.html
    ├── software_versions.yml
    └── qc_assertions.tsv
```

## Example results

<p align="center">
  <img src="examples/report_top_text.png" alt="Report header" width="100%"/>
</p>

<p align="center">
  <img src="examples/report_composition.png" alt="Lineage composition" width="100%"/>
</p>

Example `metrics.json` for a clonal and a polyclonal sample are in
[`examples/`](examples/). `bash tests/test_smoke.sh` builds equivalent reports
from mock AIRR into `results_smoke_test/`; it does not write to `examples/`.

## Build the containers

Only the Nextflow workflow uses these; the analysis layer runs from the pip
install.

```bash
make build                                       # build all 5 images locally
make push                                        # push to ghcr.io/trethewey/lymphix
make push REGISTRY=ghcr.io/myorg/lymphix         # override registry if needed
```

`make build` tags every image `$(TAG)`, default `0.1.0`, but the modules ask
for tool-versioned tags — `fastp:0.23.4`, `trust4:1.0.13`, `igblast:1.22.0`,
`fgbio:0.1.0`, `clonality:0.1.0`. Build the first three with a matching
`TAG=`, or the workflow will look for images that do not exist — another reason
to treat the Nextflow route as unvalidated.

## Tests

```bash
python -m pytest tests/ -q   # full suite: metrics maths + regression tests
bash tests/test_smoke.sh     # mock AIRR → clonality_metrics → report
make test                    # both, via make test-unit + make test-smoke
```

Neither exercises TRUST4, IgBLAST, Docker or Nextflow: the smoke test starts
from mock AIRR fixtures generated by `tests/make_mock_airr.py` and asserts the
clonal sample scores a clonality index > 0.3 and the polyclonal one < 0.1.
`make test-unit` runs `tests/test_clonality_metrics.py` alone — use `pytest
tests/` to include `tests/test_lymphix_common.py` (shared constants and the
single verdict rule) and `tests/test_regressions_2026_07.py` (one test per
defect fixed in July 2026) as well.

## Validation cohort

Eleven samples with known biology — nine hematological cell lines plus a
no-signal and a polyclonal control. Expected outcomes are encoded in
[`tests/validation_expected.json`](tests/validation_expected.json) with
Cellosaurus IDs and original references.

| Sample | Cellosaurus | Disease | Expected verdict |
|---|---|---|---|
| JURKAT       | CVCL_0065 | T-ALL                       | clonal (TRA + TRB) |
| MOLT-4       | CVCL_0013 | T-ALL                       | clonal (TRB) |
| KARPAS-299   | CVCL_1324 | ALK+ ALCL                   | clonal (TRB / TRA) |
| NAMALWA      | CVCL_0067 | Burkitt lymphoma            | clonal (IGH) |
| DAUDI        | CVCL_0008 | Burkitt lymphoma            | clonal (IGH, mutated IGHV) |
| RAJI         | CVCL_0511 | Burkitt lymphoma            | clonal (IGH, mutated IGHV) |
| OCI-LY1      | CVCL_1879 | Germinal-centre DLBCL       | clonal (IGH, typically mutated IGHV) |
| U-266        | CVCL_0566 | Multiple myeloma            | clonal (IGH) |
| MM.1S        | CVCL_8792 | Multiple myeloma            | clonal (IGH) |
| PBMC_HEALTHY | —         | Healthy donor, 3' scRNA-seq | no_signal (chemistry-correct) |
| POLYCLONAL_SIM | —       | Synthetic polyclonal        | no_clonal |

Run it end to end — this is the analysis-layer route, not Nextflow: the script
downloads FASTQ from ENA, runs TRUST4 natively, then `clonality_metrics.py`,
`generate_report.py`, `grade_validation.py` and `cohort_report.py`. It needs
`run-trust4`, `wget` and about 25 GB of free disk; runtime depends on `THREADS`
(default 8). Re-runs are idempotent — completed samples are skipped.

```bash
tests/run_validation_cohort.sh [DATA_DIR]
```

IgBLAST is not run: the TRUST4 AIRR table is passed as both inputs, so
V-identity comes from TRUST4.

The cohort overview it produces is checked in: verdict table,
lineage-composition stacked bar, and per-locus clonality-index heatmap on one
self-contained page. Two builds:

| File | Size | Use |
|---|---|---|
| [`examples/cohort_overview.html`](examples/cohort_overview.html)         | 4.9 MB | inline Plotly, offline-safe (firewalled clinical networks) |
| [`examples/cohort_overview_cdn.html`](examples/cohort_overview_cdn.html) |  39 KB | loads Plotly from CDN, needs internet |

**View rendered in your browser** (uses the CDN build):
[htmlpreview.github.io / cohort_overview_cdn.html](https://htmlpreview.github.io/?https://github.com/Trethewey/lymphix/blob/main/examples/cohort_overview_cdn.html)

For a permanent URL on this repo, enable
[GitHub Pages](https://github.com/Trethewey/lymphix/settings/pages)
(source: `main`, folder `/`) — the landing page then lives at
[trethewey.github.io/lymphix/](https://trethewey.github.io/lymphix/).

## Panel BED

`assets/regions/clonality_BCR_TCR.bed` (hg38, `chr` prefix) and
`assets/regions/clonality_BCR_TCR.no_chr.bed` (Ensembl/GIAB-style), with
`assets/regions/clonality_BCR_TCR_regions.tsv` giving the same intervals
annotated with locus, BIOMED-2 tube, tiling and role.

These are reference material for on-target QC and panel design — no Lymphix
step reads them. `params.panel_bed` in `nextflow.config` points at the
`chr`-prefixed BED but is not consumed by any process.

Note the role column: IGH V-region coverage is flagged as needed for IGHV
mutation status. A panel that anchors only on J and constant regions can still
call clonality but cannot support a CLL prognostic call.

## Citations

Please cite TRUST4 and IgBLAST when publishing results.

### Primary tools

**TRUST4** — V(D)J assembly from bulk and single-cell sequencing.
> Song L, Cohen D, Ouyang Z, Cao Y, Hu X, Liu XS.
> *TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data.*
> Nature Methods 18, 627–630 (2021).
> [doi:10.1038/s41592-021-01142-2](https://doi.org/10.1038/s41592-021-01142-2)
> · [github.com/liulab-dfci/TRUST4](https://github.com/liulab-dfci/TRUST4)

**IgBLAST** — V/D/J alignment, productivity, somatic-hypermutation (V-identity %).
> Ye J, Ma N, Madden TL, Ostell JM.
> *IgBLAST: an immunoglobulin variable domain sequence analysis tool.*
> Nucleic Acids Research 41(W1), W34–W40 (2013).
> [doi:10.1093/nar/gkt382](https://doi.org/10.1093/nar/gkt382)
> · [ncbi.nlm.nih.gov/igblast](https://www.ncbi.nlm.nih.gov/igblast/)

### Clinical interpretation

**IGHV mutation status** — the 98% V-identity cutoff for mutated vs unmutated
CLL is taken from the original prognostic studies and the iwCLL guidelines.
> Damle RN *et al.* *Ig V gene mutation status and CD38 expression as novel prognostic indicators in chronic lymphocytic leukemia.* Blood 94, 1840–1847 (1999). [doi:10.1182/blood.V94.6.1840](https://doi.org/10.1182/blood.V94.6.1840)
>
> Hamblin TJ *et al.* *Unmutated Ig V(H) genes are associated with a more aggressive form of chronic lymphocytic leukemia.* Blood 94, 1848–1854 (1999). [doi:10.1182/blood.V94.6.1848](https://doi.org/10.1182/blood.V94.6.1848)
>
> Hallek M *et al.* *iwCLL guidelines for diagnosis, indications for treatment, response assessment, and supportive management of CLL.* Blood 131, 2745–2760 (2018). [doi:10.1182/blood-2017-09-806398](https://doi.org/10.1182/blood-2017-09-806398)

**AIRR-C clonotype schema** — clonotype tables follow the Adaptive Immune
Receptor Repertoire Community (AIRR-C) v1.4 Rearrangement schema.
> Vander Heiden JA *et al.* *AIRR Community Standardized Representations for Annotated Immune Repertoires.* Frontiers in Immunology 9, 2206 (2018). [doi:10.3389/fimmu.2018.02206](https://doi.org/10.3389/fimmu.2018.02206)

### Supporting tools

| Tool | Used for | Citation |
|---|---|---|
| **fastp**    | Adapter/quality trimming | Chen et al. *Bioinformatics* 34, i884–i890 (2018). [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560) |
| **fgbio**    | UMI grouping + consensus calling | [github.com/fulcrumgenomics/fgbio](https://github.com/fulcrumgenomics/fgbio) |
| **BWA**      | Alignment for UMI-grouped reads | Li & Durbin. *Bioinformatics* 25, 1754–1760 (2009). [doi:10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324) |
| **samtools / htslib** | BAM/SAM handling | Danecek et al. *GigaScience* 10, giab008 (2021). [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| **Nextflow** | Workflow orchestration | Di Tommaso et al. *Nature Biotechnology* 35, 316–319 (2017). [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820) |
| **Plotly**   | Interactive HTML figures | [plotly.com/python](https://plotly.com/python/) |

See [`CITATION.cff`](CITATION.cff) for the machine-readable form.

## Licence

Copyright © 2026 C.S. Trethewey.

Lymphix is free software under the **GNU Affero General Public License v3.0 or
later** (AGPL-3.0-or-later) — see [`LICENSE`](LICENSE). Commercial use is
permitted. If you modify Lymphix and either distribute it or run it as a
network service, you must release your modified source under the same licence.

Lymphix was previously distributed under the MIT licence. Copies already
obtained under MIT remain usable under those terms; the AGPL applies to this
and all subsequent versions.

Third-party tools keep their own licences: TRUST4 (MIT) · IgBLAST (public
domain) · fastp (MIT) · fgbio (MIT) · BWA (GPL-3) · samtools (MIT) ·
Nextflow (Apache-2.0) · Plotly (MIT).
