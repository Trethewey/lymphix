<p align="center"><img src="assets/logo.svg" alt="Lymphix" width="420"/></p>

# Lymphix

BCR / TCR clonality and V(D)J rearrangement calling from 2×150 bp Illumina
DNA capture NGS.

Supports plain Illumina libraries, **TWIST UMI** chemistry, and **IDT xGen
UMI-UDI** duplex UMI chemistry. Runs locally (Docker), on HPC (Singularity),
and on **DNAnexus** via the native Nextflow integration.

---

## Quickstart

```bash
# 1. Smoke test — no Docker / Nextflow needed
./run.sh test

# 2. Real data — paired FASTQ
./run.sh --samplesheet samples.csv --outdir results/

# 3. DNAnexus
./run.sh dnanexus --samplesheet dx://project:/samples.csv --outdir dx://project:/results/
```

## Install

| Requirement | For what | Install |
|---|---|---|
| **Nextflow ≥ 23.10** | Pipeline orchestration | `curl -fsSL https://get.nextflow.io \| bash` |
| **Docker Desktop** *or* **Singularity** | Containers | https://docker.com / Singularity per HPC docs |
| Python 3.10+ + pandas + plotly + jinja2 | Smoke test / analysis only | `pip install pandas numpy plotly jinja2 scipy pytest` |

For DNAnexus deployment see [`docs/DNANEXUS.md`](docs/DNANEXUS.md).

## Sample sheet

`samples.csv`:

```csv
sample_id,fastq_1,fastq_2,bam,umi_preset,expected_status
PT001,/data/PT001_R1.fq.gz,/data/PT001_R2.fq.gz,,none,
PT002,/data/PT002_R1.fq.gz,/data/PT002_R2.fq.gz,,twist,
PT003,,,/data/PT003.bam,none,
NEG,,,/data/blank.bam,none,negative
```

| Column | Required | Values |
|---|---|---|
| `sample_id`         | yes | Unique per row |
| `fastq_1`,`fastq_2` | one-of | Paired Illumina reads (gz ok) **or** use `bam` |
| `bam`               | one-of | Aligned BAM (sorted + indexed), hg38/GRCh38 |
| `umi_preset`        | no | `none` (default), `twist`, `xgen_duplex`, `xgen_simplex`, `custom` |
| `expected_status`   | no | `clonal` / `polyclonal` / `negative` — drives QC assertions |

**Constraint:** `bam + umi_preset != none` is rejected — UMI processing
requires inline UMI bases in the FASTQs.

See [`assets/samplesheet_template.csv`](assets/samplesheet_template.csv) for
annotated examples.

## What it does

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
6. **Composition call** — eight read pools summing to 100% of total input:
   - Clonal B-cell (IGH-, κ-, λ-restricted)
   - Polyclonal B-cell
   - Clonal T-cell (αβ TRB, γδ TRG/TRD)
   - Polyclonal T-cell
   - Background / germline
7. **κ:λ ratio** — light-chain restriction flag (clinical range 0.5–2.5).
8. **IGHV mutation status** — CLL prognostic call (≥ 98% V-identity = unmutated).
9. **HTML report** per sample — Plotly + Jinja2, self-contained.
10. **Cohort QC** — per-sample pass/fail vs `expected_status` written to
    `pipeline_info/qc_assertions.tsv`.

## UMI presets

| Preset | Read structure (fgbio) | Consensus | Typical use |
|---|---|---|---|
| `none`         | — | none | Plain Illumina, or pre-deduplicated BAMs |
| `twist`        | `5M2S+T 5M2S+T` | molecular | TWIST UMI Adapter System |
| `xgen_duplex`  | `9M+T 9M+T`     | duplex | IDT xGen UMI-UDI (MRD-grade) |
| `xgen_simplex` | `9M+T 9M+T`     | molecular | xGen UMI without duplex calling |
| `custom`       | `--umi_read_structure ...` | molecular | Other chemistries |

UMI samples require `--bwa_index <dir>` (a BWA-indexed reference).

## Outputs

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

Two cropped views from the report Lymphix produces on a synthetic clonal sample.

**Header — banner, clinical verdict, aggregate metrics:**

<p align="center">
  <img src="examples/report_top_text.png" alt="Report header — verdict and KPIs" width="100%"/>
</p>

**Section 1 — Lineage composition (thermometer + donut + Sankey):**

<p align="center">
  <img src="examples/report_composition.png" alt="Report Section 1 — lineage composition" width="100%"/>
</p>

**Both example reports** (download the raw `.html` and open locally for the
fully-interactive Plotly version):

- **Clonal sample** — dominant IGH + TRB clones · verdict *"Bi-clonal — B-cell (IGH) and T-cell (TRB)"*
  · [raw .html](examples/clonal_sample_report.html)
  · [metrics.json](examples/clonal_sample.metrics.json)
- **Polyclonal sample** — diverse repertoire · verdict *"Polyclonal repertoire — no dominant clone detected"*
  · [raw .html](examples/polyclonal_sample_report.html)
  · [metrics.json](examples/polyclonal_sample.metrics.json)

Headline numbers from the two examples:

| Metric | Clonal sample | Polyclonal sample |
|---|---:|---:|
| Aggregate clonotypes | 50 | 2,297 |
| Aggregate reads | 21,652 | 47,776 |
| Top clone fraction | 49.8% | 0.2% |
| **Clonality index** | **0.501** | **0.006** |
| IGH clonality | 0.884 (D50 = 1) | 0.001 (D50 = 138) |
| TRB clonality | 0.728 (D50 = 1) | 0.002 (D50 = 147) |
| κ:λ ratio | 0.72 (balanced) | 1.28 (balanced) |
| **Verdict banner** | "Bi-clonal — B-cell (IGH) and T-cell (TRB)" | "Polyclonal repertoire — no dominant clone detected" |

The raw [`examples/*.metrics.json`](examples/) files are also committed so the
full numeric output is available without running the pipeline.

To regenerate these locally:

```bash
make test-smoke
# results land in results_smoke_test/<sample>/<sample>.report.html
```

## Build the containers

```bash
make build                                       # build all 5 images locally
make push                                        # push to ghcr.io/trethewey/lymphix
make push REGISTRY=ghcr.io/myorg/lymphix         # override registry if needed
```

## Tests

```bash
make test           # unit tests + end-to-end smoke
make test-unit      # pytest math tests
make test-smoke     # mock-AIRR → clonality_metrics → report
```

Smoke test runs in <10 s, requires only Python + pandas + plotly + jinja2.

## Panel BED

`regions/clonality_BCR_TCR.bed` (hg38, `chr` prefix) and
`regions/clonality_BCR_TCR.no_chr.bed` (same coords, naked contigs for
Ensembl/GIAB-style references). Used for on-target QC; not required by
TRUST4 itself.

## Citations

Lymphix is a thin orchestration layer around two established immune-repertoire
tools. **If you publish results derived from this pipeline, please cite both
primary tools as well as Lymphix.**

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

### Supporting tools

| Tool | Used for | Citation |
|---|---|---|
| **fastp**    | Adapter/quality trimming | Chen et al. *Bioinformatics* 34, i884–i890 (2018). [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560) |
| **fgbio**    | UMI grouping + consensus calling | [github.com/fulcrumgenomics/fgbio](https://github.com/fulcrumgenomics/fgbio) |
| **BWA**      | Alignment for UMI-grouped reads | Li & Durbin. *Bioinformatics* 25, 1754–1760 (2009). [doi:10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324) |
| **samtools / htslib** | BAM/SAM handling | Danecek et al. *GigaScience* 10, giab008 (2021). [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| **Nextflow** | Workflow orchestration | Di Tommaso et al. *Nature Biotechnology* 35, 316–319 (2017). [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820) |
| **Plotly**   | Interactive HTML figures | [plotly.com/python](https://plotly.com/python/) |

A machine-readable [`CITATION.cff`](CITATION.cff) is also provided so GitHub
renders a "Cite this repository" button on the repo page.

## Licence

MIT — see [`LICENSE`](LICENSE).

Third-party licences (all permissive / open):
TRUST4 (MIT) · IgBLAST (NCBI public domain) · fastp (MIT) · fgbio (MIT) ·
BWA (GPL-3) · samtools (MIT) · Nextflow (Apache-2.0) · Plotly (MIT).
