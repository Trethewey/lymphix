<p align="center"><img src="assets/logo.svg" alt="Lymphix" width="720"/></p>

# Lymphix

BCR / TCR clonality and V(D)J rearrangement calling from 2×150 bp Illumina
DNA capture NGS. Plain, TWIST UMI, and IDT xGen UMI-UDI libraries.
Local Docker, HPC Singularity, or DNAnexus.

---

## Quickstart

```bash
./run.sh test                                              # smoke test, no Docker
./run.sh --samplesheet samples.csv --outdir results/       # local run
./run.sh dnanexus --samplesheet dx://project:/samples.csv  # DNAnexus
```

## Install

- Nextflow ≥ 23.10 — `curl -fsSL https://get.nextflow.io | bash`
- Docker Desktop *or* Singularity
- Python 3.10+ with pandas, numpy, plotly, jinja2, scipy, pytest (smoke test only)

DNAnexus: [`docs/DNANEXUS.md`](docs/DNANEXUS.md).

## Sample sheet

A CSV with one row per sample. Pick the template matching your data, copy
it, replace the file names with your own.

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

<p align="center">
  <img src="examples/report_top_text.png" alt="Report header" width="100%"/>
</p>

<p align="center">
  <img src="examples/report_composition.png" alt="Lineage composition" width="100%"/>
</p>

Full HTML reports and `metrics.json` in [`examples/`](examples/).
Regenerate with `make test-smoke`.

## Build the containers

```bash
make build                                       # build all 5 images locally
make push                                        # push to ghcr.io/trethewey/lymphix
make push REGISTRY=ghcr.io/myorg/lymphix         # override registry if needed
```

## Tests

```bash
make test           # unit + smoke
make test-unit      # pytest math tests
make test-smoke     # mock-AIRR → clonality_metrics → report
```

## Panel BED

`regions/clonality_BCR_TCR.bed` (hg38, `chr` prefix) and
`regions/clonality_BCR_TCR.no_chr.bed` (Ensembl/GIAB-style). On-target QC only.

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

MIT. Third-party: TRUST4 (MIT) · IgBLAST (public domain) · fastp (MIT) ·
fgbio (MIT) · BWA (GPL-3) · samtools (MIT) · Nextflow (Apache-2.0) · Plotly (MIT).
