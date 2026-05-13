# Deploying to DNAnexus

This pipeline supports two deployment modes on DNAnexus.

## Option A — Native Nextflow Pipeline (recommended)

DNAnexus has a first-class Nextflow integration. Each Nextflow process is
scheduled as its own DNAnexus job, which gives proper per-task scaling and
billing visibility.

### Prerequisites

```bash
pip install dxpy
dx login
dx select <your-project>
```

You'll also need the `dx-toolkit` (>= v0.366.0) which includes
`dx build-nextflow-applet`.

### Build & deploy

```bash
# From the repo root
dx build-nextflow-applet . \
    --destination /applets/bcr_tcr_clonality \
    --extra-args '{"runSpec":{"timeoutPolicy":{"*":{"hours":24}}}}'
```

This compiles `main.nf` into a DNAnexus applet that runs Nextflow on a
head-node instance and dispatches each process as a sub-job.

### Run on the platform

Upload a samplesheet referencing FASTQs already on the platform:

```csv
sample_id,fastq_1,fastq_2
PT001,dx://project-xxxx:/fastq/PT001_R1.fastq.gz,dx://project-xxxx:/fastq/PT001_R2.fastq.gz
PT002,dx://project-xxxx:/fastq/PT002_R1.fastq.gz,dx://project-xxxx:/fastq/PT002_R2.fastq.gz
```

```bash
dx upload samplesheet.csv --path /samplesheets/
dx run /applets/bcr_tcr_clonality \
    -i nextflow_pipeline_params="
        --samplesheet dx://project-xxxx:/samplesheets/samplesheet.csv
        --outdir      dx://project-xxxx:/results/
        --species     human
        --igh_mutated_cutoff 98.0
    " \
    --priority normal \
    --watch
```

### Containers

The pipeline pulls five images (fgbio is only invoked when `umi_preset != none`):

```
ghcr.io/trethewey/lymphix/fastp:0.23.4
ghcr.io/trethewey/lymphix/fgbio:0.1.0      # UMI consensus (TWIST / xGen)
ghcr.io/trethewey/lymphix/trust4:1.0.13
ghcr.io/trethewey/lymphix/igblast:1.22.0
ghcr.io/trethewey/lymphix/clonality:0.1.0
```

For private registries, configure DNAnexus to authenticate (see
"Container registry credentials" in the DNAnexus Nextflow docs).
For air-gapped projects, upload the images as `.tar.gz` to your project
and override `process.container` in `conf/dnanexus.config`.

## Option B — Classic DNAnexus applet (single-job)

If your project policy disallows Nextflow on the platform, you can deploy
the whole pipeline as a single applet:

```bash
dx build . --destination /applets/bcr_tcr_clonality --overwrite
dx run /applets/bcr_tcr_clonality \
    -i samplesheet=/samplesheets/samplesheet.csv \
    --instance-type mem2_ssd1_v2_x16
```

This uses `src/code.sh` as the wrapper, which installs Nextflow inside the
worker and runs the pipeline end-to-end on a single (larger) instance.
Less efficient than Option A but simpler permissions-wise.

## Cost / instance sizing

| Step       | DX instance              | Why |
|------------|--------------------------|-----|
| fastp      | `mem1_ssd1_v2_x4`        | I/O bound |
| TRUST4     | `mem2_ssd1_v2_x16`       | CPU + ~16 GB |
| IgBLAST    | `mem2_ssd1_v2_x8`        | CPU |
| Clonality  | `mem1_ssd1_v2_x4`        | trivial |

Per sample, expect ~30–60 min wallclock and a few US$ on default instances.
