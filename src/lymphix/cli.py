"""The `lymphix` command.

Two ways of running the pipeline are wrapped here:

  * the Nextflow route (`lymphix --samplesheet ...`), which is the documented
    end-to-end path, and
  * the analysis scripts under bin/ (`lymphix metrics`, `report`, `cohort`),
    which is how the real cohorts have actually been processed — TRUST4 run
    natively, then the Python layer over its AIRR output.

Unrecognised arguments are passed straight through to the underlying tool, so
any Nextflow parameter works without this wrapper knowing about it.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__

GITHUB_PIPELINE = "Trethewey/lymphix"

# Subcommand -> script in bin/. These take their own arguments; we do not
# interpret them, we only put the interpreter and sys.path in the right place.
BIN_COMMANDS = {
    "metrics": "clonality_metrics.py",
    "report": "generate_report.py",
    "cohort": "cohort_summary.py",
    "compare": "cohort_compare.py",
    "grade": "grade_validation.py",
    "simulate": "simulate_repertoire.py",
    "merge-airr": "merge_airr.py",
    "airr-to-fasta": "airr_to_fasta.py",
}

USAGE = f"""\
Lymphix {__version__} — BCR/TCR clonality from custom-panel NGS.

Pipeline (Nextflow + containers):
  lymphix --samplesheet SAMPLES.CSV --outdir results/   Run the pipeline
  lymphix dnanexus --samplesheet dx://...               Submit to DNAnexus
  lymphix test                                          Analysis-layer smoke test

Analysis layer directly (TRUST4 already run):
  lymphix metrics --help          Clonality metrics from AIRR
  lymphix report --help           Per-sample HTML report
  lymphix cohort --help           Cohort summary across samples
  lymphix compare --help          Cohort comparison
  lymphix grade --help            Grade a validation cohort
  lymphix simulate --help         Synthetic repertoire generator
  lymphix merge-airr --help       Merge AIRR tables
  lymphix airr-to-fasta --help    AIRR to FASTA

Common pipeline options (passed through to Nextflow):
  --samplesheet PATH        CSV: sample_id,fastq_1,fastq_2,bam,umi_preset,expected_status
  --outdir PATH             Results directory (default: results/)
  --species human|mouse     Reference species (default: human)
  --umi_preset PRESET       none|twist|xgen_duplex|xgen_simplex|custom
  --total_input_reads N     For accurate background fraction in composition
  --filter_dups_in_bam      Strip flag-marked duplicates from BAM input

Set LYMPHIX_HOME to point at a pipeline checkout explicitly.
See README.md for the full sample-sheet schema.
"""


def pipeline_root() -> Path | None:
    """Locate the checkout holding main.nf.

    Checks LYMPHIX_HOME, then walks up from this file (which covers an
    editable install made from inside the clone), then from the working
    directory. Returns None when nothing looks like a checkout.
    """
    candidates: list[Path] = []

    env_home = os.environ.get("LYMPHIX_HOME")
    if env_home:
        candidates.append(Path(env_home))

    here = Path(__file__).resolve()
    candidates.extend(here.parents)

    cwd = Path.cwd().resolve()
    candidates.append(cwd)
    candidates.extend(cwd.parents)

    for path in candidates:
        if (path / "main.nf").is_file() and (path / "nextflow.config").is_file():
            return path
    return None


def _require_root(what: str) -> Path:
    root = pipeline_root()
    if root is None:
        sys.stderr.write(
            f"[lymphix] {what} needs a pipeline checkout and none was found.\n"
            f"[lymphix] Run from inside a clone, or set LYMPHIX_HOME to one.\n"
        )
        raise SystemExit(2)
    return root


def _run_nextflow(args: list[str]) -> int:
    if not shutil.which("nextflow"):
        sys.stderr.write(
            "[lymphix] Nextflow not found. Install it with:\n"
            "[lymphix]   curl -fsSL https://get.nextflow.io | bash\n"
        )
        return 1
    if not shutil.which("docker"):
        sys.stderr.write(
            "[lymphix] WARNING: Docker not found. Use -profile singularity, "
            "or install Docker Desktop.\n"
        )

    root = pipeline_root()
    if root is None:
        sys.stderr.write(
            f"[lymphix] No local checkout found; pulling {GITHUB_PIPELINE} from GitHub.\n"
        )
        target = GITHUB_PIPELINE
    else:
        target = str(root / "main.nf")

    return subprocess.call(["nextflow", "run", target, "-profile", "docker", *args])


def _run_dnanexus(args: list[str]) -> int:
    if not shutil.which("dx"):
        sys.stderr.write(
            "[lymphix] dx toolkit not found — see docs/DNANEXUS.md for prerequisites.\n"
        )
        return 1
    sys.stderr.write("[lymphix] Submitting to DNAnexus — see docs/DNANEXUS.md.\n")
    return subprocess.call(
        ["dx", "run", "/applets/lymphix", "-i", f"nextflow_pipeline_params={' '.join(args)}", "--watch"]
    )


def _run_smoke_test() -> int:
    root = _require_root("The smoke test")
    script = root / "tests" / "test_smoke.sh"
    if not script.is_file():
        sys.stderr.write(f"[lymphix] Smoke test not found at {script}\n")
        return 2
    if not shutil.which("bash"):
        sys.stderr.write("[lymphix] bash not found; the smoke test is a shell script.\n")
        return 1
    # Relative path, resolved by bash from cwd: a Windows-style absolute path
    # is not something Git Bash can open.
    return subprocess.call(["bash", "tests/test_smoke.sh"], cwd=str(root))


def _run_bin_script(script_name: str, args: list[str]) -> int:
    """Run a bin/ script with bin/ on sys.path.

    The scripts import each other by bare module name (generate_report.py does
    `import clonality_metrics`), which only resolves when bin/ is importable.
    Running them through this wrapper makes that work from any directory.
    """
    root = _require_root(f"`lymphix {script_name}`")
    script = root / "bin" / script_name
    if not script.is_file():
        sys.stderr.write(f"[lymphix] Script not found: {script}\n")
        return 2

    env = dict(os.environ)
    bin_dir = str(root / "bin")
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = bin_dir + (os.pathsep + existing if existing else "")

    return subprocess.call([sys.executable, str(script), *args], env=env)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)

    if not args or args[0] in {"-h", "--help", "help"}:
        sys.stdout.write(USAGE)
        return 0

    command = args[0]

    if command in {"-V", "--version", "version"}:
        sys.stdout.write(f"lymphix {__version__}\n")
        return 0
    if command == "test":
        return _run_smoke_test()
    if command == "dnanexus":
        return _run_dnanexus(args[1:])
    if command in BIN_COMMANDS:
        return _run_bin_script(BIN_COMMANDS[command], args[1:])

    # Anything else is a Nextflow invocation: `lymphix --samplesheet x.csv ...`
    if not command.startswith("-"):
        sys.stderr.write(
            f"[lymphix] Unknown subcommand '{command}'. Run `lymphix --help` for usage.\n"
        )
        return 2

    return _run_nextflow(args)


if __name__ == "__main__":
    raise SystemExit(main())
