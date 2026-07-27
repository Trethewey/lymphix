"""
Regression tests for the July 2026 correctness fixes.

Each test corresponds to a defect that shipped silently — the pipeline
produced a confident, well-formatted, wrong answer rather than failing. They
are grouped by the defect they pin down so a future change that reintroduces
one fails here rather than in a clinical report.

Run:    pytest tests/test_regressions_2026_07.py -v
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
sys.path.insert(0, str(BIN))

import clonality_metrics as cm  # noqa: E402


AIRR_COLUMNS = [
    "sequence_id", "v_call", "j_call", "junction", "junction_aa",
    "duplicate_count", "v_identity", "v_cigar",
]


def write_airr(path: Path, rows: list[dict]) -> None:
    """Write a minimal but valid AIRR TSV. No rows == header only."""
    frame = pd.DataFrame(rows, columns=AIRR_COLUMNS)
    frame.to_csv(path, sep="\t", index=False)


def igh_row(seq_id: str, reads: int, identity: float, junction: str) -> dict:
    return {
        "sequence_id": seq_id,
        "v_call": "IGHV3-23*01",
        "j_call": "IGHJ4*02",
        "junction": junction,
        "junction_aa": "CARDYW",
        "duplicate_count": reads,
        "v_identity": identity,
        "v_cigar": "150M",
    }


def run_metrics(tmp_path: Path, trust4_rows, igblast_rows, extra_args=()) -> dict:
    """Invoke clonality_metrics.py as the pipeline does, return metrics.json."""
    t4 = tmp_path / "trust4.airr.tsv"
    ig = tmp_path / "igblast.airr.tsv"
    write_airr(t4, trust4_rows)
    write_airr(ig, igblast_rows)
    out_metrics = tmp_path / "out.metrics.json"

    result = subprocess.run(
        [sys.executable, str(BIN / "clonality_metrics.py"),
         "--sample-id", "TEST",
         "--trust4-airr", str(t4),
         "--igblast-airr", str(ig),
         "--out-metrics", str(out_metrics),
         "--out-clonotypes", str(tmp_path / "out.clonotypes.tsv"),
         "--out-top", str(tmp_path / "out.top.tsv"),
         *extra_args],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"clonality_metrics failed:\n{result.stderr}"
    # Strict parse: this is also the metrics.json validity regression test.
    return json.loads(out_metrics.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# BLOCKER: IGHV status was a repertoire-wide read majority, not the clone's
# ---------------------------------------------------------------------------
class TestIghvIsDominantCloneNotRepertoire:
    """A mutated tumour clone on a naive-B tail must not be called unmutated.

    Naive B cells are unmutated by definition, so a polyclonal background can
    outvote the tumour clone on read count. The prognostic call inverts:
    mutated IGHV is favourable in CLL, unmutated is poor-prognosis.
    """

    @staticmethod
    def _mutated_clone_on_unmutated_tail():
        # 900-read tumour clone at 94% identity (mutated, favourable) ...
        trust4 = [igh_row("tumour", 900, 94.0, "TGTGCGAGAGATTAC" + "TGG")]
        # ... under 1200 reads of unmutated naive-B background across 400 clones
        for i in range(400):
            trust4.append(igh_row(f"naive{i}", 3, 100.0, f"TGTGCGAGA{i:04d}TGG"))
        return trust4

    def test_dominant_clone_status_follows_the_clone(self, tmp_path):
        rows = self._mutated_clone_on_unmutated_tail()
        metrics = run_metrics(tmp_path, rows, rows)
        ighv = metrics["ighv_status"]

        assert ighv["dominant_clone_status"] == "mutated", (
            "The dominant clone is at 94% V-identity, i.e. mutated. A "
            "repertoire read-majority would call this unmutated and invert "
            "the prognosis."
        )
        assert ighv["dominant_clone_reads"] == 900
        assert ighv["dominant_clone_v_identity"] == pytest.approx(94.0)

    def test_repertoire_tally_kept_but_clearly_separate(self, tmp_path):
        rows = self._mutated_clone_on_unmutated_tail()
        ighv = run_metrics(tmp_path, rows, rows)["ighv_status"]

        # The old majority is still available for description, under a name
        # that cannot be mistaken for the clinical call.
        assert "dominant_status" not in ighv, "the ambiguous key must be gone"
        assert ighv["repertoire_unmutated_read_fraction"] > 0.5
        assert ighv["dominant_clone_status"] == "mutated"

    def test_unassessed_dominant_clone_is_not_assessable(self, tmp_path):
        """An unassessed clone must not inherit the background's status."""
        # Dominant clone has no IgBLAST identity at all ...
        trust4 = [igh_row("dominant", 1500, None, "TGTGCGAGATTTTGG")]
        for i in range(20):
            trust4.append(igh_row(f"bg{i}", 10, 100.0, f"TGTGCGAGA{i:04d}TGG"))
        # ... because it is absent from the IgBLAST output entirely.
        igblast = [r for r in trust4 if r["sequence_id"] != "dominant"]

        ighv = run_metrics(tmp_path, trust4, igblast)["ighv_status"]
        assert ighv["dominant_clone_status"] == "unknown", (
            "200 background reads must not decide the status of a clone that "
            "was never assessed."
        )
        assert ighv["dominant_clone_v_identity"] is None

    def test_unknown_reads_excluded_from_the_fraction(self, tmp_path):
        """The fraction must be over assessed reads, not diluted by unknowns."""
        trust4 = [igh_row("dominant", 1500, None, "TGTGCGAGATTTTGG")]
        for i in range(20):
            trust4.append(igh_row(f"bg{i}", 10, 100.0, f"TGTGCGAGA{i:04d}TGG"))
        igblast = [r for r in trust4 if r["sequence_id"] != "dominant"]

        ighv = run_metrics(tmp_path, trust4, igblast)["ighv_status"]
        assert ighv["reads_unknown"] == 1500
        assert ighv["reads_assessed"] == 200
        # All 200 assessed reads are unmutated → 1.0, not 200/1700.
        assert ighv["repertoire_unmutated_read_fraction"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# HIGH: a missing input produced a signed "no signal" report and exit 0
# ---------------------------------------------------------------------------
class TestMissingInputIsAnError:
    def test_absent_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            cm.read_airr(tmp_path / "never_written.tsv")

    def test_zero_byte_file_raises(self, tmp_path):
        empty = tmp_path / "empty.tsv"
        empty.write_bytes(b"")
        with pytest.raises(ValueError):
            cm.read_airr(empty)

    def test_header_only_is_a_legitimate_empty_repertoire(self, tmp_path):
        """A genuinely empty result is header-only, and must still parse."""
        header_only = tmp_path / "header.tsv"
        write_airr(header_only, [])
        assert cm.read_airr(header_only).empty

    def test_pipeline_exits_nonzero_on_missing_input(self, tmp_path):
        write_airr(tmp_path / "present.tsv", [igh_row("a", 10, 99.0, "TGTGCGTGG")])
        result = subprocess.run(
            [sys.executable, str(BIN / "clonality_metrics.py"),
             "--sample-id", "TEST",
             "--trust4-airr", str(tmp_path / "present.tsv"),
             "--igblast-airr", str(tmp_path / "absent.tsv"),
             "--out-metrics", str(tmp_path / "m.json"),
             "--out-clonotypes", str(tmp_path / "c.tsv"),
             "--out-top", str(tmp_path / "t.tsv")],
            capture_output=True, text=True,
        )
        assert result.returncode != 0, (
            "A missing AIRR input must fail loudly, not produce a signed "
            "'no V(D)J signal detected' report."
        )


# ---------------------------------------------------------------------------
# HIGH/MEDIUM: metrics.json was not valid JSON for monoclonal samples
# ---------------------------------------------------------------------------
class TestMetricsJsonIsStrictlyValid:
    def test_monoclonal_sample_emits_null_not_nan(self, tmp_path):
        """One clone → clonality_index is undefined, and must serialise null."""
        rows = [igh_row("only", 500, 99.5, "TGTGCGAGATTTTGG")]
        raw = tmp_path / "out.metrics.json"
        run_metrics(tmp_path, rows, rows)  # asserts strict json.loads internally

        text = (tmp_path / "out.metrics.json").read_text(encoding="utf-8")
        assert "NaN" not in text, "bare NaN is not valid JSON"
        assert "Infinity" not in text

        metrics = json.loads(text)
        igh = metrics["per_locus"]["IGH"]
        ci = igh.get("clonality_index")
        assert ci is None or math.isfinite(ci), (
            "an undefined clonality index must be null, never 0.0, so a "
            "monoclonal sample is not rendered as perfectly polyclonal"
        )

    def test_json_safe_converts_non_finite(self):
        out = cm._json_safe({"a": float("nan"), "b": [float("inf"), 1.5], "c": "x"})
        assert out == {"a": None, "b": [None, 1.5], "c": "x"}


# ---------------------------------------------------------------------------
# HIGH: total_input_reads was fabricated, disabling the coverage warning
# ---------------------------------------------------------------------------
class TestTotalInputReadsNotFabricated:
    def test_unknown_library_size_is_null(self, tmp_path):
        rows = [igh_row(f"c{i}", 10, 99.0, f"TGTGCGAGA{i:04d}TGG") for i in range(5)]
        comp = run_metrics(tmp_path, rows, rows)["composition"]

        assert comp["total_input_reads_known"] is False
        assert comp["total_input_reads"] is None, (
            "substituting the V(D)J count makes the V(D)J fraction exactly "
            "1.0, which makes the capture-underperformance warning unreachable"
        )

    def test_supplied_library_size_is_kept(self, tmp_path):
        rows = [igh_row(f"c{i}", 10, 99.0, f"TGTGCGAGA{i:04d}TGG") for i in range(5)]
        comp = run_metrics(
            tmp_path, rows, rows, extra_args=("--total-input-reads", "1000000")
        )["composition"]

        assert comp["total_input_reads_known"] is True
        assert comp["total_input_reads"] == 1_000_000


# ---------------------------------------------------------------------------
# MEDIUM: validation returned pass=True when nothing was actually checked
# ---------------------------------------------------------------------------
class TestValidationCannotPassVacuously:
    def test_no_matching_expectation_is_not_a_pass(self, tmp_path):
        import grade_validation as gv

        # A sample that processed fine and has real signal, but for which
        # validation_expected.json carries no expectations at all.
        metrics_path = tmp_path / "m.json"
        metrics_path.write_text(json.dumps({
            "aggregate": {"n_clonotypes": 5, "top_clone_fraction": 0.4},
            "per_locus": {"IGH": {"n_clonotypes": 5, "clonality_index": 0.5,
                                  "n_reads": 500, "top_clone_fraction": 0.4}},
            "composition": {"vdj_assigned_reads": 500},
            "ighv_status": {"reads_total": 500, "dominant_clone_status": "mutated"},
        }), encoding="utf-8")

        clonotypes_path = tmp_path / "c.tsv"
        pd.DataFrame([{"locus": "IGH", "read_count": 500, "v_call": "IGHV3-23*01",
                       "j_call": "IGHJ4*02", "junction_aa": "CARDYW"}]).to_csv(
            clonotypes_path, sep="\t", index=False)

        graded = gv.grade_sample(
            sample_id="UNKNOWN_SAMPLE",
            metrics_path=metrics_path,
            clonotypes_path=clonotypes_path,
            expected={},
        )
        assert graded["pass"] is False, (
            "an empty checklist must not report a pass — it green-lights a "
            "null result as validated"
        )
        assert graded["reasons_fail"], "the reason must be stated explicitly"


# ---------------------------------------------------------------------------
# A negative call must not be made on read depth that cannot support it
# ---------------------------------------------------------------------------
class TestLowYieldNegativesStayIndeterminate:
    """You may call a positive on low input; you may not call a negative.

    The consolidated verdict rule assesses clonality before yield so that a
    low-input sample with an unambiguous dominant clone is still reported
    clonal. Applied to the negative case that same ordering produced
    "No clonal expansion" — rendered as a green pill with no caveat — for
    samples the assay had not sampled deeply enough to exclude a clone.
    """

    @staticmethod
    def _diverse_no_clonal_locus():
        # 40 clonotypes, none dominant: diverse-looking, but on what depth?
        return {"IGH": {"n_clonotypes": 40, "n_reads": 150,
                        "clonality_index": 0.05, "top_clone_fraction": 0.04}}

    def test_shallow_negative_is_indeterminate(self):
        import lymphix_common as lc
        category, loci = lc.verdict_category(
            self._diverse_no_clonal_locus(), vdj_reads=150, n_clonotypes=40)
        assert category == "indeterminate", (
            "150 V(D)J reads cannot exclude a clone; reporting 'no clonal "
            "expansion' states a confident negative the data cannot support"
        )
        assert loci == []

    def test_deep_negative_is_a_real_negative(self):
        import lymphix_common as lc
        category, _ = lc.verdict_category(
            {"IGH": {"n_clonotypes": 40, "n_reads": 50_000,
                     "clonality_index": 0.05, "top_clone_fraction": 0.04}},
            vdj_reads=50_000, n_clonotypes=40)
        assert category == "no_clonal", (
            "'no_clonal' must remain reachable — it is the only category that "
            "means 'adequately sampled and genuinely diverse'"
        )

    def test_shallow_positive_still_stands(self):
        import lymphix_common as lc
        category, loci = lc.verdict_category(
            {"IGH": {"n_clonotypes": 1, "n_reads": 150,
                     "clonality_index": 1.0, "top_clone_fraction": 1.0}},
            vdj_reads=150, n_clonotypes=1)
        assert category == "clonal", (
            "a real positive on low input must not be buried as indeterminate"
        )
        assert loci == ["IGH"]


# ---------------------------------------------------------------------------
# Report must render when the library size is unknown
# ---------------------------------------------------------------------------
class TestReportRendersWithoutTotalInputReads:
    """total_input_reads became null for unknown library sizes; the report
    header formatted it unguarded and raised TypeError on every sample that
    did not pass --total-input-reads. The smoke test missed it because both
    its samples supply the value.
    """

    @staticmethod
    def _metrics(total_input):
        return {
            "sample_id": "NULLTOTAL",
            "aggregate": {"n_clonotypes": 3, "n_reads": 300,
                          "top_clone_fraction": 0.5, "clonality_index": 0.4},
            "per_locus": {"IGH": {"n_clonotypes": 3, "n_reads": 300,
                                  "clonality_index": 0.4, "top_clone_fraction": 0.5}},
            "composition": {"vdj_assigned_reads": 300,
                            "total_input_reads": total_input,
                            "total_input_reads_known": total_input is not None,
                            "fractions": {}, "reads": {}},
            "ighv_status": None,
            "cdr3_inference": {"n_clonotypes": 3},
        }

    def test_verdict_survives_unknown_total(self):
        import generate_report as gr
        verdict = gr.compute_verdict(self._metrics(None), None)
        assert verdict["category"], "an unknown library size must not crash the verdict"

    def test_unknown_total_is_not_fabricated_in_the_warning(self):
        import generate_report as gr
        verdict = gr.compute_verdict(self._metrics(None), None)
        warnings = " ".join(verdict.get("warnings") or [])
        assert "could not be assessed" in warnings or "not supplied" in warnings, (
            "an unknown library size must be stated, not silently skipped"
        )

    def test_known_total_still_computes_the_check(self):
        import generate_report as gr
        # 300 V(D)J reads out of 10,000,000 is far below the yield floor
        verdict = gr.compute_verdict(self._metrics(10_000_000), None)
        warnings = " ".join(verdict.get("warnings") or [])
        assert "%" in warnings, "with a known total the capture check must actually fire"
