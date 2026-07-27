"""
Unit tests for bin/lymphix_common.py — the consolidated clonality rule.

Six copies of "is this locus clonal" used to live in bin/, and they had
drifted. These tests pin the one rule that replaced them and, just as
importantly, pin the property that made the drift dangerous in the first
place: every consumer must reach the same verdict from the same metrics.json.
The last class asserts exactly that, driving the real entry point in the
per-sample report, the cohort overview, the cohort summary table and the
validation grader rather than a re-implementation of any of them.

Run:    pytest tests/test_lymphix_common.py -v
"""
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import lymphix_common as lc  # noqa: E402


def locus(ci=None, top=None, clones=0, reads=0) -> dict:
    """A per-locus metrics block, shaped as clonality_metrics.py writes it."""
    return {"clonality_index": ci, "top_clone_fraction": top,
            "n_clonotypes": clones, "n_reads": reads}


def clonal(ci=None, top=None, clones=0, reads=0, **kwargs) -> bool:
    """Shorthand for the locus rule, so the assertions stay readable."""
    return lc.is_locus_clonal(clonality_index=ci, top_clone_fraction=top,
                              n_clonotypes=clones, n_reads=reads, **kwargs)


def metrics_json(sample_id: str, per_locus: dict,
                 vdj_reads: int, n_clonotypes: int) -> dict:
    """A metrics.json carrying only the fields the verdict consumers read."""
    full = {L: locus() for L in lc.LOCI}
    full.update(per_locus)
    return {
        "sample_id":   sample_id,
        "per_locus":   full,
        "aggregate":   {"n_clonotypes": n_clonotypes, "n_reads": vdj_reads,
                        "top_clone_fraction": 0.0},
        "composition": {"vdj_assigned_reads": vdj_reads, "fractions": {}},
        "ighv_status": {},
    }


# ---------------------------------------------------------------------------
# The locus-level rule
# ---------------------------------------------------------------------------
class TestIsLocusClonal:
    def test_skewed_and_dominant_is_clonal(self):
        assert clonal(ci=0.55, top=0.80, clones=40, reads=5000) is True

    def test_skewed_without_dominance_is_not_clonal(self):
        # Oligoclonal / reactive: the locus is uneven but no single clone owns
        # it. Calling that clonal is how a reactive node becomes a lymphoma.
        assert clonal(ci=0.55, top=0.10, clones=40, reads=5000) is False

    def test_dominant_without_skew_is_not_clonal(self):
        assert clonal(ci=0.05, top=0.90, clones=3, reads=5000) is False

    def test_thresholds_are_inclusive(self):
        assert clonal(ci=lc.LOCUS_CLONAL_INDEX_THRESHOLD,
                      top=lc.TOP_CLONE_FRACTION_THRESHOLD,
                      clones=10, reads=1000) is True

    def test_just_below_either_threshold_is_not_clonal(self):
        assert clonal(ci=0.2999, top=0.99, clones=10, reads=1000) is False
        assert clonal(ci=0.99, top=0.1999, clones=10, reads=1000) is False

    def test_single_clone_with_enough_reads_is_clonal(self):
        # The clonality index is undefined at n=1 — there is no distribution to
        # measure — so the multi-clone route can never fire on a genuinely
        # monoclonal locus. Without this route a cell line with one IGH clone
        # and hundreds of reads reads as 'no clonal expansion'.
        assert clonal(clones=1, reads=943) is True

    def test_single_clone_with_too_few_reads_is_not_clonal(self):
        assert clonal(clones=1, reads=lc.SINGLE_CLONE_READS_MIN - 1) is False

    def test_single_clone_read_threshold_is_inclusive(self):
        assert clonal(clones=1, reads=lc.SINGLE_CLONE_READS_MIN) is True

    def test_empty_locus_is_not_clonal(self):
        assert lc.is_locus_clonal() is False

    def test_nan_clonality_index_is_not_clonal(self):
        # pandas hands us NaN where JSON gives null. NaN >= 0.30 happens to be
        # False in Python, but leaving that implicit is how a comparison flips
        # silently on the next refactor.
        assert clonal(ci=float("nan"), top=0.9, clones=5, reads=500) is False

    def test_require_dominance_false_drops_the_top_clone_test(self):
        # The composition classifier's variant: it has already applied its own
        # per-clone dominance gate to every read it is about to bin.
        assert clonal(ci=0.55, top=0.10, clones=40, reads=5000,
                      require_dominance=False) is True

    def test_clonality_min_override(self):
        assert clonal(ci=0.40, top=0.90, clones=10, reads=500,
                      clonality_min=0.50) is False
        assert clonal(ci=0.60, top=0.90, clones=10, reads=500,
                      clonality_min=0.50) is True


# ---------------------------------------------------------------------------
# Which loci, and in what order
# ---------------------------------------------------------------------------
class TestClonalLoci:
    def test_returns_loci_in_canonical_order(self):
        per_locus = {
            "TRB": locus(ci=0.7, top=0.6, clones=30, reads=4000),
            "IGH": locus(ci=0.8, top=0.9, clones=50, reads=9000),
        }
        assert lc.clonal_loci(per_locus) == ["IGH", "TRB"]

    def test_missing_locus_block_is_absent_not_clonal(self):
        assert lc.clonal_loci({"IGH": None}) == []

    def test_empty_input(self):
        assert lc.clonal_loci({}) == []
        assert lc.clonal_loci(None) == []


# ---------------------------------------------------------------------------
# The sample-level category
# ---------------------------------------------------------------------------
class TestVerdictCategory:
    def test_no_vdj_reads_is_no_signal(self):
        assert lc.verdict_category({}, 0, 0) == ("no_signal", [])

    def test_reads_but_no_clonotypes_is_no_signal(self):
        # Everything was filtered out by --min-clone-count. There is nothing
        # left to assess, and reporting that as 'no clonal expansion' states a
        # negative the data does not support.
        assert lc.verdict_category({}, 5000, 0) == ("no_signal", [])

    def test_clonal_locus_gives_clonal(self):
        per_locus = {"IGH": locus(ci=0.8, top=0.9, clones=50, reads=9000)}
        assert lc.verdict_category(per_locus, 9000, 50) == ("clonal", ["IGH"])

    def test_few_clonotypes_and_no_clone_is_indeterminate(self):
        assert lc.verdict_category({}, 400, 3)[0] == "indeterminate"

    def test_many_clonotypes_and_no_clone_is_no_clonal(self):
        assert lc.verdict_category({}, 40000, 600)[0] == "no_clonal"

    def test_indeterminate_boundary(self):
        below = lc.INDETERMINATE_MAX_CLONOTYPES - 1
        at    = lc.INDETERMINATE_MAX_CLONOTYPES
        assert lc.verdict_category({}, 400, below)[0] == "indeterminate"
        assert lc.verdict_category({}, 400, at)[0] == "no_clonal"

    def test_low_yield_does_not_erase_a_clonal_call(self):
        # THE ONE DELIBERATE BEHAVIOUR CHANGE OF THE CONSOLIDATION.
        # cohort_summary.py used to short-circuit to 'indeterminate' for any
        # sample under LOW_VDJ_YIELD_ABSOLUTE V(D)J reads, before it looked at
        # clonality at all, while the per-sample report called the same sample
        # clonal and attached a low-yield warning. Clonality is now assessed
        # first everywhere; low yield is a warning on the report, not a verdict
        # that hides the finding.
        assert 150 < lc.LOW_VDJ_YIELD_ABSOLUTE
        per_locus = {"IGH": locus(clones=1, reads=150)}
        assert lc.verdict_category(per_locus, 150, 1) == ("clonal", ["IGH"])

    def test_biclonal_sample_lists_both_loci(self):
        per_locus = {
            "IGH": locus(ci=0.8, top=0.9, clones=50, reads=9000),
            "TRB": locus(ci=0.7, top=0.6, clones=30, reads=4000),
        }
        assert lc.verdict_category(per_locus, 13000, 80) == ("clonal", ["IGH", "TRB"])


class TestLineageVerdict:
    def test_bcr_locus_gives_clonal_b(self):
        assert lc.lineage_verdict("clonal", ["IGH"]) == "clonal_B"

    def test_tcr_locus_gives_clonal_t(self):
        assert lc.lineage_verdict("clonal", ["TRB"]) == "clonal_T"

    def test_biclonal_is_labelled_clonal_b(self):
        # A presentation compromise: the cohort strip plot has one colour per
        # sample. The loci list travels alongside, so nothing is lost.
        assert lc.lineage_verdict("clonal", ["IGH", "TRB"]) == "clonal_B"

    def test_non_clonal_categories_pass_through(self):
        for category in ("no_signal", "no_clonal", "indeterminate"):
            assert lc.lineage_verdict(category, []) == category

    def test_every_label_has_a_colour_and_a_name(self):
        for category in ("clonal", "clonal_B", "clonal_T",
                         "no_clonal", "indeterminate", "no_signal"):
            assert category in lc.VERDICT_COLORS
            assert category in lc.VERDICT_LABELS


class TestPerLocusFromFlatRow:
    def test_round_trip_through_the_flat_cohort_row(self):
        row = {f"{L}_clonality": None for L in lc.LOCI}
        row.update({f"{L}_top_fraction": None for L in lc.LOCI})
        row.update({f"{L}_n_clones": 0 for L in lc.LOCI})
        row.update({f"{L}_n_reads": 0 for L in lc.LOCI})
        row.update({"IGH_clonality": 0.8, "IGH_top_fraction": 0.9,
                    "IGH_n_clones": 50, "IGH_n_reads": 9000})
        assert lc.clonal_loci(lc.per_locus_from_flat_row(row)) == ["IGH"]


# ---------------------------------------------------------------------------
# Palette and asset consistency
# ---------------------------------------------------------------------------
class TestPalette:
    def test_every_composition_pool_has_a_label_and_a_colour(self):
        for key in lc.COMP_ORDER:
            assert key in lc.COMP_LABELS
            assert key in lc.COMP_LABELS_LONG
            assert lc.COMP_COLORS[key].startswith("#")

    def test_short_and_long_labels_cover_the_same_pools(self):
        assert set(lc.COMP_LABELS) == set(lc.COMP_LABELS_LONG)

    def test_every_locus_has_a_colour(self):
        for L in lc.LOCI:
            assert L in lc.LOCUS_COLORS

    def test_lineage_groups_partition_the_locus_list(self):
        assert set(lc.BCR_LOCI) | set(lc.TCR_LOCI) == set(lc.LOCI)
        assert not set(lc.BCR_LOCI) & set(lc.TCR_LOCI)
        assert sorted(lc.LOCUS_ORDER) == sorted(lc.LOCI)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
class TestHelpers:
    @pytest.mark.parametrize("value,expected", [
        (None, None), (float("nan"), None), ("", None), ("abc", None),
        ("0.5", 0.5), (1, 1.0), (0, 0.0),
    ])
    def test_safe_float(self, value, expected):
        assert lc.safe_float(value) == expected

    def test_safe_float_keeps_zero_distinct_from_missing(self):
        # A missing clonality index is not a clonality index of zero.
        assert lc.safe_float(0) == 0.0
        assert lc.safe_float(None) is None

    @pytest.mark.parametrize("value,expected", [
        (None, ""), (float("nan"), ""), ("IGHV3-23*01", "IGHV3-23*01"), (7, "7"),
    ])
    def test_safe_str(self, value, expected):
        assert lc.safe_str(value) == expected

    @pytest.mark.parametrize("call,expected", [
        ("IGHV3-23*01", "IGHV3-23"),
        ("IGHV3-23*01,IGHV3-23*04", "IGHV3-23"),
        ("IGHV3-23", "IGHV3-23"),
        (None, ""),
        (float("nan"), ""),
    ])
    def test_gene_name(self, call, expected):
        assert lc.gene_name(call) == expected

    def test_load_logo_svg_returns_embeddable_markup(self):
        svg = lc.load_logo_svg()
        assert 'class="logo"' in svg
        assert svg.lstrip().startswith("<")


# ---------------------------------------------------------------------------
# Every consumer must agree — the point of the whole exercise
# ---------------------------------------------------------------------------
CONSUMER_CASES = [
    # (name, per_locus, vdj_reads, n_clonotypes, expected_category, expected_loci)
    ("monoclonal_igh",
     {"IGH": locus(clones=1, reads=943)}, 943, 1, "clonal", ["IGH"]),
    ("low_yield_monoclonal_igh",
     {"IGH": locus(clones=1, reads=150)}, 150, 1, "clonal", ["IGH"]),
    ("polyclonal",
     {}, 40000, 600, "no_clonal", []),
    ("too_few_clonotypes",
     {}, 400, 3, "indeterminate", []),
    ("no_signal",
     {}, 0, 0, "no_signal", []),
    ("clonal_trb",
     {"TRB": locus(ci=0.7, top=0.6, clones=30, reads=4000)}, 4000, 30, "clonal", ["TRB"]),
    ("skewed_but_not_dominant",
     {"IGH": locus(ci=0.55, top=0.10, clones=40, reads=5000)}, 5000, 40, "no_clonal", []),
    ("biclonal",
     {"IGH": locus(ci=0.8, top=0.9, clones=50, reads=9000),
      "TRB": locus(ci=0.7, top=0.6, clones=30, reads=4000)},
     13000, 80, "clonal", ["IGH", "TRB"]),
]


class TestConsumersAgree:
    @pytest.mark.parametrize("name,per_locus,vdj,n,category,loci",
                             CONSUMER_CASES, ids=[c[0] for c in CONSUMER_CASES])
    def test_all_consumers_reach_the_same_verdict(self, tmp_path, name, per_locus,
                                                  vdj, n, category, loci):
        import cohort_report
        import cohort_summary
        import generate_report
        import grade_validation

        metrics = metrics_json(name, per_locus, vdj, n)
        path = tmp_path / name / f"{name}.metrics.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(metrics), encoding="utf-8")

        # 1. per-sample clinical report
        verdict = generate_report.compute_verdict(metrics)
        assert verdict["category"] == category, "generate_report category"
        assert verdict["clonal_loci"] == loci, "generate_report clonal loci"

        # 2. cohort overview
        assert cohort_report.derive_verdict(metrics) == (category, loci), "cohort_report"

        # 3. cohort summary table (flat row, lineage-split label)
        row = cohort_summary.summarise_sample(path)
        assert row["verdict"] == lc.lineage_verdict(category, loci), "cohort_summary"

        # 4. validation grader
        graded = grade_validation.grade_sample(name, path,
                                               tmp_path / "no_clonotypes.tsv",
                                               expected={})
        assert graded["observed_category"] == category, "grade_validation category"
        assert graded["observed_clonal_loci"] == loci, "grade_validation clonal loci"
