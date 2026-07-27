"""
Tests for opt-in clonotype collapsing in clonality_metrics.py.

TRUST4 emits one AIRR row per (assembly, CDR3 variant), so one rearrangement
arrives as a dominant row plus a tail of near-identical rows carrying a few
reads each. Counting those rows as clonotypes inflates N in Shannon, Simpson,
Gini, D50 and the clonality index at once, and can hand "dominant clone" to a
fragment.

Collapsing is therefore available, and DEFAULT OFF. These tests pin down all
three halves of that contract:

  * with the flag absent, every number is exactly what it was before;
  * with the flag on, assembly variants of one clone become one clone whose
    reads are the sum of its parts;
  * with the flag on, two genuinely distinct clones stay two clones.

Run:    pytest tests/test_clonotype_collapse.py -v
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
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

# A 36-nt IGH junction. Variants below are made by substituting single bases,
# which is how the real assembly-noise tail differs from its parent.
JUNCTION = "TGTGATTACTACTACTACTACGGTATGGACGTCTGG"


def _variant(junction: str, position: int, base: str) -> str:
    """Return `junction` with one base substituted — a Hamming-1 neighbour."""
    assert junction[position] != base
    return junction[:position] + base + junction[position + 1:]


def igh_row(seq_id: str, reads: int, junction: str,
            v_call: str = "IGHV3-15*03",
            junction_aa: str = "CDYYYYYGMDVW",
            identity: float | None = 99.0) -> dict:
    return {
        "sequence_id":     seq_id,
        "v_call":          v_call,
        "j_call":          "IGHJ6*02",
        "junction":        junction,
        "junction_aa":     junction_aa,
        "duplicate_count": reads,
        "v_identity":      identity,
        "v_cigar":         "150M",
    }


def write_airr(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows, columns=AIRR_COLUMNS).to_csv(path, sep="\t", index=False)


def run_metrics(tmp_path: Path, rows: list[dict], extra_args=()) -> dict:
    """Invoke clonality_metrics.py as the pipeline does, return metrics.json."""
    t4 = tmp_path / "trust4.airr.tsv"
    ig = tmp_path / "igblast.airr.tsv"
    write_airr(t4, rows)
    write_airr(ig, rows)
    out_metrics = tmp_path / "out.metrics.json"

    result = subprocess.run(
        [sys.executable, str(BIN / "clonality_metrics.py"),
         "--sample-id", "COLLAPSE_TEST",
         "--trust4-airr", str(t4),
         "--igblast-airr", str(ig),
         "--out-metrics", str(out_metrics),
         "--out-clonotypes", str(tmp_path / "out.clonotypes.tsv"),
         "--out-top", str(tmp_path / "out.top.tsv"),
         *extra_args],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"clonality_metrics failed:\n{result.stderr}"
    return json.loads(out_metrics.read_text(encoding="utf-8"))


def one_assembly_with_noise() -> list[dict]:
    """One clone, 1,009 reads, split by TRUST4 across three variant rows.

    Same assembly (`assemble5`), same V and J, junctions one substitution
    apart. This is the shape the real cohorts show: a dominant row plus a tail
    at 0.1-0.5% of it.
    """
    return [
        igh_row("assemble5_0", 1000, JUNCTION),
        igh_row("assemble5_1", 5,    _variant(JUNCTION, 9, "A")),
        igh_row("assemble5_2", 4,    _variant(JUNCTION, 18, "G"),
                junction_aa="CDYYYYDGMDVW"),
    ]


# ---------------------------------------------------------------------------
# The flag is off by default, and off means unchanged
# ---------------------------------------------------------------------------
class TestDefaultOffChangesNothing:
    """Collapsing rewrites n_clonotypes and therefore every diversity metric.
    Switching it on silently would change published numbers with no way to
    tell which run produced what, so the absent flag must be a no-op.
    """

    def test_rows_are_still_counted_as_clonotypes(self, tmp_path):
        metrics = run_metrics(tmp_path, one_assembly_with_noise())
        assert metrics["aggregate"]["n_clonotypes"] == 3, (
            "without the flag, the three assembly-variant rows must still be "
            "reported as three clonotypes — the old behaviour, unchanged"
        )
        assert metrics["aggregate"]["n_reads"] == 1009

    def test_every_metric_matches_the_raw_row_counts(self, tmp_path):
        """The off-path numbers must be the metrics of the raw rows exactly."""
        metrics = run_metrics(tmp_path, one_assembly_with_noise())
        expected = cm.summarise(np.array([1000, 5, 4]))
        agg = metrics["aggregate"]
        for field in ("n_clonotypes", "n_reads", "D50"):
            assert agg[field] == expected[field]
        for field in ("top_clone_fraction", "shannon_H", "clonality_index",
                      "simpson_D", "gini"):
            assert agg[field] == pytest.approx(expected[field])

    def test_the_absence_of_collapsing_is_stated_not_implied(self, tmp_path):
        metrics = run_metrics(tmp_path, one_assembly_with_noise())
        assert metrics["collapse_clonotypes"] is False
        assert metrics["collapse_key"] is None
        block = metrics["clonotype_collapse"]
        assert block["applied"] is False
        assert block["n_rows_in"] == block["n_clones_out"] == 3
        assert block["n_rows_merged"] == 0

    def test_the_flag_is_the_only_difference(self, tmp_path):
        """Same input, flag on vs off: only the collapse-dependent numbers move.

        Anything outside the clonotype counts — the read length, the germline
        filter, the recorded thresholds — must be identical, so a difference
        between two runs can always be attributed to the flag.
        """
        rows = one_assembly_with_noise()
        (tmp_path / "off").mkdir()
        (tmp_path / "on").mkdir()
        off = run_metrics(tmp_path / "off", rows)
        on = run_metrics(tmp_path / "on", rows,
                         extra_args=("--collapse-clonotypes",
                                     "--collapse-key", "locus_junction_nt_hamming1"))
        assert off["germline_rearrangement_filter"] == on["germline_rearrangement_filter"]
        assert off["read_length"] == on["read_length"]
        assert off["min_clone_count"] == on["min_clone_count"]
        assert off["aggregate"]["n_clonotypes"] != on["aggregate"]["n_clonotypes"]


# ---------------------------------------------------------------------------
# Flag on: assembly variants of one clone become one clone
# ---------------------------------------------------------------------------
class TestAssemblyVariantsCollapse:
    def test_three_variants_of_one_assembly_become_one_clone(self, tmp_path):
        metrics = run_metrics(tmp_path, one_assembly_with_noise(),
                              extra_args=("--collapse-clonotypes",
                                          "--collapse-key",
                                          "locus_junction_nt_hamming1"))
        assert metrics["aggregate"]["n_clonotypes"] == 1
        assert metrics["aggregate"]["n_reads"] == 1009, (
            "reads from variants of the same assembly are summed: TRUST4's "
            "per-variant weights partition the reads spanning that CDR3, so "
            "the sum reconstructs the assembly's support"
        )
        assert metrics["per_locus"]["IGH"]["n_clonotypes"] == 1
        assert metrics["per_locus"]["IGH"]["top_clone_fraction"] == pytest.approx(1.0)

    def test_the_collapse_is_reported_rows_in_clones_out_per_locus(self, tmp_path):
        rows = one_assembly_with_noise()
        # A second locus, untouched by collapsing, to prove the per-locus
        # figures are not just the aggregate repeated.
        rows.append({**igh_row("assemble9_0", 50, "TGTGCCAGCAGCTTCGGGACAGGGGAGCTGTTTTTT"),
                     "v_call": "TRBV20-1*01", "j_call": "TRBJ2-1*01",
                     "junction_aa": "CASSFGTGELFF"})
        metrics = run_metrics(tmp_path, rows,
                              extra_args=("--collapse-clonotypes",
                                          "--collapse-key",
                                          "locus_junction_nt_hamming1"))
        block = metrics["clonotype_collapse"]
        assert block["applied"] is True
        assert block["key"] == "locus_junction_nt_hamming1"
        assert block["n_rows_in"] == 4
        assert block["n_clones_out"] == 2
        assert block["n_rows_merged"] == 2
        assert block["per_locus"]["IGH"] == {"rows_in": 3, "clones_out": 1,
                                            "rows_merged": 2}
        assert block["per_locus"]["TRB"] == {"rows_in": 1, "clones_out": 1,
                                            "rows_merged": 0}
        assert metrics["collapse_clonotypes"] is True
        assert metrics["collapse_key"] == "locus_junction_nt_hamming1"

    def test_the_read_aggregation_rule_is_recorded(self, tmp_path):
        metrics = run_metrics(tmp_path, one_assembly_with_noise(),
                              extra_args=("--collapse-clonotypes",))
        assert metrics["clonotype_collapse"]["read_aggregation"] == \
            cm.COLLAPSE_READ_AGGREGATION

    def test_exact_duplicate_junctions_merge_without_the_hamming_stage(self):
        """One rearrangement assembled twice against paralogous V references.

        Identical junction nt, identical J, different V paralogue
        (IGKV3-15 / IGKV3D-15). This is the only exact merge that does real
        work, which is why the default key excludes the V call.
        """
        junction = "TGTCAGCAGTATAATAACTGGCCTTGGACGTTC"
        df = pd.DataFrame([
            dict(sequence_id="assemble18_0", locus="IGK", v_call="IGKV3D-15*01",
                 j_call="IGKJ1*01", junction=junction, junction_aa="CQQYNNWPWTF",
                 read_count=315, v_identity=100.0),
            dict(sequence_id="assemble32_0", locus="IGK", v_call="IGKV3-15*01",
                 j_call="IGKJ1*01", junction=junction, junction_aa="CQQYNNWPWTF",
                 read_count=170, v_identity=100.0),
        ])
        out, stats = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_EXACT)
        assert len(out) == 1
        assert stats["n_rows_merged"] == 1
        assert out.iloc[0]["v_call"] == "IGKV3D-15*01", \
            "the dominant member supplies the V call"
        assert out.iloc[0]["read_count"] == 485, (
            "reads ARE summed across assemblies. This was originally max, on "
            "the assumption that the two contigs might share reads. Measured "
            "on the real sample these numbers come from "
            "(CMDL20001026_S127_L004): the contigs share a 334 nt block and "
            "are one rearrangement, and only 3.8% of reads in "
            "_assembled_reads.fa appear against more than one contig. Taking "
            "the max discarded 170 of 485 reads — a 35% understatement of the "
            "clone, and of the sample's read total — to avoid an overcount "
            "bounded near 4%."
        )

    def test_collapsing_conserves_the_sample_read_total(self):
        """Reads must not vanish. The discarded reads did not only shrink the
        clone: they disappeared from aggregate.n_reads, the per-locus totals
        and the composition pools, so the sample lost depth it genuinely had.
        """
        junction = "TGTCAGCAGTATAATAACTGGCCTTGGACGTTC"
        df = pd.DataFrame([
            dict(sequence_id="assemble18_0", locus="IGK", v_call="IGKV3D-15*01",
                 j_call="IGKJ1*01", junction=junction, junction_aa="CQQYNNWPWTF",
                 read_count=315, v_identity=100.0),
            dict(sequence_id="assemble32_0", locus="IGK", v_call="IGKV3-15*01",
                 j_call="IGKJ1*01", junction=junction, junction_aa="CQQYNNWPWTF",
                 read_count=170, v_identity=100.0),
        ])
        out, stats = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_EXACT)
        assert stats["n_reads_in"] == 485
        assert stats["n_reads_out"] == 485
        assert stats["n_reads_delta"] == 0, (
            "a negative delta means support was discarded; a positive one "
            "means shared reads were double-counted. Either must be visible "
            "in metrics.json rather than silent."
        )
        assert int(out["read_count"].sum()) == int(df["read_count"].sum())

    def test_the_dominant_row_supplies_v_j_and_junction(self):
        df = pd.DataFrame([
            dict(sequence_id="a1_1", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=_variant(JUNCTION, 9, "A"),
                 junction_aa="CDYYYYYGMDVW", read_count=6, v_identity=97.0),
            dict(sequence_id="a1_0", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=JUNCTION,
                 junction_aa="CDYYYYYGMDVW", read_count=900, v_identity=99.0),
        ])
        out, _ = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_HAMMING1)
        assert len(out) == 1
        assert out.iloc[0]["junction"] == JUNCTION
        assert out.iloc[0]["sequence_id"] == "a1_0"
        assert out.iloc[0]["read_count"] == 906
        assert out.iloc[0]["n_collapsed_rows"] == 2


# ---------------------------------------------------------------------------
# Flag on: distinct clones stay distinct
# ---------------------------------------------------------------------------
class TestDistinctClonesAreNotMerged:
    @staticmethod
    def _two_real_clones() -> list[dict]:
        """Two co-dominant rearrangements, same V and J, three substitutions
        apart. Real in the cohort (A3523_GEO_CappSeq, 712 and 594 reads) and
        the case a distance-only rule gets wrong."""
        other = JUNCTION[:21] + "TAA" + JUNCTION[24:]
        assert sum(a != b for a, b in zip(JUNCTION, other)) == 3
        return [
            igh_row("assemble2_0", 712, JUNCTION),
            igh_row("assemble2_1", 594, other, junction_aa="CDYYYYYYMDVW"),
        ]

    def test_two_distinct_clones_survive_collapsing(self, tmp_path):
        metrics = run_metrics(tmp_path, self._two_real_clones(),
                              extra_args=("--collapse-clonotypes",
                                          "--collapse-key",
                                          "locus_junction_nt_hamming1"))
        assert metrics["aggregate"]["n_clonotypes"] == 2, (
            "two rearrangements three substitutions apart are not one clone "
            "with assembly noise"
        )
        assert metrics["aggregate"]["n_reads"] == 712 + 594

    def test_different_loci_never_merge(self):
        """Even an identical junction is a different clone at a different
        locus. The exact key is scoped to the locus for this reason."""
        df = pd.DataFrame([
            dict(sequence_id="a1_0", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=JUNCTION, junction_aa="CDYY",
                 read_count=100, v_identity=99.0),
            dict(sequence_id="a2_0", locus="TRB", v_call="TRBV20-1*01",
                 j_call="TRBJ2-1*01", junction=JUNCTION, junction_aa="CDYY",
                 read_count=90, v_identity=99.0),
        ])
        out, _ = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_EXACT)
        assert len(out) == 2

    def test_an_abundant_neighbour_is_not_absorbed(self):
        """A single substitution is not enough on its own.

        Somatic hypermutation produces genuine lineage members one nucleotide
        from the parent, and nothing in the sequence tells them apart from a
        sequencing error — only relative abundance does. A neighbour at 30% of
        the parent is far above the gate and must stand.
        """
        df = pd.DataFrame([
            dict(sequence_id="a1_0", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=JUNCTION, junction_aa="CDYY",
                 read_count=1000, v_identity=99.0),
            dict(sequence_id="a1_1", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=_variant(JUNCTION, 9, "A"),
                 junction_aa="CDYY", read_count=300, v_identity=99.0),
        ])
        out, _ = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_HAMMING1)
        assert len(out) == 2
        # ... and the gate is a parameter, not a constant: widen it and the
        # same pair does merge.
        wide, _ = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_HAMMING1,
                                             minor_fraction_max=0.5)
        assert len(wide) == 1

    def test_blank_junctions_are_never_pooled(self):
        """Out-of-frame rows with no junction must not all key to '' and
        become one clone — that merges 'unknown' with 'unknown'."""
        df = pd.DataFrame([
            dict(sequence_id="a1_8", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction="", junction_aa="CDYY",
                 read_count=8, v_identity=99.0),
            dict(sequence_id="a1_12", locus="IGH", v_call="IGHV4-34*01",
                 j_call="IGHJ4*02", junction="", junction_aa="CARG",
                 read_count=6, v_identity=99.0),
        ])
        out, _ = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_HAMMING1)
        assert len(out) == 2

    def test_no_chaining_across_two_substitutions(self):
        """Single-linkage would walk A-B-C and merge A with C, which are two
        substitutions apart. Every member must be within one substitution of
        the anchor instead."""
        a = JUNCTION
        b = _variant(a, 9, "A")
        c = _variant(b, 18, "G")
        assert sum(x != y for x, y in zip(a, c)) == 2
        df = pd.DataFrame([
            dict(sequence_id="a1_0", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=b, junction_aa="CDYY",
                 read_count=1000, v_identity=99.0),
            dict(sequence_id="a1_1", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=a, junction_aa="CDYY",
                 read_count=5, v_identity=99.0),
            dict(sequence_id="a1_2", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=c, junction_aa="CDYY",
                 read_count=4, v_identity=99.0),
        ])
        out, _ = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_HAMMING1)
        # a and c both sit one substitution from the anchor b, so all three
        # merge here — but only through b, never a→c directly.
        assert len(out) == 1
        assert out.iloc[0]["junction"] == b

    def test_no_reads_are_lost_or_invented(self):
        """Whatever the key, the collapsed table must not exceed the input
        reads, and must not silently drop a row's worth of them."""
        rows = one_assembly_with_noise()
        df = pd.DataFrame([
            dict(sequence_id=r["sequence_id"], locus="IGH", v_call=r["v_call"],
                 j_call=r["j_call"], junction=r["junction"],
                 junction_aa=r["junction_aa"], read_count=r["duplicate_count"],
                 v_identity=r["v_identity"])
            for r in rows
        ])
        for key in cm.COLLAPSE_KEYS:
            out, _ = cm.collapse_clonotype_rows(df, key=key)
            assert out["read_count"].sum() <= df["read_count"].sum()
            assert int(out["n_collapsed_rows"].sum()) == len(df)


# ---------------------------------------------------------------------------
# The IGHV call must still read the right identity after collapsing
# ---------------------------------------------------------------------------
class TestIghvIdentitySurvivesCollapse:
    def test_dominant_clone_keeps_its_own_identity(self, tmp_path):
        """The clinical call belongs to the dominant clone. A minor variant
        at a different identity must not move it."""
        rows = [
            igh_row("assemble5_0", 1000, JUNCTION, identity=94.0),
            igh_row("assemble5_1", 5, _variant(JUNCTION, 9, "A"), identity=100.0),
        ]
        metrics = run_metrics(tmp_path, rows,
                              extra_args=("--collapse-clonotypes",
                                          "--collapse-key",
                                          "locus_junction_nt_hamming1"))
        ighv = metrics["ighv_status"]
        assert ighv["dominant_clone_v_identity"] == pytest.approx(94.0), (
            "not an average of 94 and 100 — 97% is on the mutated side of the "
            "98% cutoff and belongs to neither row"
        )
        assert ighv["dominant_clone_status"] == "mutated"
        assert ighv["dominant_clone_reads"] == 1005

    def test_identity_is_borrowed_only_when_the_dominant_row_has_none(self):
        df = pd.DataFrame([
            dict(sequence_id="a1_0", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=JUNCTION, junction_aa="CDYY",
                 read_count=1000, v_identity=99.0, igblast_v_identity=np.nan),
            dict(sequence_id="a1_1", locus="IGH", v_call="IGHV3-15*03",
                 j_call="IGHJ6*02", junction=_variant(JUNCTION, 9, "A"),
                 junction_aa="CDYY", read_count=5, v_identity=99.0,
                 igblast_v_identity=94.0),
        ])
        out, _ = cm.collapse_clonotype_rows(df, key=cm.COLLAPSE_KEY_HAMMING1)
        assert out.iloc[0]["igblast_v_identity"] == pytest.approx(94.0), (
            "members of a collapsed clone differ in the junction, not in V, "
            "so a minor row's identity is a valid stand-in when the dominant "
            "row was never assessed"
        )


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
class TestCollapseHelpers:
    @pytest.mark.parametrize("seq_id,expected", [
        ("assemble5_0",   "assemble5"),
        ("assemble18_12", "assemble18"),
        ("assemble5",     "assemble5"),      # no variant suffix
        ("contig_x",      "contig_x"),       # suffix is not an index
        ("",              ""),
        (None,            ""),
    ])
    def test_assembly_of(self, seq_id, expected):
        assert cm._assembly_of(seq_id) == expected

    @pytest.mark.parametrize("a,b,expected", [
        ("AAAA", "AAAT", True),
        ("AAAA", "AAAA", False),     # identical is not distance 1
        ("AAAA", "AATT", False),
        ("AAAA", "AAA",  False),     # different lengths never compare
    ])
    def test_hamming_is_one(self, a, b, expected):
        assert cm._hamming_is_one(a, b) is expected

    def test_unknown_key_is_rejected(self):
        with pytest.raises(ValueError):
            cm.collapse_clonotype_rows(pd.DataFrame(), key="v_j_junction_aa")

    def test_empty_table_collapses_to_nothing(self):
        out, stats = cm.collapse_clonotype_rows(pd.DataFrame(),
                                                key=cm.COLLAPSE_KEY_EXACT)
        assert out.empty
        assert stats["applied"] is True
        assert stats["n_rows_in"] == 0
        assert stats["n_clones_out"] == 0
