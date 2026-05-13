"""
Unit tests for clonality_metrics.py.

Run:    pytest tests/test_clonality_metrics.py -v
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import clonality_metrics as cm  # noqa: E402


# ---------------------------------------------------------------------------
# Pure-math metric tests
# ---------------------------------------------------------------------------
class TestDiversityMath:
    def test_shannon_uniform(self):
        # 100 clones, 10 reads each → max entropy = log(100)
        counts = np.full(100, 10)
        assert math.isclose(cm.shannon(counts), math.log(100), rel_tol=1e-6)

    def test_shannon_empty(self):
        assert cm.shannon(np.array([])) == 0.0

    def test_shannon_single_clone(self):
        # Only one clone present → entropy is 0
        assert cm.shannon(np.array([1000])) == 0.0

    def test_simpson_d_uniform(self):
        # Simpson D = 1/N when uniform
        assert math.isclose(cm.simpson_d(np.full(50, 5)), 1 / 50, rel_tol=1e-6)

    def test_gini_uniform_is_zero(self):
        assert cm.gini(np.full(100, 10)) == pytest.approx(0.0, abs=1e-6)

    def test_gini_extreme_inequality_close_to_one(self):
        # One dominant clone + lots of singletons → Gini close to 1
        counts = np.concatenate(([1000], np.ones(100)))
        assert cm.gini(counts) > 0.85

    def test_d50_monoclonal(self):
        # 99% of reads in clone 0 → D50 = 1
        counts = np.concatenate(([9900], np.full(100, 1)))
        assert cm.d50(counts) == 1

    def test_d50_uniform(self):
        # 100 uniform clones → need 50 to cover 50%
        assert cm.d50(np.full(100, 10)) == 50


class TestClonalityIndex:
    def test_polyclonal_is_zero(self):
        # Uniform repertoire → clonality_index = 0
        assert cm.clonality_index(np.full(100, 10)) == pytest.approx(0.0, abs=1e-9)

    def test_monoclonal_is_near_one(self):
        # One huge clone + a tail of singletons → clonality near 1
        counts = np.concatenate(([10000], np.ones(20)))
        ci = cm.clonality_index(counts)
        assert ci > 0.7
        assert ci < 1.0

    def test_single_clone_is_nan(self):
        # Adaptive convention: undefined when N=1
        assert math.isnan(cm.clonality_index(np.array([500])))

    def test_empty_is_nan(self):
        assert math.isnan(cm.clonality_index(np.array([])))


# ---------------------------------------------------------------------------
# Composition / locus inference
# ---------------------------------------------------------------------------
class TestLocusInference:
    @pytest.mark.parametrize("v_call,expected", [
        ("IGHV1-2*02",   "IGH"),
        ("IGKV3-20*01",  "IGK"),
        ("IGLV2-14*01",  "IGL"),
        ("TRAV21*01",    "TRA"),
        ("TRBV20-1*01",  "TRB"),
        ("TRGV9*01",     "TRG"),
        ("TRDV1*01",     "TRD"),
        ("",             None),
        (None,           None),
        ("FOOBAR",       None),
    ])
    def test_infer_locus(self, v_call, expected):
        assert cm.infer_locus(v_call) == expected


class TestComposition:
    @pytest.fixture
    def clonal_igh_df(self):
        import pandas as pd
        # One dominant IGH clone (95% of locus) + small background everywhere
        rows = [
            dict(locus="IGH", read_count=9500, junction="AAA", junction_aa="C"),
            dict(locus="IGH", read_count=500,  junction="GGG", junction_aa="C"),
            dict(locus="TRB", read_count=200,  junction="TTT", junction_aa="C"),
            dict(locus="IGK", read_count=100,  junction="CCC", junction_aa="C"),
            dict(locus="IGL", read_count=80,   junction="AAA", junction_aa="C"),
        ]
        return pd.DataFrame(rows)

    def test_clonal_igh_dominant(self, clonal_igh_df):
        out = cm.compute_composition(clonal_igh_df,
                                      total_input_reads=20000,
                                      clonal_threshold=0.05)
        # Both IGH clones meet the >=5% locus-fraction bar (9500/10000=95%, 500/10000=5%)
        assert out["reads"]["clonal_IGH"] == 10000
        # The lone TRB clone is 100% of TRB locus, so it's clonal_TRB (not polyclonal)
        assert out["reads"]["clonal_TRB"] == 200
        # Lone IGK and IGL clones are also 100% of their loci, so each is clonal-restricted
        assert out["reads"]["clonal_IGK_kappa"] == 100
        assert out["reads"]["clonal_IGL_lambda"] == 80
        assert out["reads"]["polyclonal_B"] == 0
        assert out["reads"]["polyclonal_T"] == 0
        assert out["reads"]["background"] == 20000 - (10000 + 200 + 100 + 80)
        assert sum(out["fractions"].values()) == pytest.approx(1.0, abs=1e-6)

    def test_background_when_total_unknown(self, clonal_igh_df):
        out = cm.compute_composition(clonal_igh_df,
                                      total_input_reads=None,
                                      clonal_threshold=0.05)
        # Without a denominator, background = 0
        assert out["reads"]["background"] == 0
        assert out["total_input_reads_known"] is False

    def test_kappa_lambda_balanced(self):
        import pandas as pd
        df = pd.DataFrame([
            dict(locus="IGK", read_count=1000, junction="A", junction_aa="C"),
            dict(locus="IGL", read_count=1000, junction="A", junction_aa="C"),
        ])
        out = cm.compute_composition(df, 5000, 0.05)
        assert out["kappa_lambda_ratio"] == 1.0
        assert out["kappa_lambda_call"] == "balanced"

    def test_kappa_restricted(self):
        import pandas as pd
        df = pd.DataFrame([
            dict(locus="IGK", read_count=4000, junction="A", junction_aa="C"),
            dict(locus="IGL", read_count=100,  junction="A", junction_aa="C"),
        ])
        out = cm.compute_composition(df, 5000, 0.05)
        assert out["kappa_lambda_ratio"] == 40.0
        assert out["kappa_lambda_call"] == "kappa_restricted"

    def test_lambda_restricted(self):
        import pandas as pd
        df = pd.DataFrame([
            dict(locus="IGK", read_count=100,  junction="A", junction_aa="C"),
            dict(locus="IGL", read_count=4000, junction="A", junction_aa="C"),
        ])
        out = cm.compute_composition(df, 5000, 0.05)
        assert out["kappa_lambda_ratio"] == pytest.approx(0.025)
        assert out["kappa_lambda_call"] == "lambda_restricted"


# ---------------------------------------------------------------------------
# Summary roll-up
# ---------------------------------------------------------------------------
def test_summarise_polyclonal_repertoire():
    counts = np.full(50, 20)
    s = cm.summarise(counts)
    assert s["n_clonotypes"] == 50
    assert s["n_reads"] == 1000
    assert s["top_clone_fraction"] == pytest.approx(0.02)
    assert s["clonality_index"] == pytest.approx(0.0, abs=1e-9)
    assert s["D50"] == 25
