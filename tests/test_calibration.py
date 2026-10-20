import numpy as np
import pandas as pd
import pytest

from plate_truth.calibration import (
    bin_labels,
    bin_rates,
    bootstrap_groups,
    calibration_by_bin,
    identity_bin,
    interval,
    logo_bin_predictions,
    performance,
    smoothed_rate,
    weighted_recovery,
    wilson_interval,
)

EDGES = [90.0, 94.5, 97.0, 99.0]


def test_bin_labels_and_assignment():
    assert bin_labels([90, 97]) == ["none", "<90", "90-97", ">=97"]
    b = identity_bin(pd.Series([np.nan, 85.0, 94.5, 96.99, 97.0, 99.0, 100.0]), EDGES)
    assert list(b) == ["none", "<90", "94.5-97", "94.5-97", "97-99", ">=99", ">=99"]


def test_smoothed_rate():
    assert smoothed_rate(0, 0) == 0.5
    assert smoothed_rate(1, 3) == pytest.approx(0.375)


def test_bin_rates_counts():
    df = pd.DataFrame({"b": pd.Categorical(["x", "x", "y"], categories=["x", "y", "z"]),
                       "g": [True, False, True]})
    r = bin_rates(df, "b", "g")
    assert r.loc["x", "n"] == 2 and r.loc["x", "k"] == 1 and r.loc["x", "rate"] == 0.5
    assert r.loc["z", "n"] == 0 and r.loc["z", "rate_smoothed"] == 0.5


def test_logo_uses_only_other_groups():
    # Group A: bin x always grows. Group B: bin x never grows.
    df = pd.DataFrame({"grp": ["A", "A", "B", "B"],
                       "b": pd.Categorical(["x"] * 4), "g": [True, True, False, False]})
    p = logo_bin_predictions(df, "grp", "b", "g")
    # A's prediction comes from B only: (0 + 0.5) / (2 + 1)
    assert p[0] == pytest.approx(0.5 / 3)
    assert p[2] == pytest.approx(2.5 / 3)


def test_wilson_interval_known_value():
    lo, hi = wilson_interval(5, 10)
    assert lo == pytest.approx(0.2366, abs=1e-3) and hi == pytest.approx(0.7634, abs=1e-3)


def test_performance_perfect_ranking():
    out = performance(pd.Series([0, 0, 1, 1]), pd.Series([0.1, 0.2, 0.8, 0.9]))
    assert out["auc"] == 1.0
    assert out["brier"] == pytest.approx((0.01 + 0.04 + 0.04 + 0.01) / 4)


def test_calibration_by_bin():
    df = pd.DataFrame({"b": ["x", "x", "y"], "g": [1, 0, 1], "p": [0.4, 0.6, 0.9]})
    c = calibration_by_bin(df, "b", "g", "p")
    assert c.loc["x", "predicted"] == 0.5 and c.loc["x", "observed"] == 0.5


def test_bootstrap_groups_constant_stat():
    df = pd.DataFrame({"grp": ["A", "B", "C"], "v": [1.0, 1.0, 1.0]})
    vals = bootstrap_groups(df, "grp", lambda d: d["v"].mean(), n_boot=20, seed=1)
    assert np.all(vals == 1.0)
    assert interval(np.array([1.0, 2.0, 3.0]), 0.5) == (1.5, 2.5)


def test_bootstrap_relabels_duplicates():
    df = pd.DataFrame({"grp": ["A", "B"], "v": [1.0, 2.0]})
    vals = bootstrap_groups(df, "grp", lambda d: d["grp"].nunique(), n_boot=10, seed=0)
    assert np.all(vals == 2)


def test_weighted_recovery():
    rel = pd.Series([0.5, 0.3, 0.2])
    assert weighted_recovery(rel, pd.Series([1, 0, 1])) == pytest.approx(0.7)
    assert weighted_recovery(rel, pd.Series([0.5, 0.5, 0.5])) == pytest.approx(0.5)
