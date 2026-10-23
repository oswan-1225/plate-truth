import numpy as np
import pandas as pd
import pytest
from scipy.special import expit

from plate_truth.phase3 import (
    criteria_checks,
    averaged_curve,
    bootstrap_mean,
    calibration_metrics,
    condition_recovery,
    effect_code,
    group_metrics,
    group_weights,
    keep_levels,
    known_levels,
    local_recalibration,
    marginal_probability,
    pilot_split,
    ratio_with_interval,
    recalibrate_intercept,
    recalibrate_logistic,
    recalibration_by_group,
    required_groups,
    source_recovery,
    temperature_band,
)
from plate_truth.rfit import fixed_formula, spline_term

CFG = {"identity_spline": {"boundary_knots": [75.0, 100.0],
                           "interior_knots": [90.0, 94.5, 97.0, 99.0]}}


def test_temperature_band():
    t = pd.Series([7.0, 20.0, 25.0, 33.0, 45.0, np.nan])
    assert temperature_band(t, 20.0, 33.0).tolist() == [
        "cold", "cold", "moderate", "warm", "warm", "unknown"]


def test_keep_levels_and_known_levels():
    v = keep_levels(pd.Series(["a", "b", None, "z"]), ["a", "b"], "other")
    assert v.tolist() == ["a", "b", "other", "other"]
    assert known_levels(pd.Series(["b", "unknown", "a", "b"])) == ["a", "b"]


def test_effect_code_sum_to_zero_and_unknown_is_average():
    v = pd.Series(["a", "b", "c", "unknown", "d"])
    x = effect_code(v, ["a", "b", "c"], "f")
    assert list(x.columns) == ["f_a", "f_b"]
    assert x.loc[0].tolist() == [1, 0] and x.loc[1].tolist() == [0, 1]
    assert x.loc[2].tolist() == [-1, -1]                  # last level
    assert x.loc[3].tolist() == [0, 0] and x.loc[4].tolist() == [0, 0]   # unknown, unseen
    # the three known levels' codes sum to zero, so all-zero = their average
    assert x.loc[0:2].sum().tolist() == [0, 0]
    assert effect_code(v, ["a"], "f").shape[1] == 0


def test_group_weights_equal_totals():
    g = pd.Series(["x", "x", "x", "y"])
    w = group_weights(g)
    assert w.tolist() == pytest.approx([2 / 3, 2 / 3, 2 / 3, 2.0])
    assert w.groupby(g).sum().tolist() == pytest.approx([2.0, 2.0])


def test_marginal_probability():
    assert marginal_probability(np.array([0.3]), 0.0, 20)[0] == pytest.approx(expit(0.3))
    # symmetric around 0 -> exactly 0.5
    assert marginal_probability(np.array([0.0]), 2.0, 20)[0] == pytest.approx(0.5)
    # a rare outcome: averaging over study effects raises the expected rate
    eta = np.array([-5.0])
    mc = expit(eta[0] + 1.5 * np.random.default_rng(1).standard_normal(400_000)).mean()
    assert marginal_probability(eta, 1.5, 40)[0] == pytest.approx(mc, rel=0.01)
    assert marginal_probability(eta, 1.5, 40)[0] > expit(-5.0)


def test_calibration_metrics_perfect_and_too_extreme():
    rng = np.random.default_rng(0)
    p = rng.uniform(0.01, 0.6, 50_000)
    y = rng.random(50_000) < p
    m = calibration_metrics(pd.Series(y), pd.Series(p))
    assert m["oe"] == pytest.approx(1.0, abs=0.03)
    assert m["slope"] == pytest.approx(1.0, abs=0.05)
    # predictions doubled on the logit scale are too extreme -> slope about 0.5
    p2 = expit(2 * np.log(p / (1 - p)))
    assert calibration_metrics(pd.Series(y), pd.Series(p2))["slope"] == pytest.approx(0.5, abs=0.05)
    # hand check of O/E: 2 events, predictions sum to 4
    assert calibration_metrics(pd.Series([1, 1, 0, 0]), pd.Series([1, 1, 1, 1.0]))["oe"] == \
        pytest.approx(0.5)


def test_group_metrics_auc():
    df = pd.DataFrame({"study_group": ["g"] * 4, "y": [0, 0, 1, 1],
                       "p": [0.1, 0.2, 0.3, 0.4]})
    m = group_metrics(df, "p").iloc[0]
    assert m["auc"] == 1.0 and m["events"] == 2
    assert m["brier"] == pytest.approx((0.01 + 0.04 + 0.49 + 0.36) / 4)


def test_bootstrap_mean_constant():
    assert np.all(bootstrap_mean(np.array([0.7, 0.7, np.nan]), 50, 1) == 0.7)


def test_source_and_condition_recovery():
    rows = pd.DataFrame({
        "study_id": ["s"] * 4, "study_group": ["s"] * 4, "source_key": ["s|1"] * 4,
        "condition_id": ["c1", "c1", "c2", "c2"], "asv_id": ["a", "b", "a", "b"],
        "env_rel": [0.6, 0.4, 0.6, 0.4], "identity": [99.0, 90.0, 99.0, 90.0],
        "y": [0, 1, 0, 0], "p": [0.5, 0.0, 0.5, 0.0]})
    s = source_recovery(rows, 94.5).iloc[0]
    assert s["proxy_has_relative"] == pytest.approx(0.6)
    assert s["observed_any_condition"] == pytest.approx(0.4)
    c = condition_recovery(rows, "p").set_index("condition_id")
    assert c.loc["c1", "observed"] == pytest.approx(0.4)
    assert c.loc["c2", "predicted"] == pytest.approx(0.3)


def test_ratio_with_interval():
    est, lo, hi = ratio_with_interval(np.array([0.5, 0.5]), np.array([0.1, 0.1]), 100, 0)
    assert est == pytest.approx(5.0) and lo == pytest.approx(5.0) and hi == pytest.approx(5.0)
    assert np.isnan(ratio_with_interval(np.array([0.5]), np.array([0.1]), 100, 0)[1])


def test_fixed_formula():
    cols = ["identity", "no_hit", "log_env", "host", "medium_a"]
    s = spline_term(CFG)
    assert s == "ns(identity, knots = c(90, 94.5, 97, 99), Boundary.knots = c(75, 100))"
    assert fixed_formula("primary", cols, CFG) == f"{s} + no_hit + log_env + medium_a"
    assert fixed_formula("identity_only", cols, CFG) == f"{s} + no_hit"
    assert fixed_formula("no_identity", cols, CFG) == "log_env + medium_a"
    assert fixed_formula("habitat", cols, CFG).endswith(f"host + host:{s} + host:no_hit")


def test_averaged_curve():
    # one row at 0, identity contribution log(3), no random effects -> 3/(1+3) = 0.75
    assert averaged_curve(np.array([[0.0]]), np.array([[np.log(3)]]), 0.0, 20)[0, 0] ==         pytest.approx(0.75)
    # rows at -1 and +1 average to exactly 0.5 by symmetry
    assert averaged_curve(np.array([[-1.0], [1.0]]), np.array([[0.0]]), 0.0, 20)[0, 0] ==         pytest.approx(0.5)
    # shapes: (rows, draws) and (grid, draws) -> (grid, draws)
    assert averaged_curve(np.zeros((3, 2)), np.zeros((4, 2)), 1.0, 20).shape == (4, 2)


def test_criteria_checks_hand_values():
    crit = {"calibration_slope_range": [0.8, 1.25], "oe_range": [0.67, 1.5],
            "auc_supported_groups": 8, "calibrated_groups": 7, "proxy_over_groups": 8,
            "reference_groups": 10}
    auc = np.full(10, 0.8)
    oe = np.array([1.0] * 7 + [0.3] * 3)
    slope = np.array([1.0] * 6 + [2.0] * 4)          # 6 groups in both ranges
    ratio = np.array([5.0] * 9 + [0.5])
    c = criteria_checks(auc, oe, slope, ratio, crit, 200, 1)
    assert c["groups_auc_above_0.5"] == 10 and c["H1"]
    assert c["groups_calibrated"] == 6 and not c["H2"]
    assert c["groups_ratio_above_1"] == 9 and c["median_ratio"] == 5.0 and c["H3"]
    # 7 groups: 8 of 10 scales to ceil(5.6) = 6; all 7 above 1 passes, 5 of 7 AUCs does not
    c7 = criteria_checks(np.array([0.8] * 5 + [0.45] * 2), oe[:7], slope[:7],
                         np.full(7, 5.0), crit, 200, 1)
    assert c7["groups_needed_h1"] == 6 and not c7["H1"]
    assert c7["groups_needed_h3"] == 6 and c7["H3"]


def test_required_groups_scales_the_share():
    assert required_groups(8, 10, 10) == 8      # unchanged at full size
    assert required_groups(7, 10, 10) == 7      # no float rounding up (0.7 * 10 = 7.000...01)
    assert required_groups(8, 7, 10) == 6       # 5.6 -> 6
    assert required_groups(7, 7, 10) == 5       # 4.9 -> 5
    assert required_groups(8, 9, 10) == 8       # 7.2 -> 8


def test_pilot_split_sizes():
    rng = np.random.default_rng(0)
    assert len(pilot_split(np.arange(10), 0.25, rng)) == 3      # ceil(2.5)
    assert len(pilot_split(np.arange(4), 0.1, rng)) == 1        # at least one
    assert len(pilot_split(np.arange(4), 1.0, rng)) == 3        # at least one left to score


def test_recalibrate_intercept_matches_pilot_events():
    # pilot: 10 rows at p = 0.1 (E = 1), 2 events: the shift makes each 0.2 (E = 2)
    y, p = np.array([1, 1] + [0] * 8), np.full(10, 0.1)
    assert recalibrate_intercept(y, p, np.array([0.1]), 0.5)[0] == pytest.approx(0.2)
    # no events: target 0.5 events, so each row becomes 0.05, not 0
    assert recalibrate_intercept(np.zeros(10), p, np.array([0.1]), 0.5)[0] == \
        pytest.approx(0.05)


def test_recalibrate_logistic_reproduces_two_level_rates():
    # two prediction levels with observed rates 1/4 and 2/4: a two-parameter fit is saturated
    p = np.array([expit(-2.0)] * 4 + [0.5] * 4)
    y = np.array([1, 0, 0, 0, 1, 1, 0, 0])
    new = recalibrate_logistic(y, p, np.array([expit(-2.0), 0.5]))
    assert new == pytest.approx([0.25, 0.5], abs=1e-6)


def test_local_recalibration_scores_only_held_out_conditions():
    rng = np.random.default_rng(1)
    n = 4000
    df = pd.DataFrame({"study_group": "g", "condition_id": rng.integers(0, 8, n).astype(str),
                       "pred": np.full(n, 0.05)})
    df["y"] = (rng.random(n) < 0.2).astype(int)                # true rate 4x the prediction
    s = local_recalibration(df, "pred", 0.25, 5, 10, 0.5, 2)
    assert set(s["method"]) == {"none", "intercept", "logistic"} and len(s) == 15
    assert (s["pilot_conditions"] == 2).all()                   # ceil(0.25 * 8)
    none = s[s["method"] == "none"]["oe"]
    shifted = s[s["method"] == "intercept"]["oe"]
    assert (none > 3).all() and (abs(shifted - 1) < 0.2).all()  # the shift fixes the level


def test_recalibration_by_group_medians():
    splits = pd.DataFrame({"study_group": "g", "method": "none", "split": [0, 1, 2],
                           "oe": [0.5, 1.0, 1.2], "slope": [1.0, 1.0, np.nan]})
    r = recalibration_by_group(splits, [0.8, 1.25], [0.67, 1.5]).iloc[0]
    assert r["median_oe"] == 1.0 and r["median_slope"] == 1.0 and r["splits_scored"] == 3
    assert r["share_splits_in_range"] == pytest.approx(1 / 3) and r["calibrated"]
