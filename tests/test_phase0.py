import numpy as np
import pandas as pd
import pytest

from plate_truth.phase0 import attach_identity, growth_inputs, label_growth, recovery_table


def _sheet():
    return pd.DataFrame({
        "run_accession": ["e1", "e2", "p1", "p2"],
        "kind": ["environment", "environment", "plate", "plate"],
        "source_sample": ["S"] * 4,
        "condition_id": [None, None, "S|m|ae", "S|m|ae"],
        "medium": [None, None, "m", "m"], "oxygen": [None, None, "ae", "ae"]})


def test_growth_inputs_and_labels():
    counts = pd.DataFrame({"run": ["e1", "e2", "e1", "p1", "p2", "p1"],
                           "asv_id": ["a", "a", "b", "a", "a", "b"],
                           "count": [30, 30, 40, 5, 5, 90]})
    env, plates, conds = growth_inputs(counts, _sheet(), "environment")
    e = env.set_index("asv_id")
    assert e.loc["a", "reads"] == 60 and e.loc["a", "rel"] == pytest.approx(0.6)
    t = label_growth(env, plates, conds, env_min_reads=2, grew_min_reads=5, min_fold=2.0,
                     detected_min_reads=2).set_index("asv_id")
    # a: 10% on plates vs 60% in soil -> detected, not enriched; b: 90% vs 40% -> fold 2.25
    assert bool(t.loc["a", "detected_on_plates"]) and not t.loc["a", "enriched_on_plates"]
    assert bool(t.loc["b", "enriched_on_plates"])


def test_growth_inputs_excludes_asvs_before_normalising():
    counts = pd.DataFrame({"run": ["e1", "e1"], "asv_id": ["a", "c"], "count": [50, 50]})
    env, _, _ = growth_inputs(counts, _sheet(), "environment", exclude_asvs={"c"})
    assert env.set_index("asv_id").loc["a", "rel"] == 1.0


def test_attach_identity_bins():
    t = pd.DataFrame({"asv_id": ["a", "b"]})
    nr = pd.DataFrame({"asv_id": ["a", "b"], "reference_set": ["all", "all"],
                       "identity_id1": [98.0, np.nan]})
    out = attach_identity(t, nr, "all", "id1", [90.0, 94.5, 97.0, 99.0])
    assert list(out["identity_bin"]) == ["97-99", "none"]


def test_recovery_table():
    d = pd.DataFrame({"source_sample": ["S"] * 3, "env_rel": [0.5, 0.3, 0.2],
                      "y": [True, False, False], "identity": [99.0, 95.0, 80.0],
                      "pred_bin": [0.6, 0.2, 0.0]})
    r = recovery_table(d, 94.5).loc["S"]
    assert r["observed_grew"] == pytest.approx(0.5)
    assert r["proxy_has_relative"] == pytest.approx(0.8)
    assert r["predicted_calibrated"] == pytest.approx(0.36)
    assert r["proxy_over_observed"] == pytest.approx(1.6)
