import numpy as np
import pandas as pd

from plate_truth.growth import (
    any_condition,
    asv_condition_table,
    fold_enrichment,
    is_detected,
    is_enriched,
    plate_only_asvs,
    pooled_abundance,
)


def test_pooled_abundance_sums_replicates():
    counts = pd.DataFrame({"g": ["s", "s", "s"], "asv_id": ["a", "a", "b"], "count": [3, 1, 4]})
    out = pooled_abundance(counts, "g").set_index("asv_id")
    assert out.loc["a", "reads"] == 4 and out.loc["a", "rel"] == 0.5


def test_pooled_abundance_copy_number_weight():
    counts = pd.DataFrame({"g": ["s", "s"], "asv_id": ["a", "b"], "count": [40, 10]})
    # a has 4 copies per genome, b has 1: 10 'cells' each, so 50/50 after correction
    out = pooled_abundance(counts, "g", weight=pd.Series({"a": 4.0, "b": 1.0})).set_index("asv_id")
    assert out.loc["a", "rel"] == 0.5 and out.loc["a", "reads"] == 40


def test_fold_and_definitions():
    fold = fold_enrichment(pd.Series([0.02, 0.01, 0.0]), pd.Series([0.01, 0.01, 0.0]))
    assert fold[0] == 2.0 and fold[1] == 1.0 and np.isnan(fold[2])
    reads = pd.Series([5, 100, 4])
    f = pd.Series([2.0, 1.9, 10.0])
    assert is_enriched(reads, f, 5, 2.0).tolist() == [True, False, False]
    assert is_detected(pd.Series([1, 2]), 2).tolist() == [False, True]


def _toy():
    env = pd.DataFrame({"source_sample": ["S"] * 3, "asv_id": ["a", "b", "c"],
                        "reads": [50, 49, 1], "rel": [0.5, 0.49, 0.01]})
    plates = pd.DataFrame({"condition_id": ["S|m|ae"] * 3, "asv_id": ["a", "c", "d"],
                           "reads": [10, 80, 10], "rel": [0.1, 0.8, 0.1]})
    conds = pd.DataFrame({"condition_id": ["S|m|ae"], "source_sample": ["S"],
                          "medium": ["m"], "oxygen": ["ae"]})
    return env, plates, conds


def test_asv_condition_table_restricts_and_fills():
    env, plates, conds = _toy()
    t = asv_condition_table(env, plates, conds, env_min_reads=2).set_index("asv_id")
    assert set(t.index) == {"a", "b"}  # c has 1 env read: not present; d not in env
    assert t.loc["b", "plate_reads"] == 0
    assert t.loc["a", "fold"] == 0.1 / 0.5


def test_plate_only_asvs():
    env, plates, conds = _toy()
    out = plate_only_asvs(env, plates, conds, env_min_reads=2, plate_min_reads=5)
    assert sorted(out["asv_id"]) == ["c", "d"]  # c below env presence, d absent


def test_any_condition():
    t = pd.DataFrame({"source_sample": ["S", "S"], "asv_id": ["a", "a"],
                      "condition_id": ["x", "y"], "env_reads": [5, 5], "env_rel": [0.1, 0.1],
                      "g": [False, True]})
    out = any_condition(t, "g")
    assert out["g"].tolist() == [True] and out["n_conditions"].tolist() == [2]
