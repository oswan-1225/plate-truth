import pandas as pd

from plate_truth.qc import (
    abundance_distortion,
    bray_curtis,
    pool_depth,
    read_retention,
    usable_conditions,
    usable_summary,
)


def test_abundance_distortion_random_and_biased_loss():
    # a and b each lose half their reads: no distortion. c loses everything it had beyond 10.
    kept = pd.Series({"a": 50, "b": 30, "c": 20})
    lost = pd.Series({"a": 50, "b": 30, "c": 80})
    d = abundance_distortion(kept, lost, min_reads=1)
    # total 260; kept 100. a: (50/100)/(100/260) = 1.3; c: (20/100)/(100/260) = 0.52
    assert round(d.loc["a", "distortion"], 3) == 1.3
    assert round(d.loc["c", "distortion"], 3) == 0.52
    assert len(abundance_distortion(kept, lost, min_reads=101)) == 0


def test_bray_curtis():
    a = pd.Series({"x": 2, "y": 2})
    assert bray_curtis(a, a * 10) == 0.0
    assert bray_curtis(a, pd.Series({"z": 5})) == 1.0
    assert bray_curtis(a, pd.Series({"x": 1, "y": 3})) == 0.25


def test_pool_depth_sums_runs():
    counts = pd.DataFrame({"study_id": ["A", "A", "A"], "condition_id": ["c1", "c1", "c2"],
                           "count": [300, 700, 50]})
    d = pool_depth(counts, ["study_id", "condition_id"])
    assert d[("A", "c1")] == 1000 and d[("A", "c2")] == 50


def _flags():
    conditions = pd.DataFrame({"study_id": ["A", "A", "A"], "source_sample": ["S", "S", "T"],
                               "condition_id": ["c1", "c2", "c3"]})
    plate = pd.Series({("A", "c1"): 1000, ("A", "c2"): 999, ("A", "c3"): 5000})
    env = pd.Series({("A", "S"): 2000})  # T has no environmental pool
    return usable_conditions(conditions, plate, env, 1000)


def test_usable_conditions_needs_both_pools():
    f = _flags().set_index("condition_id")
    # c1: both pools reach 1000; c2: plate one read short; c3: no environment pool (depth 0)
    assert f["usable"].tolist() == [True, False, False]
    assert f.loc["c3", "env_depth"] == 0


def test_read_retention():
    qc = pd.DataFrame({"study_id": ["A", "A"], "reads_in": [100, 300], "reads_window": [90, 270],
                       "reads_filtered": [80, 240], "reads_mapped": [50, 150]})
    r = read_retention(qc).iloc[0]
    assert r["runs"] == 2 and r["reads_mapped"] == 200 and r["frac_mapped"] == 0.5
    assert read_retention(qc.assign(read_count=[300, 500])).iloc[0]["read_count"] == 800


def test_usable_summary_counts_only_usable_conditions():
    rows = pd.DataFrame({"study_id": ["A"] * 4, "condition_id": ["c1", "c1", "c2", "c2"],
                         "source_sample": ["S"] * 4, "asv_id": ["a", "b", "a", "b"],
                         "enriched_on_plates": [True, False, True, True],
                         "detected_on_plates": [True, True, True, True]})
    s = usable_summary(_flags(), rows).iloc[0]
    assert (s["conditions"], s["conditions_usable"], s["sources"], s["sources_usable"]) == \
        (3, 1, 2, 1)
    # only c1 survives: 2 rows, 1 enriched, 2 detected
    assert (s["asvs_usable"], s["rows_usable"], s["enriched_rows"], s["detected_rows"]) == \
        (2, 2, 1, 2)
