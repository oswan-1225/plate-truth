import pandas as pd
import pytest

from plate_truth.database import write_tables
from plate_truth.queries import connect, environment_counts, plate_counts, plating_conditions


@pytest.fixture
def con(tmp_path):
    db = tmp_path / "t.sqlite"
    samples = pd.DataFrame({
        "run_accession": ["e1", "e2", "p1", "p2", "p3"],
        "study_id": ["s"] * 5,
        "kind": ["environment", "environment", "plate", "plate", "plate"],
        "source_sample": ["A", "A", "A", "A", "B"],
        "condition_id": [None, None, "c1", "c1", "c2"]})
    abund = pd.DataFrame({"run_accession": ["e1", "e2", "e1", "p1", "p2", "p1", "p3"],
                          "asv_id": ["x", "x", "chl", "x", "x", "y", "y"],
                          "count": [3, 4, 50, 10, 5, 2, 9]})
    tax = pd.DataFrame({"asv_id": ["x", "y", "chl"], "is_organelle": [0, 0, 1]})
    conds = pd.DataFrame({"condition_id": ["c1", "c2"], "study_id": ["s", "s"],
                          "source_sample": ["A", "B"], "medium": ["TSA", "TSA"]})
    write_tables(db, {"samples": samples, "abundances": abund, "taxonomy": tax,
                      "conditions": conds, "studies": pd.DataFrame({"study_id": ["s"]})})
    c = connect(db)
    yield c
    c.close()


def test_environment_counts_pool_runs_and_drop_organelles(con):
    env = environment_counts(con)
    assert env.to_dict("records") == [{"study_id": "s", "source_sample": "A", "asv_id": "x",
                                       "count": 7}]
    assert len(environment_counts(con, prokaryotes_only=False)) == 2


def test_plate_counts_pool_replicate_plates(con):
    p = plate_counts(con).set_index(["condition_id", "asv_id"])["count"]
    assert p[("c1", "x")] == 15 and p[("c1", "y")] == 2 and p[("c2", "y")] == 9


def test_plating_conditions_linked_only(con):
    assert plating_conditions(con)["condition_id"].tolist() == ["c1"]   # B has no environment
    assert len(plating_conditions(con, linked_only=False)) == 2
