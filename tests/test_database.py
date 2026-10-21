import sqlite3

import numpy as np
import pandas as pd

from plate_truth.database import condition_ids, conditions_table, is_organelle, write_tables


def sheet():
    return pd.DataFrame({
        "study_id": ["s", "s", "s", "s"],
        "run_accession": ["r1", "r2", "r3", "r4"],
        "source_sample": ["A", "A", "A", "A"],
        "kind": ["plate", "plate", "plate", "environment"],
        "medium": ["TSA", "TSA", "R2A", np.nan],
        "oxygen": ["aerobic", "aerobic", np.nan, np.nan],
        "temperature_c": [30.0, 30.0, 7.0, np.nan],
        "incubation_days": [7, 7, 7, 7],
        "condition_detail": [np.nan, np.nan, "reference plate", np.nan],
    })


def test_condition_ids_pool_replicate_plates_and_skip_environment():
    cid = condition_ids(sheet())
    assert cid[0] == cid[1] == "s|A|TSA|aerobic|30|NA"     # replicate plates share an id
    assert cid[2] == "s|A|R2A|NA|7|reference plate"
    assert pd.isna(cid[3])


def test_conditions_table_counts_plates():
    s = sheet()
    s["condition_id"] = condition_ids(s)
    t = conditions_table(s).set_index("condition_id")
    assert t.loc["s|A|TSA|aerobic|30|NA", "n_plates"] == 2
    assert len(t) == 2


def test_is_organelle():
    tax = pd.DataFrame({"domain": ["Bacteria", "Bacteria", "Bacteria", "Eukaryota", np.nan],
                        "order": ["Chloroplast", "Rickettsiales", "Bacillales", np.nan, np.nan],
                        "family": [np.nan, "Mitochondria", "Bacillaceae", np.nan, np.nan]})
    assert is_organelle(tax).tolist() == [True, True, False, True, True]


def test_write_tables_keeps_schema_columns_only(tmp_path):
    db = tmp_path / "t.sqlite"
    write_tables(db, {"medium_class": pd.DataFrame({"medium_class": ["x"], "definition": ["d"],
                                                    "extra": [1]})})
    with sqlite3.connect(db) as con:
        assert con.execute("SELECT medium_class, definition FROM medium_class").fetchall() == [
            ("x", "d")]
