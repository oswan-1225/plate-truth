import numpy as np
import pandas as pd
import pytest

from plate_truth.app_model import (
    basis_from_grid,
    community_summary,
    family_composition,
    phylum_effects,
    prepare_table,
)
from plate_truth.phase5 import recovery_shares, taxon_probabilities


def test_basis_from_grid_interpolates_and_clips():
    grid_id = np.array([75.0, 80.0, 100.0])
    grid_b = np.array([[0.0, 1.0], [1.0, 1.0], [5.0, 0.0]])
    out = basis_from_grid(np.array([77.5, 90.0, 60.0, 101.0]), grid_id, grid_b)
    assert out.tolist() == [[0.5, 1.0], [3.0, 0.5], [0.0, 1.0], [5.0, 0.0]]


def test_prepare_table_filters_and_normalises():
    df = pd.DataFrame({"ASV_ID": ["a", "b", "c", "d"], "Count": [6, 2, 1, 2],
                       "identity": [0.99, 0.70, 0.98, None]})
    t, notes = prepare_table(df, present_min_reads=2, no_hit_below=75.0)
    assert t["asv_id"].tolist() == ["a", "b", "d"]          # c has 1 read
    assert t["a"].tolist() == pytest.approx([6 / 11, 2 / 11, 2 / 11])  # of all reads
    assert t["identity"].iloc[0] == pytest.approx(99.0)      # fraction read as percent
    assert t["no_hit"].tolist() == [0.0, 1.0, 1.0]           # 70% -> no hit; empty -> no hit
    assert t["identity_model"].tolist() == pytest.approx([99.0, 75.0, 75.0])
    assert any("fractions" in n for n in notes) and any("no hit" in n for n in notes)
    with pytest.raises(ValueError, match="count"):
        prepare_table(pd.DataFrame({"asv_id": ["a"]}), 2, 75.0)


def test_phylum_effects():
    eff = pd.Series({"Bacillota": 1.0, "other": -2.0})
    phyla = pd.Series(["Bacillota", "Nitrospirota", None])
    assert phylum_effects(phyla, eff, ["Bacillota"], True).tolist() == [1.0, -2.0, -2.0]
    assert phylum_effects(phyla, eff, ["Bacillota"], False).tolist() == [0.0, 0.0, 0.0]


def test_taxon_probabilities_match_recovery_shares():
    eta = np.array([-1.0, 0.5, 2.0])
    a = np.array([0.2, 0.3, 0.5])
    pts = np.array([-0.5, 0.0, 0.7])
    p1, pk = taxon_probabilities(eta, pts, 0.3, 0.4, 8, 10)
    r1, rk = recovery_shares(eta, a, pts, 0.3, 0.4, 8, 10)
    assert a @ p1 == pytest.approx(r1)
    assert a @ pk == pytest.approx(rk)


def test_community_summary_and_composition():
    pred = pd.DataFrame({"a": [0.5, 0.25, 0.25], "identity": [99.0, 90.0, np.nan],
                         "family": ["F1", "F2", None], "phylum": ["P", "P", "Q"],
                         "predicted_one_plate_p10": [0.1, 0.0, 0.0],
                         "predicted_one_plate_p50": [0.2, 0.0, 0.0],
                         "predicted_one_plate_p90": [0.4, 0.2, 0.0],
                         "predicted_any_of_8_p90": [0.8, 0.4, 0.0]})
    s = community_summary(pred, 8, 94.5)
    assert s["shortcut"] == pytest.approx(0.5)
    assert s["one_plate_p90"] == pytest.approx(0.25)
    assert s["any_of_8_p90"] == pytest.approx(0.5)
    comp = family_composition(pred, pd.Series({"F1": "under",
                                               "unassigned family (Q)": "over"}))
    assert comp.to_dict() == {"under": 0.5, "neither": 0.0, "over": 0.25,
                              "not evaluated": 0.25}
