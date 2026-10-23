import numpy as np
import pandas as pd
import pytest
from scipy.special import expit
from scipy.stats import norm

from plate_truth.phase5 import (
    best_reference,
    choose_runs,
    composition_shares,
    condition_weights,
    gh_nodes,
    group_share,
    habitat_of,
    has_word,
    merge_possible,
    mixture_points,
    pair_has_primers,
    primer_mismatches,
    primer_pair_shares,
    pseudo_taxon,
    recovery_shares,
    replicate_means,
    required_share,
    screening_order,
    shortcut_share,
    shortcut_verdict,
    study_effect,
)

HABITATS = {
    "soil": {"lineages": ["root:Environmental:Terrestrial:Soil"], "exclude": []},
    "marine_water": {"lineages": ["root:Environmental:Aquatic:Marine"],
                     "exclude": ["Sediment"]},
    "marine_sediment": {"lineages": ["root:Environmental:Aquatic:Marine:Sediment"],
                        "exclude": []},
}


def test_habitat_of_prefix_and_exclusion():
    assert habitat_of("root:Environmental:Terrestrial:Soil:Forest soil", HABITATS) == "soil"
    assert habitat_of("root:Environmental:Terrestrial:Soil", HABITATS) == "soil"
    assert habitat_of("root:Environmental:Terrestrial:Soilless", HABITATS) is None
    assert habitat_of("root:Environmental:Aquatic:Marine:Coastal", HABITATS) == "marine_water"
    # excluded from marine water, caught by the sediment habitat listed later
    assert habitat_of("root:Environmental:Aquatic:Marine:Sediment", HABITATS) == \
        "marine_sediment"
    assert habitat_of(None, HABITATS) is None


def test_has_word_whole_words_only():
    assert has_word("Soil isolates from Chile", ["isolates"])
    assert has_word("MOCK community", ["mock"])
    assert not has_word("Plateau soils", ["plate"])
    assert not has_word("", ["plate"])


def test_screening_order_reproducible_and_complete():
    a = screening_order(["c", "a", "b", "a"], 1)
    assert sorted(a) == ["a", "b", "c"]
    assert a == screening_order(["b", "c", "a"], 1)


def test_choose_runs_one_per_sample():
    runs = pd.DataFrame({"run_accession": ["R3", "R1", "R2", "R4"],
                         "sample_accession": ["S1", "S1", "S2", "S3"]})
    out = choose_runs(runs, 5, 0)
    assert out["run_accession"].tolist() == ["R1", "R2", "R4"]
    assert len(choose_runs(runs, 2, 0)) == 2


def test_primer_mismatches_iupac():
    assert primer_mismatches("GTGCCAGCAGCCGCGGTAA", "GTGYCAGCMGCCGCGGTAA") == 0
    assert primer_mismatches("GTGACAGCAGCCGCGGTAA", "GTGYCAGCMGCCGCGGTAA") == 1   # A not in Y
    assert primer_mismatches("GTG", "GTGY") == 1                                    # too short


def test_pair_has_primers_either_orientation_and_offset():
    fwd, rev = ["GTGYCAGCMGCCGCGGTAA"], ["GGACTACNVGGGTWTCTAAT"]
    r1 = "ACGT" + "GTGCCAGCAGCCGCGGTAA" + "TTTT"
    r2 = "GGACTACAAGGGTATCTAAT" + "CCCC"
    assert pair_has_primers(r1, r2, fwd, rev, 12, 2)
    assert pair_has_primers(r2, r1, fwd, rev, 12, 2)
    assert not pair_has_primers(r1, r1, fwd, rev, 12, 2)
    assert not pair_has_primers("A" * 13 + r1[4:], r2, fwd, rev, 12, 2)          # offset 13


def test_primer_pair_shares_and_best_reference():
    refs = {"V4": {"fwd": ["GTGYCAGCMGCCGCGGTAA"], "rev": ["GGACTACNVGGGTWTCTAAT"]},
            "V3V4": {"fwd": ["CCTACGGGNGGCWGCAG"], "rev": ["GACTACHVGGGTATCTAATCC"]}}
    good = ("GTGCCAGCAGCCGCGGTAAAA", "GGACTACAAGGGTATCTAATAA")
    bad = ("A" * 30, "C" * 30)
    s = primer_pair_shares([good, good, good, bad], refs, 12, 2)
    assert s == {"V4": 0.75, "V3V4": 0.0}
    assert best_reference(s, 0.5) == "V4"
    assert best_reference({"V4": 0.4}, 0.5) is None


def test_merge_possible():
    assert merge_possible(150, 150, 292, 8)
    assert not merge_possible(150, 150, 292, 20)
    assert merge_possible(250, 250, 465, 20)


def test_condition_weights_equal_group_totals():
    w = condition_weights(pd.Series(["a", "a", "a", "b"]))
    assert w.tolist() == pytest.approx([1 / 6, 1 / 6, 1 / 6, 1 / 2])
    assert w.sum() == pytest.approx(1.0)


def test_mixture_points_single_normal_and_point_masses():
    p = mixture_points(np.array([1.0]), np.array([1.0]), 2.0, 4)
    assert p == pytest.approx(1.0 + 2.0 * norm.ppf([0.125, 0.375, 0.625, 0.875]), abs=1e-6)
    # sigma 0: weighted quantiles of the point masses
    p0 = mixture_points(np.array([0.0, 10.0]), np.array([0.75, 0.25]), 0.0, 4)
    assert p0.tolist() == [0.0, 0.0, 0.0, 10.0]
    # two far-apart equal components: half the points near each
    p2 = mixture_points(np.array([-50.0, 50.0]), np.array([0.5, 0.5]), 1.0, 2)
    assert p2 == pytest.approx([-50.0, 50.0], abs=1e-3)


def test_gh_nodes_moments():
    x, w = gh_nodes(2.0, 20)
    assert w.sum() == pytest.approx(1.0)
    assert (w * x ** 2).sum() == pytest.approx(4.0)
    assert gh_nodes(0.0, 20)[0].tolist() == [0.0]


def test_recovery_shares_hand_checked():
    # no source effect, one condition point: p = expit(eta + u)
    eta = np.array([0.0, np.log(3.0)])            # p = 0.5 and 0.75
    a = np.array([1.0, 3.0])
    r1, r2 = recovery_shares(eta, a, np.array([0.0]), 0.0, 0.0, 2, 20)
    assert r1 == pytest.approx((0.5 + 3 * 0.75) / 4)
    assert r2 == pytest.approx((0.75 + 3 * (1 - 0.25 ** 2)) / 4)
    # study effect shifts the logit
    r1u, _ = recovery_shares(np.array([0.0]), np.array([1.0]), np.array([0.0]), np.log(3.0),
                             0.0, 1, 20)
    assert r1u == pytest.approx(0.75)
    # two condition points are averaged
    r1c, _ = recovery_shares(np.array([0.0]), np.array([1.0]), np.array([-np.log(3), np.log(3)]),
                             0.0, 0.0, 1, 20)
    assert r1c == pytest.approx(0.5)


def test_recovery_shares_source_effect_matches_quadrature():
    s, w = gh_nodes(1.5, 40)
    expect = (expit(-2.0 + s) * w).sum()
    r1, rk = recovery_shares(np.array([-2.0]), np.array([1.0]), np.array([0.0]), 0.0, 1.5, 1, 40)
    assert r1 == pytest.approx(expect)
    assert rk == pytest.approx(expect)             # any of 1 plate = one plate


def test_study_effect_quantiles():
    assert study_effect(0.5, 1.0) == pytest.approx(0.0)
    assert study_effect(0.9, 2.0) == pytest.approx(2.0 * 1.2815516, abs=1e-6)


def test_shortcut_share_no_hit_counts_below():
    s = shortcut_share(np.array([99.0, 94.5, 90.0, np.nan]), np.array([1.0, 1.0, 1.0, 1.0]), 94.5)
    assert s == 0.5


def test_pseudo_taxon_and_group_share():
    tax = pd.DataFrame({"domain": ["Bacteria", "Bacteria", "Bacteria"],
                        "phylum": ["P1", "P1", None],
                        "class": ["C1", "C1", None], "order": ["O1", "O1", None],
                        "family": ["F1", "F1", None], "genus": ["G1", None, None]})
    k = pseudo_taxon(tax, "genus")
    assert k.tolist() == ["d:Bacteria;p:P1;c:C1;o:O1;f:F1;g:G1", "d:Bacteria;p:P1;c:C1;o:O1;f:F1",
                          "d:Bacteria"]
    assert group_share(pd.Series(["x", "x", "y", "x"])).tolist() == [0.75, 0.75, 0.25, 0.75]


def test_replicate_means_single_study_uses_draws():
    v = np.array([[1.0, 2.0, 3.0]])
    assert replicate_means(v, 6, 0).tolist() == [1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
    two = replicate_means(np.array([[0.0, 0.0], [1.0, 1.0]]), 2000, 0)
    assert set(np.unique(two)) <= {0.0, 0.5, 1.0}
    assert two.mean() == pytest.approx(0.5, abs=0.05)


def test_required_share_two_thirds():
    assert [required_share(n, 2, 3) for n in (1, 2, 3, 4, 5)] == [1, 2, 2, 3, 4]


def test_shortcut_verdict():
    v = shortcut_verdict(np.array([0.3, 0.2, -0.1]), np.linspace(0.05, 0.4, 101), 2, 3)
    assert v["studies_positive"] == 2 and v["studies_needed"] == 2
    assert v["shortcut_exceeds"]
    v2 = shortcut_verdict(np.array([0.3, -0.2, -0.1]), np.linspace(0.05, 0.4, 101), 2, 3)
    assert not v2["shortcut_exceeds"]
    v3 = shortcut_verdict(np.array([0.3, 0.2]), np.linspace(-0.1, 0.4, 101), 2, 3)
    assert not v3["shortcut_exceeds"]


def test_composition_shares():
    lab = pd.Series(["under", "over", "under", "neither"])
    s = composition_shares(np.array([1.0, 1.0, 2.0, 4.0]), lab)
    assert s["under"] == 3 / 8 and s["over"] == 1 / 8 and s["neither"] == 0.5


def test_spline_basis_zero_at_left_boundary_and_matches_fit():
    from plate_truth.config import data_path, load_config
    from plate_truth.phase5 import spline_basis
    cfg = load_config("phase3")
    b = spline_basis(np.array([75.0, 100.0]), cfg)
    assert b.shape == (2, 5)
    assert b[0] == pytest.approx(np.zeros(5))
    fit = data_path("interim") / "phase3" / "fits_rev1" / "primary" / "primary" / "all"
    if not (fit / "grid_X.csv").exists():
        pytest.skip("fitted model not available")
    gx = pd.read_csv(fit / "grid_X.csv")
    grid = pd.read_csv(fit / "grid.csv")
    cols = [c for c in gx.columns if c.startswith("ns(")]
    assert spline_basis(grid["identity"].to_numpy(), cfg) == pytest.approx(gx[cols].to_numpy())
