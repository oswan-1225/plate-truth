import numpy as np
import pandas as pd
import pytest

from plate_truth.phase4 import (
    asv_strains,
    boundary_lrt_p,
    classify_collection_media,
    classify_taxa,
    combine,
    culture_media_records,
    discordant_pairs,
    genome_base,
    has_organic_carbon,
    ingredient_roles,
    mantel_haenszel_or,
    match_term,
    matched_condition_pairs,
    mean_log_loss,
    medium_facts,
    member_strains,
    pair_counts,
    recorded_by_study,
    recorded_shares,
    reportable_taxa,
    strain_classes,
    taxon_label,
    taxon_oe,
)

RULES = {"antibacterial": ["Ampicillin"], "habitat_material": ["Soil extract"],
         "not_complex": ["Agar", "Starch"], "organic_excluded": ["Agar", "HEPES"]}


def test_culture_media_records_parses_mediadive_links():
    rec = {"Culture and growth conditions": {"culture medium": [
        {"name": "NUTRIENT AGAR (DSMZ Medium 1)", "growth": "yes",
         "link": "https://mediadive.dsmz.de/medium/1"},
        {"name": "CIP Medium 6", "growth": "no",
         "link": "https://catalogue-crbip.pasteur.fr/fiche_milieu.xhtml?crbip=6"}]}}
    out = culture_media_records(7, rec)
    assert [r["mediadive_id"] for r in out] == ["1", None]
    assert [r["growth"] for r in out] == ["yes", "no"]
    single = {"Culture and growth conditions": {"culture medium": {
        "name": "x", "growth": "yes", "link": "https://mediadive.dsmz.de/medium/514a"}}}
    assert culture_media_records(1, single)[0]["mediadive_id"] == "514a"


@pytest.mark.parametrize("formula,expected", [
    ("C6H12O6", True), ("NaCl", False), ("CaCl2", False), ("NaHCO3", False),
    ("Na2CO3", False), ("CO2", False), ("KCN", False), ("C2H3NaO2", True), (None, False),
    ("CuSO4", False)])
def test_has_organic_carbon(formula, expected):
    assert has_organic_carbon(formula) is expected


def ingredients():
    return pd.DataFrame({
        "id": [1, 2, 3, 4, 5, 6, 7, 8, 9],
        "name": ["Peptone", "Agar", "Glucose", "NaCl", "Ampicillin", "Soil extract",
                 "R2A agar", "HEPES", "Starch"],
        "complex_compound": [1, 1, 0, 0, 0, 1, 1, 0, 1],
        "formula": [None, "(C12H18O9)n", "C6H12O6", "NaCl", "C16H19N3O4S", None, None,
                    "C8H18N2O4S", None]})


def test_ingredient_roles():
    r = ingredient_roles(ingredients(), RULES).set_index("name")
    assert r.loc["Peptone", "complex_source"] and r.loc["Peptone", "organic"]
    assert not r.loc["Agar", "complex_source"] and not r.loc["Agar", "organic"]
    assert r.loc["Glucose", "organic"] and not r.loc["Glucose", "complex_source"]
    assert not r.loc["NaCl", "organic"]
    assert r.loc["Ampicillin", "antibacterial"]
    assert r.loc["Soil extract", "habitat_material"]
    assert r.loc["R2A agar", "composite_agar"] and not r.loc["Peptone", "composite_agar"]
    assert not r.loc["HEPES", "organic"]
    # starch: flagged complex by MediaDive, not a complex source, but organic matter
    assert r.loc["Starch", "organic"] and not r.loc["Starch", "complex_source"]


def test_medium_facts_and_classes():
    roles = ingredient_roles(ingredients(), RULES)
    comp = pd.DataFrame({
        # m1 nutrient-agar-like: peptone 5 + glucose 3 = 8 g/L organic, agar ignored -> rich
        # m2 dilute: peptone 1 + glucose 1 = 2 g/L -> dilute
        # m3 defined: glucose 10 -> defined_minimal
        # m4 R2A agar 18 g/L as one ingredient: 18 - 15 = 3 g/L -> dilute
        # m5 peptone with unknown amount, glucose 1 -> unresolved; m6 same + glucose 5 -> rich
        # m7 ampicillin -> selective; m8 soil extract -> environment_extract
        # m9 optional peptone only -> defined_minimal (optional ingredients ignored)
        "medium_id": ["m1", "m1", "m1", "m2", "m2", "m3", "m4", "m5", "m5", "m6", "m6", "m7",
                      "m7", "m8", "m9", "m9"],
        "ingredient_id": [1, 3, 2, 1, 3, 3, 7, 1, 3, 1, 3, 1, 5, 6, 1, 3],
        "g_l": [5, 3, 15, 1, 1, 10, 18, None, 1, None, 5, 5, 0.1, None, 5, 1],
        "optional": [0] * 14 + [1, 0]})
    facts = medium_facts(comp, roles, 15.0).set_index("medium_id")
    assert facts.loc["m1", "organic_known_g"] == pytest.approx(8.0)
    assert facts.loc["m4", "organic_known_g"] == pytest.approx(3.0)
    assert facts.loc["m5", "organic_unknown"]
    cls = classify_collection_media(facts.reset_index()).set_index("medium_id")["medium_class"]
    assert cls.to_dict() == {"m1": "rich_complex", "m2": "dilute_oligotrophic",
                             "m3": "defined_minimal", "m4": "dilute_oligotrophic",
                             "m5": "unresolved", "m6": "rich_complex", "m7": "selective",
                             "m8": "environment_extract", "m9": "defined_minimal"}


def test_strain_classes_keeps_positive_classed_records():
    rec = pd.DataFrame({"bacdive_id": [1, 1, 2, 3], "growth": ["yes", "no", None, "yes"],
                        "mediadive_id": ["1", "2", "1", None]})
    mc = pd.DataFrame({"medium_id": ["1", "2"], "medium_class": ["rich_complex",
                                                                 "dilute_oligotrophic"]})
    out = strain_classes(rec, mc, ["rich_complex", "dilute_oligotrophic"])
    assert sorted(map(tuple, out.to_numpy())) == [(1, "rich_complex"), (2, "rich_complex")]


def test_genome_base_and_member_strains():
    assert genome_base("RS_GCF_000005845.2") == "GCA_000005845"
    assert genome_base("GCA_016900935") == "GCA_016900935"
    members = pd.DataFrame({
        "reference": ["V4"] * 4, "ref_id": ["r1", "r2", "r3", "r3"],
        "member_id": ["bacdive_strain|9", "gtdb_isolate|RS_GCF_000005845.2|1",
                      "refseq_type|NR_1", "gtdb_isolate|GB_GCA_999.1|1"]})
    genomes = pd.DataFrame({"bacdive_id": [50, 51], "assembly": ["GCA_000005845", "GCA_1"]})
    out = member_strains(members, genomes)
    assert sorted(map(tuple, out.to_numpy())) == [("V4", "r1", 9), ("V4", "r2", 50)]


def test_asv_strains_through_tied_hits():
    nr = pd.DataFrame({"asv_id": ["a", "b"], "reference": ["V4", "V4"],
                       "targets": ["r1;r2", None]})
    links = pd.DataFrame({"reference": ["V4", "V4", "V4"], "ref_id": ["r1", "r2", "r2"],
                          "bacdive_id": [9, 9, 50]})
    out = asv_strains(nr, links)
    assert sorted(out["bacdive_id"]) == [9, 50] and set(out["asv_id"]) == {"a"}


def test_recorded_shares_and_match_term():
    links = pd.DataFrame({"asv_id": ["a", "a", "a", "b"], "bacdive_id": [1, 2, 3, 4]})
    sc = pd.DataFrame({"bacdive_id": [1, 1, 2, 4], "medium_class": [
        "rich_complex", "dilute_oligotrophic", "rich_complex", "dilute_oligotrophic"]})
    rec = recorded_shares(links, sc, ["rich_complex", "dilute_oligotrophic"]).set_index("asv_id")
    # strain 3 has no classed record, so a has 2 strains: both rich, one dilute
    assert rec.loc["a", "n_strains"] == 2
    assert rec.loc["a", "rec_rich_complex"] == 1.0
    assert rec.loc["a", "rec_dilute_oligotrophic"] == 0.5
    assert rec.loc["b", "rec_rich_complex"] == 0.0 and rec.loc["b", "any_dilute_oligotrophic"]
    rows = pd.DataFrame({"asv_id": ["a", "a", "a"],
                         "medium": ["rich_complex", "dilute_oligotrophic", "selective"]})
    rows = rows.merge(rec.reset_index(), on="asv_id")
    assert match_term(rows, ["rich_complex", "dilute_oligotrophic"]).tolist() == [1.0, 0.5, 0.0]


def test_taxon_label_and_combine():
    tax = pd.DataFrame({"phylum": ["Bacillota", "Bacillota", None],
                        "family": ["Bacillaceae", None, None]})
    assert taxon_label(tax, "family").tolist() == [
        "Bacillaceae", "unassigned family (Bacillota)", "unassigned"]
    assert combine(pd.Series(["f"]), pd.Series(["s"])).tolist() == ["f:s"]


def test_boundary_lrt_p():
    # stat = 2 * (ll1 - ll0) = 3.84 -> chi2_1 tail 0.05 -> mixture p 0.025
    assert boundary_lrt_p(-100.0, -98.08) == pytest.approx(0.025, abs=2e-4)
    assert boundary_lrt_p(-100.0, -100.0) == 1.0


def test_mean_log_loss():
    assert mean_log_loss([1, 0], [0.5, 0.5]) == pytest.approx(np.log(2))
    assert mean_log_loss([1], [1.0]) < 1e-9


def test_reportable_taxa_uses_row_counts_only():
    rows = pd.DataFrame({"fam": ["A"] * 3 + ["B"] * 3, "y": [0, 0, 0, 1, 1, 1],
                         "medium": ["x", "y", "y", "x", "x", "x"]})
    assert reportable_taxa(rows, "fam", 3, 2) == ["A"]


def test_matched_pairs_and_discordant_cells():
    conds = pd.DataFrame({
        "condition_id": ["c1", "c2", "c3", "c4"], "medium": ["rich", "dilute", "rich", "dilute"],
        "study_id": ["s"] * 4, "source_key": ["s|1"] * 4,
        "oxygen": ["aerobic", "aerobic", "anaerobic", "aerobic"],
        "temperature": ["moderate"] * 4, "incubation_days": [7.0, 7.0, 7.0, 14.0]})
    pairs = matched_condition_pairs(conds)
    # only c1-c2 match on every key; orientation: dilute < rich alphabetically
    assert pairs[["condition_id_x", "condition_id_y", "class_x", "class_y"]].values.tolist() \
        == [["c2", "c1", "dilute", "rich"]]
    rows = pd.DataFrame({
        "condition_id": ["c1", "c2"] * 4, "asv_id": ["a", "a", "b", "b", "c", "c", "d", "d"],
        "study_group": ["g"] * 8, "y": [1, 0, 0, 1, 1, 1, 1, 0]})
    recorded = pd.DataFrame({"asv_id": ["a", "b", "c", "d"],
                             "any_rich": [1.0, 1.0, 1.0, 1.0],
                             "any_dilute": [0.0, 0.0, 0.0, 1.0]})
    disc = discordant_pairs(rows, pairs, recorded).set_index("asv_id")
    # a: recorded rich (y side) only, grew on rich (y) -> d; b: recorded y, grew x -> c;
    # c: not discordant; d: recorded on both -> dropped
    assert disc["cell"].to_dict() == {"a": "d", "b": "c"}


def test_mantel_haenszel_or():
    # one stratum a=10 b=5 c=5 d=10: OR = (10*10)/(5*5) = 4
    t = pd.DataFrame({"a": [10.0], "b": [5.0], "c": [5.0], "d": [10.0]})
    assert mantel_haenszel_or(t) == pytest.approx(4.0)
    # two identical strata keep OR 4; a stratum with no b or c gives no denominator
    assert mantel_haenszel_or(pd.concat([t, t])) == pytest.approx(4.0)


def test_pair_counts_with_weights():
    disc = pd.DataFrame({"study_id": ["s", "s", "s"], "study_group": ["g"] * 3,
                         "asv_id": ["a", "a", "b"], "stratum": ["x|y"] * 3,
                         "cell": ["a", "d", "b"]})
    t = pair_counts(disc).iloc[0]
    assert (t["a"], t["b"], t["c"], t["d"]) == (1.0, 1.0, 0.0, 1.0)
    tw = pair_counts(disc, pd.Series({"a": 2.0, "b": 0.0})).iloc[0]
    assert (tw["a"], tw["b"], tw["d"]) == (2.0, 0.0, 2.0)


def test_taxon_oe_and_classify_taxa():
    df = pd.DataFrame({"fam": ["A"] * 6, "study_group": ["g1", "g1", "g2", "g2", "g3", "g3"],
                       "y": [0, 0, 0, 1, 0, 0], "p": [0.5] * 6})
    t = taxon_oe(df, "fam", "p")
    assert t["expected"].tolist() == [1.0, 1.0, 1.0] and t["observed"].tolist() == [0, 1, 0]
    # A: O = 1, E = 3 in 3 groups, two of them below 1; every bootstrap ratio is at most 1
    big = t.assign(observed=t["observed"] * 10, expected=t["expected"] * 10)
    out = classify_taxa(big, "fam", 5.0, 3, 0.6667, 200, 1).iloc[0]
    assert out["oe"] == pytest.approx(1 / 3)
    assert out["share_groups_below_1"] == pytest.approx(2 / 3)
    # the bootstrap can draw g2 three times (O/E 1), so the upper bound reaches 1: not 'under'
    assert out["oe_hi"] == pytest.approx(1.0) and out["call"] == "neither"
    under = big.assign(observed=[0, 5, 0])
    assert classify_taxa(under, "fam", 5.0, 3, 0.6667, 200, 1).iloc[0]["call"] == "under"
    # fewer than min_groups groups with enough expected events: not evaluated
    assert classify_taxa(t, "fam", 5.0, 3, 0.6667, 200, 1).empty


def test_recorded_by_study():
    rows = pd.DataFrame({"study_id": ["s1", "s2", "s1"], "rec_rich": [0.5, 1.0, 0.0]})
    out = recorded_by_study(rows, ["rich"])
    assert out["rec_rich__s1"].tolist() == [0.5, 0.0, 0.0]
    assert out["rec_rich__s2"].tolist() == [0.0, 1.0, 0.0]
