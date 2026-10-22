import pandas as pd

from plate_truth.config import load_config
from plate_truth.samples import (
    PARSERS, condition_id, link_summary, kruger_compartment, parse_hinsu2021,
    parse_yu2020, parse_verhoeven2023, parse_villanelo2026, parse_ordonez2025,
    parse_ortizseverin2024, parse_indraningrat2019, parse_kruger2019, parse_li2023,
    parse_zehavi2018, parse_zheng2024, parse_zheng_library, li2023_incubation_days)


def test_parse_zheng_library_kinds():
    names = pd.Series(["TA_3_ORI_ORI", "TA_7_FIL_FIL", "TA_12_AN_TSA", "GD_1_MI_inSITU",
                       "XM2_AE_31_FIL_FIL"])
    out = parse_zheng_library(names)
    assert out["kind"].tolist()[:4] == ["environment", "inoculum", "plate", "plate"]
    assert out.loc[2, "medium"] == "1% TSA" and out.loc[2, "oxygen"] == "anaerobic"
    assert out.loc[3, "medium"] == "soil extract" and out.loc[3, "oxygen"] == "microaerophilic"
    assert pd.isna(out.loc[0, "medium"])  # environment rows carry no medium
    assert out.loc[2, "replicate"] == 12
    assert out.iloc[4].isna().all()  # sludge-style name does not fit the soil pattern


def test_condition_id_only_for_plates():
    df = parse_zheng_library(pd.Series(["SM_1_AE_TSA", "SM_2_ORI_ORI"]))
    cid = condition_id(df)
    assert cid[0] == "SM|1% TSA|aerobic"
    assert pd.isna(cid[1])


# --- Phase 2 parsers: each test uses names copied from the real ENA records ---------------------

def meta(**cols):
    """A small ENA-like metadata frame; missing columns are filled with None."""
    n = len(next(iter(cols.values())))
    base = {c: [None] * n for c in ("library_name", "sample_title", "sample_alias",
                                     "description", "isolation_source")}
    return pd.DataFrame({**base, **cols})


def test_zheng2024_sludge_merges_tank_labels_and_reads_r2a():
    out = parse_zheng2024(meta(library_name=[
        "XM2_AN_4_AN_R2A", "XM2_AE_26_ORI_ORI", "SZ_24_AN_R2A", "EBPR_14_AE_R2A", "TA_3_ORI_ORI",
        "XM4_2_AE_Vb12"]), "sludge")
    assert out.loc[0, ["source_sample", "kind", "medium", "oxygen"]].tolist() == [
        "XM2", "plate", "R2A", "anaerobic"]
    assert out.loc[0, "condition_detail"] == "AN" and out.loc[0, "replicate"] == 4
    assert out.loc[1, ["source_sample", "kind"]].tolist() == ["XM2", "environment"]
    assert out.loc[2, "medium"] == "R2A"                    # a digit inside the medium code
    assert out.loc[3:5, "kind"].tolist() == ["other"] * 3   # EBPR, soil, unlisted XM4


def test_zheng2024_soil_matches_phase0_parser():
    out = parse_zheng2024(meta(library_name=["TA_12_AN_TSA", "GD_1_MI_inSITU", "SZ_24_AN_R2A"]),
                          "soil")
    assert out.loc[0, ["medium", "oxygen"]].tolist() == ["1% TSA", "anaerobic"]
    assert out.loc[1, "medium"] == "soil extract"
    assert out.loc[2, "kind"] == "other"


def test_kruger_compartments_and_sample_types():
    assert kruger_compartment("Biopore soil") == "biopore"
    assert kruger_compartment("matrix 80-120 cmbs extracted soil community") == "matrix_80_120"
    out = parse_kruger2019(meta(
        description=["type of sample:soil sample_plow layer", "type of sample:Nycodenz extraction",
                     "type of sample:reference plate", "type of sample:PSM_incubation:NH",
                     "type of sample:PSM_incubation:15mm"],
        isolation_source=["plow layer soil", "extracted plow layer soil community",
                          "extracted red fracture soil community grown on solid medium",
                          "biopore soil community grown onto a hydrated surface",
                          "biopore soil community grown onto a hydrated surface"]))
    assert out["kind"].tolist() == ["environment", "inoculum", "plate", "plate", "other"]
    assert out["source_sample"].tolist()[:4] == ["plow_layer", "plow_layer", "red_fracture",
                                                  "biopore"]
    assert out.loc[3, "condition_detail"] == "PSM full plate"


def test_li2023_temperature_and_medium_from_name():
    out = parse_li2023(meta(library_name=["BS", "RH41", "BL32", "XX"]))
    assert out.loc[0, ["source_sample", "kind"]].tolist() == ["bulk", "environment"]
    assert out.loc[1, ["source_sample", "medium", "temperature_c"]].tolist() == [
        "rhizosphere", "M41 (TP)", 45]
    assert out.loc[2, "temperature_c"] == 15
    assert out.loc[3, "kind"] == "other"
    assert parse_li2023(meta(library_name=["BM11", "RL60"]))["medium"].tolist() == [
        "M11 (HA)", "M60 (0.01xMA)"]                  # HA is M11 in the paper's Fig. 7 caption


def test_verhoeven2023_batches_are_separate_sources():
    out = parse_verhoeven2023(meta(
        library_name=["MQ210611-179", "MQ210611-183", "MQ210611-155", "MQ210129-1",
                      "MQ210129-2", "MQ210129-173", "MQ210611-9"],
        sample_title=["MV-PHN-AS-A", "MV-PHN-AS-str-A", "MV-PHN-PL7-Anaer-C", "MV-PHN-PL-AS2",
                      "MV-PHN-PL-NEG", "MV-PHN-PL-101", "MV-PHN-PL1-Med-C-1"]))
    assert out["kind"].tolist() == ["environment", "inoculum", "plate", "environment",
                                    "control", "plate", "other"]
    assert out["source_sample"].tolist()[:4] == ["AS_exp2", "AS_exp2", "AS_exp2", "AS_exp1"]
    assert out.loc[2, ["medium", "oxygen", "replicate"]].tolist() == ["PL7", "anaerobic", 3]
    assert out.loc[5, ["medium", "condition_detail"]].tolist() == ["unresolved", "plate 101"]


def test_villanelo2026_antibiotics_controls_and_raw():
    out = parse_villanelo2026(meta(library_name=["RAW_2", "CTRL_1", "CIP_HC_1", "PCRNEG1",
                                               "isolate_7"]))
    assert out["kind"].tolist() == ["environment", "plate", "plate", "control", "other"]
    assert out.loc[2, ["medium", "condition_detail", "replicate"]].tolist() == [
        "ASF agarose + CIP", "HC", 1]


def test_indraningrat2019_links_by_sponge_number():
    out = parse_indraningrat2019(meta(
        library_name=["Indra_Curacao_AS6", "Indra_Cultivation_Mu.AS11.67",
                      "Indra_Cultivation_XM10.16", "Indra_Curacao_SW7"],
        sample_title=["SCR50", "Mu.AS11", "MA.XM10", "W2"],
        description=["sponges", "scraping colonies on plates", "sponge_inoculum",
                     "filter_water"]))
    assert out["kind"].tolist() == ["environment", "plate", "inoculum", "other"]
    assert out["source_sample"].tolist()[:3] == ["AS6", "AS11", "XM10"]
    assert out.loc[1, "medium"] == "mucin agar"


def test_indraningrat2019_cryo_libraries_are_not_plates():
    out = parse_indraningrat2019(meta(library_name=["Indra_Cultivation_Co.AS11.49"],
                                      sample_title=["Co.AS11"],
                                      description=["scraping colonies on plates"]))
    assert out.loc[0, "kind"] == "other"


def test_ordonez2025_keeps_16s_only():
    out = parse_ordonez2025(meta(
        sample_alias=["16S_Ale_Root_3", "ITS_Ale_Stem_3", "16s-isolation_10", "its-isolation_10"],
        description=["Rhizosphere microbiota", "Stem microbiota", "leaf isolates",
                     "leaf isolates"]))
    assert out["kind"].tolist() == ["environment", "other", "plate", "other"]
    assert out.loc[0, "source_sample"] == "rhizosphere" and out.loc[2, "source_sample"] == "leaf"


def test_ortizseverin2024_only_direct_pools_are_plates():
    out = parse_ortizseverin2024(meta(sample_alias=["PSh-A", "AFS-B", "PSh-D", "PSh-metaG"],
                                 library_name=["PSh_A1", "AFS_B10", "PSh_D1", "PSh"]))
    assert out["kind"].tolist() == ["environment", "plate", "other", "other"]
    assert out.loc[1, ["source_sample", "replicate"]].tolist() == ["AFS", 10]


def test_yu2020_carcass_area_temperature():
    out = parse_yu2020(meta(library_name=["Back_D_SH", "Breast_D_BA7", "Breast_D_TSA30"]))
    assert out.loc[0, ["source_sample", "kind"]].tolist() == ["D_back", "environment"]
    assert out.loc[1, ["source_sample", "medium", "temperature_c"]].tolist() == [
        "D_breast", "blood agar", 7]
    assert out.loc[2, "temperature_c"] == 30


def test_hinsu2021_links_farm_replicates_only():
    out = parse_hinsu2021(meta(
        library_name=["R31", "R35", "R13", "R3-SEA"],
        sample_title=["Rhizosphere soil sample R31", "Rhizosphere soil sample R35",
                      "Rhizosphere soil sample R13",
                      "Cultured diversity of Sample R3 on media Soil Extract agar"]))
    assert out["kind"].tolist() == ["environment", "environment", "other", "plate"]
    assert out.loc[3, ["source_sample", "medium"]].tolist() == ["R3", "Soil Extract agar"]


def test_zehavi2018_links_plates_to_the_plated_cow_only():
    out = parse_zehavi2018(meta(sample_title=[
        "Rumen inoculum source Dilution factor 0", "Rumen inoculum Dilution factor 10-3",
        "Agar plate media M10 with rumen fluid Dilution factor 10-4 rep1",
        "Agar plate media M10 without rumen fluid Dilution factor 10-6 rep2",
        "Bovine rumen sample 35"]))
    assert out["kind"].tolist() == ["environment", "inoculum", "plate", "plate", "other"]
    assert out.loc[2, ["medium", "condition_detail", "replicate"]].tolist() == [
        "M10 + clarified rumen fluid", "dilution 10-4", 1]
    assert out.loc[3, "medium"] == "M10 defined"
    assert out["source_sample"].tolist()[:4] == ["plated_cow"] * 4


def test_every_included_study_has_a_parser():
    included = {s["id"] for s in load_config("census_studies")["studies"]
                if s["decision"] == "include"}
    assert included <= set(PARSERS)


def test_link_summary_counts_linked_plates():
    sheet = pd.DataFrame({
        "study_id": ["s"] * 5,
        "source_sample": ["a", "a", "a", "b", "b"],
        "kind": ["environment", "plate", "plate", "plate", "inoculum"]})
    row = link_summary(sheet).iloc[0]
    assert (row["environment_runs"], row["plate_runs"], row["linked_plate_runs"],
            row["linked_sources"]) == (1, 3, 2, 1)


def test_li2023_incubation_days_reads_table9_rows():
    text = "\n".join([
        "Medi Abbr. 15 C Incubation",
        "M1     BEN       BL1   8           BM1   3                 BH1   6           RL1   9"
        "           RM1   3                RH1                      3",
        "M11    HA        BL11  8   BM11  3   BH11  6   RL11  13  RM11  4   RH11  13"])
    d = li2023_incubation_days(text)
    assert len(d) == 12
    assert d["BL1"] == 8 and d["RH1"] == 3 and d["RL11"] == 13 and d["BM11"] == 3
