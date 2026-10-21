import math

import pandas as pd

from plate_truth.media import classify_medium, media_table


def test_rule_order_selective_first():
    # an antibiotic plate made on sludge filtrate is selective, not environment_extract
    assert classify_medium(True, True, True, 2.0) == ("selective", "R1")
    assert classify_medium(False, True, True, 2.0) == ("environment_extract", "R2")


def test_defined_minimal_ignores_organic_content():
    assert classify_medium(False, False, False, 30.0) == ("defined_minimal", "R3")


def test_dilute_threshold_puts_r2a_in_dilute():
    # R2A: ~3.05 g/L organic (Wang et al. Suppl. Table 8 recipe); 1x nutrient agar: 13 g/L
    assert classify_medium(False, False, True, 3.05)[0] == "dilute_oligotrophic"
    assert classify_medium(False, False, True, 4.0)[0] == "dilute_oligotrophic"
    assert classify_medium(False, False, True, 13.0)[0] == "rich_complex"


def test_media_table_override_facts_and_pending():
    entries = [
        {"study_id": "s", "medium": "pool", "medium_class_override": "other",
         "override_reason": "pooled"},
        {"study_id": "s", "medium": "M3", "antibacterial": True, "habitat_material": False,
         "complex_source": True, "organic_g_per_l": 11.1},
        {"study_id": "s", "medium": "PCA"},
    ]
    t = media_table(entries).set_index("medium")
    assert t.loc["pool", ["medium_class", "class_rule"]].tolist() == ["other", "override: pooled"]
    assert t.loc["M3", "medium_class"] == "selective"
    assert pd.isna(t.loc["PCA", "medium_class"]) and t.loc["PCA", "class_rule"].startswith("pending")


def test_unknown_organic_content_is_unresolved():
    assert classify_medium(False, False, True, None)[0] == "unresolved"
    assert classify_medium(False, False, True, math.nan)[0] == "unresolved"
