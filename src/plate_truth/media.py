"""Medium classes from recipe facts (D-025).

Each medium in config/media.yaml carries facts transcribed from its recipe: whether it holds
an antibacterial agent, whether it is made from habitat material, whether it holds any complex
nutrient source, and its total organic content in g/L. The class is computed from those facts
by one ordered rule, so a reviewer can check a class by checking the recipe, not a judgement.
"""

import math

import pandas as pd

# Organic matter (g/L) at or below which a complex medium counts as dilute. R2A, the reference
# low-nutrient medium, has about 3 g/L; full-strength general media have 6 g/L or more.
DILUTE_MAX_G_PER_L = 4.0


def classify_medium(antibacterial: bool, habitat_material: bool, complex_source: bool,
                    organic_g_per_l: float | None) -> tuple[str, str]:
    """(medium_class, rule) for one medium; rules are applied in this order:

    R1 selective            an antibiotic or other antibacterial agent is added
    R2 environment_extract  habitat material (soil extract, wastewater, sludge filtrate,
                            seawater base, plant extract) is an ingredient
    R3 defined_minimal      no complex nutrient source (peptone, extract, hydrolysate,
                            infusion, milk, blood)
    R4 dilute_oligotrophic  complex, total organic matter <= DILUTE_MAX_G_PER_L
    R5 rich_complex         complex, total organic matter above that
    'unresolved' when organic content is needed but unknown.
    """
    if antibacterial:
        return "selective", "R1"
    if habitat_material:
        return "environment_extract", "R2"
    if not complex_source:
        return "defined_minimal", "R3"
    if organic_g_per_l is None or (isinstance(organic_g_per_l, float)
                                   and math.isnan(organic_g_per_l)):
        return "unresolved", "R4/R5 needs organic content"
    if organic_g_per_l <= DILUTE_MAX_G_PER_L:
        return "dilute_oligotrophic", "R4"
    return "rich_complex", "R5"


def media_table(entries: list[dict]) -> pd.DataFrame:
    """config/media.yaml entries -> one row per (study, medium) with the computed class.

    An explicit medium_class_override wins (rule 'override: <reason>'); entries without recipe
    facts get class NULL and rule 'pending recipe source'.
    """
    rows = []
    for e in entries:
        if e.get("medium_class_override"):
            cls, rule = e["medium_class_override"], f"override: {e.get('override_reason')}"
        elif "complex_source" in e:
            cls, rule = classify_medium(bool(e.get("antibacterial")),
                                        bool(e.get("habitat_material")),
                                        bool(e["complex_source"]), e.get("organic_g_per_l"))
        else:
            cls, rule = None, "pending recipe source"
        rows.append({"study_id": e["study_id"], "medium": e["medium"],
                     "medium_name": e.get("medium_name"), "medium_class": cls,
                     "class_rule": rule})
    return pd.DataFrame(rows, columns=["study_id", "medium", "medium_name", "medium_class",
                                       "class_rule"])
