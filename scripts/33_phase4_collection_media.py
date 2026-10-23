"""Phase 4, step 1: medium classes of culture-collection records, linked to every ASV.

1. Every MediaDive medium is classed by the plate-media rule from its MediaDive recipe
   (plate_truth.phase4, rules in config/phase4.yaml).
2. Every culture-medium entry in the cached BacDive records is read; growth-positive entries
   that link a MediaDive medium give each strain its recorded medium classes.
3. Each ASV's tied best hits (primary identity, all-sources reference) are linked to BacDive
   strains (a strain's own 16S, or a GTDB isolate genome that the strain's record lists), and
   the share of those strains recorded on each class is computed. A second table uses the
   strains' own 16S links only (sensitivity).

Inputs: data/raw/mediadive/ (script 32), data/raw/bacdive/, data/plate_truth.sqlite,
data/interim/reference/bacdive_genomes.csv (script 03), Phase 3 rows (for the analysed ASVs).
Outputs:
  results/phase4/collection_media_classes.csv   one row per MediaDive medium: facts and class
  results/phase4/collection_coverage.csv        per study group: analysed ASVs, share linked
  results/phase4/collection_class_counts.csv    recorded strains per class
  data/interim/phase4/bacdive_culture_media.parquet, strain_classes.parquet,
  asv_recorded.parquet, asv_recorded_strain16s.parquet

Usage: python scripts/33_phase4_collection_media.py
"""

import json
import sqlite3

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase4 import (
    GROUP,
    NO_RECIPE,
    asv_strains,
    classify_collection_media,
    culture_media_records,
    ingredient_roles,
    medium_facts,
    member_strains,
    recorded_shares,
    strain_classes,
)


def load_mediadive(raw) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(media, ingredients, flattened compositions) from the script-32 cache."""
    media = pd.concat([pd.DataFrame(json.loads(p.read_text(encoding="utf-8"))["data"])
                       for p in sorted(raw.glob("media_*.json"))], ignore_index=True)
    media["id"] = media["id"].astype(str)
    ingr = pd.DataFrame([json.loads(p.read_text(encoding="utf-8"))["data"]
                         for p in sorted((raw / "ingredient").glob("*.json"))])
    comp = []
    for mid in media["id"]:
        body = json.loads((raw / "composition" / f"{mid}.json").read_text(encoding="utf-8"))
        for r in body.get("data") or []:
            comp.append({"medium_id": mid, "ingredient_id": r["id"], "g_l": r.get("g_l"),
                         "optional": r.get("optional", 0)})
    return media, ingr, pd.DataFrame(comp)


def load_bacdive_media(raw) -> pd.DataFrame:
    """Culture-medium entries of every cached BacDive record."""
    rows = []
    for path in sorted((raw / "bacdive").rglob("batch_*.json")):
        results = json.loads(path.read_text(encoding="utf-8")).get("results") or {}
        if isinstance(results, dict):
            for bid, rec in results.items():
                rows += culture_media_records(int(bid), rec)
    return pd.DataFrame(rows)


def main() -> None:
    cfg = load_config("phase4")
    classes = cfg["medium_classes"]
    raw = data_path("raw")
    out = data_path("interim") / "phase4"
    res = data_path("results") / "phase4"
    out.mkdir(parents=True, exist_ok=True)
    res.mkdir(parents=True, exist_ok=True)

    media, ingr, comp = load_mediadive(raw / "mediadive")
    roles = ingredient_roles(ingr, cfg["collection_media"])
    facts = medium_facts(comp, roles, cfg["collection_media"]["composite_agar_g_per_l"])
    mc = classify_collection_media(facts)
    mc = media[["id", "name", "source", "complex_medium"]].rename(columns={"id": "medium_id"}) \
        .merge(mc, on="medium_id", how="left")
    mc["medium_class"] = mc["medium_class"].fillna(NO_RECIPE)
    mc.to_csv(res / "collection_media_classes.csv", index=False)
    print(mc["medium_class"].value_counts().to_string(), flush=True)

    rec = load_bacdive_media(raw)
    rec.to_parquet(out / "bacdive_culture_media.parquet")
    sc = strain_classes(rec, mc, classes)
    sc.to_parquet(out / "strain_classes.parquet")
    used = rec.dropna(subset=["mediadive_id"]).merge(
        mc[["medium_id", "medium_class"]], left_on="mediadive_id", right_on="medium_id",
        how="left")
    counts = pd.DataFrame({
        "records_positive": used[used["growth"].astype(str).str.lower() != "no"]
        .groupby("medium_class").size(),
        "strains": sc.groupby("medium_class")["bacdive_id"].nunique()}).fillna(0)
    counts.rename_axis("medium_class").reset_index().to_csv(
        res / "collection_class_counts.csv", index=False)
    print(counts.to_string(), flush=True)

    con = sqlite3.connect(data_path("database"))
    members = pd.read_sql("SELECT reference, ref_id, member_id FROM reference_strains", con)
    nr = pd.read_sql("SELECT asv_id, reference, targets_id1 AS targets FROM nearest_relative "
                     "WHERE reference_set = 'all'", con)
    genomes = pd.read_csv(data_path("interim") / "reference" / "bacdive_genomes.csv")
    analysed = pd.read_parquet(data_path("interim") / "phase3" / "rows" / "primary.parquet",
                               columns=["asv_id", GROUP]).drop_duplicates("asv_id")
    nr = nr[nr["asv_id"].isin(analysed["asv_id"])]

    links = member_strains(members, genomes)
    for name, lk in [("asv_recorded", links),
                     ("asv_recorded_strain16s",
                      member_strains(members[members["member_id"].str.startswith(
                          "bacdive_strain|")], genomes))]:
        a2s = asv_strains(nr, lk)
        rs = recorded_shares(a2s, sc, classes)
        rs.to_parquet(out / f"{name}.parquet")
        if name == "asv_recorded":
            cov = analysed.assign(linked=analysed["asv_id"].isin(a2s["asv_id"]),
                                  classed=analysed["asv_id"].isin(rs["asv_id"]))
            cov = cov.groupby(GROUP).agg(asvs=("asv_id", "size"), linked=("linked", "mean"),
                                         with_classed_record=("classed", "mean"))
            cov.loc["all"] = [len(analysed), analysed["asv_id"].isin(a2s["asv_id"]).mean(),
                              analysed["asv_id"].isin(rs["asv_id"]).mean()]
            cov.reset_index().to_csv(res / "collection_coverage.csv", index=False)
            print(cov.round(3).to_string(), flush=True)
        print(f"{name}: {len(rs)} ASVs with a classed record", flush=True)


if __name__ == "__main__":
    main()
