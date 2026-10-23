"""Phase 4, step 2: analysis rows for Q5 (medium class by taxon) and Q6 (curation hypothesis).

Starts from the Phase 3 primary rows (same definitions, outcome and features) and:
- recomputes each condition's medium class from config/media.yaml, so the four Hinsu et al.
  media classified for Phase 4 are included (Phase 3 keeps them unknown, as registered);
- adds the SILVA family and full phylum of each ASV and the interaction keys for the random
  effects (taxon x study, taxon x medium class);
- keeps the study groups with at least min_known_classes known medium classes, and only their
  rows with a known class (Q5);
- for Q6, keeps the Q5 rows whose ASV has a classed culture-collection record (script 33) and
  adds the recorded-class shares.
Row counts are printed; outcomes are not summarised here.

Outputs: data/interim/phase4/rows/q5.parquet, q6.parquet, q6_strain16s.parquet
Usage: python scripts/34_phase4_rows.py
"""

import sqlite3

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.media import media_table
from plate_truth.phase3 import keep_levels
from plate_truth.phase4 import GROUP, combine, taxon_label
from plate_truth.queries import asv_taxonomy


def main() -> None:
    cfg = load_config("phase4")
    classes = cfg["medium_classes"]
    out = data_path("interim") / "phase4" / "rows"
    out.mkdir(parents=True, exist_ok=True)

    rows = pd.read_parquet(data_path("interim") / "phase3" / "rows" / "primary.parquet")
    mt = media_table(load_config("media")["media"])
    cls = mt.set_index(["study_id", "medium"])["medium_class"]
    keys = pd.MultiIndex.from_arrays([rows["study_id"], rows["medium_label"]])
    rows["medium_class"] = cls.reindex(keys).to_numpy()
    rows["medium"] = keep_levels(rows["medium_class"], classes, "unknown")

    tax = asv_taxonomy(sqlite3.connect(data_path("database"))).set_index("asv_id")
    rows["family"] = rows["asv_id"].map(taxon_label(tax, "family"))
    rows["phylum_full"] = rows["asv_id"].map(tax["phylum"].fillna("unassigned"))

    known = rows[rows["medium"] != "unknown"]
    n_classes = known.groupby(GROUP)["medium"].nunique()
    groups = sorted(n_classes[n_classes >= cfg["min_known_classes"]].index)
    q5 = known[known[GROUP].isin(groups)].copy()
    for taxon in ["family", "phylum_full"]:
        q5[f"{taxon}_study"] = combine(q5[taxon], q5["study_id"])
        q5[f"{taxon}_medium"] = combine(q5[taxon], q5["medium"])
    q5 = q5.reset_index(drop=True)
    q5.to_parquet(out / "q5.parquet")
    print(f"Q5 groups: {groups}", flush=True)
    print(q5.groupby([GROUP, "medium"]).agg(rows=("asv_id", "size"),
                                             asvs=("asv_id", "nunique"),
                                             conditions=("condition_id", "nunique")).to_string())
    print(f"Q5: {len(q5)} rows, {q5['family'].nunique()} families, "
          f"{q5['phylum_full'].nunique()} phyla", flush=True)

    for name in ["q6", "q6_strain16s"]:
        tag = "asv_recorded" if name == "q6" else "asv_recorded_strain16s"
        rec = pd.read_parquet(data_path("interim") / "phase4" / f"{tag}.parquet")
        q6 = q5.merge(rec, on="asv_id").reset_index(drop=True)
        q6.to_parquet(out / f"{name}.parquet")
        print(f"{name}: {len(q6)} rows, {q6['asv_id'].nunique()} ASVs "
              f"({q6['asv_id'].nunique() / q5['asv_id'].nunique():.1%} of Q5 ASVs)", flush=True)
        print((q6.drop_duplicates("asv_id")[[f"any_{c}" for c in classes]] > 0).mean()
              .round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
