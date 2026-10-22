"""Phase 2 sensitivity: do the Zheng et al. growth calls change when every read is counted?

Builds the (ASV, plating condition) growth table of each Zheng et al. study twice, from the
paired-read counts (script 16, the primary table) and from the every-read counts (script 23),
with the same definitions, thresholds and depth rule (D-027), and compares them: agreement of
the enriched and detected calls on rows present in both, and the in-sample share of rows that
grew per identity band. In-sample and descriptive: this compares two count tables, it does not
measure prediction performance.

Output: results/phase2/zheng_count_sensitivity.csv (per study, counts and identity band),
        results/phase2/zheng_count_agreement.csv (per study)
"""

import pandas as pd

from plate_truth.calibration import identity_bin
from plate_truth.config import data_path, load_config, project_root
from plate_truth.database import condition_ids, is_organelle
from plate_truth.growth import pooled_abundance
from plate_truth.phase0 import label_growth
from plate_truth.qc import pool_depth, usable_conditions

STUDIES = ["zheng2024_soil", "zheng2024_sludge"]


def growth_table(counts: pd.DataFrame, sheet: pd.DataFrame, organelles: set, th: dict
                 ) -> pd.DataFrame:
    """(ASV, condition) rows with outcomes for one study's count table, usable conditions only."""
    c = counts[~counts["asv_id"].isin(organelles)].merge(
        sheet[["run_accession", "kind", "source_sample", "condition_id"]],
        left_on="run", right_on="run_accession")
    env_c, plate_c = c[c["kind"] == "environment"], c[c["kind"] == "plate"]
    conds = sheet.loc[sheet["kind"] == "plate", ["condition_id", "source_sample", "medium",
                                                 "oxygen"]].drop_duplicates()
    keys = conds.assign(study_id="s")
    flags = usable_conditions(keys, pool_depth(plate_c.assign(study_id="s"),
                                               ["study_id", "condition_id"]),
                              pool_depth(env_c.assign(study_id="s"),
                                         ["study_id", "source_sample"]),
                              th["depth"]["min_pool_reads"])
    conds = conds[conds["condition_id"].isin(flags.loc[flags["usable"], "condition_id"])]
    env = pooled_abundance(env_c, "source_sample")
    plates = pooled_abundance(plate_c, "condition_id")
    return label_growth(env, plates, conds, th["environment"]["min_reads"],
                        th["grew"]["min_reads"], th["grew"]["min_fold_enrichment"],
                        th["detected"]["min_reads"])


def main() -> None:
    th = load_config("thresholds")
    edges = th["identity"]["bin_edges_percent"]
    sheet_all = pd.read_csv(data_path("processed") / "sample_sheet.csv")
    sheet_all["condition_id"] = condition_ids(sheet_all)
    bands, agree = [], []
    for sid in STUDIES:
        d = data_path("interim") / "phase2" / sid
        sheet = sheet_all[sheet_all["study_id"] == sid]
        tax = pd.read_csv(d / "taxonomy.csv")
        organelles = set(tax.loc[is_organelle(tax), "asv_id"])
        nr = pd.read_csv(d / "nearest_relative.csv")
        ident = nr.loc[nr["reference_set"] == "all"].set_index("asv_id")["identity_id1"]
        tables = {name: growth_table(pd.read_csv(d / f), sheet, organelles, th)
                  for name, f in (("paired", "counts.csv"), ("every_read", "counts_r1all.csv"))}
        for name, t in tables.items():
            t["identity_bin"] = identity_bin(t["asv_id"].map(ident), edges)
            g = t.groupby("identity_bin", observed=False).agg(
                rows=("asv_id", "size"), enriched=("enriched_on_plates", "mean"),
                detected=("detected_on_plates", "mean"))
            bands.append(g.reset_index().assign(study_id=sid, counts=name))
            bands.append(pd.DataFrame([{"identity_bin": "all", "rows": len(t),
                                        "enriched": t["enriched_on_plates"].mean(),
                                        "detected": t["detected_on_plates"].mean(),
                                        "study_id": sid, "counts": name}]))
        key = ["condition_id", "asv_id"]
        both = tables["paired"].merge(tables["every_read"], on=key, suffixes=("_p", "_e"))
        row = {"study_id": sid, "rows_paired": len(tables["paired"]),
               "rows_every_read": len(tables["every_read"]), "rows_both": len(both)}
        for out in ("enriched_on_plates", "detected_on_plates"):
            p, e = both[f"{out}_p"], both[f"{out}_e"]
            row[f"{out}_agree"] = (p == e).mean()
            row[f"{out}_only_paired"] = int((p & ~e).sum())
            row[f"{out}_only_every_read"] = int((~p & e).sum())
            row[f"{out}_both"] = int((p & e).sum())
        agree.append(row)
        print(pd.Series(row).to_string(), flush=True)
    out = project_root() / "results" / "phase2"
    cols = ["study_id", "counts", "identity_bin", "rows", "enriched", "detected"]
    pd.concat(bands)[cols].round(5).to_csv(out / "zheng_count_sensitivity.csv", index=False)
    pd.DataFrame(agree).round(4).to_csv(out / "zheng_count_agreement.csv", index=False)


if __name__ == "__main__":
    main()
