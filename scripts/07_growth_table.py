"""Observed growth table for the Zheng et al. soils (primary settings).

Primary settings: environment = original soil (ORI) pooled over PCR replicates; ASV present
in a soil if >= environment.min_reads reads; plates pooled per condition; thresholds from
config/thresholds.yaml. ASVs matching the paper's candidate-contaminant list (Supplementary
Data 2) are kept and flagged (is_contaminant_candidate); excluding them is a sensitivity run.

Outputs (data/processed/):
  zheng2024_contaminant_asvs.csv     our ASVs matching the paper's contaminant list
  zheng2024_asv_condition.csv        (ASV, plating condition) with enriched/detected
  zheng2024_asv_source.csv           (ASV, soil): outcome in any condition
  zheng2024_plate_only.csv           on plates but not present in that soil
"""

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.growth import plate_only_asvs
from plate_truth.identity import read_hits, search
from plate_truth.phase0 import growth_inputs, label_growth, per_source
from plate_truth.reference import write_fasta

SUPP = "41467_2024_53446_MOESM3_ESM.xlsx"


def contaminant_asvs(asvs, supp_xlsx, work) -> pd.DataFrame:
    """ASVs within 1% (terminal gaps ignored) of a sequence on the paper's contaminant list.

    The list includes the spiked internal standard (strain NS6-2) and ASVs whose abundance
    tracked its dilution series in extraction blanks."""
    raw = pd.read_excel(supp_xlsx, sheet_name="Supplementary Data 2", header=1)
    raw = raw.rename(columns={raw.columns[0]: "contaminant_id", raw.columns[2]: "seq"})
    raw = raw.dropna(subset=["seq"])
    db = work / "contaminants.fasta"
    write_fasta(((f"c{i}", s.strip().upper()) for i, s in enumerate(raw["seq"])), db)
    search(asvs, db, work / "contaminant_hits.tsv", iddef=2, min_id=0.99, maxaccepts=1,
           maxrejects=0)
    hits = read_hits(work / "contaminant_hits.tsv")
    ids = dict(zip([f"c{i}" for i in range(len(raw))], raw["contaminant_id"], strict=True))
    hits["contaminant_id"] = hits["target"].map(ids)
    return hits.rename(columns={"query": "asv_id"})[["asv_id", "contaminant_id", "id2"]]


def main() -> None:
    th = load_config("thresholds")
    z = data_path("interim") / "zheng2024"
    work = data_path("interim") / "growth"
    work.mkdir(parents=True, exist_ok=True)
    proc = data_path("processed")
    sheet = pd.read_csv(proc / "zheng2024_samples.csv")
    counts = pd.read_csv(z / "counts.csv")

    cont = contaminant_asvs(z / "asv_seqs.fasta", data_path("raw") / "zheng2024" / "supp" / SUPP,
                            work)
    cont.to_csv(proc / "zheng2024_contaminant_asvs.csv", index=False)
    share = counts[counts["asv_id"].isin(cont["asv_id"])]["count"].sum() / counts["count"].sum()
    print(f"contaminant-list ASVs: {len(cont)} ({share:.2%} of all reads)")

    env, plates, conds = growth_inputs(counts, sheet, "environment")
    table = label_growth(env, plates, conds, th["environment"]["min_reads"],
                         th["grew"]["min_reads"], th["grew"]["min_fold_enrichment"],
                         th["detected"]["min_reads"])
    table["is_contaminant_candidate"] = table["asv_id"].isin(cont["asv_id"])
    table.to_csv(proc / "zheng2024_asv_condition.csv", index=False)
    src = per_source(table)
    src.to_csv(proc / "zheng2024_asv_source.csv", index=False)
    only = plate_only_asvs(env, plates, conds, th["environment"]["min_reads"],
                           th["grew"]["min_reads"])
    only.to_csv(proc / "zheng2024_plate_only.csv", index=False)

    print(f"(ASV, condition) rows: {len(table)}; enriched {table['enriched_on_plates'].mean():.3f}"
          f"; detected {table['detected_on_plates'].mean():.3f}")
    print(f"(ASV, soil) rows: {len(src)}; enriched in any condition "
          f"{src['enriched_on_plates'].mean():.3f}; "
          f"detected {src['detected_on_plates'].mean():.3f}")
    print(f"plate-only (ASV, condition) rows with >= {th['grew']['min_reads']} reads: {len(only)}; "
          f"distinct ASVs {only['asv_id'].nunique()}")
    print(src.groupby("source_sample")[["enriched_on_plates", "detected_on_plates"]].mean())


if __name__ == "__main__":
    main()
