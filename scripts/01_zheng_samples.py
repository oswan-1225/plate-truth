"""Build the Zheng et al. 2024 soil sample sheet from the ENA run table.

Source: ENA Portal API filereport for PRJNA1029248 (read_run), cached in data/raw/zheng2024/.
Output: data/processed/zheng2024_samples.csv, one row per soil run (environment, inoculum, plate).
"""

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.http import download
from plate_truth.samples import condition_id, parse_zheng_library

FIELDS = ("run_accession,sample_accession,library_name,library_strategy,instrument_model,"
          "read_count,fastq_ftp,fastq_md5")


def main() -> None:
    study = next(s for s in load_config("studies")["studies"] if s["id"] == "zheng2024_soil")
    raw = data_path("raw") / "zheng2024" / "ena_runs.tsv"
    download("https://www.ebi.ac.uk/ena/portal/api/filereport", raw,
             params={"accession": study["accession"], "result": "read_run", "format": "tsv",
                     "fields": FIELDS})
    runs = pd.read_csv(raw, sep="\t")
    runs = runs[runs["library_strategy"] == "AMPLICON"]
    parsed = parse_zheng_library(runs["library_name"])
    sheet = pd.concat([runs, parsed], axis=1)
    sheet = sheet[sheet["source_sample"].isin(study["source_samples"])].copy()
    sheet["condition_id"] = condition_id(sheet)
    sheet["study_id"] = study["id"]
    cols = ["study_id", "run_accession", "sample_accession", "library_name", "source_sample",
            "replicate", "kind", "medium", "oxygen", "condition_id", "read_count",
            "fastq_ftp", "fastq_md5"]
    sheet = sheet[cols].sort_values(["source_sample", "kind", "condition_id", "replicate"])
    out = data_path("processed") / "zheng2024_samples.csv"
    sheet.to_csv(out, index=False)
    print(f"{len(sheet)} soil runs -> {out}")
    print(sheet.groupby(["source_sample", "kind"]).size().unstack(fill_value=0))
    print(sheet.dropna(subset=["condition_id"]).groupby("condition_id").size().describe())


if __name__ == "__main__":
    main()
