"""Phase 2, step 1: one sample sheet for every included study.

For each included study in config/census_studies.yaml, joins the cached ENA run and sample
tables, applies that study's parser from plate_truth.samples (which source sample each run
belongs to, and whether it is environment, inoculum, plate, control or other), and fills
study-wide fields (temperature, oxygen, incubation) from the config where the study used a
single value. Only amplicon runs are kept. Wang et al. 2023 plates get their own incubation
time from Supplementary Table 9 (text extracted from page 26 of the supplement PDF with
pdftotext -layout into data/raw/census/supp/PMC10516943/table9.txt).

Outputs: data/processed/sample_sheet.csv, data/processed/sample_links.csv (per-study link counts)
"""

import pandas as pd

from plate_truth.census import RUN_FIELDS, SAMPLE_FIELDS, ena_table
from plate_truth.config import data_path, load_config
from plate_truth.samples import PARSERS, link_summary, wang_incubation_days

# Per-plate incubation times published as a table (study id -> function of the raw data dir).
WANG_TABLE9 = ("census", "supp", "PMC10516943", "table9.txt")

KEEP = ["study_id", "study_group", "run_accession", "sample_accession", "library_name",
        "instrument_platform", "read_count", "source_sample", "kind", "replicate", "medium",
        "oxygen", "temperature_c", "incubation_days", "condition_detail", "region",
        "primer_fwd", "primer_rev"]


def study_metadata(study: dict) -> pd.DataFrame:
    """The study's amplicon runs joined to their ENA sample records (cached tables)."""
    cache = data_path("raw") / "census" / "ena"
    acc = study["accession"]
    runs = ena_table(acc, "read_run", RUN_FIELDS, cache / f"{acc}.tsv")
    samples = ena_table(acc, "sample", SAMPLE_FIELDS, cache / f"{acc}_samples.tsv")
    extra = [c for c in samples.columns if c not in runs.columns or c == "sample_accession"]
    meta = runs.merge(samples[extra], on="sample_accession", how="left")
    return meta[meta["library_strategy"] == "AMPLICON"].reset_index(drop=True)


def fill_study_fields(sheet: pd.DataFrame, study: dict) -> pd.DataFrame:
    """Fill temperature/oxygen from the config when the study used a single value."""
    out = sheet.copy()
    for col, key in (("temperature_c", "temperature_c"), ("oxygen", "oxygen")):
        values = study.get(key) or []
        if len(values) == 1:
            out[col] = out[col].fillna(values[0]).where(out["kind"] == "plate")
    out["incubation_days"] = study.get("incubation_days")
    for key in ("study_group", "region", "primer_fwd", "primer_rev"):
        out[key] = study.get(key)
    out["study_id"] = study["id"]
    return out


def main() -> None:
    sheets = []
    for study in load_config("census_studies")["studies"]:
        if study["decision"] != "include":
            continue
        meta = study_metadata(study)
        parsed = PARSERS[study["id"]](meta)
        sheet = fill_study_fields(pd.concat([meta, parsed], axis=1), study)
        if study["id"] == "wang2023_desertsoil":
            days = wang_incubation_days(data_path("raw").joinpath(*WANG_TABLE9).read_text())
            plate = sheet["kind"] == "plate"
            sheet.loc[plate, "incubation_days"] = sheet.loc[plate, "library_name"].map(days)
            assert sheet.loc[plate, "incubation_days"].notna().all(), "Table 9 is incomplete"
        sheets.append(sheet[sheet["kind"] != "other"][KEEP])
        print(f"{study['id']}: {len(meta)} amplicon runs -> "
              f"{sheet['kind'].value_counts().to_dict()}")
    sheet = pd.concat(sheets, ignore_index=True)
    out = data_path("processed")
    sheet.to_csv(out / "sample_sheet.csv", index=False)
    links = link_summary(sheet)
    links.to_csv(out / "sample_links.csv", index=False)
    print()
    print(links.to_string(index=False))
    unlinked = sheet[(sheet["kind"] == "plate")
                     & ~sheet.set_index(["study_id", "source_sample"]).index.isin(
                         sheet[sheet["kind"] == "environment"]
                         .set_index(["study_id", "source_sample"]).index)]
    if len(unlinked):
        print("\nplate runs with no environment run for their source sample:")
        print(unlinked.groupby(["study_id", "source_sample"], dropna=False).size().to_string())


if __name__ == "__main__":
    main()
