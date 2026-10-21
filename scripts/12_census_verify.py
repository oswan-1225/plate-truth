"""Phase 1 census, step 2: verify each curated candidate against its ENA records.

Reads config/census_studies.yaml (judgement calls, made by reading the papers), fetches and
caches each study's ENA run and sample tables, labels every run as environment / plate / other
using the per-study patterns, and writes the census with counts computed from the real tables.
Nothing countable is typed by hand.

Outputs: data/processed/census.csv, data/processed/census_runs.csv (per-run labels)
"""

import pandas as pd

from plate_truth.census import RUN_FIELDS, SAMPLE_FIELDS, ena_table
from plate_truth.config import data_path, load_config


def label_runs(runs: pd.DataFrame, samples: pd.DataFrame, study: dict) -> pd.DataFrame:
    """Add a 'role' column (environment / plate / other) using the study's patterns.

    The match field may live in the run table or only in the sample table, so the two are
    joined on sample_accession first.
    """
    df = runs.copy()
    field = study["match_field"]
    if field not in df.columns and field in samples.columns:
        df = df.merge(samples[["sample_accession", field]], on="sample_accession", how="left")
    df[field] = df.get(field, pd.Series("", index=df.index)).fillna("")
    role = pd.Series("other", index=df.index)
    if study.get("plate_pattern"):
        role[df[field].str.contains(study["plate_pattern"], regex=True, na=False)] = "plate"
    if study.get("env_pattern"):
        role[df[field].str.contains(study["env_pattern"], regex=True, na=False)] = "environment"
    df["role"] = role
    df["study_id"] = study["id"]
    return df


def main() -> None:
    cache = data_path("raw") / "census" / "ena"
    cache.mkdir(parents=True, exist_ok=True)
    rows, per_run = [], []
    for study in load_config("census_studies")["studies"]:
        acc = study["accession"]
        try:
            runs = ena_table(acc, "read_run", RUN_FIELDS, cache / f"{acc}.tsv")
            samples = ena_table(acc, "sample", SAMPLE_FIELDS, cache / f"{acc}_samples.tsv")
        except Exception as exc:                                   # noqa: BLE001
            rows.append({**summary_fields(study), "n_runs": 0, "note": f"ENA error: {exc}"})
            continue
        labelled = label_runs(runs, samples, study)
        amp = labelled[labelled["library_strategy"] == "AMPLICON"]
        counts = amp["role"].value_counts()
        per_run.append(labelled[labelled["role"] != "other"])
        rows.append({
            **summary_fields(study),
            "n_runs": len(runs),
            "n_amplicon_runs": len(amp),
            "n_environment_runs": int(counts.get("environment", 0)),
            "n_plate_runs": int(counts.get("plate", 0)),
            "n_other_runs": int(counts.get("other", 0)),
            "platform_observed": ";".join(sorted(runs["instrument_platform"].unique())),
            "note": "",
        })
    census = pd.DataFrame(rows)
    out = data_path("processed")
    census.to_csv(out / "census.csv", index=False)
    if per_run:
        pd.concat(per_run).to_csv(out / "census_runs.csv", index=False)
    show = ["study_id", "accession", "habitat", "host_associated", "decision",
            "n_amplicon_runs", "n_environment_runs", "n_plate_runs", "platform_observed"]
    print(census[show].to_string(index=False))
    usable = census[(census["decision"] == "include") & (census["n_environment_runs"] > 0)
                    & (census["n_plate_runs"] > 0)]
    print(f"\nincluded with both sides present: {len(usable)}")
    print("habitats (non-host):",
          sorted(usable.loc[~usable["host_associated"], "habitat"].unique()))
    print("habitats (host-associated):",
          sorted(usable.loc[usable["host_associated"], "habitat"].unique()))


def summary_fields(study: dict) -> dict:
    """The curated fields copied through to the census table."""
    keys = ["id", "accession", "doi", "citation", "habitat", "host_associated", "media",
            "temperature_c", "oxygen", "incubation_days", "region", "primers", "decision",
            "reason"]
    out = {k if k != "id" else "study_id": study.get(k) for k in keys}
    for k in ("media", "temperature_c", "oxygen"):
        if isinstance(out[k], list):
            out[k] = ";".join(str(v) for v in out[k])
    return out


if __name__ == "__main__":
    main()
