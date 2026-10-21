"""Phase 1 census, step 3: screen every candidate accession against its ENA records.

For each project accession mentioned in a shortlisted paper, fetch the run table (cached) and
score the evidence that the project contains BOTH plate-grown community libraries and
culture-independent libraries of the same material. The score only orders manual screening;
inclusion still comes from reading the methods (DECISIONS.md D-015).

Output: data/interim/census/accession_evidence.csv
"""

import re
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from plate_truth.census import RUN_FIELDS, SAMPLE_FIELDS, ena_table
from plate_truth.config import data_path

# Words in a library/sample label that indicate a plate-grown community, or a direct
# (culture-independent) sample of the source material.
PLATE_WORDS = re.compile(
    r"plate|agar|cultur|colon|isolat|medium|media|TSA|R2A|LB\b|BHI|MRS|PCA\b|MacConkey|"
    r"blood agar|scrap|wash|sweep|grown", re.I)
ENV_WORDS = re.compile(
    r"original|direct|uncultur|culture.independ|environment|source|inocul|total community|"
    r"raw |bulk|_ORI|community from|metagenom|swab|influent|sludge sample|soil sample", re.I)
ISOLATE_WORDS = re.compile(r"^isolate|single ?colony|strain|genome|MAG", re.I)
MIN_AMPLICON_RUNS = 6


def labels(runs: pd.DataFrame, samples: pd.DataFrame | None) -> pd.Series:
    """One text label per run, combining run and sample metadata."""
    cols = [c for c in ("library_name", "sample_title", "sample_alias") if c in runs.columns]
    text = runs[cols].astype(str).agg(" ".join, axis=1)
    if samples is not None and "description" in samples.columns:
        # map rather than merge: a sample accession can appear in several runs, and a merge
        # would change the row count and reindex the text
        extra = (samples["description"].fillna("") + " "
                 + samples.get("isolation_source", "").fillna(""))
        lookup = dict(zip(samples["sample_accession"], extra, strict=True))
        text = text + " " + runs["sample_accession"].map(lookup).fillna("")
    return text.fillna("").astype(str)


def screen(accession: str, cache) -> dict:
    """Evidence counts for one project accession."""
    row = {"accession": accession}
    try:
        runs = ena_table(accession, "read_run", RUN_FIELDS, cache / f"{accession}.tsv")
    except Exception as exc:                                        # noqa: BLE001
        return {**row, "status": f"run table error: {str(exc)[:60]}"}
    if len(runs) == 0:
        return {**row, "status": "no runs in ENA"}
    amp = runs[runs["library_strategy"] == "AMPLICON"]
    row.update(status="ok", n_runs=len(runs), n_amplicon=len(amp),
               platforms=";".join(sorted(runs["instrument_platform"].unique())))
    if len(amp) < MIN_AMPLICON_RUNS:
        row["status"] = "too few amplicon runs"
        return row
    samples = None
    try:
        samples = ena_table(accession, "sample", SAMPLE_FIELDS,
                            cache / f"{accession}_samples.tsv")
    except Exception:                                               # noqa: BLE001
        pass
    text = labels(amp, samples)
    plate = text.str.contains(PLATE_WORDS)
    env = text.str.contains(ENV_WORDS)
    isolate = text.str.contains(ISOLATE_WORDS)
    row.update(n_plate_like=int((plate & ~isolate).sum()), n_env_like=int(env.sum()),
               n_isolate_like=int(isolate.sum()),
               distinct_labels=int(text.map(lambda s: re.sub(r"\d+", "#", s)).nunique()),
               example_labels=" || ".join(
                   sorted({re.sub(r"\d+", "#", s)[:60] for s in text})[:6]))
    row["both_sides"] = row["n_plate_like"] > 0 and row["n_env_like"] > 0
    return row


def main() -> None:
    cache = data_path("raw") / "census" / "ena"
    cache.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(data_path("interim") / "census" / "candidates.csv")
    cand = cand[(cand["ft_plate_terms"] >= 1) & cand["accessions"].notna()]
    owner = {}
    for _, r in cand.iterrows():
        for acc in str(r["accessions"]).split(";"):
            owner.setdefault(acc, []).append(r["pmcid"])
    accessions = sorted(owner)
    print(f"screening {len(accessions)} accessions from {len(cand)} papers", flush=True)
    with ThreadPoolExecutor(6) as pool:
        rows = list(pool.map(lambda a: screen(a, cache), accessions))
    out = pd.DataFrame(rows)
    out["papers"] = out["accession"].map(lambda a: ";".join(owner[a][:3]))
    out = out.sort_values(["both_sides", "n_plate_like"], ascending=False)
    dest = data_path("interim") / "census" / "accession_evidence.csv"
    out.to_csv(dest, index=False)
    print(out["status"].value_counts().to_string())
    print(f"both sides present: {int(out['both_sides'].fillna(False).sum())} -> {dest}")


if __name__ == "__main__":
    main()
