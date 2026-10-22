"""Phase 3, step 1: build the (ASV, plating condition) analysis rows for the primary analysis
and every sensitivity analysis (ANALYSIS_PLAN.md sections 3 and 7).

Reads data/plate_truth.sqlite. One parquet file per variant in data/interim/phase3/rows/
(skipped if present; delete a file to rebuild it). Prints and writes a per-study row and event
summary for each variant (results/phase3/rows_summary.csv): descriptive counts, not performance.

Usage: python scripts/26_phase3_rows.py [variant ...]   (default: all variants)
"""

import sys

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import VARIANTS, analysis_rows
from plate_truth.queries import connect


def main(names: list[str]) -> None:
    th, cfg = load_config("thresholds"), load_config("phase3")
    out = data_path("interim") / "phase3" / "rows"
    out.mkdir(parents=True, exist_ok=True)
    res = data_path("results") / "phase3"
    res.mkdir(parents=True, exist_ok=True)
    summaries = []
    for name in names:
        path = out / f"{name}.parquet"
        tag = VARIANTS[name].get("db_tag")
        db = data_path("interim") / f"plate_truth__{tag}.sqlite" if tag else None
        if path.exists():
            rows = pd.read_parquet(path)
        elif db is not None and not db.exists():
            print(f"{name}: {db.name} not built yet, skipped", flush=True)
            continue
        else:
            rows = analysis_rows(connect(db), th, cfg, VARIANTS[name])
            rows.to_parquet(path, index=False)
        s = rows.groupby(["study_group", "study_id"]).agg(
            conditions=("condition_id", "nunique"), asvs=("asv_id", "nunique"),
            rows=("y", "size"), events=("y", "sum"), no_hit=("no_hit", "sum"))
        summaries.append(s.reset_index().assign(variant=name))
        print(f"{name}: {len(rows):,} rows, {int(rows['y'].sum()):,} events", flush=True)
    summ = pd.concat(summaries, ignore_index=True)
    old = res / "rows_summary.csv"
    if old.exists():
        prev = pd.read_csv(old)
        summ = pd.concat([prev[~prev["variant"].isin(names)], summ], ignore_index=True)
    summ.to_csv(old, index=False)


if __name__ == "__main__":
    main(sys.argv[1:] or list(VARIANTS))
