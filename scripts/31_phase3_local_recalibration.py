"""Phase 3, exploratory (not registered): does knowing a study's level fix calibration?

H2 (calibrated across studies) failed: a new study's overall growth rate is the main unknown.
This asks how much a small local pilot would help. Each held-out study group's predictions
(primary model, held out; script 28) are split by plating condition: a random pilot share of
conditions recalibrates them, either by one intercept shift (the study's level) or by logistic
recalibration (level and spread), and only the remaining conditions are scored with the H2
ranges. Repeated over random splits.

Outputs: results/phase3/exploratory/local_recalibration_splits.csv (every split),
local_recalibration_by_group.csv and local_recalibration_summary.csv (groups calibrated per
pilot share and method).

Usage: python scripts/31_phase3_local_recalibration.py
"""

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import GROUP, local_recalibration, recalibration_by_group


def summarise(by_group: pd.DataFrame) -> pd.DataFrame:
    """Per pilot share and method: groups whose median O/E, median slope, and both, are in range."""
    g = by_group.groupby(["pilot_fraction", "method"])
    return pd.DataFrame({"groups": g[GROUP].count(),
                         "groups_oe_in_range": g["oe_in_range"].sum(),
                         "groups_slope_in_range": g["slope_in_range"].sum(),
                         "groups_calibrated": g["calibrated"].sum()}).reset_index()


def main() -> None:
    cfg = load_config("phase3")
    crit, lr = cfg["criteria"], cfg["local_recalibration"]
    seed = load_config("thresholds")["resampling"]["seed"]
    df = pd.read_parquet(data_path("interim") / "phase3" / "heldout" / "primary.parquet",
                         columns=[GROUP, "condition_id", "y", "pred_primary"])
    splits, groups = [], []
    for f in lr["pilot_fractions"]:
        s = local_recalibration(df, "pred_primary", f, lr["n_splits"], lr["min_pilot_events"],
                                lr["continuity"], seed).assign(pilot_fraction=f)
        b = recalibration_by_group(s, crit["calibration_slope_range"], crit["oe_range"])
        b["oe_in_range"] = b["median_oe"].between(*crit["oe_range"])
        b["slope_in_range"] = b["median_slope"].between(*crit["calibration_slope_range"])
        splits.append(s)
        groups.append(b.assign(pilot_fraction=f))
        print(f"pilot share {f}: done", flush=True)
    out = data_path("results") / "phase3" / "exploratory"
    out.mkdir(parents=True, exist_ok=True)
    pd.concat(splits).to_csv(out / "local_recalibration_splits.csv", index=False)
    by_group = pd.concat(groups)
    by_group.to_csv(out / "local_recalibration_by_group.csv", index=False)
    summary = summarise(by_group)
    summary.to_csv(out / "local_recalibration_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
