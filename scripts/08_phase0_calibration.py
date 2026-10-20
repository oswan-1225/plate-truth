"""Phase 0 core result: does nearest-relative identity predict observed growth on plates,
in soils held out from fitting?

Units: (ASV, soil; outcome in any condition) [comparable with the pilot] and
(ASV, plating condition). Outcomes: enriched_on_plates (primary) and detected_on_plates.
Every performance number comes from leave-one-soil-out predictions. Intervals: cluster
bootstrap over soils (held-out predictions fixed, soils resampled).

Outputs (results/phase0/): performance.csv, per_soil_performance.csv, calibration_bins.csv,
recovery_by_soil.csv, recovery_summary.csv, heldout_predictions.csv.gz
"""

import numpy as np
import pandas as pd

from plate_truth.calibration import bootstrap_groups, interval, performance
from plate_truth.config import data_path, load_config
from plate_truth.phase0 import GROUP, attach_identity, evaluate

OUTCOMES = ["enriched_on_plates", "detected_on_plates"]


def load_units(cfg: dict) -> dict[str, pd.DataFrame]:
    """Unit tables with primary identity (all references, id1) attached."""
    proc = data_path("processed")
    nr = pd.read_csv(proc / "zheng2024_nearest_relative.csv")
    edges = cfg["identity"]["bin_edges_percent"]
    units = {"asv_soil": pd.read_csv(proc / "zheng2024_asv_source.csv"),
             "asv_condition": pd.read_csv(proc / "zheng2024_asv_condition.csv")}
    return {k: attach_identity(v, nr, "all", "id1", edges) for k, v in units.items()}


def bootstrap_auc_brier(rows: pd.DataFrame, pred: str, n_boot: int, seed: int) -> dict:
    """95% cluster-bootstrap intervals for held-out AUC and Brier score."""
    auc = bootstrap_groups(rows, GROUP, lambda d: performance(d["y"], d[pred])["auc"],
                           n_boot, seed)
    brier = bootstrap_groups(rows, GROUP, lambda d: performance(d["y"], d[pred])["brier"],
                             n_boot, seed)
    (alo, ahi), (blo, bhi) = interval(auc), interval(brier)
    return {"auc_lo": alo, "auc_hi": ahi, "brier_lo": blo, "brier_hi": bhi}


def bootstrap_bins(rows: pd.DataFrame, n_boot: int, seed: int) -> pd.DataFrame:
    """Cluster-bootstrap 95% interval of the observed rate in each identity bin."""
    cats = rows["identity_bin"].cat.categories
    def rates(d):
        return d.groupby("identity_bin", observed=False)["y"].mean().reindex(cats).to_numpy()
    rng = np.random.default_rng(seed)
    groups = rows[GROUP].unique()
    parts = {g: rows[rows[GROUP] == g] for g in groups}
    draws = np.array([rates(pd.concat([parts[g] for g in rng.choice(groups, len(groups))]))
                      for _ in range(n_boot)])
    lo, hi = np.nanquantile(draws, [0.025, 0.975], axis=0)
    return pd.DataFrame({"boot_lo": lo, "boot_hi": hi}, index=cats)


def recovery_summary(rec: pd.DataFrame, n_boot: int, seed: int) -> dict:
    """Mean across soils of observed / proxy / calibrated recovery, and the proxy
    overestimate factor (mean proxy / mean observed), with soil-bootstrap intervals."""
    r = rec.reset_index().rename(columns={"index": GROUP})
    out = {}
    for col in ["observed_grew", "proxy_has_relative", "predicted_calibrated"]:
        vals = bootstrap_groups(r, GROUP, lambda d, c=col: d[c].mean(), n_boot, seed)
        out[col] = r[col].mean()
        out[f"{col}_lo"], out[f"{col}_hi"] = interval(vals)
    ratio = bootstrap_groups(r, GROUP, lambda d: d["proxy_has_relative"].mean()
                             / d["observed_grew"].mean(), n_boot, seed)
    out["proxy_over_observed"] = r["proxy_has_relative"].mean() / r["observed_grew"].mean()
    out["proxy_over_observed_lo"], out["proxy_over_observed_hi"] = interval(ratio)
    return out


def main() -> None:
    cfg = load_config("thresholds")
    nb, seed = cfg["resampling"]["n_bootstrap"], cfg["resampling"]["seed"]
    proxy = cfg["identity"]["database_proxy_percent"]
    out = data_path("results") / "phase0"
    out.mkdir(parents=True, exist_ok=True)
    perf, per_soil, calib, rec_rows, rec_sum, preds = [], [], [], [], [], []
    for unit, table in load_units(cfg).items():
        for outcome in OUTCOMES:
            res = evaluate(table, outcome, proxy, with_models=True)
            rows = res["rows"]
            tag = {"unit": unit, "outcome": outcome}
            for model, p in res["performance"].items():
                pred = {"bins": "pred_bin", "spline": "pred_spline",
                        "spline_abundance": "pred_spline_abund",
                        "abundance_only": "pred_abund_only"}[model]
                perf.append({**tag, "model": model, **p,
                             **bootstrap_auc_brier(rows, pred, nb, seed)})
            per_soil.append(res["per_source"].assign(**tag))
            c = res["calibration"].join(bootstrap_bins(rows, nb, seed))
            calib.append(c.reset_index().assign(**tag))
            rec_rows.append(res["recovery"].assign(**tag))
            if unit == "asv_soil":
                rec_sum.append({**tag, **recovery_summary(res["recovery"], nb, seed)})
            keep = [GROUP, "asv_id", "identity", "identity_bin", "env_rel", "y", "pred_bin",
                    "pred_spline", "pred_spline_abund"]
            if "condition_id" in rows:
                keep.insert(1, "condition_id")
            preds.append(rows[keep].assign(**tag))
            print(f"{unit} / {outcome}: " + "; ".join(
                f"{m} AUC {p['auc']:.3f} Brier {p['brier']:.4f}"
                for m, p in res["performance"].items()))
    pd.DataFrame(perf).to_csv(out / "performance.csv", index=False)
    pd.concat(per_soil).rename_axis(GROUP).to_csv(out / "per_soil_performance.csv")
    pd.concat(calib).to_csv(out / "calibration_bins.csv", index=False)
    pd.concat(rec_rows).rename_axis(GROUP).to_csv(out / "recovery_by_soil.csv")
    pd.DataFrame(rec_sum).to_csv(out / "recovery_summary.csv", index=False)
    pd.concat(preds).to_csv(out / "heldout_predictions.csv.gz", index=False)
    print(pd.concat(calib)[["unit", "outcome", "identity_bin", "n", "observed", "predicted"]]
          .to_string(index=False))
    print(pd.DataFrame(rec_sum).T.to_string())


if __name__ == "__main__":
    main()
