"""Phase 3, step 3: held-out evaluation, the identity-growth curve and recovery shares
(ANALYSIS_PLAN.md sections 5, 6 and 8).

For every variant whose fits are complete (script 27): held-out predictions for a new study
(fixed effects plus the phylum effect, averaged over the study, source and plate-pool random
effects; ANALYSIS_PLAN.md Deviation 1), per-group AUC, Brier, O/E and
calibration slope, calibration tables, the curve with 95% intervals and the new-study spread,
fixed-effect estimates, recovery shares and the registered hypothesis criteria. The binned-rate
and gradient-boosting companions are fitted here (Python). Only held-out predictions enter any
performance number.

Outputs: results/phase3/<variant>/ (tables) and results/phase3/summary_variants.csv;
held-out predictions in data/interim/phase3/heldout/<variant>.parquet.

With --registered, the fits of the model exactly as registered (data/interim/phase3/fits/,
phylum as fixed effects) are evaluated into results/phase3/as_registered/. Their held-out
predictions are fixed effects only, so the same new-study averaging applies; no curve is drawn,
because the registered "average conditions" curve is the artefact Deviation 1 describes.

Usage: python scripts/28_phase3_evaluate.py [variant ...] [--registered]
"""

import json
import sys

import lightgbm as lgb
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score

from plate_truth.calibration import calibration_by_bin, identity_bin, interval, logo_bin_predictions
from plate_truth.config import data_path, load_config
from plate_truth.phase3 import (
    CONTINUOUS,
    FACTORS,
    GROUP,
    VARIANTS,
    averaged_curve,
    bootstrap_mean,
    bootstrap_median,
    calibration_by_decile,
    condition_recovery,
    design_columns,
    group_metrics,
    marginal_probability,
    ratio_with_interval,
    source_recovery,
)
from plate_truth.rfit import RANDOM

GLMM_LOSO = ["primary", "identity_only", "no_identity"]
LGB_FEATURES = ["identity", "log_env", "log_depth", "log2_days", *FACTORS]


def new_study_sigma(info: dict) -> float:
    """SD of the summed random effects that are unknown for a new study (study, source sample,
    plate pool); phylum is known for every row and is not averaged over."""
    return float(np.sqrt(sum(info["sd"][g] ** 2 for g in RANDOM)))


def fit_info(path) -> dict:
    """fit.json of one fit, or None if the fit is missing."""
    f = path / "fit.json"
    return json.loads(f.read_text()) if f.exists() else None


def glmm_heldout(rows: pd.DataFrame, fit_dir, cfg: dict) -> tuple[pd.Series, pd.DataFrame]:
    """Held-out new-study predictions from every fold of one model, and a fold status table.
    Returns (None, status) if any fold is missing or failed."""
    pred = pd.Series(np.nan, index=rows.index)
    status = []
    for g in sorted(rows[GROUP].unique()):
        info = fit_info(fit_dir / g)
        status.append({"fold": g, **({k: info.get(k) for k in
                                      ("status", "converged", "pdHess", "seconds")}
                                     if info else {"status": "missing"})})
        if not info or info["status"] != "ok":
            continue
        eta = pd.read_csv(fit_dir / g / "test_eta.csv")
        sigma = new_study_sigma(info)
        pred[eta["row_id"].to_numpy()] = marginal_probability(
            eta["eta"].to_numpy(), sigma, cfg["marginal"]["gauss_hermite_nodes"])
    status = pd.DataFrame(status)
    return (None if pred.isna().any() else pred), status


def lightgbm_heldout(rows: pd.DataFrame, cfg: dict) -> pd.Series:
    """Leave-one-group-out gradient boosting with the registered fixed settings."""
    x = rows[LGB_FEATURES].copy()
    for f in FACTORS:
        x[f] = x[f].astype("category")
    p = cfg["lightgbm"]
    pred = pd.Series(np.nan, index=rows.index)
    for g in rows[GROUP].unique():
        held = rows[GROUP] == g
        m = lgb.LGBMClassifier(n_estimators=p["n_estimators"], learning_rate=p["learning_rate"],
                               num_leaves=p["num_leaves"],
                               min_child_samples=p["min_child_samples"],
                               subsample=p["subsample"], subsample_freq=p["subsample_freq"],
                               colsample_bytree=p["colsample_bytree"], random_state=p["seed"],
                               verbose=-1)
        m.fit(x[~held], rows.loc[~held, "y"], sample_weight=rows.loc[~held, "w"])
        pred[held] = m.predict_proba(x[held])[:, 1]
    return pred


def summarize(per_group: pd.DataFrame, rows: pd.DataFrame, pred: pd.Series, crit: dict,
              n_boot: int, seed: int) -> dict:
    """Headline numbers for one model: mean per-group AUC with a group-bootstrap interval,
    pooled AUC, and counts of groups meeting the registered criteria."""
    auc = per_group["auc"].to_numpy(float)
    lo, hi = interval(bootstrap_mean(auc, n_boot, seed))
    s_lo, s_hi = crit["calibration_slope_range"]
    o_lo, o_hi = crit["oe_range"]
    calibrated = per_group["slope"].between(s_lo, s_hi) & per_group["oe"].between(o_lo, o_hi)
    return {"groups": len(per_group), "mean_auc": np.nanmean(auc), "mean_auc_lo": lo,
            "mean_auc_hi": hi, "pooled_auc": roc_auc_score(rows["y"], pred),
            "groups_auc_above_0.5": int((auc > 0.5).sum()),
            "median_oe": per_group["oe"].median(), "median_slope": per_group["slope"].median(),
            "groups_calibrated": int(calibrated.sum()),
            "mean_brier": per_group["brier"].mean()}


def curve(fit_dir, rows: pd.DataFrame, centers: dict, cfg: dict, seed: int,
          host: float | None = None) -> pd.DataFrame:
    """The identity curve for a new study, averaged over the analysed rows (Deviation 1): at
    each grid identity every sampled row gets that identity and the new-study predictions are
    averaged. Point estimate over cfg curve.point_rows rows; 95% interval over fixed-effect
    draws (curve.draw_rows rows); 10th-90th percentile spread over study effects. host: for the
    habitat model, the grid half (+1 host-associated, -1 environmental) and rows to use."""
    info = fit_info(fit_dir)
    sd = info["sd"]
    total = new_study_sigma(info)
    within = float(np.sqrt(total ** 2 - sd["study_id"] ** 2))
    c = cfg["curve"]
    grid = pd.read_csv(fit_dir / "grid.csv")
    gx = pd.read_csv(fit_dir / "grid_X.csv")
    # factor coding from all rows, exactly as in the fit; rows are sampled afterwards
    x_all = design_columns(rows, pd.Series(True, index=rows.index), centers)
    x_all["(Intercept)"] = 1.0
    if host is not None:
        keep = (grid["host"] == host).to_numpy()
        grid, gx = grid[keep].reset_index(drop=True), gx[keep].reset_index(drop=True)
        rows = rows[np.where(rows["host_associated"].astype(bool), 1.0, -1.0) == host]
    beta = pd.read_csv(fit_dir / "coef.csv").set_index("term")["estimate"]
    draws = pd.read_csv(fit_dir / "draws.csv")
    ident = [t for t in beta.index if "identity" in t or "no_hit" in t]
    other = [t for t in beta.index if t not in ident]
    phylum = pd.read_csv(fit_dir / "ranef_phylum.csv").set_index("level")["effect"]

    def row_design(n: int, offset: int):
        sample = rows.sample(min(n, len(rows)), random_state=seed + offset)
        return (x_all.loc[sample.index, other].to_numpy(),
                sample["phylum"].map(phylum).fillna(0).to_numpy())

    n = cfg["marginal"]["gauss_hermite_nodes"]
    x, u = row_design(c["point_rows"], 0)
    row_eta = (x @ beta[other].to_numpy() + u)[:, None]
    grid_eta = (gx[ident].to_numpy() @ beta[ident].to_numpy())[:, None]
    out = grid[[k for k in ("identity", "no_hit", "host") if k in grid]].copy()
    out["p"] = averaged_curve(row_eta, grid_eta, total, n)[:, 0]
    z = norm.ppf([0.1, 0.9]) * sd["study_id"]
    out["study_p10"] = averaged_curve(row_eta + z[0], grid_eta, within, n)[:, 0]
    out["study_p90"] = averaged_curve(row_eta + z[1], grid_eta, within, n)[:, 0]
    xd, ud = row_design(c["draw_rows"], 1)
    row_d = xd @ draws[other].to_numpy().T + ud[:, None]
    grid_d = gx[ident].to_numpy() @ draws[ident].to_numpy().T
    pdraw = averaged_curve(row_d, grid_d, total, c["draw_nodes"])
    out["p_lo"], out["p_hi"] = np.quantile(pdraw, [0.025, 0.975], axis=1)
    out["rows_averaged"] = len(x)
    return out


def coefficients(fit_dir) -> pd.DataFrame:
    """Fixed effects as odds ratios with 95% Wald intervals."""
    c = pd.read_csv(fit_dir / "coef.csv")
    c["odds_ratio"] = np.exp(c["estimate"])
    c["or_lo"] = np.exp(c["estimate"] - 1.96 * c["se"])
    c["or_hi"] = np.exp(c["estimate"] + 1.96 * c["se"])
    return c


def level_support(rows: pd.DataFrame) -> pd.DataFrame:
    """Which studies contribute rows to each condition and phylum level."""
    parts = []
    for f in FACTORS:
        g = rows.groupby(f).agg(rows=("y", "size"), events=("y", "sum"),
                                studies=("study_id", lambda s: ";".join(sorted(s.unique()))))
        parts.append(g.reset_index().rename(columns={f: "level"}).assign(factor=f))
    return pd.concat(parts, ignore_index=True)


def recovery(rows: pd.DataFrame, pred: pd.Series, proxy: float, n_boot: int, seed: int
             ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per group: shortcut vs observed (any condition) shares with the ratio and its interval
    over source samples; per condition: observed vs predicted shares."""
    src = source_recovery(rows, proxy)
    out = []
    for g, s in src.groupby(GROUP):
        est, lo, hi = ratio_with_interval(s["proxy_has_relative"].to_numpy(),
                                          s["observed_any_condition"].to_numpy(), n_boot, seed)
        out.append({GROUP: g, "sources": len(s),
                    "proxy_has_relative": s["proxy_has_relative"].mean(),
                    "observed_any_condition": s["observed_any_condition"].mean(),
                    "proxy_over_observed": est, "ratio_lo": lo, "ratio_hi": hi})
    cond = condition_recovery(rows.assign(pred=pred), "pred")
    per_cond = cond.groupby(GROUP).agg(conditions=("observed", "size"),
                                       observed_median=("observed", "median"),
                                       predicted_median=("predicted", "median"))
    return pd.DataFrame(out).merge(per_cond.reset_index(), on=GROUP), cond


def results_root(registered: bool):
    """results/phase3, or results/phase3/as_registered for the model as registered."""
    base = data_path("results") / "phase3"
    return base / "as_registered" if registered else base


def evaluate_variant(name: str, cfg: dict, th: dict, registered: bool = False) -> dict | None:
    """All outputs for one variant; None if its primary-model folds are not complete."""
    root = data_path("interim") / "phase3"
    rows = pd.read_parquet(root / "rows" / f"{name}.parquet")
    fits = root / ("fits" if registered else cfg["fits_dir"]) / name
    n_boot, seed = th["resampling"]["n_bootstrap"], th["resampling"]["seed"]
    out = results_root(registered) / name
    out.mkdir(parents=True, exist_ok=True)

    preds, status = {}, []
    for model in GLMM_LOSO:
        p, st = glmm_heldout(rows, fits / model, cfg)
        status.append(st.assign(model=model))
        if p is not None:
            preds[model] = p
    pd.concat(status).to_csv(out / "fit_status.csv", index=False)
    if "primary" not in preds:
        print(f"{name}: primary model folds incomplete, skipped", flush=True)
        return None
    rows["identity_bin"] = identity_bin(rows["identity"], th["identity"]["bin_edges_percent"])
    preds["bins"] = logo_bin_predictions(rows, GROUP, "identity_bin", "y")
    if name == "primary":
        preds["lightgbm"] = lightgbm_heldout(rows, cfg)

    crit = cfg["criteria"]
    per_group, summary, calib = [], [], []
    for model, p in preds.items():
        rows[f"pred_{model}"] = p
        pg = group_metrics(rows, f"pred_{model}").assign(model=model)
        per_group.append(pg)
        summary.append({"model": model, **summarize(pg, rows, p, crit, n_boot, seed)})
        calib.append(calibration_by_bin(rows, "identity_bin", "y", f"pred_{model}")
                     .reset_index().assign(model=model))
    per_group = pd.concat(per_group, ignore_index=True)
    per_group.to_csv(out / "heldout_by_group.csv", index=False)
    summary = pd.DataFrame(summary)
    summary.to_csv(out / "heldout_summary.csv", index=False)
    pd.concat(calib).to_csv(out / "calibration_by_identity_bin.csv", index=False)
    pd.concat([calibration_by_bin(s, "identity_bin", "y", "pred_primary").reset_index()
               .assign(**{GROUP: g}) for g, s in rows.groupby(GROUP)]).to_csv(
        out / "calibration_by_identity_bin_by_group.csv", index=False)
    calibration_by_decile(rows, "pred_primary").to_csv(out / "calibration_by_decile.csv",
                                                       index=False)
    (root / "heldout").mkdir(exist_ok=True)
    rows.filter(regex="^(study_id|study_group|condition_id|asv_id|y|pred_.*)$").to_parquet(
        root / "heldout" / f"{name}{'__registered' if registered else ''}.parquet", index=False)

    rec, cond = recovery(rows, preds["primary"], th["identity"]["database_proxy_percent"],
                         n_boot, seed)
    rec.to_csv(out / "recovery_by_group.csv", index=False)
    cond.to_csv(out / "recovery_by_condition.csv", index=False)

    centers = {c: float(rows[c].median()) for c in CONTINUOUS}    # as in script 27
    full = fits / "primary" / "all"
    if fit_info(full) and fit_info(full)["status"] == "ok":
        coefficients(full).to_csv(out / "coefficients.csv", index=False)
        if not registered:
            curve(full, rows, centers, cfg, seed).to_csv(out / "curve.csv", index=False)
            pd.read_csv(full / "ranef_phylum.csv").to_csv(out / "phylum_effects.csv",
                                                          index=False)
    if name == "primary":
        level_support(rows).to_csv(out / "level_support.csv", index=False)
        hab = fits / "habitat" / "all"
        if fit_info(hab) and fit_info(hab)["status"] == "ok" and not registered:
            pd.concat([curve(hab, rows, centers, cfg, seed, h) for h in (-1.0, 1.0)]).to_csv(
                out / "curve_by_host.csv", index=False)
            coefficients(hab).to_csv(out / "coefficients_habitat.csv", index=False)
        hyp = hypotheses(per_group, rec, crit, n_boot, seed)
        hyp.to_csv(out / "hypotheses.csv", index=False)
        print(hyp.to_string(index=False))

    head = summary.set_index("model").loc["primary"].to_dict()
    head.update({"variant": name, "rows": len(rows), "events": int(rows["y"].sum()),
                 "median_proxy_over_observed": rec["proxy_over_observed"].median()})
    print(f"{name}:\n{summary.to_string(index=False)}", flush=True)
    return head


def hypotheses(per_group: pd.DataFrame, rec: pd.DataFrame, crit: dict, n_boot: int,
               seed: int) -> pd.DataFrame:
    """The registered criteria (ANALYSIS_PLAN.md section 8), evaluated as written."""
    pg = per_group.set_index(["model", GROUP])
    io = pg.loc["identity_only"]
    lo1, hi1 = interval(bootstrap_mean(io["auc"].to_numpy(float), n_boot, seed))
    n1 = int((io["auc"] > 0.5).sum())
    pr = pg.loc["primary"]
    s_lo, s_hi = crit["calibration_slope_range"]
    o_lo, o_hi = crit["oe_range"]
    n2 = int((pr["slope"].between(s_lo, s_hi) & pr["oe"].between(o_lo, o_hi)).sum())
    ratio = rec.set_index(GROUP)["proxy_over_observed"]
    lo3, hi3 = interval(bootstrap_median(ratio.to_numpy(float), n_boot, seed))
    n3 = int((ratio > 1).sum())
    rows = [
        {"hypothesis": "H1 identity ranks growth", "statistic": f"identity-only AUC > 0.5 in {n1} "
         f"groups; mean AUC {io['auc'].mean():.3f} ({lo1:.3f}-{hi1:.3f})",
         "supported": n1 >= crit["auc_supported_groups"] and lo1 > 0.5},
        {"hypothesis": "H2 calibrated across studies", "statistic": f"{n2} groups with slope "
         f"{s_lo}-{s_hi} and O/E {o_lo}-{o_hi}", "supported": n2 >= crit["calibrated_groups"]},
        {"hypothesis": "H3 shortcut overstates recovery", "statistic": f"ratio > 1 in {n3} "
         f"groups; median {ratio.median():.2f} ({lo3:.2f}-{hi3:.2f})",
         "supported": n3 >= crit["proxy_over_groups"] and lo3 > 1},
    ]
    if "no_identity" in pg.index.get_level_values(0):
        d = (pr["auc"] - pg.loc["no_identity"]["auc"]).to_numpy(float)
        lo4, hi4 = interval(bootstrap_mean(d, n_boot, seed))
        rows.append({"hypothesis": "H4 identity adds to covariates", "statistic":
                     f"mean AUC gain {np.nanmean(d):.3f} ({lo4:.3f}-{hi4:.3f})",
                     "supported": np.nanmean(d) > 0 and lo4 > 0})
    return pd.DataFrame(rows)


def main(names: list[str], registered: bool) -> None:
    cfg, th = load_config("phase3"), load_config("thresholds")
    heads = [h for n in names if (h := evaluate_variant(n, cfg, th, registered)) is not None]
    if heads:
        path = results_root(registered) / "summary_variants.csv"
        new = pd.DataFrame(heads)
        if path.exists():
            old = pd.read_csv(path)
            new = pd.concat([old[~old["variant"].isin(new["variant"])], new], ignore_index=True)
        new.to_csv(path, index=False)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args or list(VARIANTS), "--registered" in sys.argv)
