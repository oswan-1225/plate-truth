"""Phase 3: the cross-study calibration model (ANALYSIS_PLAN.md).

Builds the (ASV, plating condition) analysis rows from the database, codes the model features,
turns fixed-effect predictions into predictions for an unseen study, and computes held-out
metrics. The mixed model itself is fitted in R (src/R/fit_glmm.R, run by plate_truth.rfit);
everything that transforms data or computes a metric lives here so that it is tested in Python.
"""

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit, logit
from sklearn.metrics import brier_score_loss, roc_auc_score

from plate_truth.growth import pooled_abundance
from plate_truth.phase0 import label_growth
from plate_truth.qc import pool_depth, usable_conditions
from plate_truth.queries import (
    asv_taxonomy,
    environment_counts,
    nearest_relative,
    plate_counts,
    plating_conditions,
)

UNKNOWN = "unknown"
GROUP = "study_group"
OXYGEN_LEVELS = ["aerobic", "anaerobic", "microaerophilic"]
FACTORS = ["medium", "temperature", "oxygen", "phylum"]
CODED = ["medium", "temperature", "oxygen"]        # sum-to-zero fixed effects; phylum is a
                                                    # random intercept (Deviation 1)
CONTINUOUS = ["log_env", "log_depth", "log2_days"]

# The primary analysis; each sensitivity analysis overrides some of these (ANALYSIS_PLAN.md s7).
PRIMARY_VARIANT = {
    "outcome": "enriched_on_plates",
    "grew_min_reads": None,          # None = value from config/thresholds.yaml
    "min_fold": None,
    "identity_def": "id1",
    "reference_set": "all",
    "min_pool_reads": None,
    "every_read_studies": ["zheng2024_soil", "zheng2024_sludge"],
    "copy_number": False,
    "equal_group_weight": False,
    "drop_tiers": [],
    "db_tag": None,                  # a tagged database (scripts 16-19 --tag) instead of the primary
}

# Sensitivity analyses in the order of ANALYSIS_PLAN.md section 7 (fits run in this order).
# verhoeven_maxee2 reads data/interim/plate_truth__maxee2.sqlite, built by scripts 16-19 with
# --tag=maxee2 (verhoeven2023_sludge reprocessed with --max-ee=2).
VARIANTS = {
    "primary": {},
    "detected": {"outcome": "detected_on_plates"},
    "reads2": {"grew_min_reads": 2},
    "reads10": {"grew_min_reads": 10},
    "fold1.5": {"min_fold": 1.5},
    "fold5": {"min_fold": 5.0},
    "iddef2": {"identity_def": "id2"},
    "idmid": {"identity_def": "mid"},
    "ref_bacdive": {"reference_set": "bacdive_strain"},
    "ref_refseq": {"reference_set": "refseq_type"},
    "depth5000": {"min_pool_reads": 5000},
    "depth7500": {"min_pool_reads": 7500},
    "zheng_paired": {"every_read_studies": []},
    "verhoeven_maxee2": {"db_tag": "maxee2"},
    "copy_number": {"copy_number": True},
    "equal_weight": {"equal_group_weight": True},
    "no_tier_c": {"drop_tiers": ["C"]},
}


# --- feature coding --------------------------------------------------------------------------

def temperature_band(temp: pd.Series, cold_max: float, warm_min: float) -> pd.Series:
    """'cold' (<= cold_max), 'warm' (>= warm_min), 'moderate' between, 'unknown' if missing."""
    out = pd.Series("moderate", index=temp.index, dtype=object)
    out[temp <= cold_max] = "cold"
    out[temp >= warm_min] = "warm"
    out[temp.isna()] = UNKNOWN
    return out


def keep_levels(values: pd.Series, levels: list[str], other: str) -> pd.Series:
    """Values in levels are kept; everything else (including missing) becomes other."""
    return values.where(values.isin(levels), other).astype(object)


def known_levels(values: pd.Series) -> list[str]:
    """Sorted levels present in values, excluding 'unknown'."""
    return sorted(v for v in values.dropna().unique() if v != UNKNOWN)


def effect_code(values: pd.Series, levels: list[str], prefix: str) -> pd.DataFrame:
    """Sum-to-zero (deviation) coding: one column per level except the last.

    Level j (not last) -> 1 in column j; the last level -> -1 in every column; any value not in
    levels ('unknown', or a level unseen in training) -> 0 everywhere, i.e. the unweighted mean
    of the level effects. One level or none -> no columns.
    """
    cols = {}
    for lev in levels[:-1]:
        col = (values == lev).astype(float)
        col[values == levels[-1]] = -1.0
        cols[f"{prefix}_{_safe(lev)}"] = col
    return pd.DataFrame(cols, index=values.index)


def _safe(level: str) -> str:
    """Column-name-safe version of a level label."""
    return "".join(c if c.isalnum() else "_" for c in str(level))


def design_columns(rows: pd.DataFrame, train: pd.Series, centers: dict) -> pd.DataFrame:
    """Numeric model inputs: identity, no_hit, centred continuous covariates and effect-coded
    condition factors, with factor levels taken from the training rows only."""
    out = pd.DataFrame({"identity": rows["identity_model"], "no_hit": rows["no_hit"]},
                       index=rows.index)
    for c in CONTINUOUS:
        out[c] = rows[c] - centers[c]
    out["host"] = np.where(rows["host_associated"].astype(bool), 1.0, -1.0)
    for f in CODED:
        out = out.join(effect_code(rows[f], known_levels(rows.loc[train, f]), f))
    return out


# --- analysis rows -----------------------------------------------------------------------------

def mixed_counts(query, con, every_read_studies: list[str]) -> pd.DataFrame:
    """Counts from the paired table, with the listed studies taken from the every-read table."""
    base = query(con)
    if not every_read_studies:
        return base
    alt = query(con, table="abundances_every_read")
    return pd.concat([base[~base["study_id"].isin(every_read_studies)],
                      alt[alt["study_id"].isin(every_read_studies)]], ignore_index=True)


def analysis_rows(con, th: dict, cfg: dict, variant: dict) -> pd.DataFrame:
    """One row per (environment-present ASV, usable plating condition) with outcome y and every
    model feature. variant: PRIMARY_VARIANT keys (None = take the threshold from th)."""
    v = {**PRIMARY_VARIANT, **variant}
    grew_reads = v["grew_min_reads"] or th["grew"]["min_reads"]
    fold = v["min_fold"] or th["grew"]["min_fold_enrichment"]
    depth_min = v["min_pool_reads"] or th["depth"]["min_pool_reads"]

    env = mixed_counts(environment_counts, con, v["every_read_studies"])
    plates = mixed_counts(plate_counts, con, v["every_read_studies"])
    conds = plating_conditions(con, linked_only=True)
    conds = conds[~conds["tier"].isin(v["drop_tiers"])]
    flags = usable_conditions(conds[["study_id", "source_sample", "condition_id"]],
                              pool_depth(plates, ["study_id", "condition_id"]),
                              pool_depth(env, ["study_id", "source_sample"]), depth_min)
    conds = conds.merge(flags[["study_id", "condition_id", "plate_depth", "usable"]],
                        on=["study_id", "condition_id"])
    conds = conds[conds["usable"]]

    tax = asv_taxonomy(con).set_index("asv_id")
    weight = tax["copy_number"] if v["copy_number"] else None
    parts = []
    for sid, c in conds.groupby("study_id"):
        e = pooled_abundance(env[env["study_id"] == sid], "source_sample", weight)
        p = pooled_abundance(plates[plates["study_id"] == sid], "condition_id", weight)
        t = label_growth(e, p, c[["condition_id", "source_sample", "medium", "oxygen"]],
                         th["environment"]["min_reads"], grew_reads, fold,
                         th["detected"]["min_reads"])
        parts.append(t.drop(columns=["medium", "oxygen"]).assign(study_id=sid))
    rows = pd.concat(parts, ignore_index=True)
    rows = rows.merge(conds.drop(columns=["source_sample", "usable"]),
                      on=["study_id", "condition_id"])
    return add_features(rows, nearest_relative(con, v["reference_set"]), tax, cfg, v)


def add_features(rows: pd.DataFrame, nr: pd.DataFrame, tax: pd.DataFrame, cfg: dict,
                 v: dict) -> pd.DataFrame:
    """Attach identity, taxonomy and condition features to labelled rows."""
    ident = nr.set_index("asv_id")[f"identity_{v['identity_def']}"]
    out = rows.rename(columns={"medium": "medium_label"})
    out["identity"] = out["asv_id"].map(ident)
    out["no_hit"] = out["identity"].isna().astype(float)
    out["identity_model"] = out["identity"].fillna(cfg["identity_spline"]["no_hit_value"])
    out["log_env"] = np.log10(out["env_rel"])
    out["log_depth"] = np.log10(out["plate_depth"].astype(float))
    out["log2_days"] = np.log2(out["incubation_days"].astype(float))
    tb = cfg["temperature_bands"]
    out["temperature"] = temperature_band(out["temperature_c"], tb["cold_max"], tb["warm_min"])
    out["oxygen"] = keep_levels(out["oxygen"], OXYGEN_LEVELS, UNKNOWN)
    out["medium"] = keep_levels(out["medium_class"], cfg["medium_classes"], UNKNOWN)
    out["phylum"] = keep_levels(out["asv_id"].map(tax["phylum"]), cfg["phyla"], "other")
    out["source_key"] = out["study_id"] + "|" + out["source_sample"].astype(str)
    out["y"] = out[v["outcome"]].astype(int)
    out["w"] = group_weights(out[GROUP]) if v["equal_group_weight"] else 1.0
    return out.reset_index(drop=True)


def group_weights(groups: pd.Series) -> pd.Series:
    """Row weights giving every group the same total weight, with mean weight 1."""
    n = groups.map(groups.value_counts())
    return (len(groups) / groups.nunique()) / n


# --- predictions for an unseen study -------------------------------------------------------------

def marginal_probability(eta: np.ndarray, sigma: float, n_nodes: int) -> np.ndarray:
    """E[expit(eta + b)] for b ~ Normal(0, sigma^2), by Gauss-Hermite quadrature.

    The expected rate in a new study whose random effects are unknown; sigma is the square root
    of the summed random-effect variances. sigma = 0 gives expit(eta)."""
    eta = np.asarray(eta, dtype=float)
    if sigma == 0:
        return expit(eta)
    x, w = np.polynomial.hermite_e.hermegauss(n_nodes)   # weight exp(-x^2/2)
    w = w / w.sum()
    return expit(eta[..., None] + sigma * x) @ w


def averaged_curve(row_eta: np.ndarray, grid_eta: np.ndarray, sigma: float, n_nodes: int
                   ) -> np.ndarray:
    """Expected enrichment rate at each grid point, averaged over a set of rows.

    row_eta: (rows, draws) linear predictor of each row WITHOUT its identity terms (intercept,
    covariates, phylum effect). grid_eta: (grid, draws) contribution of the identity terms at
    each grid point. Returns (grid, draws): the mean over rows of the new-study prediction
    marginal_probability(row_eta + grid_eta, sigma) (ANALYSIS_PLAN.md Deviation 1)."""
    row_eta = np.atleast_2d(np.asarray(row_eta, dtype=float))
    grid_eta = np.atleast_2d(np.asarray(grid_eta, dtype=float))
    out = np.empty(grid_eta.shape)
    for g in range(grid_eta.shape[0]):
        out[g] = marginal_probability(row_eta + grid_eta[g], sigma, n_nodes).mean(axis=0)
    return out


# --- metrics ---------------------------------------------------------------------------------

def calibration_metrics(y: pd.Series, p: pd.Series) -> dict:
    """Observed/expected ratio and calibration slope (logistic regression of y on logit(p)).

    O/E = sum(y) / sum(p): 1 means the overall level is right. Slope 1 means the spread of the
    predictions is right; below 1 they are too extreme."""
    y = np.asarray(y, dtype=float)
    p = np.clip(np.asarray(p, dtype=float), 1e-9, 1 - 1e-9)
    oe = y.sum() / p.sum()
    if y.min() == y.max() or np.ptp(p) == 0:     # slope undefined
        return {"oe": oe, "slope": np.nan, "slope_lo": np.nan, "slope_hi": np.nan}
    fit = sm.GLM(y, sm.add_constant(logit(p), has_constant="add"),
                 family=sm.families.Binomial()).fit()
    lo, hi = fit.conf_int()[1]
    return {"oe": oe, "slope": fit.params[1], "slope_lo": lo, "slope_hi": hi}


def group_metrics(df: pd.DataFrame, pred: str, group: str = GROUP) -> pd.DataFrame:
    """Per held-out group: rows, events, AUC, Brier, O/E and calibration slope."""
    out = {}
    for g, s in df.groupby(group):
        y = s["y"].astype(int)
        out[g] = {"n": len(s), "events": int(y.sum()),
                  "auc": roc_auc_score(y, s[pred]) if y.nunique() == 2 else np.nan,
                  "brier": brier_score_loss(y, s[pred]),
                  **calibration_metrics(y, s[pred])}
    return pd.DataFrame(out).T.rename_axis(group).reset_index()


def bootstrap_mean(values: np.ndarray, n_boot: int, seed: int) -> np.ndarray:
    """Bootstrap distribution of the mean of per-group values (groups resampled with
    replacement; NaN values ignored)."""
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    rng = np.random.default_rng(seed)
    return v[rng.integers(0, len(v), size=(n_boot, len(v)))].mean(axis=1)


def bootstrap_median(values: np.ndarray, n_boot: int, seed: int) -> np.ndarray:
    """As bootstrap_mean, for the median."""
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    rng = np.random.default_rng(seed)
    return np.median(v[rng.integers(0, len(v), size=(n_boot, len(v)))], axis=1)


def calibration_by_decile(df: pd.DataFrame, pred: str, n_bins: int = 10) -> pd.DataFrame:
    """Mean predicted vs observed rate within quantile bins of the prediction."""
    q = pd.qcut(df[pred].rank(method="first"), n_bins, labels=False)
    return df.groupby(q).agg(n=("y", "size"), observed=("y", "mean"),
                             predicted=(pred, "mean")).rename_axis("decile").reset_index()


# --- recovery shares -----------------------------------------------------------------------------

def source_recovery(rows: pd.DataFrame, proxy_percent: float) -> pd.DataFrame:
    """Per (study, source sample): read-weighted shares of the environment that have a cultured
    relative (proxy_has_relative) and that were enriched in at least one condition
    (observed_any_condition)."""
    a = rows.groupby(["study_id", GROUP, "source_key", "asv_id"], as_index=False).agg(
        env_rel=("env_rel", "first"), identity=("identity", "first"), y=("y", "max"))
    a["proxy"] = (a["identity"] >= proxy_percent).astype(float)
    a["w_proxy"] = a["env_rel"] * a["proxy"]
    a["w_obs"] = a["env_rel"] * a["y"]
    g = a.groupby(["study_id", GROUP, "source_key"])[["env_rel", "w_proxy", "w_obs"]].sum()
    return pd.DataFrame({"proxy_has_relative": g["w_proxy"] / g["env_rel"],
                         "observed_any_condition": g["w_obs"] / g["env_rel"]}).reset_index()


def condition_recovery(rows: pd.DataFrame, pred: str) -> pd.DataFrame:
    """Per plating condition: read-weighted observed and predicted shares of the environment."""
    r = rows.assign(w_obs=rows["env_rel"] * rows["y"], w_pred=rows["env_rel"] * rows[pred])
    g = r.groupby(["study_id", GROUP, "condition_id"])[["env_rel", "w_obs", "w_pred"]].sum()
    return pd.DataFrame({"observed": g["w_obs"] / g["env_rel"],
                         "predicted": g["w_pred"] / g["env_rel"]}).reset_index()


def ratio_with_interval(proxy: np.ndarray, observed: np.ndarray, n_boot: int, seed: int
                        ) -> tuple[float, float, float]:
    """mean(proxy) / mean(observed) over source samples, with a percentile bootstrap interval
    over source samples (NaN interval for a single source)."""
    proxy, observed = np.asarray(proxy, float), np.asarray(observed, float)
    est = proxy.mean() / observed.mean()
    if len(proxy) < 2:
        return est, np.nan, np.nan
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(proxy), size=(n_boot, len(proxy)))
    with np.errstate(divide="ignore", invalid="ignore"):
        boots = proxy[idx].mean(axis=1) / observed[idx].mean(axis=1)
    lo, hi = np.nanquantile(boots, [0.025, 0.975])
    return est, lo, hi


# --- registered criteria -------------------------------------------------------------------------

def criteria_checks(auc: np.ndarray, oe: np.ndarray, slope: np.ndarray, ratio: np.ndarray,
                    crit: dict, n_boot: int, seed: int) -> dict:
    """The registered pass/fail rules (ANALYSIS_PLAN.md section 8) from per-group values.

    auc: held-out AUC per group (identity-only model for H1 as registered; the primary model in
    sensitivity analyses, which refit only the primary model). oe, slope: primary model.
    ratio: shortcut / observed recovery per group. Group-bootstrap intervals as in the plan."""
    auc, oe, slope, ratio = (np.asarray(a, dtype=float) for a in (auc, oe, slope, ratio))
    auc_lo = np.nanquantile(bootstrap_mean(auc, n_boot, seed), 0.025)
    s_lo, s_hi = crit["calibration_slope_range"]
    o_lo, o_hi = crit["oe_range"]
    n_cal = int(((slope >= s_lo) & (slope <= s_hi) & (oe >= o_lo) & (oe <= o_hi)).sum())
    ratio_lo = np.nanquantile(bootstrap_median(ratio, n_boot, seed), 0.025)
    n_auc, n_ratio = int((auc > 0.5).sum()), int((ratio > 1).sum())
    return {"groups_auc_above_0.5": n_auc, "mean_auc_lo": auc_lo,
            "H1": bool(n_auc >= crit["auc_supported_groups"] and auc_lo > 0.5),
            "groups_calibrated": n_cal, "H2": bool(n_cal >= crit["calibrated_groups"]),
            "groups_ratio_above_1": n_ratio, "median_ratio": float(np.nanmedian(ratio)),
            "median_ratio_lo": ratio_lo,
            "H3": bool(n_ratio >= crit["proxy_over_groups"] and ratio_lo > 1)}

