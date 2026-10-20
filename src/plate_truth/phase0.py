"""Phase 0 analysis steps shared by the main run and the sensitivity grid.

Each function takes explicit parameters (never reads thresholds itself) so a sensitivity run
is the same code with different arguments.
"""

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler

from plate_truth.calibration import (
    calibration_by_bin,
    identity_bin,
    logo_bin_predictions,
    logo_model_predictions,
    performance,
    weighted_recovery,
)
from plate_truth.growth import (
    any_condition,
    asv_condition_table,
    is_detected,
    is_enriched,
    pooled_abundance,
)

GROUP = "source_sample"


def growth_inputs(counts: pd.DataFrame, sheet: pd.DataFrame, env_kind: str,
                  weight: pd.Series | None = None, exclude_asvs: set | None = None
                  ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Pooled environment and plate abundances plus the condition table.

    counts: run, asv_id, count. sheet: run_accession, kind, source_sample, condition_id, ...
    env_kind: 'environment' (original soil) or 'inoculum' (plated filtrate).
    exclude_asvs are removed before relative abundances are computed.
    """
    c = counts.merge(sheet, left_on="run", right_on="run_accession")
    if exclude_asvs:
        c = c[~c["asv_id"].isin(exclude_asvs)]
    env = pooled_abundance(c[c["kind"] == env_kind], GROUP, weight)
    plates = pooled_abundance(c[c["kind"] == "plate"], "condition_id", weight)
    conds = sheet.loc[sheet["kind"] == "plate",
                      ["condition_id", GROUP, "medium", "oxygen"]].drop_duplicates()
    return env, plates, conds


def label_growth(env, plates, conds, env_min_reads: int, grew_min_reads: int,
                 min_fold: float, detected_min_reads: int) -> pd.DataFrame:
    """(ASV, condition) table with enriched_on_plates and detected_on_plates columns."""
    t = asv_condition_table(env, plates, conds, env_min_reads)
    t["enriched_on_plates"] = is_enriched(t["plate_reads"], t["fold"], grew_min_reads, min_fold)
    t["detected_on_plates"] = is_detected(t["plate_reads"], detected_min_reads)
    return t


def per_source(table: pd.DataFrame) -> pd.DataFrame:
    """(ASV, source sample) table: each outcome true if true in any condition."""
    a = any_condition(table, "enriched_on_plates")
    b = any_condition(table, "detected_on_plates")[[GROUP, "asv_id", "detected_on_plates"]]
    return a.merge(b, on=[GROUP, "asv_id"])


def attach_identity(table: pd.DataFrame, nr: pd.DataFrame, reference_set: str, id_col: str,
                    edges: list[float]) -> pd.DataFrame:
    """Add identity (percent, NaN = no hit) and its bin for one reference set/definition."""
    sub = nr.loc[nr["reference_set"] == reference_set, ["asv_id", f"identity_{id_col}"]]
    out = table.merge(sub.rename(columns={f"identity_{id_col}": "identity"}), on="asv_id",
                      how="left")
    out["identity_bin"] = identity_bin(out["identity"], edges)
    return out


def spline_logit(extra: list[str] | None = None):
    """fit_predict for logo_model_predictions: logistic regression on a cubic spline of
    identity (no-hit ASVs get 75, the search floor) plus extra features entering linearly
    after standardising."""
    extra = extra or []

    def fit_predict(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
        def design(d):
            x = d[["identity", *extra]].copy()
            x["identity"] = x["identity"].fillna(75.0)
            return x
        ct = ColumnTransformer([("spline", SplineTransformer(n_knots=5, degree=3), ["identity"])],
                               remainder=StandardScaler())
        model = make_pipeline(ct, LogisticRegression(max_iter=1000))
        model.fit(design(train), train["y"].astype(int))
        return model.predict_proba(design(test))[:, 1]
    return fit_predict


def evaluate(df: pd.DataFrame, outcome: str, proxy_percent: float,
             with_models: bool = True) -> dict:
    """Leave-one-source-out evaluation of P(outcome | identity) on one unit table.

    Returns pooled held-out performance, per-source performance, the held-out calibration
    table and read-weighted recovery per source. Every prediction is held out.
    """
    d = df.copy()
    d["y"] = d[outcome].astype(bool)
    d["pred_bin"] = logo_bin_predictions(d, GROUP, "identity_bin", "y")
    res = {"performance": {"bins": performance(d["y"], d["pred_bin"])}}
    d["log10_env_rel"] = np.log10(d["env_rel"])
    if with_models:
        d["pred_spline"] = logo_model_predictions(d, GROUP, "y", spline_logit())
        d["pred_spline_abund"] = logo_model_predictions(
            d, GROUP, "y", spline_logit(["log10_env_rel"]))
        res["performance"]["spline"] = performance(d["y"], d["pred_spline"])
        res["performance"]["spline_abundance"] = performance(d["y"], d["pred_spline_abund"])
        d["pred_abund_only"] = logo_model_predictions(
            d, GROUP, "y", _abundance_only)
        res["performance"]["abundance_only"] = performance(d["y"], d["pred_abund_only"])
    res["per_source"] = pd.DataFrame(
        {g: performance(s["y"], s["pred_bin"]) for g, s in d.groupby(GROUP)}).T
    res["calibration"] = calibration_by_bin(d, "identity_bin", "y", "pred_bin")
    res["recovery"] = recovery_table(d, proxy_percent)
    res["rows"] = d
    return res


def _abundance_only(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
    """Logistic regression on log10 environmental abundance alone (baseline)."""
    m = LogisticRegression(max_iter=1000).fit(train[["log10_env_rel"]], train["y"].astype(int))
    return m.predict_proba(test[["log10_env_rel"]])[:, 1]


def recovery_table(d: pd.DataFrame, proxy_percent: float) -> pd.DataFrame:
    """Per source sample: read-weighted share of the environment that is
    observed_grew (outcome), proxy_has_relative (identity >= proxy), predicted_calibrated
    (held-out bin probability). Needs columns env_rel, y, identity, pred_bin."""
    rows = {}
    for g, s in d.groupby(GROUP):
        rows[g] = {
            "observed_grew": weighted_recovery(s["env_rel"], s["y"].astype(float)),
            "proxy_has_relative": weighted_recovery(
                s["env_rel"], (s["identity"] >= proxy_percent).astype(float)),
            "predicted_calibrated": weighted_recovery(s["env_rel"], s["pred_bin"]),
            "n_asvs": len(s),
        }
    t = pd.DataFrame(rows).T
    t["proxy_over_observed"] = t["proxy_has_relative"] / t["observed_grew"]
    return t
