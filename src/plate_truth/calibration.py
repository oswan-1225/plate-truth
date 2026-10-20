"""Calibrating P(grew | nearest-relative identity) with leave-one-group-out validation.

A 'group' is the unit held out (a soil in Phase 0, a study in Phase 3). Every prediction for
a held-out group is computed only from the other groups, so performance numbers never use
the data a rate was estimated from.
"""

from collections.abc import Callable

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score


def bin_labels(edges: list[float]) -> list[str]:
    """Ordered bin labels for identity edges, e.g. [90, 97] -> ['<90', '90-97', '>=97'],
    preceded by 'none' (no hit in the reference)."""
    labs = [f"<{edges[0]:g}"]
    labs += [f"{lo:g}-{hi:g}" for lo, hi in zip(edges[:-1], edges[1:], strict=False)]
    labs.append(f">={edges[-1]:g}")
    return ["none"] + labs


def identity_bin(identity: pd.Series, edges: list[float]) -> pd.Series:
    """Assign each identity (percent; NaN = no hit) to a bin. Left-closed: 97.0 -> '97-99'."""
    labs = bin_labels(edges)
    cut = pd.cut(identity, [-np.inf, *edges, np.inf], right=False, labels=labs[1:])
    out = cut.astype(object).where(identity.notna(), "none")
    return pd.Categorical(out, categories=labs, ordered=True)


def smoothed_rate(k: float, n: float) -> float:
    """Jeffreys-prior estimate of a proportion, (k + 0.5) / (n + 1).

    Never exactly 0 or 1, so a bin with few or no observations still gives a usable,
    appropriately uncertain prediction (0.5 when n = 0)."""
    return (k + 0.5) / (n + 1.0)


def bin_rates(df: pd.DataFrame, bin_col: str, outcome: str) -> pd.DataFrame:
    """n, k (successes), raw rate and smoothed rate per bin (all categories kept)."""
    g = df.groupby(bin_col, observed=False)[outcome].agg(n="size", k="sum")
    g["rate"] = g["k"] / g["n"].where(g["n"] > 0)
    g["rate_smoothed"] = [smoothed_rate(k, n) for k, n in zip(g["k"], g["n"], strict=False)]
    return g


def wilson_interval(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score 95% interval for a proportion k/n."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (centre - half, centre + half)


def logo_bin_predictions(df: pd.DataFrame, group_col: str, bin_col: str, outcome: str
                         ) -> pd.Series:
    """For each row, the smoothed bin rate estimated from all OTHER groups."""
    pred = pd.Series(np.nan, index=df.index)
    for g in df[group_col].unique():
        held = df[group_col] == g
        rates = bin_rates(df[~held], bin_col, outcome)["rate_smoothed"]
        pred[held] = df.loc[held, bin_col].map(rates).astype(float)
    return pred


def logo_model_predictions(df: pd.DataFrame, group_col: str, outcome: str,
                           fit_predict: Callable[[pd.DataFrame, pd.DataFrame], np.ndarray]
                           ) -> pd.Series:
    """Generic leave-one-group-out: fit_predict(train, test) returns probabilities for test."""
    pred = pd.Series(np.nan, index=df.index)
    for g in df[group_col].unique():
        held = df[group_col] == g
        pred[held] = fit_predict(df[~held], df[held])
    return pred


def performance(y: pd.Series, p: pd.Series) -> dict:
    """AUC (discrimination) and Brier score (mean squared error of probabilities)."""
    y = y.astype(int)
    auc = roc_auc_score(y, p) if y.nunique() == 2 else np.nan
    return {"n": len(y), "events": int(y.sum()), "auc": auc, "brier": brier_score_loss(y, p)}


def calibration_by_bin(df: pd.DataFrame, bin_col: str, outcome: str, pred_col: str
                       ) -> pd.DataFrame:
    """Held-out mean predicted vs observed rate per bin, with Wilson intervals."""
    g = df.groupby(bin_col, observed=True).agg(n=(outcome, "size"), k=(outcome, "sum"),
                                                predicted=(pred_col, "mean"))
    g["observed"] = g["k"] / g["n"]
    ci = [wilson_interval(int(k), int(n)) for k, n in zip(g["k"], g["n"], strict=False)]
    g["observed_lo"], g["observed_hi"] = zip(*ci, strict=False)
    return g


def bootstrap_groups(df: pd.DataFrame, group_col: str, stat: Callable[[pd.DataFrame], float],
                     n_boot: int, seed: int) -> np.ndarray:
    """Cluster bootstrap: resample whole groups with replacement and recompute stat.

    Resampled copies of the same group are relabelled so they count as distinct groups.
    """
    rng = np.random.default_rng(seed)
    groups = df[group_col].unique()
    parts = {g: df[df[group_col] == g] for g in groups}
    out = []
    for _ in range(n_boot):
        pick = rng.choice(groups, size=len(groups), replace=True)
        sample = pd.concat([parts[g].assign(**{group_col: f"{g}#{i}"})
                            for i, g in enumerate(pick)], ignore_index=True)
        out.append(stat(sample))
    return np.asarray(out, dtype=float)


def interval(values: np.ndarray, level: float = 0.95) -> tuple[float, float]:
    """Percentile interval of bootstrap values, ignoring NaN."""
    a = (1 - level) / 2
    return tuple(np.nanquantile(values, [a, 1 - a]))


def weighted_recovery(rel: pd.Series, weight: pd.Series) -> float:
    """Read-weighted fraction: sum(rel * weight) / sum(rel). weight may be 0/1 (observed or
    proxy) or a probability (calibrated prediction)."""
    return float((rel * weight).sum() / rel.sum())
