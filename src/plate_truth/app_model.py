"""The calibration model as used by the Streamlit app (app/), without R or the pipeline data.

The app loads a small bundle written by scripts/48_app_bundle.py (app/model/): the fixed
effects and their draws, the random-effect SDs, phylum effects, the training condition mix as
quantile points, and the identity spline basis tabulated on a fine grid (so R is not needed to
evaluate it). Predictions follow ANALYSIS_PLAN.md section 22 exactly as in script 43, through
plate_truth.phase5.taxon_probabilities. Everything here is a prediction for plating conditions
like those of the training studies, never an observation.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from plate_truth.phase4 import taxon_label
from plate_truth.phase5 import (
    sequence_eta,
    shortcut_share,
    study_effect,
    taxon_probabilities,
)


def load_bundle(folder: Path) -> dict:
    """Read the model bundle from app/model/."""
    meta = json.loads((folder / "model.json").read_text(encoding="utf-8"))
    grid = pd.read_csv(folder / "spline_grid.csv")
    return {
        "beta": pd.Series(meta["beta"]),
        "draws": pd.read_csv(folder / "draws.csv", index_col=0),
        "sd": meta["sd"], "env_center": meta["env_center"], "settings": meta["settings"],
        "phylum": pd.Series(meta["phylum_effects"]),
        "points": np.asarray(meta["condition_points"], float),
        "points_draws": np.loadtxt(folder / "condition_points_draws.csv", delimiter=",",
                                   ndmin=2),
        "grid_identity": grid["identity"].to_numpy(float),
        "grid_basis": grid.drop(columns="identity").to_numpy(float),
        "family_calls": pd.read_csv(folder / "family_calls.csv").set_index("family")["call"],
    }


def basis_from_grid(identity: np.ndarray, grid_identity: np.ndarray,
                    grid_basis: np.ndarray) -> np.ndarray:
    """Spline basis at each identity by linear interpolation of the tabulated basis (one
    column at a time); values outside the grid are clipped to its ends."""
    x = np.clip(np.asarray(identity, float), grid_identity[0], grid_identity[-1])
    return np.column_stack([np.interp(x, grid_identity, grid_basis[:, j])
                            for j in range(grid_basis.shape[1])])


def prepare_table(df: pd.DataFrame, present_min_reads: int, no_hit_below: float
                  ) -> tuple[pd.DataFrame, list[str]]:
    """Check and normalise an uploaded table (one row per ASV).

    Required: asv_id, count (every prokaryotic ASV of one sample). Optional: identity (percent;
    empty = no cultured relative), phylum, family. Returns the rows present in the sample
    (count >= present_min_reads) with relative abundance a among ALL the sample's reads (as
    in the surveys), no_hit and identity_model, and a list of notes for the user.
    Identities given as fractions (all <= 1) are converted to percent; identities below
    no_hit_below count as no hit, as in the training data."""
    notes = []
    cols = {c.lower().strip(): c for c in df.columns}
    missing = [c for c in ("asv_id", "count") if c not in cols]
    if missing:
        raise ValueError(f"missing column(s): {', '.join(missing)}")
    t = df.rename(columns={v: k for k, v in cols.items()}).copy()
    t["count"] = pd.to_numeric(t["count"], errors="coerce").fillna(0)
    if "identity" not in t:
        t["identity"] = np.nan
        notes.append("No identity column: every ASV is treated as having no cultured relative.")
    t["identity"] = pd.to_numeric(t["identity"], errors="coerce")
    if t["identity"].notna().any() and t["identity"].max() <= 1:
        t["identity"] = 100 * t["identity"]
        notes.append("Identities were between 0 and 1 and were read as fractions (x 100).")
    low = t["identity"] < no_hit_below
    if low.any():
        notes.append(f"{int(low.sum())} ASVs below {no_hit_below:g}% identity count as no hit.")
        t.loc[low, "identity"] = np.nan
    t.attrs["phylum_given"] = "phylum" in t
    for c in ("phylum", "family"):
        if c not in t:
            t[c] = np.nan
    if not t.attrs["phylum_given"]:
        notes.append("No phylum column: the phylum effect is set to its average (0), which "
                     "widens the uncertainty that cannot be shown here.")
    total = t["count"].sum()
    t["a"] = t["count"] / total if total > 0 else 0.0
    given = t.attrs["phylum_given"]
    t = t[t["count"] >= present_min_reads].copy()
    t.attrs["phylum_given"] = given
    if t.empty:
        raise ValueError(f"no ASV has at least {present_min_reads} reads")
    t["no_hit"] = t["identity"].isna().astype(float)
    t["identity_model"] = t["identity"].fillna(no_hit_below)
    t = t.reset_index(drop=True)
    t.attrs["phylum_given"] = given
    return t, notes


def phylum_effects(phylum: pd.Series, effects: pd.Series, levels: list[str],
                   given: bool) -> np.ndarray:
    """Phylum random effect per ASV: the model's level for the modelled phyla and 'other' for
    any other phylum or an unassigned one (as in training). If no phylum column was given at
    all, every ASV gets 0, the average phylum."""
    if not given:
        return np.zeros(len(phylum))
    lev = phylum.where(phylum.isin(levels), "other")
    return np.array(lev.map(effects), dtype=float)


def predict(t: pd.DataFrame, b: dict, n_draws: int = 100) -> pd.DataFrame:
    """Per-ASV predictions: probability of enrichment on one plate in a low-, median- and
    high-recovery study (study effect at its 10th, 50th, 90th percentile), on any of K plates
    in a high-recovery study, and a 95% interval for the median-study value from fixed-effect
    draws. Adds the shortcut flag (identity >= threshold)."""
    s = b["settings"]
    basis = basis_from_grid(t["identity_model"], b["grid_identity"], b["grid_basis"])
    phy = phylum_effects(t["phylum"], b["phylum"], s["phyla"], t.attrs.get("phylum_given", True))
    args = (basis, t["no_hit"].to_numpy(), np.log10(t["a"].to_numpy()), phy)
    eta = sequence_eta(*args, b["beta"], b["env_center"])
    sd, k = b["sd"], s["conditions_k"]
    out = t[["asv_id", "count", "a", "identity", "phylum", "family"]].copy()
    out["has_cultured_relative"] = t["identity"].fillna(0).to_numpy() >= s["shortcut_percent"]
    for q in (0.1, 0.5, 0.9):
        p1, pk = taxon_probabilities(eta, b["points"], study_effect(q, sd["study_id"]),
                                     sd["source_key"], k, s["source_nodes"])
        out[f"predicted_one_plate_p{int(q * 100)}"] = p1
        if q == 0.9:
            out[f"predicted_any_of_{k}_p90"] = pk
    draws = b["draws"].iloc[:, :n_draws]
    eta_d = sequence_eta(*args, draws, b["env_center"])
    u = study_effect(0.5, sd["study_id"])
    pd_ = np.column_stack([taxon_probabilities(eta_d[:, d], b["points_draws"][d], u,
                                               sd["source_key"], k, s["source_nodes_draws"])[0]
                           for d in range(eta_d.shape[1])])
    out["predicted_one_plate_p50_lo"] = np.quantile(pd_, 0.025, axis=1)
    out["predicted_one_plate_p50_hi"] = np.quantile(pd_, 0.975, axis=1)
    return out


def community_summary(pred: pd.DataFrame, k: int, threshold: float) -> dict:
    """Read-weighted shares for the whole sample: the shortcut (identity >= threshold), and
    the predicted share in taxa enriched on one plate (low-, median- and high-recovery study)
    and on any of k plates (high-recovery study)."""
    a = pred["a"].to_numpy()

    def share(col: str) -> float:
        return float(a @ pred[col].to_numpy() / a.sum())

    return {"shortcut": shortcut_share(pred["identity"].to_numpy(), a, threshold),
            "one_plate_p10": share("predicted_one_plate_p10"),
            "one_plate_p50": share("predicted_one_plate_p50"),
            "one_plate_p90": share("predicted_one_plate_p90"),
            f"any_of_{k}_p90": share(f"predicted_any_of_{k}_p90")}


def family_composition(pred: pd.DataFrame, calls: pd.Series) -> pd.Series:
    """Read share in families called under-recovered, over-recovered or neither in the
    training studies, and in families not evaluated. Families are labelled as in Phase 4, so
    an ASV without a family belongs to 'unassigned family (<phylum>)'."""
    call = taxon_label(pred, "family").map(calls).fillna("not evaluated")
    return (pred.groupby(call)["a"].sum() / pred["a"].sum()).reindex(
        ["under", "neither", "over", "not evaluated"], fill_value=0.0)
