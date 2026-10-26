"""Phase 5, step 4: the Phase 5 procedure applied to the training studies (ANALYSIS_PLAN.md
section 23).

Each training study group is treated as a survey and predicted by the Phase 3 primary model
fitted without it (leave-one-group-out fits): its source samples' taxa (abundance, identity,
phylum) as the survey, the other nine groups' plating conditions as the condition mix (each
group weighted equally), the fold's own random-effect SDs. The predicted range of recovery
(one plate; any of 8 plates; any of the group's own number of plates) is compared with what
grew on the group's plates.

Output: results/phase5/training_check.csv (one row per group) and
results/phase5/training_check_criteria.csv.
"""

import numpy as np
import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import CONTINUOUS, GROUP, condition_recovery, design_columns, \
    source_recovery
from plate_truth.phase5 import (
    condition_columns,
    condition_effects,
    condition_weights,
    mixture_points,
    recovery_shares,
    sequence_eta,
    spline_basis,
    study_effect,
)


def load_fit(fit_dir) -> dict:
    """Coefficients, phylum effects and random-effect SDs of one fit."""
    import json
    return {"beta": pd.read_csv(fit_dir / "coef.csv").set_index("term")["estimate"],
            "phylum": pd.read_csv(fit_dir / "ranef_phylum.csv").set_index("level")["effect"],
            "sd": json.loads((fit_dir / "fit.json").read_text())["sd"]}


def condition_points(rows: pd.DataFrame, train: pd.Series, centers: dict, fit: dict,
                     n_points: int) -> np.ndarray:
    """Equally weighted points for (condition contribution + condition effect) over the
    training conditions, each study group weighted equally."""
    cond = rows[train].drop_duplicates("condition_id")
    x = design_columns(cond, pd.Series(True, index=cond.index), centers)
    cols = condition_columns(x.columns)
    missing = [c for c in cols if c not in fit["beta"].index]
    if missing:
        raise ValueError(f"design columns without a coefficient: {missing}")
    g = condition_effects(x[cols], fit["beta"])
    return mixture_points(g, condition_weights(cond[GROUP]).to_numpy(),
                          fit["sd"]["condition_id"], n_points)


def group_check(rows: pd.DataFrame, group: str, fit: dict, centers: dict, basis: dict,
                cfg: dict, proxy: float) -> dict:
    m = cfg["model"]
    held = rows[GROUP] == group
    pts = condition_points(rows, ~held, centers, fit, m["condition_points"])
    r = rows[held]
    k_own = int(np.floor(r.groupby("source_key")["condition_id"].nunique().median()))
    seq = r.drop_duplicates(["source_key", "asv_id"]).copy()
    eta = sequence_eta(np.vstack(seq["identity_model"].map(basis).to_numpy()),
                       seq["no_hit"].to_numpy(), seq["log_env"].to_numpy(),
                       seq["phylum"].map(fit["phylum"]).fillna(0.0).to_numpy(),
                       fit["beta"], centers["log_env"])
    seq["eta"] = eta
    out = {GROUP: group, "habitat": " + ".join(sorted(r["habitat"].unique())),
           "host_associated": bool(r["host_associated"].iloc[0]),
           "sources": r["source_key"].nunique(), "conditions": r["condition_id"].nunique(),
           "k_own": k_own}
    per_q = {}
    for q in m["study_quantiles"]:
        u = study_effect(q, fit["sd"]["study_id"])
        vals = []
        for _, s in seq.groupby("source_key"):
            a, e = s["env_rel"].to_numpy(), s["eta"].to_numpy()
            r1, r8 = recovery_shares(e, a, pts, u, fit["sd"]["source_key"], m["conditions_k"],
                                     m["source_nodes"])
            _, rown = recovery_shares(e, a, pts, u, fit["sd"]["source_key"], k_own,
                                      m["source_nodes"])
            vals.append((r1, r8, rown))
        per_q[q] = np.mean(vals, axis=0)
        tag = f"p{int(round(q * 100))}"
        out[f"r1_{tag}"], out[f"r8_{tag}"], out[f"rown_{tag}"] = per_q[q]
    cond = condition_recovery(r.assign(pred=0.0), "pred")
    cond["source_key"] = cond["condition_id"].map(r.drop_duplicates("condition_id")
                                                  .set_index("condition_id")["source_key"])
    out["observed_per_condition"] = cond.groupby("source_key")["observed"].median().mean()
    src = source_recovery(r, proxy)
    out["observed_any_condition"] = src["observed_any_condition"].mean()
    out["shortcut"] = src["proxy_has_relative"].mean()
    return out


def main() -> None:
    cfg = load_config("phase5")
    cfg3 = load_config("phase3")
    proxy = load_config("thresholds")["identity"]["database_proxy_percent"]
    rows = pd.read_parquet(data_path("interim") / "phase3" / "rows" / "primary.parquet")
    centers = {c: float(rows[c].median()) for c in CONTINUOUS}
    ids = np.unique(rows["identity_model"].to_numpy())
    basis = dict(zip(ids, spline_basis(ids, cfg3)))
    fits = data_path("interim") / "phase3" / cfg3["fits_dir"] / "primary" / "primary"
    out = []
    for g in sorted(rows[GROUP].unique()):
        res = group_check(rows, g, load_fit(fits / g), centers, basis, cfg, proxy)
        out.append(res)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items()},
              flush=True)
    t = pd.DataFrame(out)
    t["inside_range"] = t["observed_per_condition"].between(t["r1_p10"], t["r1_p90"])
    t["below_upper"] = t["observed_per_condition"] <= t["r1_p90"]
    t["any_below_upper_own_k"] = t["observed_any_condition"] <= t["rown_p90"]
    t["any_below_upper_k8"] = t["observed_any_condition"] <= t["r8_p90"]
    t["shortcut_above_upper"] = t["shortcut"] > t["r1_p90"]
    t["shortcut_above_upper_k8"] = t["shortcut"] > t["r8_p90"]
    res = data_path("results") / "phase5"
    res.mkdir(parents=True, exist_ok=True)
    t.to_csv(res / "training_check.csv", index=False)
    c = cfg["criteria"]
    crit = pd.DataFrame([
        {"check": "range describes the spread (observed per-condition inside R1 range)",
         "groups": int(t["inside_range"].sum()), "needed": c["check_coverage_groups"]},
        {"check": "upper end (observed per-condition <= R1 90th)",
         "groups": int(t["below_upper"].sum()), "needed": c["check_upper_groups"]},
        {"check": "generous upper end (observed any-condition <= any-of-own-K 90th)",
         "groups": int(t["any_below_upper_own_k"].sum()), "needed": c["check_upper_groups"]},
    ])
    crit["of"] = len(t)
    crit["passed"] = crit["groups"] >= crit["needed"]
    crit.to_csv(res / "training_check_criteria.csv", index=False)
    print(t.round(4).to_string(index=False))
    print(crit.to_string(index=False))


if __name__ == "__main__":
    main()
