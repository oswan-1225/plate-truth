"""Phase 3, step 5: apply the registered pass/fail criteria to every analysis variant
(ANALYSIS_PLAN.md sections 7 and 8).

Reads results/phase3/<variant>/heldout_by_group.csv and recovery_by_group.csv (script 28) and
writes results/phase3/sensitivity_criteria.csv: per variant, the H1-H3 statistics and verdicts,
and whether any verdict differs from the primary analysis ("changes the conclusion" in the
plan's sense). Sensitivity analyses refit only the primary model, so H1 there uses the primary
model's held-out AUC; for the primary variant both that and the registered identity-only H1 are
given. H4 needs the no-identity model and is tested in the primary analysis only.
"""

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import GROUP, VARIANTS, criteria_checks


def variant_checks(name: str, crit: dict, n_boot: int, seed: int) -> list[dict]:
    """Criteria rows for one evaluated variant (two rows for primary: both H1 models)."""
    d = data_path("results") / "phase3" / name
    if not (d / "heldout_by_group.csv").exists():
        return []
    g = pd.read_csv(d / "heldout_by_group.csv")
    rec = pd.read_csv(d / "recovery_by_group.csv").set_index(GROUP)["proxy_over_observed"]
    prim = g[g["model"] == "primary"].set_index(GROUP)
    ratio = rec.reindex(prim.index).to_numpy()
    models = ["primary", "identity_only"] if name == "primary" else ["primary"]
    out = []
    for m in models:
        auc = g[g["model"] == m].set_index(GROUP)["auc"].reindex(prim.index).to_numpy()
        out.append({"variant": name, "h1_model": m,
                    **criteria_checks(auc, prim["oe"], prim["slope"], ratio, crit, n_boot,
                                      seed)})
    return out


def main() -> None:
    crit = load_config("phase3")["criteria"]
    th = load_config("thresholds")["resampling"]
    rows = [r for v in VARIANTS for r in variant_checks(v, crit, th["n_bootstrap"], th["seed"])]
    t = pd.DataFrame(rows)
    ref = t[(t["variant"] == "primary") & (t["h1_model"] == "primary")].iloc[0]
    t["changes_conclusion"] = [any(r[h] != ref[h] for h in ("H1", "H2", "H3"))
                               for _, r in t.iterrows()]
    t.to_csv(data_path("results") / "phase3" / "sensitivity_criteria.csv", index=False)
    pd.set_option("display.width", 220)
    print(t.round(3).to_string(index=False))
    missing = [v for v in VARIANTS if v not in set(t["variant"])]
    if missing:
        print("not evaluated:", ", ".join(missing))


if __name__ == "__main__":
    main()
