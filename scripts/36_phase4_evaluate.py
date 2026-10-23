"""Phase 4, step 4: evaluate the fits (ANALYSIS_PLAN.md sections 13-16).

Q5: likelihood-ratio test of family x medium class (H5a), held-out log loss with and without it
    (H5b), and the family- and phylum-by-class effect tables.
Q6: odds ratio for the match term, likelihood-ratio test, held-out log loss with and without
    match (H6), the sensitivity analyses and the matched-plate companion statistic.
Q7: per-taxon observed/expected against the relatedness model's level-shifted held-out
    predictions, with the pre-set call; sensitivity analyses.

Outputs in results/phase4/. Usage: python scripts/36_phase4_evaluate.py
"""

import json
import sqlite3

import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.metrics import roc_auc_score

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import bootstrap_mean, marginal_probability, recalibrate_intercept
from plate_truth.phase4 import (
    GROUP,
    PAIR_KEYS,
    boundary_lrt_p,
    bootstrap_pair_or,
    classify_taxa,
    discordant_pairs,
    mantel_haenszel_or,
    matched_condition_pairs,
    mean_log_loss,
    pair_counts,
    reportable_taxa,
    taxon_label,
    taxon_oe,
)
from plate_truth.queries import asv_taxonomy

ROOT = data_path("interim") / "phase4"
FITS = ROOT / "fits"
RES = data_path("results") / "phase4"


def fit_info(model: str, fold: str) -> dict:
    """fit.json of one fit; an empty dict if the fit is missing."""
    path = FITS / model / fold / "fit.json"
    return json.loads(path.read_text()) if path.exists() else {}


def spec(model: str) -> dict:
    """Random and known groups of a model, from the spec.json of any of its fits."""
    return json.loads(next((FITS / model).glob("*/spec.json")).read_text())


def heldout(model: str, rows: pd.DataFrame, nodes: int) -> pd.Series:
    """New-study predictions for every row from the leave-one-group-out fits: fixed effects plus
    known random effects, averaged over the unknown random effects."""
    s = spec(model)
    unknown = [g for g in s["random"] if g not in s["known"]]
    pred = pd.Series(np.nan, index=rows.index)
    for g in sorted(rows[GROUP].unique()):
        info = fit_info(model, g)
        if info.get("status") != "ok":
            raise RuntimeError(f"{model}/{g}: fit missing or failed")
        eta = pd.read_csv(FITS / model / g / "test_eta.csv")
        sigma = float(np.sqrt(sum(info["sd"][u] ** 2 for u in unknown)))
        pred[eta["row_id"].to_numpy()] = marginal_probability(eta["eta"].to_numpy(), sigma, nodes)
    return pred


def compare_heldout(rows: pd.DataFrame, p_with: pd.Series, p_without: pd.Series, need: int,
                    n_boot: int, seed: int) -> tuple[pd.DataFrame, dict]:
    """Per group: log loss and AUC with and without the term; and the transfer verdict."""
    out = []
    for g, s in rows.groupby(GROUP):
        y = s["y"].to_numpy()
        out.append({GROUP: g, "host_associated": int(s["host_associated"].iloc[0]),
                    "n": len(s), "events": int(y.sum()),
                    "logloss_with": mean_log_loss(y, p_with[s.index]),
                    "logloss_without": mean_log_loss(y, p_without[s.index]),
                    "auc_with": roc_auc_score(y, p_with[s.index]),
                    "auc_without": roc_auc_score(y, p_without[s.index])})
    t = pd.DataFrame(out)
    t["logloss_diff"] = t["logloss_with"] - t["logloss_without"]
    t["auc_diff"] = t["auc_with"] - t["auc_without"]
    boots = bootstrap_mean(t["logloss_diff"].to_numpy(), n_boot, seed)
    lo, hi = np.quantile(boots, [0.025, 0.975])
    improved = int((t["logloss_diff"] < 0).sum())
    summary = {"groups": len(t), "groups_improved": improved, "groups_needed": need,
               "mean_logloss_diff": t["logloss_diff"].mean(), "mean_logloss_diff_lo": lo,
               "mean_logloss_diff_hi": hi, "mean_auc_diff": t["auc_diff"].mean()}
    return t, summary


def lrt(model: str, null: str) -> dict:
    """Log-likelihoods and likelihood-ratio statistic of a full-data fit against its null
    model; the caller turns the statistic into a p-value (boundary or not)."""
    a, b = fit_info(model, "all"), fit_info(null, "all")
    stat = 2 * (a["logLik"] - b["logLik"])
    return {"loglik": a["logLik"], "loglik_null": b["logLik"], "lr_stat": stat,
            "converged": a.get("converged"), "pdHess": a.get("pdHess"),
            "null_converged": b.get("converged"), "null_pdHess": b.get("pdHess")}


def taxon_class_table(model: str, rows: pd.DataFrame, taxon: str, cfg: dict) -> pd.DataFrame:
    """Taxon x medium-class effects (conditional modes, 95% intervals) for reportable taxa, with
    rows and enriched rows per cell."""
    re_ = pd.read_csv(FITS / model / "all" / f"ranef_{taxon}_medium.csv")
    split = re_["level"].str.rsplit(":", n=1)
    re_[taxon], re_["medium"] = split.str[0], split.str[1]
    re_["lo"], re_["hi"] = re_["effect"] - 1.96 * re_["se"], re_["effect"] + 1.96 * re_["se"]
    re_["interval_excludes_0"] = (re_["lo"] > 0) | (re_["hi"] < 0)
    cells = rows.groupby([taxon, "medium"]).agg(rows=("y", "size"), enriched=("y", "sum"),
                                               groups=(GROUP, "nunique")).reset_index()
    keep = reportable_taxa(rows, taxon, cfg["q5"]["report_min_rows"],
                           cfg["q5"]["report_min_classes"])
    t = re_.merge(cells, on=[taxon, "medium"])
    return t[t[taxon].isin(keep)].drop(columns="level").sort_values([taxon, "medium"])


def coef_row(model: str, term: str) -> dict:
    """Estimate, SE, odds ratio and 95% Wald interval of one fixed effect in a full fit."""
    c = pd.read_csv(FITS / model / "all" / "coef.csv").set_index("term").loc[term]
    lo, hi = c["estimate"] - 1.96 * c["se"], c["estimate"] + 1.96 * c["se"]
    return {"estimate": c["estimate"], "se": c["se"], "odds_ratio": np.exp(c["estimate"]),
            "or_lo": np.exp(lo), "or_hi": np.exp(hi)}


def q5(cfg: dict, p3: dict, seed: int) -> None:
    rows = pd.read_parquet(ROOT / "rows" / "q5.parquet")
    nodes = p3["marginal"]["gauss_hermite_nodes"]
    tests = []
    for taxon, model in [("family", "q5_family"), ("phylum_full", "q5_phylum")]:
        t = lrt(model, f"{model}_null")
        t.update({"taxon": taxon, "p_value": boundary_lrt_p(t["loglik_null"], t["loglik"]),
                  "sd_taxon_medium": fit_info(model, "all")["sd"][f"{taxon}_medium"]})
        t["H5a"] = bool(t["p_value"] < cfg["q5"]["alpha"])
        tests.append(t)
        taxon_class_table(model, rows, taxon, cfg).to_csv(
            RES / f"q5_{'family' if taxon == 'family' else 'phylum'}_medium.csv", index=False)
    by_group, s = compare_heldout(rows, heldout("q5_family", rows, nodes),
                                  heldout("q5_family_null", rows, nodes),
                                  cfg["q5"]["transfer_groups_required"],
                                  cfg["q6"]["n_boot"], seed)
    by_group.to_csv(RES / "q5_heldout_by_group.csv", index=False)
    s["H5b"] = bool(s["groups_improved"] >= s["groups_needed"] and s["mean_logloss_diff_hi"] < 0)
    tests[0].update(s)
    pd.DataFrame(tests).to_csv(RES / "q5_tests.csv", index=False)
    print(pd.DataFrame(tests).T.to_string(), flush=True)
    print(by_group.round(4).to_string(), flush=True)


def q6(cfg: dict, p3: dict, seed: int) -> None:
    rows = pd.read_parquet(ROOT / "rows" / "q6.parquet")
    nodes = p3["marginal"]["gauss_hermite_nodes"]
    summary = []
    for label, m in [("primary", "q6"), ("identity_ge_97", "q6close"),
                     ("indicators", "q6binary"), ("strain_16s_links", "q6s16")]:
        if not (FITS / f"{m}_match" / "all" / "fit.json").exists():
            summary.append({"analysis": label, "status": "not run"})
            continue
        t = lrt(f"{m}_match", f"{m}_nomatch")
        t["p_value"] = float(chi2.sf(max(t["lr_stat"], 0.0), 1))   # one fixed effect
        summary.append({"analysis": label, **coef_row(f"{m}_match", "match"), **t,
                        "n": fit_info(f"{m}_match", "all")["n"],
                        "events": fit_info(f"{m}_match", "all")["events"]})
    by_group, s = compare_heldout(rows, heldout("q6_match", rows, nodes),
                                  heldout("q6_nomatch", rows, nodes),
                                  cfg["q6"]["transfer_groups_required"], cfg["q6"]["n_boot"],
                                  seed)
    by_group.to_csv(RES / "q6_heldout_by_group.csv", index=False)
    prim = summary[0]
    prim.update(s)
    prim["H6"] = bool(prim["or_lo"] > 1 and s["groups_improved"] >= s["groups_needed"])

    classes = cfg["medium_classes"]
    conds = rows[["condition_id", "medium", GROUP, *PAIR_KEYS]].drop_duplicates("condition_id")
    pairs = matched_condition_pairs(conds)
    rec = rows.drop_duplicates("asv_id")[["asv_id", *[f"any_{c}" for c in classes]]]
    disc = discordant_pairs(rows, pairs, rec)
    counts = pair_counts(disc)
    counts.to_csv(RES / "q6_matched_pairs.csv", index=False)
    boots = bootstrap_pair_or(disc, cfg["q6"]["n_boot"], seed)
    lo, hi = np.nanquantile(boots, [0.025, 0.975])
    prim.update({"pairs_or": mantel_haenszel_or(counts), "pairs_or_lo": lo, "pairs_or_hi": hi,
                 "pairs_discordant": len(disc), "pairs_asvs": disc["asv_id"].nunique()})
    pd.DataFrame(summary).to_csv(RES / "q6_summary.csv", index=False)
    print(pd.DataFrame(summary).T.to_string(), flush=True)
    print(by_group.round(4).to_string(), flush=True)


def q7(cfg: dict, p3: dict, seed: int) -> None:
    rows = pd.read_parquet(data_path("interim") / "phase3" / "rows" / "primary.parquet")
    tax = asv_taxonomy(sqlite3.connect(data_path("database"))).set_index("asv_id")
    rows["family"] = rows["asv_id"].map(taxon_label(tax, "family"))
    rows["phylum_full"] = rows["asv_id"].map(tax["phylum"].fillna("unassigned"))
    p3h = pd.read_parquet(data_path("interim") / "phase3" / "heldout" / "primary.parquet")
    assert (p3h["asv_id"].to_numpy() == rows["asv_id"].to_numpy()).all()
    preds = {"notaxon": heldout("q7_notaxon", rows, p3["marginal"]["gauss_hermite_nodes"]),
             "identity_only": p3h["pred_identity_only"]}
    c = cfg["q7"]
    cont = p3["local_recalibration"]["continuity"]
    variants = []
    for name, p in preds.items():
        shifted = pd.Series(np.nan, index=rows.index)
        for _, s in rows.groupby(GROUP):
            shifted[s.index] = recalibrate_intercept(s["y"].to_numpy(), p[s.index].to_numpy(),
                                                     p[s.index].to_numpy(), cont)
        variants.append((name, shifted))
        if name == "notaxon":
            variants.append(("notaxon_unshifted", p))
    for name, p in variants:
        df = rows.assign(pred=p.to_numpy())
        for taxon in ["family", "phylum_full"]:
            per = taxon_oe(df, taxon, "pred")
            calls = classify_taxa(per, taxon, c["min_expected_per_group"], c["min_groups"],
                                  c["consistency_share"], c["n_boot"], seed)
            tag = "family" if taxon == "family" else "phylum"
            suffix = "" if name == "notaxon" else f"_{name}"
            calls.sort_values("oe").to_csv(RES / f"q7_{tag}{suffix}.csv", index=False)
            if name == "notaxon":
                per.to_csv(RES / f"q7_{tag}_by_group.csv", index=False)
                print(tag, calls["call"].value_counts().to_dict(), flush=True)


def fit_status() -> None:
    """One row per fit: status, convergence, positive-definite Hessian, warnings."""
    out = []
    for path in sorted(FITS.glob("*/*/fit.json")):
        info = json.loads(path.read_text())
        w = info.get("warnings") or []
        out.append({"model": path.parent.parent.name, "fold": path.parent.name,
                    "status": info.get("status"), "converged": info.get("converged"),
                    "pdHess": info.get("pdHess"), "n": info.get("n"),
                    "events": info.get("events"),
                    "warnings": " | ".join(w if isinstance(w, list) else [w])})
    pd.DataFrame(out).to_csv(RES / "fit_status.csv", index=False)


def main() -> None:
    cfg = load_config("phase4")
    p3 = load_config("phase3")
    seed = load_config("thresholds")["resampling"]["seed"]
    RES.mkdir(parents=True, exist_ok=True)
    fit_status()
    q5(cfg, p3, seed)
    q6(cfg, p3, seed)
    q7(cfg, p3, seed)


if __name__ == "__main__":
    main()
