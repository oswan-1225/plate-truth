"""Phase 5, step 5: the Phase 3 model applied to the survey samples (ANALYSIS_PLAN.md
sections 22 and 24).

For every survey sample: the shortcut share (identity >= 94.5%, and >= 97%), the predicted
recovery share for one plate (R1) and for any of 8 plates (R8) with the study effect at its
10th, 50th and 90th percentiles (point estimates, and 200 fixed-effect draws at the 10th and
90th), and composition by the Phase 4 Q7 calls. Studies are the mean of their samples;
habitats the mean of their studies, with intervals from replicates that resample studies and
take one draw each. The shortcut verdict per habitat compares the shortcut share with the upper
end of the range (R1, primary; R8, generous).

Predictions average over the 687 usable training plating conditions, each study group weighted
equally (sensitivity: each condition weighted equally).

Outputs: results/phase5/<arm>/{samples,studies,habitats,top_under_families}.csv,
results/phase5/amplicon/run_qc.csv (primer-window share per run) and
results/phase5/sensitivity.csv; per-sample draws in data/interim/phase5/pred/.
Usage: python scripts/43_phase5_apply.py [amplicon|shotgun ...]
"""

import json
import sys

import numpy as np
import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import CONTINUOUS, GROUP, design_columns, keep_levels
from plate_truth.phase4 import taxon_label
from plate_truth.phase5 import (
    amplicon_samples,
    composition_shares,
    condition_columns,
    condition_effects,
    condition_weights,
    mixture_points,
    recovery_shares,
    replicate_means,
    run_window_qc,
    sequence_eta,
    shortcut_share,
    shortcut_verdict,
    shotgun_samples,
    spline_basis,
    study_effect,
)

VARIANTS = {
    "primary": {},
    "iddef2": {"arms": ["amplicon"], "id_col": "identity_id2"},
    "gap_blind": {"arms": ["amplicon"], "id_col": "identity_mid"},
    "conditions_equal": {"weights": "condition"},
    "depth5000": {"min_reads": 5000},
    "v4_only": {"arms": ["amplicon"], "v4_only": True},
    "no_lab_reactor": {"arms": ["amplicon"], "drop_flag": "lab_reactor"},
    "abundance_median": {"arms": ["shotgun"], "abundance": "median"},
    "overlap150": {"arms": ["shotgun"], "min_overlap": 150},
}
CALLS = ["under", "over", "neither", "not evaluated"]


# --- model --------------------------------------------------------------------------------------

def load_model(cfg3: dict, cfg5: dict) -> dict:
    """Full-data primary fit, the training condition designs and centres."""
    fit = data_path("interim") / "phase3" / cfg3["fits_dir"] / "primary" / "primary" / "all"
    rows = pd.read_parquet(data_path("interim") / "phase3" / "rows" / "primary.parquet")
    centers = {c: float(rows[c].median()) for c in CONTINUOUS}
    cond = rows.drop_duplicates("condition_id")
    x = design_columns(cond, pd.Series(True, index=cond.index), centers)
    beta = pd.read_csv(fit / "coef.csv").set_index("term")["estimate"]
    draws = pd.read_csv(fit / "draws.csv").iloc[:cfg5["model"]["n_draws"]].T
    cols = condition_columns(x.columns)
    return {"beta": beta, "draws": draws, "x_cond": x[cols], "groups": cond[GROUP],
            "phylum": pd.read_csv(fit / "ranef_phylum.csv").set_index("level")["effect"],
            "sd": json.loads((fit / "fit.json").read_text())["sd"], "centers": centers,
            "env_median": float(rows.drop_duplicates(["source_key", "asv_id"])["log_env"]
                                .median())}


def condition_mix(model: dict, weights: str, cfg5: dict) -> tuple[np.ndarray, list]:
    """Points for the point estimate and for each draw."""
    m = cfg5["model"]
    w = (condition_weights(model["groups"]) if weights == "group"
         else pd.Series(1.0 / len(model["groups"]), index=model["groups"].index)).to_numpy()
    sd_c = model["sd"]["condition_id"]
    point = mixture_points(condition_effects(model["x_cond"], model["beta"]), w, sd_c,
                           m["condition_points"])
    g = condition_effects(model["x_cond"], model["draws"])
    draws = [mixture_points(g[:, d], w, sd_c, m["condition_points_draws"])
             for d in range(g.shape[1])]
    return point, draws


# --- survey sequences -----------------------------------------------------------------------------

def q7_calls() -> tuple[pd.Series, pd.Series]:
    res = data_path("results") / "phase4"
    fam = pd.read_csv(res / "q7_family.csv").set_index("family")["call"]
    phy = pd.read_csv(res / "q7_phylum.csv").set_index("phylum_full")["call"]
    return fam, phy


def load_amplicon(v: dict, cfg5: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    runs = pd.read_csv(data_path("processed") / "survey_runs_amplicon.csv", dtype=str)
    seqs, depth = [], []
    for (habitat, study, ref), _ in runs.groupby(["habitat", "study", "reference"], sort=False):
        if v.get("v4_only") and not ref.startswith("V4"):
            continue
        if v.get("drop_flag") and                 cfg5.get("study_notes", {}).get(study, {}).get("flag") == v["drop_flag"]:
            continue
        d = data_path("interim") / "phase5" / "amplicon" / study
        if not (d / "taxonomy.csv").exists():
            continue
        counts = pd.read_csv(d / "counts.csv")
        if counts.empty:                     # no run produced usable reads (see run_qc.csv)
            continue
        nr = pd.read_csv(d / "nearest_relative.csv").set_index("asv_id")
        tax = pd.read_csv(d / "taxonomy.csv").set_index("asv_id")
        s, dp = amplicon_samples(counts, nr[v.get("id_col", "identity_id1")], tax,
                                 cfg5["sample"]["present_min_reads"],
                                 v.get("min_reads", cfg5["sample"]["min_prokaryotic_reads"]))
        s["log_env"] = np.log10(s["a"])
        meta = {"habitat": habitat, "study": study, "reference": ref}
        seqs.append(s.assign(**meta))
        depth.append(dp.assign(**meta))
    return pd.concat(seqs, ignore_index=True), pd.concat(depth, ignore_index=True)


def write_run_qc(cfg5: dict) -> None:
    """Amplicon runs' merged and primer-window shares (results/phase5/amplicon/run_qc.csv);
    runs below the census screen's window threshold are flagged."""
    runs = pd.read_csv(data_path("processed") / "survey_runs_amplicon.csv", dtype=str)
    min_share = cfg5["selection"]["amplicon"]["min_primer_pair_share"]
    pairs = runs.set_index("run_accession")["read_count"].astype(float)
    out = []
    for (habitat, study), _ in runs.groupby(["habitat", "study"], sort=False):
        f = data_path("interim") / "phase5" / "amplicon" / study / "read_qc.csv"
        if f.exists():
            out.append(run_window_qc(pd.read_csv(f), min_share, pairs).assign(habitat=habitat,
                                                                       study=study))
    if out:
        t = pd.concat(out, ignore_index=True)
        d = data_path("results") / "phase5" / "amplicon"
        d.mkdir(parents=True, exist_ok=True)
        t.to_csv(d / "run_qc.csv", index=False)
        print(f"run QC: {int(t['low_window_share'].sum())} of {len(t)} runs below "
              f"{min_share:g} window share", flush=True)


def load_shotgun(v: dict, cfg5: dict, model: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    runs = pd.read_csv(data_path("processed") / "survey_runs_shotgun.csv", dtype=str)
    w = cfg5["shotgun_window"]
    seqs, depth = [], []
    for (habitat, study), _ in runs.groupby(["habitat", "study"], sort=False):
        f = data_path("interim") / "phase5" / "shotgun" / study / "fragments.csv"
        if not f.exists():
            continue
        frags = pd.read_csv(f)
        frags = frags[frags["overlap"] >= v.get("min_overlap", w["min_overlap"])]
        s, dp = shotgun_samples(frags, w["pseudo_taxon_rank"],
                                v.get("min_reads", cfg5["sample"]["min_window_fragments"]))
        s["log_env"] = (np.full(len(s), model["env_median"]) if v.get("abundance") == "median"
                        else np.log10(s["abundance"]))
        meta = {"habitat": habitat, "study": study, "reference": w["reference"]}
        seqs.append(s.assign(**meta))
        depth.append(dp.assign(**meta))
    return pd.concat(seqs, ignore_index=True), pd.concat(depth, ignore_index=True)


def add_model_columns(s: pd.DataFrame, model: dict, cfg3: dict) -> pd.DataFrame:
    s = s.copy()
    s["no_hit"] = s["identity"].isna().astype(float)
    s["identity_model"] = s["identity"].fillna(cfg3["identity_spline"]["no_hit_value"])
    s["phylum_level"] = keep_levels(s["phylum"], cfg3["phyla"], "other")
    s["phylum_effect"] = s["phylum_level"].map(model["phylum"]).fillna(0.0)
    return s


# --- predictions per sample -----------------------------------------------------------------------

def predict_samples(s: pd.DataFrame, model: dict, mix, cfg3: dict, cfg5: dict
                    ) -> tuple[pd.DataFrame, dict]:
    m = cfg5["model"]
    sd = model["sd"]
    point_pts, draw_pts = mix
    ids = np.unique(s["identity_model"].to_numpy())
    basis = dict(zip(ids, spline_basis(ids, cfg3)))
    k = m["conditions_k"]
    qs = m["study_quantiles"]
    out, draws = [], {}
    for (habitat, study, run), g in s.groupby(["habitat", "study", "run"], sort=False):
        b = np.vstack(g["identity_model"].map(basis).to_numpy())
        args = (b, g["no_hit"].to_numpy(), g["log_env"].to_numpy(),
                g["phylum_effect"].to_numpy())
        a = g["a"].to_numpy()
        eta = sequence_eta(*args, model["beta"], model["centers"]["log_env"])
        row = {"habitat": habitat, "study": study, "run": run, "sequences": len(g)}
        for q in qs:
            u = study_effect(q, sd["study_id"])
            tag = f"p{int(round(q * 100))}"
            row[f"r1_{tag}"], row[f"r8_{tag}"] = recovery_shares(
                eta, a, point_pts, u, sd["source_key"], k, m["source_nodes"])
        eta_d = sequence_eta(*args, model["draws"], model["centers"]["log_env"])
        arr = np.empty((4, eta_d.shape[1]))
        for d in range(eta_d.shape[1]):
            for j, q in enumerate((qs[0], qs[-1])):
                arr[2 * j:2 * j + 2, d] = recovery_shares(
                    eta_d[:, d], a, draw_pts[d], study_effect(q, sd["study_id"]),
                    sd["source_key"], k, m["source_nodes_draws"])
        draws[run] = arr          # rows: r1_p10, r8_p10, r1_p90, r8_p90
        out.append(row)
    return pd.DataFrame(out), draws


def composition(s: pd.DataFrame, proxy: float) -> pd.DataFrame:
    fam_call, phy_call = q7_calls()
    s = s.copy()
    s["family_label"] = taxon_label(s, "family")
    s["family_call"] = s["family_label"].map(fam_call).fillna("not evaluated")
    s["phylum_call"] = s["phylum"].fillna("unassigned").map(phy_call).fillna("not evaluated")
    rows = []
    for (habitat, study, run), g in s.groupby(["habitat", "study", "run"], sort=False):
        a = g["a"].to_numpy()
        fam = composition_shares(a, g["family_call"]).reindex(CALLS, fill_value=0.0)
        phy = composition_shares(a, g["phylum_call"]).reindex(CALLS, fill_value=0.0)
        sc = g[g["identity"] >= proxy]
        rows.append({"habitat": habitat, "study": study, "run": run,
                     "shortcut": shortcut_share(g["identity"], a, proxy),
                     "shortcut_97": shortcut_share(g["identity"], a, 97.0),
                     "no_hit": float(a @ g["no_hit"].to_numpy() / a.sum()),
                     **{f"family_{c.replace(' ', '_')}": fam[c] for c in CALLS},
                     **{f"phylum_{c.replace(' ', '_')}": phy[c] for c in CALLS},
                     "shortcut_in_under_family": (sc["a"][sc["family_call"] == "under"].sum()
                                                  / sc["a"].sum() if len(sc) else np.nan)})
    return pd.DataFrame(rows), s


def top_under(s: pd.DataFrame) -> pd.DataFrame:
    """Per habitat: mean read share (studies weighted equally) of each under-recovered family.
    'unassigned family (<phylum>)' groups are modelled in Phase 4 but not interpreted, so they
    are not listed (they still count in the composition shares)."""
    u = s[(s["family_call"] == "under") & ~s["family_label"].str.startswith("unassigned")]
    per_run = u.groupby(["habitat", "study", "run", "family_label"])["a"].sum()
    runs = s.groupby(["habitat", "study"])["run"].nunique()
    per_study = per_run.groupby(["habitat", "study", "family_label"]).sum()
    per_study = per_study / per_study.index.droplevel(2).map(runs).to_numpy()
    studies = s.groupby("habitat")["study"].nunique()
    hab = per_study.groupby(["habitat", "family_label"]).sum()
    hab = hab / hab.index.get_level_values(0).map(studies).to_numpy()
    hab = hab.rename("mean_read_share").reset_index()
    hab["rank"] = hab.groupby("habitat")["mean_read_share"].rank(ascending=False, method="first")
    return hab[hab["rank"] <= 5].sort_values(["habitat", "rank"])


# --- aggregation ------------------------------------------------------------------------------

def interval(x: np.ndarray) -> tuple[float, float]:
    lo, hi = np.quantile(x, [0.025, 0.975])
    return float(lo), float(hi)


def summarise(samples: pd.DataFrame, draws: dict, cfg5: dict, seed: int) -> tuple:
    """Study means and the habitat table with intervals and the shortcut verdicts."""
    m, crit = cfg5["model"], cfg5["criteria"]
    num, den = crit["shortcut_studies_share"]
    value_cols = [c for c in samples.columns if c not in ("habitat", "study", "run", "sequences")]
    studies = samples.groupby(["habitat", "study"], sort=False)[value_cols].mean()
    studies["samples"] = samples.groupby(["habitat", "study"], sort=False).size()
    studies = studies.reset_index()
    rows = []
    for habitat, st in studies.groupby("habitat", sort=False):
        # per study, per draw: mean over the study's samples (rows r1_p10, r8_p10, r1_p90, r8_p90)
        sd_arr = np.stack([np.mean([draws[r] for r in samples.loc[
            (samples["habitat"] == habitat) & (samples["study"] == s), "run"]], axis=0)
            for s in st["study"]])                                # (studies, 4, draws)
        h = cfg5["habitats"][habitat]
        row = {"habitat": habitat, "training_habitat": h["training"], "report": h["report"],
               "studies": len(st), "samples": int(st["samples"].sum())}
        reps = {}
        for j, name in enumerate(["r1_p10", "r8_p10", "r1_p90", "r8_p90"]):
            reps[name] = replicate_means(sd_arr[:, j, :], m["n_replicates"], seed)
            row[name] = st[name].mean()
            row[f"{name}_lo"], row[f"{name}_hi"] = interval(reps[name])
        s_reps = replicate_means(st[["shortcut"]].to_numpy(), m["n_replicates"], seed)
        row["shortcut"] = st["shortcut"].mean()
        row["shortcut_lo"], row["shortcut_hi"] = interval(s_reps)
        row["shortcut_97"] = st["shortcut_97"].mean()
        for which, j in (("r1", 2), ("r8", 3)):
            d_draws = st["shortcut"].to_numpy()[:, None] - sd_arr[:, j, :]
            d_reps = replicate_means(d_draws, m["n_replicates"], seed)
            v = shortcut_verdict(st["shortcut"].to_numpy() - st[f"{which}_p90"].to_numpy(),
                                 d_reps, num, den)
            row.update({f"{k}_{which}": val for k, val in v.items() if k != "studies"})
        for c in [c for c in st.columns if c.startswith(("family_", "phylum_", "no_hit",
                                                          "shortcut_in"))]:
            vals = st[c].to_numpy(float)
            row[c] = np.nanmean(vals)
            if np.isfinite(vals).all():
                row[f"{c}_lo"], row[f"{c}_hi"] = interval(
                    replicate_means(vals[:, None], m["n_replicates"], seed))
        rows.append(row)
    return studies, pd.DataFrame(rows)


def run_arm(arm: str, vname: str, v: dict, model: dict, cfg3: dict, cfg5: dict, proxy: float,
            seed: int) -> pd.DataFrame | None:
    if arm not in v.get("arms", [arm]):
        return None
    s, depth = (load_amplicon(v, cfg5) if arm == "amplicon" else load_shotgun(v, cfg5, model))
    s = add_model_columns(s, model, cfg3)
    mix = condition_mix(model, v.get("weights", "group"), cfg5)
    pred, draws = predict_samples(s, model, mix, cfg3, cfg5)
    comp, s_lab = composition(s, proxy)
    samples = pred.merge(comp, on=["habitat", "study", "run"]).merge(
        depth.rename(columns={"run": "run"})[["run", "prokaryotic_reads"]], on="run")
    studies, habitats = summarise(samples, draws, cfg5, seed)
    out = data_path("results") / "phase5" / arm
    pdir = data_path("interim") / "phase5" / "pred"
    pdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(pdir / f"{arm}_{vname}.npz", **draws)
    if vname == "primary":
        out.mkdir(parents=True, exist_ok=True)
        samples.to_csv(out / "samples.csv", index=False)
        studies.to_csv(out / "studies.csv", index=False)
        habitats.to_csv(out / "habitats.csv", index=False)
        depth.to_csv(out / "sample_depth.csv", index=False)
        top_under(s_lab).to_csv(out / "top_under_families.csv", index=False)
    print(f"{arm} {vname}:", flush=True)
    print(habitats[["habitat", "studies", "samples", "shortcut", "r1_p10", "r1_p90", "r8_p90",
                    "d_mean_r1", "d_lo_r1", "shortcut_exceeds_r1", "shortcut_exceeds_r8"]]
          .round(3).to_string(index=False), flush=True)
    return habitats.assign(arm=arm, variant=vname)


def main(arms: list[str], variants: list[str]) -> None:
    cfg3, cfg5 = load_config("phase3"), load_config("phase5")
    th = load_config("thresholds")
    proxy, seed = th["identity"]["database_proxy_percent"], cfg5["selection"]["seed"]
    model = load_model(cfg3, cfg5)
    if "amplicon" in arms:
        write_run_qc(cfg5)
    res = []
    for vname in variants:
        for arm in arms:
            h = run_arm(arm, vname, VARIANTS[vname], model, cfg3, cfg5, proxy, seed)
            if h is not None:
                res.append(h)
    sens = data_path("results") / "phase5" / "sensitivity.csv"
    new = pd.concat(res, ignore_index=True)
    if sens.exists():
        old = pd.read_csv(sens)
        keep = ~old.set_index(["arm", "variant"]).index.isin(
            new.set_index(["arm", "variant"]).index)
        new = pd.concat([old[keep], new], ignore_index=True)
    new.to_csv(sens, index=False)


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if x in ("amplicon", "shotgun")] or ["amplicon", "shotgun"]
    vs = [x for x in sys.argv[1:] if x in VARIANTS] or list(VARIANTS)
    main(a, vs)
