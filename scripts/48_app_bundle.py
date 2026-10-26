"""Phase 6: the model bundle for the Streamlit app (app/model/), and a check that the app's
code reproduces the registered survey predictions.

The bundle holds only what predictions need, in small text files: the sequence fixed effects
(intercept, identity spline, no-hit, source abundance) and their first draws, the random-effect
SDs, phylum effects, the training plating-condition mix as quantile points (point estimate and
per draw, exactly as script 43 builds them), the identity spline basis tabulated every 0.01%
(so the app does not need R), the Phase 4 family calls, and the Phase 5 habitat tables.

Check: for the first runs of the survey samples, plate_truth.app_model.predict must reproduce
the one-plate and any-of-K shares in results/phase5/amplicon/samples.csv; the script stops if
the largest difference exceeds 1e-4. The spline interpolation error is printed as well.
Run after scripts 43 and 46.
"""

import importlib.util
import json
import shutil
import sys

import numpy as np
import pandas as pd

from plate_truth.app_model import community_summary, load_bundle, predict
from plate_truth.config import data_path, load_config, project_root
from plate_truth.phase5 import spline_basis, spline_terms

N_DRAWS = 100          # fixed-effect draws shipped with the app (intervals per taxon)
GRID_STEP = 0.01
CHECK_RUNS = 6
TOLERANCE = 1e-4
EXAMPLE = ("MGYS00006176", "ERR11027514")   # first soil survey sample (results/phase5)


def script43():
    path = project_root() / "scripts" / "43_phase5_apply.py"
    spec = importlib.util.spec_from_file_location("phase5_apply", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def example(cfg5: dict, folder) -> None:
    """One public survey sample in the app's upload format (every prokaryotic ASV of the run,
    with its count, identity, phylum and family), for the app's example button."""
    from plate_truth.phase5 import is_prokaryote
    study, run = EXAMPLE
    d = data_path("interim") / "phase5" / "amplicon" / study
    counts = pd.read_csv(d / "counts.csv")
    tax = pd.read_csv(d / "taxonomy.csv").set_index("asv_id")
    nr = pd.read_csv(d / "nearest_relative.csv").set_index("asv_id")["identity_id1"]
    c = counts[counts["run"] == run].join(tax, on="asv_id")
    c = c[is_prokaryote(c.fillna("")).to_numpy()]
    c["identity"] = c["asv_id"].map(nr)
    c[["asv_id", "count", "identity", "phylum", "family"]].sort_values(
        "count", ascending=False).to_csv(folder / "example_soil_sample.csv", index=False)
    print(f"example: {study} {run}, {len(c)} prokaryotic ASVs, {c['count'].sum()} reads")


def main() -> None:
    cfg3, cfg5 = load_config("phase3"), load_config("phase5")
    th = load_config("thresholds")
    s43 = script43()
    model = s43.load_model(cfg3, cfg5)
    point, draws = s43.condition_mix(model, "group", cfg5)
    out = project_root() / "app" / "model"
    out.mkdir(parents=True, exist_ok=True)

    seq_terms = ["(Intercept)", *spline_terms(model["beta"].index), "no_hit", "log_env"]
    m = cfg5["model"]
    meta = {
        "beta": model["beta"].loc[seq_terms].to_dict(),
        "sd": model["sd"], "env_center": model["centers"]["log_env"],
        "phylum_effects": model["phylum"].to_dict(),
        "condition_points": point.tolist(),
        "settings": {"conditions_k": m["conditions_k"], "source_nodes": m["source_nodes"],
                     "source_nodes_draws": m["source_nodes_draws"],
                     "study_quantiles": m["study_quantiles"], "phyla": cfg3["phyla"],
                     "shortcut_percent": th["identity"]["database_proxy_percent"],
                     "no_hit_value": cfg3["identity_spline"]["no_hit_value"],
                     "present_min_reads": cfg5["sample"]["present_min_reads"],
                     "min_prokaryotic_reads": cfg5["sample"]["min_prokaryotic_reads"]},
        "source": "Phase 3 primary model (revised, Deviation 1), fitted on all ten study "
                  "groups; condition mix: the 687 usable training conditions, each study "
                  "group weighted equally.",
    }
    (out / "model.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    model["draws"].loc[seq_terms].iloc[:, :N_DRAWS].to_csv(out / "draws.csv")
    np.savetxt(out / "condition_points_draws.csv", np.vstack(draws[:N_DRAWS]), delimiter=",",
               fmt="%.10g")

    lo, hi = cfg3["identity_spline"]["boundary_knots"]
    grid = np.round(np.arange(lo, hi + GRID_STEP / 2, GRID_STEP), 4)
    basis = spline_basis(grid, cfg3)
    pd.DataFrame(basis, columns=[f"b{j + 1}" for j in range(basis.shape[1])]).assign(
        identity=grid)[["identity", *[f"b{j + 1}" for j in range(basis.shape[1])]]].to_csv(
        out / "spline_grid.csv", index=False, float_format="%.10g")

    res = data_path("results")
    fam = pd.read_csv(res / "phase4" / "q7_family.csv")[["family", "call", "oe", "oe_lo",
                                                         "oe_hi"]]
    fam.to_csv(out / "family_calls.csv", index=False)
    data = project_root() / "app" / "data"
    data.mkdir(exist_ok=True)
    for arm in ("amplicon", "shotgun"):
        shutil.copyfile(res / "phase5" / arm / "habitats.csv", data / f"habitats_{arm}.csv")
    shutil.copyfile(res / "phase3" / "primary" / "curve.csv", data / "curve.csv")

    example(cfg5, data)

    # --- check against the registered survey predictions
    b = load_bundle(out)
    test_ids = np.random.default_rng(1).uniform(lo, hi, 500)
    from plate_truth.app_model import basis_from_grid
    err = np.abs(basis_from_grid(test_ids, b["grid_identity"], b["grid_basis"])
                 - spline_basis(test_ids, cfg3)).max()
    print(f"spline grid: largest interpolation error {err:.2e}")
    seqs, _ = s43.load_amplicon({}, cfg5)
    seqs = s43.add_model_columns(seqs, model, cfg3)
    ref = pd.read_csv(res / "phase5" / "amplicon" / "samples.csv").set_index("run")
    k = m["conditions_k"]
    worst = 0.0
    for run in list(dict.fromkeys(seqs["run"]))[:CHECK_RUNS]:
        g = seqs[seqs["run"] == run].reset_index(drop=True)
        t = g[["asv_id", "a", "identity", "no_hit", "identity_model", "phylum",
               "family"]].assign(count=0)
        t.attrs["phylum_given"] = True
        summ = community_summary(predict(t, b, n_draws=2), k, meta["settings"]
                                 ["shortcut_percent"])
        diffs = [abs(summ["one_plate_p10"] - ref.loc[run, "r1_p10"]),
                 abs(summ["one_plate_p90"] - ref.loc[run, "r1_p90"]),
                 abs(summ[f"any_of_{k}_p90"] - ref.loc[run, "r8_p90"]),
                 abs(summ["shortcut"] - ref.loc[run, "shortcut"])]
        worst = max(worst, *diffs)
        print(f"{run}: one plate p90 {summ['one_plate_p90']:.5f} vs {ref.loc[run, 'r1_p90']:.5f}"
              f"; shortcut {summ['shortcut']:.4f} vs {ref.loc[run, 'shortcut']:.4f}")
    print(f"largest difference from results/phase5: {worst:.2e}")
    if worst > TOLERANCE:
        sys.exit("app predictions do not reproduce the registered survey numbers")
    size = sum(f.stat().st_size for f in out.iterdir()) / 1e6
    print(f"wrote app/model/ ({size:.1f} MB) and app/data/")


if __name__ == "__main__":
    main()
