"""Phase 3, step 2: fit the mixed models (ANALYSIS_PLAN.md sections 4, 5 and 7).

Every fit is one model on one training set: a leave-one-study-group-out fold (held-out group
predicted from fixed effects) or 'all' (full data; curve grid predicted with fixed-effect draws).
Fits run in the plan's priority order: the primary model's folds and full fit, the identity-only
and no-identity models' folds, the habitat model's full fit, then each sensitivity analysis
(primary model folds and full fit) in plan order.

Resumable and stoppable: each fit writes data/interim/phase3/<fits_dir>/<variant>/<model>/<fold>/
(fits_dir from config/phase3.yaml: the revised model of ANALYSIS_PLAN.md Deviation 1) and
is skipped when its fit.json exists. Creating the file data/interim/phase3/STOP stops new fits
from starting (running fits finish); delete it to continue.

Usage: python scripts/27_phase3_fit.py [variant ...]   (default: all variants built by script 26)
"""

import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import CONTINUOUS, GROUP, VARIANTS, design_columns
from plate_truth.rfit import fixed_formula, is_done, random_effects, run_fit

PRIMARY_MODELS = [("primary", True, True), ("identity_only", True, False),
                  ("no_identity", True, False), ("habitat", False, True)]
SENSITIVITY_MODELS = [("primary", True, True)]


def jobs_for(variant: str, groups: list[str]) -> list[tuple[str, str, str]]:
    """(variant, model, fold) in priority order; fold 'all' is the full-data fit."""
    models = PRIMARY_MODELS if variant == "primary" else SENSITIVITY_MODELS
    out = []
    for model, loso, full in models:
        if loso:
            out += [(variant, model, g) for g in groups]
        if full:
            out.append((variant, model, "all"))
    return out


def curve_grid(columns: list[str], cfg: dict, host: bool) -> pd.DataFrame:
    """Identity grid at average conditions (all other inputs 0 = centre or average level),
    plus one no-hit row; duplicated for host = -1 and +1 for the habitat model."""
    s = cfg["identity_spline"]
    ident = np.arange(s["boundary_knots"][0], s["boundary_knots"][1] + 1e-9,
                      cfg["curve"]["grid_step"])
    g = pd.DataFrame(0.0, index=range(len(ident) + 1), columns=columns)
    g["identity"] = [*ident, s["no_hit_value"]]
    g["no_hit"] = [0.0] * len(ident) + [1.0]
    if not host:
        return g
    return pd.concat([g.assign(host=-1.0), g.assign(host=1.0)], ignore_index=True)


def run_job(job, rows: pd.DataFrame, cfg: dict, centers: dict, out_root, seed: int) -> str:
    """Build one fit's inputs and run it in R."""
    variant, model, fold = job
    outdir = out_root / variant / model / fold
    if is_done(outdir):
        return f"skip {variant}/{model}/{fold}"
    train = rows[GROUP] != fold if fold != "all" else pd.Series(True, index=rows.index)
    x = design_columns(rows, train, centers)
    random, known = random_effects(model)
    data = pd.concat([rows[["y", "w", *random]].reset_index(drop=True),
                      x.reset_index(drop=True)], axis=1)
    data["row_id"] = rows.index.to_numpy()
    data["train"] = train.astype(int).to_numpy()
    grid = curve_grid(list(x.columns), cfg, model == "habitat") if fold == "all" else None
    t0 = time.time()
    info = run_fit(outdir, data, fixed_formula(model, list(x.columns), cfg), grid,
                   cfg["curve"]["n_draws"], seed, random, known)
    return (f"done {variant}/{model}/{fold} in {(time.time() - t0) / 60:.1f} min: "
            f"{info['status']}, converged={info.get('converged')}, pdHess={info.get('pdHess')}")


def main(variants: list[str]) -> None:
    cfg = load_config("phase3")
    seed = load_config("thresholds")["resampling"]["seed"]
    root = data_path("interim") / "phase3"
    fits = root / cfg["fits_dir"]
    stop = root / "STOP"
    jobs = []
    for v in variants:
        groups = pd.read_parquet(root / "rows" / f"{v}.parquet", columns=[GROUP])[GROUP]
        jobs += [j for j in jobs_for(v, sorted(groups.unique()))
                 if not is_done(fits / j[0] / j[1] / j[2])]
    print(f"{len(jobs)} fits queued", flush=True)

    cache, lock = {}, threading.Lock()

    def load(v: str):
        """Rows and centres of one variant; at most two variants are held in memory (jobs run
        in variant order, so at most two are in use at once)."""
        with lock:
            if v not in cache:
                while len(cache) >= 2:
                    cache.pop(next(iter(cache)))
                rows = pd.read_parquet(root / "rows" / f"{v}.parquet")
                cache[v] = (rows, {c: float(rows[c].median()) for c in CONTINUOUS})
            return cache[v]

    def work(job):
        if stop.exists():
            return f"stopped before {'/'.join(job)}"
        try:
            rows, centers = load(job[0])
            return run_job(job, rows, cfg, centers, fits, seed)
        except Exception as e:   # one failed fit must not stop the queue; it is reported
            return f"FAILED {'/'.join(job)}: {e}"

    with ThreadPoolExecutor(max_workers=cfg["compute"]["parallel_fits"]) as pool:
        for fut in as_completed([pool.submit(work, j) for j in jobs]):
            print(time.strftime("%H:%M:%S"), fut.result(), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or [v for v in VARIANTS
                          if (data_path("interim") / "phase3" / "rows" / f"{v}.parquet").exists()])
