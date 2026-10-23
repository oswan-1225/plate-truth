"""Phase 4, step 3: fit the mixed models (ANALYSIS_PLAN.md sections 13-15).

Models (fixed effects as the Phase 3 primary model unless stated):
  q5_family        + random family, family x study, family x medium class
  q5_family_null   the same without family x medium class (H5a likelihood-ratio test)
  q5_phylum(_null) the same with phylum in place of family
  q6_match         Q6 rows; + recorded-class shares rec_<class> and match (the share recorded
                   on the row's own class); random phylum as in Phase 3
  q6_nomatch       the same without match
  q6close_*, q6binary_*, q6s16_*   sensitivity analyses (identity >= 97%; indicators instead
                   of shares; strains' own 16S links only)
  q7_notaxon       Phase 3 rows; the primary model without phylum (Q7)
Each fit is a fold: a held-out study group, or 'all' (full data). Fits are resumable (a folder
with fit.json is skipped) and stoppable (create data/interim/phase4/STOP).

Usage: python scripts/35_phase4_fit.py [model ...]   (default: every model, in priority order)
"""

import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import CONTINUOUS, GROUP, design_columns
from plate_truth.phase4 import match_term
from plate_truth.rfit import RANDOM, is_done, run_fit, spline_term

Q5 = ["study_id", "source_key", "condition_id"]


def q5_spec(taxon: str, interaction: bool) -> dict:
    """Random effects of a Q5 model: taxon, taxon x study and (optionally) taxon x medium."""
    random = Q5 + [taxon, f"{taxon}_study"] + ([f"{taxon}_medium"] if interaction else [])
    known = [taxon] + ([f"{taxon}_medium"] if interaction else [])
    return {"rows": "q5", "random": random, "known": known, "extra": "none"}


def q6_spec(rows: str, extra: str, match: bool) -> dict:
    """A Q6 model: recorded-class terms (shares 'rec' or indicators 'any'), with or without
    the match term; random effects as the Phase 3 primary model."""
    return {"rows": rows, "random": RANDOM + ["phylum"], "known": ["phylum"],
            "extra": extra, "match": match}


# name -> spec; folds: 'loso' (each group held out) and/or 'all'
MODELS = {
    "q5_family": {**q5_spec("family", True), "folds": ["all", "loso"]},
    "q5_family_null": {**q5_spec("family", False), "folds": ["all", "loso"]},
    "q6_match": {**q6_spec("q6", "rec", True), "folds": ["all", "loso"]},
    "q6_nomatch": {**q6_spec("q6", "rec", False), "folds": ["all", "loso"]},
    "q7_notaxon": {"rows": "q7", "random": RANDOM, "known": [], "extra": "none",
                   "folds": ["loso"]},
    "q6close_match": {**q6_spec("q6close", "rec", True), "folds": ["all"]},
    "q6close_nomatch": {**q6_spec("q6close", "rec", False), "folds": ["all"]},
    "q6binary_match": {**q6_spec("q6", "any", True), "folds": ["all"]},
    "q6binary_nomatch": {**q6_spec("q6", "any", False), "folds": ["all"]},
    "q6s16_match": {**q6_spec("q6_strain16s", "rec", True), "folds": ["all"]},
    "q6s16_nomatch": {**q6_spec("q6_strain16s", "rec", False), "folds": ["all"]},
    "q5_phylum": {**q5_spec("phylum_full", True), "folds": ["all"]},
    "q5_phylum_null": {**q5_spec("phylum_full", False), "folds": ["all"]},
}


def load_rows(name: str, cfg: dict) -> pd.DataFrame:
    """Analysis rows of a row set (q5, q6, q6close, q6_strain16s, q7)."""
    root = data_path("interim")
    if name == "q7":
        return pd.read_parquet(root / "phase3" / "rows" / "primary.parquet")
    if name == "q6close":
        r = pd.read_parquet(root / "phase4" / "rows" / "q6.parquet")
        return r[r["identity"] >= cfg["q6"]["close_identity"]].reset_index(drop=True)
    return pd.read_parquet(root / "phase4" / "rows" / f"{name}.parquet")


def recorded_classes(rows: pd.DataFrame, classes: list[str], min_share: float) -> list[str]:
    """Classes whose recorded-class term is fitted: at least min_share of the rows' ASVs have
    a linked strain recorded on the class (a property of the records, not of outcomes)."""
    a = rows.drop_duplicates("asv_id")
    return [c for c in classes if (a[f"any_{c}"] > 0).mean() >= min_share]


def extra_columns(rows: pd.DataFrame, spec: dict, kept: list[str]) -> pd.DataFrame:
    """Recorded-class columns (and the match term) of a Q6 model."""
    if spec["extra"] == "none":
        return pd.DataFrame(index=rows.index)
    p = spec["extra"]
    x = pd.DataFrame({f"{p}_{c}": rows[f"{p}_{c}"] for c in kept}, index=rows.index)
    if spec["match"]:
        x["match"] = match_term(rows, kept, prefix=p)
    return x


def jobs_for(models: list[str], groups: dict) -> list[tuple[str, str]]:
    """(model, fold) in priority order: every model's full fit first, then the folds."""
    full = [(m, "all") for m in models if "all" in MODELS[m]["folds"]]
    loso = [(m, g) for m in models if "loso" in MODELS[m]["folds"]
            for g in groups[MODELS[m]["rows"]]]
    return full + loso


def run_job(job, rows: pd.DataFrame, centers: dict, p3: dict, cfg: dict, out_root) -> str:
    """Build one fit's inputs and run it in R."""
    model, fold = job
    spec = MODELS[model]
    outdir = out_root / model / fold
    if is_done(outdir):
        return f"skip {model}/{fold}"
    train = rows[GROUP] != fold if fold != "all" else pd.Series(True, index=rows.index)
    x = design_columns(rows, train, centers)
    kept = (recorded_classes(rows, cfg["medium_classes"], cfg["q6"]["min_share_linked"])
            if spec["extra"] != "none" else [])
    extra = extra_columns(rows, spec, kept)
    covars = [c for c in x.columns if c not in ("identity", "no_hit", "host")]
    fixed = " + ".join([spline_term(p3), "no_hit", *covars, *extra.columns])
    data = pd.concat([rows[["y", "w", *spec["random"]]], x, extra], axis=1).reset_index(
        drop=True)
    data["row_id"] = rows.index.to_numpy()
    data["train"] = train.astype(int).to_numpy()
    t0 = time.time()
    info = run_fit(outdir, data, fixed, None, 0, 0, spec["random"], spec["known"])
    return (f"done {model}/{fold} in {(time.time() - t0) / 60:.1f} min: {info['status']}, "
            f"converged={info.get('converged')}, pdHess={info.get('pdHess')}")


def main(models: list[str]) -> None:
    cfg = load_config("phase4")
    p3 = load_config("phase3")
    root = data_path("interim") / "phase4"
    fits, stop = root / "fits", root / "STOP"
    row_sets = sorted({MODELS[m]["rows"] for m in models})
    cache, lock = {}, threading.Lock()

    def load(name: str):
        with lock:
            if name not in cache:
                r = load_rows(name, cfg)
                cache[name] = (r, {c: float(r[c].median()) for c in CONTINUOUS})
            return cache[name]

    groups = {n: sorted(load(n)[0][GROUP].unique()) for n in row_sets}
    jobs = [j for j in jobs_for(models, groups) if not is_done(fits / j[0] / j[1])]
    print(f"{len(jobs)} fits queued", flush=True)

    def work(job):
        if stop.exists():
            return f"stopped before {'/'.join(job)}"
        try:
            rows, centers = load(MODELS[job[0]]["rows"])
            return run_job(job, rows, centers, p3, cfg, fits)
        except Exception as e:   # one failed fit must not stop the queue; it is reported
            return f"FAILED {'/'.join(job)}: {e}"

    with ThreadPoolExecutor(max_workers=cfg["compute"]["parallel_fits"]) as pool:
        for fut in as_completed([pool.submit(work, j) for j in jobs]):
            print(time.strftime("%H:%M:%S"), fut.result(), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:] or list(MODELS))
