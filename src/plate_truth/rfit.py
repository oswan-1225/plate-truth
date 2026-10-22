"""Run the Phase 3 mixed-model fits in R (src/R/fit_glmm.R).

Each fit is one folder holding spec.json, its outputs and, last, fit.json. A folder with
fit.json is complete and is skipped on a rerun, so a stopped run resumes where it left off.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from plate_truth.config import project_root

RANDOM = ["study_id", "source_key", "condition_id"]
FIT_SCRIPT = project_root() / "src" / "R" / "fit_glmm.R"


def spline_term(cfg: dict) -> str:
    """R natural-spline term for identity with the knots fixed in config/phase3.yaml."""
    s = cfg["identity_spline"]
    knots = ", ".join(f"{k:g}" for k in s["interior_knots"])
    lo, hi = s["boundary_knots"]
    return f"ns(identity, knots = c({knots}), Boundary.knots = c({lo:g}, {hi:g}))"


def fixed_formula(model: str, columns: list[str], cfg: dict) -> str:
    """Right-hand side of the fixed effects for a named model (ANALYSIS_PLAN.md s4).

    columns: the design columns available (effect-coded factor columns included)."""
    ident = [spline_term(cfg), "no_hit"]
    covars = [c for c in columns if c not in ("identity", "no_hit", "host")]
    if model == "primary":
        terms = ident + covars
    elif model == "identity_only":
        terms = ident
    elif model == "no_identity":
        terms = covars
    elif model == "habitat":
        terms = ident + covars + ["host", f"host:{spline_term(cfg)}", "host:no_hit"]
    else:
        raise ValueError(f"unknown model {model!r}")
    return " + ".join(terms)


def rscript_exe() -> str:
    """Rscript from the active conda environment (Windows layout), else from PATH."""
    cand = Path(sys.executable).parent / "Lib" / "R" / "bin" / "x64" / "Rscript.exe"
    if cand.exists():
        return str(cand)
    found = shutil.which("Rscript")
    if found:
        return found
    raise FileNotFoundError("Rscript not found in the environment")


def r_env() -> dict:
    """Environment for R with the conda DLL folders on PATH (needed on Windows)."""
    env = os.environ.copy()
    root = Path(sys.executable).parent
    extra = [root, root / "Library" / "mingw-w64" / "bin", root / "Library" / "bin",
             root / "Lib" / "R" / "bin" / "x64"]
    env["PATH"] = os.pathsep.join(str(p) for p in extra) + os.pathsep + env.get("PATH", "")
    return env


def is_done(outdir: Path) -> bool:
    """A fit is complete when its fit.json exists."""
    return (outdir / "fit.json").exists()


def run_fit(outdir: Path, data: pd.DataFrame, fixed: str, grid: pd.DataFrame | None = None,
            n_draws: int = 0, seed: int = 0) -> dict:
    """Write the inputs, run R, delete the (large) input CSV and return fit.json."""
    outdir.mkdir(parents=True, exist_ok=True)
    inp = outdir / "input.csv"
    data.to_csv(inp, index=False)
    spec = {"input": str(inp), "outdir": str(outdir), "fixed": fixed, "random": RANDOM,
            "n_draws": n_draws, "seed": seed}
    if grid is not None:
        grid.to_csv(outdir / "grid.csv", index=False)
        spec["grid"] = str(outdir / "grid.csv")
    (outdir / "spec.json").write_text(json.dumps(spec, indent=1))
    res = subprocess.run([rscript_exe(), str(FIT_SCRIPT), str(outdir / "spec.json")],
                         env=r_env(), capture_output=True, text=True)
    (outdir / "r_log.txt").write_text(res.stdout + res.stderr)
    inp.unlink(missing_ok=True)
    if not is_done(outdir):
        raise RuntimeError(f"R fit failed in {outdir}; see r_log.txt")
    return json.loads((outdir / "fit.json").read_text())
