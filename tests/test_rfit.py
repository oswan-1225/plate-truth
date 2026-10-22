import json

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit

from plate_truth.phase3 import design_columns
from plate_truth.rfit import fixed_formula, rscript_exe, run_fit

CFG = {"identity_spline": {"boundary_knots": [75.0, 100.0],
                           "interior_knots": [90.0, 94.5, 97.0, 99.0]}}

try:
    rscript_exe()
    HAVE_R = True
except FileNotFoundError:
    HAVE_R = False


def _rows(n: int, rng) -> pd.DataFrame:
    study = rng.integers(0, 6, n)
    ident = np.clip(75 + rng.exponential(12, n), 75, 100)
    medium = rng.choice(["a", "b", "c"], n)
    source = np.char.add(study.astype(str), rng.integers(0, 4, n).astype(str))
    cond = np.char.add(source, rng.integers(0, 4, n).astype(str))
    eff = {"a": 0.5, "b": 0.0, "c": -0.5}
    eta = -6 + 0.25 * (ident - 75) + np.vectorize(eff.get)(medium) + rng.normal(0, 0.3, 6)[study]
    return pd.DataFrame({
        "y": (rng.random(n) < expit(eta)).astype(int), "w": 1.0, "eta_true": eta,
        "study_id": study.astype(str), "source_key": source,
        "condition_id": cond, "study_group": study.astype(str),
        "identity_model": ident, "no_hit": 0.0, "log_env": 0.0, "log_depth": 0.0,
        "log2_days": 0.0, "host_associated": 0, "medium": medium, "temperature": "moderate",
        "oxygen": "aerobic", "phylum": "x"})


@pytest.mark.skipif(not HAVE_R, reason="R not available")
def test_r_fit_recovers_truth_and_predicts_held_out(tmp_path):
    rng = np.random.default_rng(3)
    rows = _rows(6000, rng)
    train = rows["study_group"] != "0"
    x = design_columns(rows, train, {"log_env": 0, "log_depth": 0, "log2_days": 0})
    data = pd.concat([rows[["y", "w", "study_id", "source_key", "condition_id"]], x], axis=1)
    data["row_id"] = rows.index
    data["train"] = train.astype(int)
    grid = pd.DataFrame(0.0, index=range(3), columns=x.columns)
    grid["identity"] = [80.0, 95.0, 100.0]
    info = run_fit(tmp_path / "fit", data, fixed_formula("primary", list(x.columns), CFG),
                   grid, n_draws=50, seed=1)
    assert info["status"] == "ok" and info["converged"]
    eta = pd.read_csv(tmp_path / "fit" / "test_eta.csv")
    assert len(eta) == (~train).sum()
    truth = rows.loc[eta["row_id"], "eta_true"].to_numpy()
    assert np.corrcoef(eta["eta"], truth)[0, 1] > 0.95
    coef = pd.read_csv(tmp_path / "fit" / "coef.csv").set_index("term")["estimate"]
    assert coef["medium_a"] == pytest.approx(0.5, abs=0.25)   # deviation coding: a vs average
    draws = pd.read_csv(tmp_path / "fit" / "grid_draws.csv")
    assert draws.shape == (3, 50)
    assert json.loads((tmp_path / "fit" / "fit.json").read_text())["n_test"] == (~train).sum()
    assert not (tmp_path / "fit" / "input.csv").exists()
