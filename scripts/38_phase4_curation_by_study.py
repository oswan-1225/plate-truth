"""Phase 4, exploratory (not registered): is the curation match effect within-study?

The registered Q6 model pools each recorded-class effect (rec_X) over studies. Most dilute-plate
rows come from one study group, so a study-specific benefit of relatives kept on dilute media
could be absorbed by the match term. Two checks, each fitted with and without match:
  q6x_bystudy   rec_X effects estimated separately per study (rec_X x study)
  q6x_bystudy_no_sponges  the same without the sponge study, whose few events (79) separate
                one of its rec_X x study terms
  q6x_<group>   the registered Q6 model within each study group alone (per-study match)
Outputs: results/phase4/exploratory/curation_by_study.csv
Usage: python scripts/38_phase4_curation_by_study.py
"""

import json

import numpy as np
import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase3 import CONTINUOUS, GROUP, design_columns
from plate_truth.phase4 import match_term, recorded_by_study
from plate_truth.rfit import RANDOM, is_done, run_fit, spline_term


def loglik(info: dict) -> float:
    """Log-likelihood of a fit; NaN when R could not compute it."""
    v = info.get("logLik")
    return float(v) if isinstance(v, (int, float)) else np.nan


def fit(rows: pd.DataFrame, name: str, with_match: bool, by_study: bool, kept: list[str],
        p3: dict, out_root) -> dict:
    """One full-data fit of a Q6 variant; returns fit.json plus the match coefficient."""
    outdir = out_root / name / ("match" if with_match else "nomatch")
    if not is_done(outdir):
        centers = {c: float(rows[c].median()) for c in CONTINUOUS}
        train = pd.Series(True, index=rows.index)
        x = design_columns(rows, train, centers)
        covars = [c for c in x.columns if c not in ("identity", "no_hit", "host")]
        extra = (recorded_by_study(rows, kept) if by_study else
                 pd.DataFrame({f"rec_{c}": rows[f"rec_{c}"] for c in kept}, index=rows.index))
        if with_match:
            extra["match"] = match_term(rows, kept)
        random = RANDOM + ["phylum"] if rows["study_id"].nunique() > 1 else \
            ["source_key", "condition_id", "phylum"]
        data = pd.concat([rows[["y", "w", *random]], x, extra], axis=1).reset_index(drop=True)
        data["row_id"] = np.arange(len(data))
        data["train"] = 1
        fixed = " + ".join([spline_term(p3), "no_hit", *covars, *extra.columns])
        run_fit(outdir, data, fixed, None, 0, 0, random, ["phylum"])
    info = json.loads((outdir / "fit.json").read_text())
    if with_match and info["status"] == "ok":
        c = pd.read_csv(outdir / "coef.csv").set_index("term").loc["match"]
        info.update({"estimate": c["estimate"], "se": c["se"]})
    return info


def main() -> None:
    cfg, p3 = load_config("phase4"), load_config("phase3")
    out_root = data_path("interim") / "phase4" / "fits_exploratory"
    res = data_path("results") / "phase4" / "exploratory"
    res.mkdir(parents=True, exist_ok=True)
    rows = pd.read_parquet(data_path("interim") / "phase4" / "rows" / "q6.parquet")
    a = rows.drop_duplicates("asv_id")
    kept = [c for c in cfg["medium_classes"]
            if (a[f"any_{c}"] > 0).mean() >= cfg["q6"]["min_share_linked"]]
    no_sponges = rows[rows[GROUP] != "indraningrat2019"].reset_index(drop=True)
    jobs = [("bystudy", rows, True), ("bystudy_no_sponges", no_sponges, True)] + [(g, s.reset_index(drop=True), False)
                                        for g, s in rows.groupby(GROUP)]
    out = []
    for name, r, by_study in jobs:
        m = fit(r, f"q6x_{name}", True, by_study, kept, p3, out_root)
        n = fit(r, f"q6x_{name}", False, by_study, kept, p3, out_root)
        est, se = m.get("estimate", np.nan), m.get("se", np.nan)
        out.append({"analysis": name, "rows": len(r), "events": int(r["y"].sum()),
                    "match_rows_nonzero": int((match_term(r, kept) > 0).sum()),
                    "odds_ratio": np.exp(est), "or_lo": np.exp(est - 1.96 * se),
                    "or_hi": np.exp(est + 1.96 * se),
                    "lr_stat": 2 * (loglik(m) - loglik(n)),
                    "converged": m.get("converged"), "pdHess": m.get("pdHess"),
                    "nomatch_converged": n.get("converged"), "nomatch_pdHess": n.get("pdHess")})
        print(out[-1], flush=True)
    pd.DataFrame(out).to_csv(res / "curation_by_study.csv", index=False)


if __name__ == "__main__":
    main()
