"""Phase 0 sensitivity grid. Every scenario re-runs the same code as the primary analysis
with one choice changed, at the (ASV, soil) unit, held-out identity-bin model.

Scenarios:
  thresholds    enriched: min reads {2,5,10} x fold {1.5,2,5}; detected: min reads {2,5,10}
  identity      id1 (primary) / id2 (terminal gaps ignored) / mid (all gaps ignored)
  reference     all sources (contaminated assemblies removed) / all sources unfiltered /
                BacDive strains only (pilot) / RefSeq type material only
  environment   original soil (primary) / plated inoculum (filtrate)
  copy_number   off (primary) / relative abundances divided by rrnDB 16S copy number
  contaminants  kept and flagged (primary) / paper's candidate-contaminant ASVs excluded

Output: results/phase0/sensitivity.csv (one row per scenario: held-out AUC and Brier,
observed rate per identity bin, proxy over-estimate factor).
"""

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.phase0 import attach_identity, evaluate, growth_inputs, label_growth, per_source


def run(counts, sheet, nr, th, *, outcome="enriched_on_plates", grew_reads=None, fold=None,
        det_reads=None, id_col="id1", ref_set="all", env_kind="environment", weight=None,
        exclude=None) -> dict:
    """One scenario; arguments left as None take the primary value from the config."""
    grew_reads = grew_reads or th["grew"]["min_reads"]
    fold = fold or th["grew"]["min_fold_enrichment"]
    det_reads = det_reads or th["detected"]["min_reads"]
    env, plates, conds = growth_inputs(counts, sheet, env_kind, weight, exclude)
    table = label_growth(env, plates, conds, th["environment"]["min_reads"], grew_reads, fold,
                         det_reads)
    unit = attach_identity(per_source(table), nr, ref_set, id_col,
                           th["identity"]["bin_edges_percent"])
    res = evaluate(unit, outcome, th["identity"]["database_proxy_percent"], with_models=False)
    perf = res["performance"]["bins"]
    rec = res["recovery"]
    row = {"outcome": outcome, "grew_min_reads": grew_reads, "min_fold": fold,
           "detected_min_reads": det_reads, "identity_def": id_col, "reference_set": ref_set,
           "environment": env_kind, "copy_number": weight is not None,
           "contaminants_excluded": bool(exclude), "n": perf["n"], "events": perf["events"],
           "heldout_auc": perf["auc"], "heldout_brier": perf["brier"],
           "observed_overall": perf["events"] / perf["n"],
           "proxy_over_observed": rec["proxy_has_relative"].mean() / rec["observed_grew"].mean()}
    for b, r in res["calibration"]["observed"].items():
        row[f"rate_{b}"] = r
    return row


def main() -> None:
    th = load_config("thresholds")
    proc = data_path("processed")
    sheet = pd.read_csv(proc / "zheng2024_samples.csv")
    counts = pd.read_csv(data_path("interim") / "zheng2024" / "counts.csv")
    nr = pd.read_csv(proc / "zheng2024_nearest_relative.csv")
    cont = set(pd.read_csv(proc / "zheng2024_contaminant_asvs.csv")["asv_id"])
    tax = pd.read_csv(proc / "zheng2024_taxonomy.csv", index_col="asv_id")
    cn = tax["copy_number"]
    base = {"counts": counts, "sheet": sheet, "nr": nr, "th": th}

    rows = [dict(run(**base), scenario="primary")]
    for r in th["sensitivity"]["min_reads"]:
        for f in th["sensitivity"]["min_fold_enrichment"]:
            rows.append(dict(run(**base, grew_reads=r, fold=f), scenario="thresholds"))
        rows.append(dict(run(**base, outcome="detected_on_plates", det_reads=r),
                         scenario="thresholds"))
    for id_col in ["id2", "mid"]:
        rows.append(dict(run(**base, id_col=id_col), scenario="identity"))
    for ref in ["all_unfiltered", "bacdive_strain", "refseq_type"]:
        rows.append(dict(run(**base, ref_set=ref), scenario="reference"))
        rows.append(dict(run(**base, ref_set=ref, outcome="detected_on_plates"),
                         scenario="reference"))
    rows.append(dict(run(**base, env_kind="inoculum"), scenario="environment"))
    rows.append(dict(run(**base, weight=cn), scenario="copy_number"))
    rows.append(dict(run(**base, exclude=cont), scenario="contaminants"))
    out = pd.DataFrame(rows)
    out.to_csv(data_path("results") / "phase0" / "sensitivity.csv", index=False)
    show = ["scenario", "outcome", "grew_min_reads", "min_fold", "detected_min_reads",
            "identity_def", "reference_set", "environment", "copy_number",
            "contaminants_excluded", "heldout_auc", "observed_overall", "proxy_over_observed"]
    print(out[show + [c for c in out if c.startswith("rate_")]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
