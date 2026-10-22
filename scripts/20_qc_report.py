"""Phase 2, step 7: per-study QC report and the clean-data check before Phase 3.

Reads data/plate_truth.sqlite (script 19). For every study: reads retained through the
pipeline, ASVs (all and prokaryotic), environment-plate links, control libraries, and what
survives the minimum-depth rule (D-027) at the primary threshold and each sensitivity
threshold. Outcome counts (enriched/detected rows) are descriptive, in-sample, and are
reported only to show how many events each study contributes; they are not performance.

Outputs (results/phase2/):
  qc_read_retention.csv   per-study read totals at each step
  qc_asvs_links.csv       per-study ASVs, source samples, conditions and link counts
  qc_controls.csv         one row per control library
  qc_pool_depths.csv      one row per linked plating condition with both pool depths
  usable_data.csv         per study x depth threshold: what survives the rule
"""

import pandas as pd

from plate_truth.config import load_config, project_root
from plate_truth.growth import pooled_abundance
from plate_truth.phase0 import label_growth
from plate_truth.qc import pool_depth, read_retention, usable_conditions, usable_summary
from plate_truth.queries import connect, environment_counts, plate_counts, plating_conditions


def growth_rows(env: pd.DataFrame, plates: pd.DataFrame, conds: pd.DataFrame,
                th: dict) -> pd.DataFrame:
    """(ASV, condition) rows with outcomes, built one study at a time (source-sample names
    are only unique within a study)."""
    parts = []
    for sid, c in conds.groupby("study_id"):
        e = pooled_abundance(env[env["study_id"] == sid], "source_sample")
        p = pooled_abundance(plates[plates["study_id"] == sid], "condition_id")
        t = label_growth(e, p, c[["condition_id", "source_sample", "medium", "oxygen"]],
                         th["environment"]["min_reads"], th["grew"]["min_reads"],
                         th["grew"]["min_fold_enrichment"], th["detected"]["min_reads"])
        parts.append(t.assign(study_id=sid))
    return pd.concat(parts, ignore_index=True)


def asv_link_table(con, conds: pd.DataFrame) -> pd.DataFrame:
    """Per-study ASV counts and environment-plate link counts."""
    asvs = pd.read_sql("""SELECT a.study_id, COUNT(*) AS asvs,
                          SUM(t.is_organelle = 0) AS asvs_prokaryotic
                          FROM asvs a JOIN taxonomy t USING (asv_id) GROUP BY a.study_id""", con)
    runs = pd.read_sql("""SELECT study_id, kind, COUNT(*) AS n FROM samples
                          GROUP BY study_id, kind""", con)
    runs = runs.pivot(index="study_id", columns="kind", values="n").fillna(0).astype(int)
    runs.columns = [f"runs_{k}" for k in runs.columns]
    allc = plating_conditions(con, linked_only=False)
    links = pd.DataFrame({
        "conditions": allc.groupby("study_id").size(),
        "conditions_linked": conds.groupby("study_id").size(),
        "sources_linked": conds.groupby("study_id")["source_sample"].nunique(),
    })
    return asvs.set_index("study_id").join(runs).join(links).fillna(0).reset_index()


def controls_table(con) -> pd.DataFrame:
    """Control libraries: reads mapped, prokaryotic reads and ASVs with any read."""
    return pd.read_sql("""
        SELECT s.study_id, s.run_accession, s.library_name, s.reads_in, s.reads_mapped,
               SUM(CASE WHEN t.is_organelle = 0 THEN a.count ELSE 0 END) AS reads_prokaryotic,
               COUNT(a.asv_id) AS asvs
        FROM samples s LEFT JOIN abundances a USING (run_accession)
        LEFT JOIN taxonomy t USING (asv_id)
        WHERE s.kind = 'control' GROUP BY s.run_accession ORDER BY s.study_id""", con)


def main() -> None:
    th = load_config("thresholds")
    out = project_root() / "results" / "phase2"
    out.mkdir(parents=True, exist_ok=True)
    con = connect()

    qc = pd.read_sql("SELECT * FROM samples WHERE reads_in IS NOT NULL", con)
    read_retention(qc).to_csv(out / "qc_read_retention.csv", index=False)
    conds = plating_conditions(con, linked_only=True)
    asv_link_table(con, conds).to_csv(out / "qc_asvs_links.csv", index=False)
    controls_table(con).to_csv(out / "qc_controls.csv", index=False)

    env, plates = environment_counts(con), plate_counts(con)
    env_depth = pool_depth(env, ["study_id", "source_sample"])
    plate_depth = pool_depth(plates, ["study_id", "condition_id"])
    rows = growth_rows(env, plates, conds, th)

    keys = conds[["study_id", "source_sample", "condition_id"]]
    usable_conditions(keys, plate_depth, env_depth, th["depth"]["min_pool_reads"]).to_csv(
        out / "qc_pool_depths.csv", index=False)
    summaries = []
    for n in [0, th["depth"]["min_pool_reads"], *th["depth"]["sensitivity"]]:
        flags = usable_conditions(keys, plate_depth, env_depth, n)
        summaries.append(usable_summary(flags, rows).assign(
            min_pool_reads=n, primary=n == th["depth"]["min_pool_reads"]))
    usable = pd.concat(summaries, ignore_index=True)
    usable.to_csv(out / "usable_data.csv", index=False)
    pd.set_option("display.width", 200)
    print(usable.to_string(index=False))


if __name__ == "__main__":
    main()
