"""Observed growth on plates: the 'enriched' (primary) and 'detected' (loose) definitions.

Unit: (ASV, plating condition), restricted to ASVs present in the matching environmental
sample. Plates in one condition are pooled by summing reads; environmental PCR replicates of
one source sample are pooled the same way.

  detected_on_plates : plate_reads >= detected_min_reads
  enriched_on_plates : plate_reads >= grew_min_reads AND plate_rel / env_rel >= min_fold
"""

import numpy as np
import pandas as pd


def pooled_abundance(counts: pd.DataFrame, group_col: str, weight: pd.Series | None = None
                     ) -> pd.DataFrame:
    """Sum reads per (group, asv) and convert to relative abundance within each group.

    counts has columns [group_col, 'asv_id', 'count']. If weight (indexed by asv_id) is
    given, relative abundance uses count / weight, e.g. 16S copy number, while 'reads'
    stays the raw read count.
    """
    g = counts.groupby([group_col, "asv_id"], as_index=False)["count"].sum()
    g = g.rename(columns={"count": "reads"})
    signal = g["reads"] / g["asv_id"].map(weight).fillna(1.0) if weight is not None \
        else g["reads"].astype(float)
    g["rel"] = signal / signal.groupby(g[group_col]).transform("sum")
    return g


def fold_enrichment(plate_rel: pd.Series, env_rel: pd.Series) -> pd.Series:
    """plate_rel / env_rel; infinite when the ASV is absent from the environment."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return plate_rel / env_rel


def is_detected(plate_reads: pd.Series, min_reads: int) -> pd.Series:
    """Loose definition: at least min_reads reads in the plating condition."""
    return plate_reads >= min_reads


def is_enriched(plate_reads: pd.Series, fold: pd.Series, min_reads: int,
                min_fold: float) -> pd.Series:
    """Primary definition: enough reads AND relative abundance raised by >= min_fold."""
    return (plate_reads >= min_reads) & (fold >= min_fold)


def asv_condition_table(env: pd.DataFrame, plates: pd.DataFrame, conditions: pd.DataFrame,
                        env_min_reads: int) -> pd.DataFrame:
    """One row per (environment-present ASV, plating condition of that source sample).

    env:        source_sample, asv_id, reads, rel     (pooled environment)
    plates:     condition_id, asv_id, reads, rel      (pooled plates per condition)
    conditions: condition_id, source_sample, medium, oxygen
    ASVs with fewer than env_min_reads environment reads are not 'present' and are left out;
    ASVs absent from a condition's plates get plate_reads 0.
    """
    present = env[env["reads"] >= env_min_reads].rename(
        columns={"reads": "env_reads", "rel": "env_rel"})
    grid = present.merge(conditions, on="source_sample")
    p = plates.rename(columns={"reads": "plate_reads", "rel": "plate_rel"})
    grid = grid.merge(p, on=["condition_id", "asv_id"], how="left")
    grid[["plate_reads", "plate_rel"]] = grid[["plate_reads", "plate_rel"]].fillna(0)
    grid["plate_reads"] = grid["plate_reads"].astype(int)
    grid["fold"] = fold_enrichment(grid["plate_rel"], grid["env_rel"])
    return grid


def plate_only_asvs(env: pd.DataFrame, plates: pd.DataFrame, conditions: pd.DataFrame,
                    env_min_reads: int, plate_min_reads: int) -> pd.DataFrame:
    """ASVs with >= plate_min_reads on a condition's plates but not present in that source
    sample's environment (rare-biosphere growth, or contamination). Reported separately."""
    p = plates.merge(conditions[["condition_id", "source_sample"]], on="condition_id")
    p = p[p["reads"] >= plate_min_reads]
    e = env[["source_sample", "asv_id", "reads"]].rename(columns={"reads": "env_reads"})
    p = p.merge(e, on=["source_sample", "asv_id"], how="left")
    p["env_reads"] = p["env_reads"].fillna(0).astype(int)
    return p[p["env_reads"] < env_min_reads].reset_index(drop=True)


def any_condition(table: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """Collapse (ASV, condition) rows to (ASV, source sample): outcome true in any condition."""
    keys = ["source_sample", "asv_id"]
    agg = table.groupby(keys, as_index=False).agg(
        env_reads=("env_reads", "first"), env_rel=("env_rel", "first"),
        n_conditions=("condition_id", "nunique"), outcome=(outcome, "any"))
    return agg.rename(columns={"outcome": outcome})
