"""Per-study quality control and the minimum-depth rule for usable data (D-027).

The enrichment ratio compares two pools: a plating condition's plates (reads summed over
replicate plates) and its source sample's environment (reads summed over environmental runs).
A pool with very few reads gives a meaningless relative abundance (5 reads of 70 is 7%), so
both pools must reach a minimum number of prokaryotic reads. The rule is applied to pools,
not to single libraries, because pools are what the grew definition compares.
"""

import pandas as pd


def pool_depth(counts: pd.DataFrame, group_cols: list[str]) -> pd.Series:
    """Total reads per pool. counts has group_cols plus 'count'; returns a Series 'depth'
    indexed by group_cols."""
    return counts.groupby(group_cols)["count"].sum().rename("depth")


def usable_conditions(conditions: pd.DataFrame, plate_depth: pd.Series, env_depth: pd.Series,
                      min_reads: int) -> pd.DataFrame:
    """Flag each plating condition as usable when its plate pool and its source sample's
    environmental pool both have at least min_reads reads.

    conditions: study_id, source_sample, condition_id
    plate_depth: indexed by (study_id, condition_id); env_depth: by (study_id, source_sample).
    A missing pool counts as depth 0. Adds plate_depth, env_depth and usable columns.
    """
    out = conditions.copy()
    out["plate_depth"] = [plate_depth.get((s, c), 0) for s, c in
                          zip(out["study_id"], out["condition_id"])]
    out["env_depth"] = [env_depth.get((s, src), 0) for s, src in
                        zip(out["study_id"], out["source_sample"])]
    out["usable"] = (out["plate_depth"] >= min_reads) & (out["env_depth"] >= min_reads)
    return out


def read_retention(qc: pd.DataFrame) -> pd.DataFrame:
    """Per-study read totals at each pipeline step and the fraction of input reads mapped.

    qc: study_id, reads_in, reads_window, reads_filtered, reads_mapped (one row per run), and
    optionally read_count (ENA's count for the run). reads_in is what entered primer cutting or
    orientation: merged pairs for paired reads, reads otherwise; read_count is spots (pairs) for
    SRA submissions but can be single reads for ENA submissions, so the two are not a ratio.
    """
    steps = ["reads_in", "reads_window", "reads_filtered", "reads_mapped"]
    if "read_count" in qc.columns:
        steps = ["read_count", *steps]
    out = qc.groupby("study_id")[steps].sum()
    out.insert(0, "runs", qc.groupby("study_id").size())
    out["frac_mapped"] = out["reads_mapped"] / out["reads_in"]
    return out.reset_index()


def abundance_distortion(kept: pd.Series, lost: pd.Series, min_reads: int) -> pd.DataFrame:
    """How read loss changes each ASV's relative abundance.

    kept, lost: reads per ASV that reached the analysis and that were lost (e.g. pairs that
    did not merge). For ASVs with at least min_reads in total, distortion = share among kept
    reads / share among all reads: 1 means loss was unrelated to the ASV, 0.5 means its
    relative abundance is halved by the loss.
    """
    t = pd.DataFrame({"kept": kept, "lost": lost}).fillna(0)
    t["total"] = t["kept"] + t["lost"]
    share_total = t["total"] / t["total"].sum()
    share_kept = t["kept"] / t["kept"].sum()
    t["distortion"] = share_kept / share_total
    return t[t["total"] >= min_reads]


def bray_curtis(a: pd.Series, b: pd.Series) -> float:
    """Bray-Curtis dissimilarity of two count profiles after scaling each to proportions
    (0 = identical composition, 1 = nothing shared)."""
    p = (a / a.sum()).reindex(a.index.union(b.index), fill_value=0)
    q = (b / b.sum()).reindex(p.index, fill_value=0)
    return float((p - q).abs().sum() / 2)


def usable_summary(flags: pd.DataFrame, rows: pd.DataFrame) -> pd.DataFrame:
    """Per-study counts of what survives the depth rule.

    flags: output of usable_conditions (one row per linked condition).
    rows: (ASV, condition) table with study_id, condition_id, asv_id, source_sample and the
          boolean outcome columns enriched_on_plates and detected_on_plates.
    Returns conditions and source samples (total and usable), and for usable conditions only:
    distinct environment-present ASVs, (ASV, condition) rows and outcome counts.
    """
    g = flags.groupby("study_id")
    out = pd.DataFrame({
        "conditions": g.size(),
        "conditions_usable": g["usable"].sum(),
        "sources": g["source_sample"].nunique(),
        "sources_usable": flags[flags["usable"]].groupby("study_id")["source_sample"].nunique(),
    })
    keep = rows.merge(flags.loc[flags["usable"], ["study_id", "condition_id"]],
                      on=["study_id", "condition_id"])
    k = keep.groupby("study_id")
    out["asvs_usable"] = k["asv_id"].nunique()
    out["rows_usable"] = k.size()
    out["enriched_rows"] = k["enriched_on_plates"].sum()
    out["detected_rows"] = k["detected_on_plates"].sum()
    return out.fillna(0).astype(int).reset_index()
