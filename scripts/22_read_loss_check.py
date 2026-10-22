"""Phase 2 QC: is the read loss before the ASV table random with respect to taxa?

Paired-read studies lose pairs that do not merge (short or low-quality overlaps), and the two
Zheng et al. studies also have many reads deposited without their mate (an unpaired file per
run; most of those reads start at the reverse primer). If the lost reads come
disproportionately from some taxa, relative abundances are distorted, and with them the
enrichment ratio.

For a few runs per study (two environmental, two plate, closest to the study's median depth):
download all read files and merge pairs exactly as script 16 does. Kept reads: merged pairs
taken through the study's read path to the forward-oriented window. Lost reads: forward reads
of unmerged pairs plus unpaired reads. Reads are compared over the SAME bases: a lost read
that starts at the forward primer is compared on the window's first 150 bases, one that
starts at the reverse primer (reverse-complemented) on the last 150, and kept reads are cut
both ways. Each part is assigned to the group of ASVs sharing those bases (>= 97% identity).
Per group, within each run and each end: share among kept reads / share among kept + lost
(the 'distortion'; 1 = loss unrelated to the taxon). The same summary is computed for a
random-loss baseline (same reads and loss rate, lost reads drawn at random), so sampling noise
in small runs is not mistaken for bias.

Output: results/phase2/read_loss_check.csv (per study), read_loss_groups.csv.gz (per run, end
and group); per-run counts cached in data/interim/phase2/read_loss/
"""

import gzip
import shutil
import sys

import numpy as np
import pandas as pd

from plate_truth.amplicon import (
    asv_end_groups,
    cut_window,
    orient_reads,
    paired_fastq_urls,
    reads_by_orientation,
    reverse_complement,
)
from plate_truth.config import data_path, load_config, project_root
from plate_truth.http import download
from plate_truth.qc import abundance_distortion, bray_curtis
from plate_truth.reference import read_fasta, write_fasta
from plate_truth.tools import vsearch

PART = 150
RUNS_PER_KIND = 2
MIN_READS = 50          # groups with fewer reads (kept + lost) in a run are not summarised


def pick_runs(sheet: pd.DataFrame, n: int) -> pd.DataFrame:
    """n environment and n plate runs closest to the study's median ENA read count."""
    out = []
    for kind in ("environment", "plate"):
        k = sheet[sheet["kind"] == kind]
        med = k["read_count"].median()
        out.append(k.assign(d=(k["read_count"] - med).abs()).nsmallest(n, "d"))
    return pd.concat(out)


def fastq_seqs(path):
    """Sequences of a FASTQ file (plain text)."""
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip()


def assign(parts: list[str], groups_fa, tmp, tag: str) -> pd.Series:
    """Number of parts assigned to each ASV group (>= 97% identity, plus strand)."""
    if not parts:
        return pd.Series(dtype=int)
    fa, hits = tmp / f"{tag}.fa", tmp / f"{tag}.tsv"
    write_fasta(((f"{tag}{i}", s) for i, s in enumerate(parts)), fa)
    vsearch("--usearch_global", fa, "--db", groups_fa, "--id", 0.97, "--strand", "plus",
            "--maxaccepts", 1, "--maxrejects", 32, "--userfields", "query+target",
            "--userout", hits)
    h = pd.read_csv(hits, sep="\t", names=["query", "group_id"])
    return h["group_id"].value_counts()


def lost_parts(lost_fq, study: dict, ref_fa, tmp) -> tuple[list[str], list[str]]:
    """Lost reads as (window-start parts, window-end parts in forward orientation)."""
    if study["read_prep"] == "primers_in_reads":
        fwd, rev = reads_by_orientation(lost_fq, study["primer_fwd"], study["primer_rev"], 0.1,
                                        tmp)
        start = [s[:PART] for s in fastq_seqs(fwd) if len(s) >= PART]
        end = [reverse_complement(s[:PART]) for s in fastq_seqs(rev) if len(s) >= PART]
        return start, end
    oriented, tab = tmp / "lost_oriented.fq", tmp / "orient.tsv"     # primers already removed
    vsearch("--orient", lost_fq, "--db", ref_fa, "--fastqout", oriented, "--tabbedout", tab)
    strand = pd.read_csv(tab, sep="\t", header=None, usecols=[1])[1]
    seqs = list(fastq_seqs(oriented))
    plus = strand[strand.isin(["+", "-"])].tolist()
    start = [s[:PART] for s, st in zip(seqs, plus) if st == "+" and len(s) >= PART]
    end = [s[-PART:] for s, st in zip(seqs, plus) if st == "-" and len(s) >= PART]
    return start, end


def run_counts(row, study: dict, groups: dict, ref_fa, tmp) -> pd.DataFrame:
    """Kept and lost reads per (end, ASV group) for one run."""
    files = [f if f.startswith("http") else "https://" + f
             for f in str(row["fastq_ftp"]).split(";")]
    pair = paired_fastq_urls(row["fastq_ftp"])
    raw = [download(u, tmp / u.rsplit("/", 1)[1]) for u in pair + [f for f in files
                                                                    if f not in pair]]
    merged, lost = tmp / "m.fq", tmp / "lost.fq"
    vsearch("--fastq_mergepairs", raw[0], "--reverse", raw[1], "--fastq_maxdiffs", 10,
            "--fastq_allowmergestagger", "--fastqout", merged,
            "--fastqout_notmerged_fwd", lost)
    with open(lost, "ab") as out:                       # unpaired reads are lost too
        for s in raw[2:]:
            out.write(gzip.decompress(s.read_bytes()))
    window = tmp / "window.fq"
    if study["read_prep"] == "primers_in_reads":
        lo, hi = load_config("thresholds")["reference"]["region_length"][study["region"]]
        cut_window(merged, window, study["primer_fwd"], study["primer_rev"], 0.1, lo, hi,
                   tmp / "cut.json")
    else:
        orient_reads(merged, ref_fa, window)
    kept = [s for s in fastq_seqs(window) if len(s) >= PART]
    l_start, l_end = lost_parts(lost, study, ref_fa, tmp)
    rows = []
    for end, k_parts, l_parts in (("start", [s[:PART] for s in kept], l_start),
                                  ("end", [s[-PART:] for s in kept], l_end)):
        k = assign(k_parts, groups[end], tmp, f"kept_{end}")
        lo_ = assign(l_parts, groups[end], tmp, f"lost_{end}")
        t = pd.DataFrame({"kept": k, "lost": lo_}).fillna(0).astype(int)
        rows.append(t.rename_axis("group_id").reset_index().assign(end=end))
    return pd.concat(rows, ignore_index=True)


def run_counts_cached(row, study, groups, ref_fa) -> pd.DataFrame:
    """run_counts, cached as a small CSV per run (the downloads are deleted)."""
    work = data_path("interim") / "phase2" / "read_loss"
    cache = work / f"{row['run_accession']}.csv"
    if not cache.exists():
        tmp = work / row["run_accession"]
        tmp.mkdir(parents=True, exist_ok=True)
        try:
            run_counts(row, study, groups, ref_fa, tmp).to_csv(cache, index=False)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    return pd.read_csv(cache)


def study_groups(sid: str, work) -> dict:
    """ASV groups by their first and by their last PART bases, written as FASTA."""
    d = data_path("interim") / "phase2" / sid
    asvs = [(h.split(";")[0], s) for h, s in read_fasta(d / "asv_seqs.fasta")]
    out = {}
    for end in ("start", "end"):
        g = asv_end_groups(asvs, PART, end).drop_duplicates("group_id")
        out[end] = work / f"{sid}_groups_{end}.fasta"
        write_fasta(g[["group_id", "part"]].itertuples(index=False, name=None), out[end])
    return out


def summarise(sid: str, counts: list[pd.DataFrame], kinds: list[str]) -> tuple[dict, list]:
    """Per-study summary over runs and ends, and the per-group rows."""
    rng = np.random.default_rng(load_config("thresholds")["resampling"]["seed"])
    parts, nulls, bc = [], [], []
    for c, kind in zip(counts, kinds):
        for end, ce in c.groupby("end"):
            ce = ce.set_index("group_id")
            if ce["lost"].sum() == 0:
                continue
            d = abundance_distortion(ce["kept"], ce["lost"], MIN_READS)
            parts.append(d.assign(run=c.attrs["run"], kind=kind, end=end))
            bc.append(bray_curtis(ce["kept"], ce["lost"]))
            # random-loss baseline: same reads and loss rate, lost reads drawn at random
            total = ce["kept"] + ce["lost"]
            k = pd.Series(rng.binomial(total, ce["kept"].sum() / total.sum()), index=ce.index)
            nulls.append(abundance_distortion(k, total - k, MIN_READS))
    d, dn = pd.concat(parts), pd.concat(nulls)
    w, wn = d["total"] / d["total"].sum(), dn["total"] / dn["total"].sum()
    lost = sum(c["lost"].sum() for c in counts)
    reads = lost + sum(c.loc[c["end"] == "start", "kept"].sum() for c in counts)
    return ({"study_id": sid, "runs": len(counts), "lost_share": lost / reads,
             "bray_curtis_kept_vs_lost_median": float(pd.Series(bc).median()),
             "group_rows_summarised": len(d),
             "distortion_p05": d["distortion"].quantile(0.05),
             "distortion_p95": d["distortion"].quantile(0.95),
             "share_reads_distortion_within_1.25": w[(d["distortion"] >= 0.8) &
                                                     (d["distortion"] <= 1.25)].sum(),
             "random_loss_within_1.25": wn[(dn["distortion"] >= 0.8) &
                                           (dn["distortion"] <= 1.25)].sum(),
             "share_reads_distortion_beyond_2": w[(d["distortion"] < 0.5) |
                                                  (d["distortion"] > 2)].sum(),
             "random_loss_beyond_2": wn[(dn["distortion"] < 0.5) | (dn["distortion"] > 2)].sum()},
            d.assign(study_id=sid).reset_index(names="group_id"))


def main() -> None:
    studies = {s["id"]: s for s in load_config("census_studies")["studies"]
               if s["decision"] == "include" and s.get("read_prep") in
               ("primers_in_reads", "primers_removed")}
    sets = pd.read_csv(data_path("processed") / "primer_sets.csv")
    regions = data_path("interim") / "reference" / "regions"
    work = data_path("interim") / "phase2" / "read_loss"
    work.mkdir(parents=True, exist_ok=True)
    sheet = pd.read_csv(data_path("processed") / "sample_sheet.csv")
    summary, per_group = [], []
    for sid in [s for s in (sys.argv[1:] or list(studies)) if s in studies]:
        st = studies[sid]
        ref = next(r.ref_name for r in sets.itertuples() if sid in r.studies.split(";"))
        groups = study_groups(sid, work)
        urls = pd.read_csv(data_path("raw") / "census" / "ena" / f"{st['accession']}_fastq.tsv",
                           sep="\t")
        runs = pick_runs(sheet[sheet["study_id"] == sid], RUNS_PER_KIND).merge(
            urls, on="run_accession")
        counts = []
        for row in runs.to_dict("records"):
            # per-run read_prep: oak plates have primers removed, its environment does not
            st_run = {**st, "read_prep": st["read_prep_plate"] if row["kind"] == "plate"
                      and st.get("read_prep_plate") else st["read_prep"]}
            c = run_counts_cached(row, st_run, groups, regions / f"cultured_{ref}.fasta")
            c.attrs["run"] = row["run_accession"]
            counts.append(c)
        s, d = summarise(sid, counts, runs["kind"].tolist())
        summary.append(s)
        per_group.append(d)
        print(pd.Series(s).to_string(), flush=True)
    out = project_root() / "results" / "phase2"
    pd.DataFrame(summary).round(4).to_csv(out / "read_loss_check.csv", index=False)
    pd.concat(per_group).to_csv(out / "read_loss_groups.csv.gz", index=False)


if __name__ == "__main__":
    main()
