"""Phase 2 sensitivity: Zheng et al. counts from every read, not only complete pairs.

In both Zheng et al. 2024 studies a large share of each run's reads were deposited without
their mate (an unpaired file next to _1/_2), so the main pipeline, which needs pairs to rebuild
the full window, never sees them. Every spot has its first read in the deposited files, so
counting first reads (complete pairs' R1 plus the unpaired reads) counts every spot once.

The libraries were sequenced in both orientations: a read starts at the forward or at the
reverse primer. Per run: each read is sorted by the primer it starts with and the primer
removed (cutadapt; reads with neither are dropped), truncated to 200 bases, expected-error
filtered (VSEARCH, maxEE 1, no Ns) and dereplicated; reads from the reverse primer are
reverse-complemented, so they show the window's last 200 bases. Reads are mapped at 97% to the
study's ASVs cut to their first (forward reads) or last (reverse reads) 200 bases. ASVs that
share those bases are one group; a group's count is shared among its ASVs in proportion to
their paired-read counts in the same run, and the two ends are summed (D-030).

Usage: python scripts/23_zheng_forward_reads.py [STUDY_ID ...]   (default: both Zheng studies)
Outputs (data/interim/phase2/{study}/): counts_r1all.csv, read_qc_r1all.csv, r1all/ (per run)
"""

import gzip
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from plate_truth.amplicon import (
    asv_end_groups,
    count_fastq_records,
    fasta_size_total,
    map_reads,
    otutab_to_long,
    paired_fastq_urls,
    reads_by_orientation,
    split_group_counts,
    window_end_reads,
)
from plate_truth.config import data_path, load_config
from plate_truth.http import download
from plate_truth.reference import read_fasta, write_fasta

LENGTH = 200
WORKERS = 4
STUDIES = ["zheng2024_soil", "zheng2024_sludge"]


def first_read_files(fastq_ftp: str) -> list[str]:
    """https URLs of every spot's first read: the _1 file and any unpaired file (never _2)."""
    files = [f if f.startswith("http") else "https://" + f for f in str(fastq_ftp).split(";")]
    pair = paired_fastq_urls(fastq_ftp)
    return [u for u in files if u not in pair[1:]]


def process_run(row, study: dict, groups_fa: dict, amp: dict, work) -> dict:
    """First reads of one run -> group counts per end (cached CSV) and read counts."""
    run = row["run_accession"]
    cache, qc_file = work / f"{run}.csv", work / f"{run}_qc.json"
    if cache.exists() and qc_file.exists():
        return pd.read_json(qc_file, typ="series").to_dict()
    tmp = work / run
    tmp.mkdir(parents=True, exist_ok=True)
    try:
        raw = [download(u, tmp / u.rsplit("/", 1)[1])
               for u in first_read_files(row["fastq_ftp"])]
        joined = tmp / "first_reads.fastq"
        with open(joined, "wb") as out:
            for f in raw:
                out.write(gzip.decompress(f.read_bytes()))
        fwd, rev = reads_by_orientation(joined, study["primer_fwd"], study["primer_rev"],
                                        amp["primer_error_rate"], tmp)
        qc = {"run": run, "status": "ok", "first_reads": count_fastq_records(joined),
              "start_fwd": count_fastq_records(fwd), "start_rev": count_fastq_records(rev)}
        parts = []
        for end, reads, rc in (("start", fwd, False), ("end", rev, True)):
            derep, tab = tmp / f"{end}.fa", tmp / f"{end}.tsv"
            window_end_reads(reads, LENGTH, amp["max_ee"], derep, tmp, revcomp=rc)
            qc[f"filtered_{end}"] = fasta_size_total(derep)
            map_reads(derep, groups_fa[end], tab, amp["map_id"], threads=1)
            parts.append(otutab_to_long(tab).rename(columns={"asv_id": "group_id"})
                         .assign(end=end))
        counts = pd.concat(parts, ignore_index=True)
        counts["run"] = run
        counts.to_csv(cache, index=False)
        qc["mapped"] = int(counts["count"].sum())
        pd.Series(qc).to_json(qc_file)
        return qc
    except Exception as exc:                                   # noqa: BLE001
        return {"run": run, "status": f"error: {exc}"}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def process_study(study: dict) -> None:
    """All runs of one study, then group counts shared out to ASVs."""
    amp = load_config("thresholds")["amplicon"]
    d = data_path("interim") / "phase2" / study["id"]
    work = d / "r1all"
    work.mkdir(exist_ok=True)
    asvs = [(h.split(";")[0], s) for h, s in read_fasta(d / "asv_seqs.fasta")]
    groups, groups_fa = {}, {}
    for end in ("start", "end"):
        groups[end] = asv_end_groups(asvs, LENGTH, end)
        groups_fa[end] = work / f"groups_{end}.fasta"
        write_fasta(groups[end].drop_duplicates("group_id")[["group_id", "part"]].itertuples(
            index=False, name=None), groups_fa[end])
    print(f"{study['id']}: {len(asvs)} ASVs; {LENGTH}-bp groups: "
          f"{groups['start']['group_id'].nunique()} by start, "
          f"{groups['end']['group_id'].nunique()} by end", flush=True)
    sheet = pd.read_csv(data_path("processed") / "sample_sheet.csv")
    sheet = sheet[sheet["study_id"] == study["id"]]
    urls = pd.read_csv(data_path("raw") / "census" / "ena" / f"{study['accession']}_fastq.tsv",
                       sep="\t")
    runs = sheet.merge(urls, on="run_accession")
    results = []
    with ThreadPoolExecutor(WORKERS) as pool:
        futs = [pool.submit(process_run, r, study, groups_fa, amp, work)
                for r in runs.to_dict("records")]
        for i, f in enumerate(as_completed(futs), 1):
            results.append(f.result())
            if i % 25 == 0:
                print(f"  {i}/{len(futs)} runs", flush=True)
    qc = pd.DataFrame(results)
    bad = qc[qc["status"] != "ok"]
    if len(bad):
        raise SystemExit(f"{len(bad)} runs failed:\n{bad.to_string()}")
    gc = pd.concat([pd.read_csv(work / f"{r}.csv") for r in runs["run_accession"]])
    if "end" not in gc.columns or gc["end"].isna().any():
        raise SystemExit("a cached run lacks the 'end' column (old format): delete r1all/")
    missing = set(runs["run_accession"]) - set(gc["run"])
    if missing:
        raise SystemExit(f"{len(missing)} runs have no counts: {sorted(missing)[:5]}")
    paired = pd.read_csv(d / "counts.csv")
    split = pd.concat([split_group_counts(gc[gc["end"] == end], groups[end], paired)
                       for end in ("start", "end")])
    split.groupby(["asv_id", "run"], as_index=False)["count"].sum().to_csv(
        d / "counts_r1all.csv", index=False)
    qc.to_csv(d / "read_qc_r1all.csv", index=False)
    print(f"{study['id']}: first reads {qc['first_reads'].sum():,}; mapped "
          f"{qc['mapped'].sum():,}", flush=True)


def main() -> None:
    studies = {s["id"]: s for s in load_config("census_studies")["studies"]}
    for sid in sys.argv[1:] or STUDIES:
        process_study(studies[sid])


if __name__ == "__main__":
    main()
