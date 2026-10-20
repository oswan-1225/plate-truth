"""Zheng et al. 2024 soil runs: raw reads -> ASV table.

Per run (parallel, cached): download the paired FASTQs from ENA, trim primers (cutadapt),
merge, filter (maxEE), dereplicate. Raw FASTQs are deleted after each run; the per-run
dereplicated FASTA is kept in data/interim/zheng2024/derep/.
Pooled: UNOISE3 denoising + uchime3 de novo chimera removal over all soil runs together,
then every run's reads are mapped back to the ASVs.

Outputs (data/interim/zheng2024/): asv_seqs.fasta, counts.csv (run, asv_id, count),
read_qc.csv (reads in, after primer trimming, after merging+filtering, mapped).
"""

import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from plate_truth.amplicon import (
    count_fastq_records,
    denoise,
    fasta_size_total,
    map_reads,
    merge_filter_derep,
    otutab_to_long,
    paired_fastq_urls,
    trim_primers,
)
from plate_truth.config import data_path, load_config
from plate_truth.http import download

STUDY = "zheng2024_soil"
WORKERS = 4


def process_run(row: pd.Series, primers: dict, amp: dict, out_dir, tmp_root) -> dict:
    """Download, trim, merge, filter and dereplicate one run; return read counts."""
    run = row["run_accession"]
    derep = out_dir / "derep" / f"{run}.fasta"
    qc_file = out_dir / "qc" / f"{run}.json"
    if derep.exists() and qc_file.exists():
        return pd.read_json(qc_file, typ="series").to_dict()
    urls = paired_fastq_urls(row["fastq_ftp"])
    if not urls:
        return {"run": run, "status": "no_paired_fastq"}
    tmp = tmp_root / run
    tmp.mkdir(parents=True, exist_ok=True)
    try:
        raw = [download(u, tmp / f"raw_{i}.fastq.gz") for i, u in enumerate(urls, 1)]
        t1, t2 = tmp / "t1.fq", tmp / "t2.fq"
        rep = trim_primers(raw[0], raw[1], t1, t2, primers["forward"], primers["reverse"],
                           amp["primer_error_rate"], tmp / "cutadapt.json")
        merge_filter_derep(t1, t2, tmp, run, derep, amp["max_ee"], amp["min_len"])
        qc = {"run": run, "status": "ok",
              "pairs_in": rep["read_counts"]["input"],
              "pairs_primer_trimmed": count_fastq_records(t1),
              "reads_filtered": fasta_size_total(derep)}
        pd.Series(qc).to_json(qc_file)
        return qc
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> None:
    study = next(s for s in load_config("studies")["studies"] if s["id"] == STUDY)
    amp = load_config("thresholds")["amplicon"]
    sheet = pd.read_csv(data_path("processed") / "zheng2024_samples.csv")
    out_dir = data_path("interim") / "zheng2024"
    tmp_root = out_dir / "tmp"
    for sub in ("derep", "qc"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    results = []
    with ThreadPoolExecutor(WORKERS) as pool:
        futures = {pool.submit(process_run, r, study["primers"], amp, out_dir, tmp_root):
                   r["run_accession"] for _, r in sheet.iterrows()}
        for i, fut in enumerate(as_completed(futures), 1):
            results.append(fut.result())
            if i % 20 == 0:
                print(f"  {i}/{len(futures)} runs processed", flush=True)
    qc = pd.DataFrame(results)
    if (qc["status"] != "ok").any():
        print(qc[qc["status"] != "ok"])
    if "--per-run-only" in sys.argv:
        return

    ok = qc.loc[qc["status"] == "ok", "run"]
    pooled = out_dir / "all_derep.fasta"
    with pooled.open("w") as handle:
        for run in ok:
            handle.write((out_dir / "derep" / f"{run}.fasta").read_text())
    work = out_dir / "work"
    work.mkdir(exist_ok=True)
    asvs = out_dir / "asv_seqs.fasta"
    denoise(pooled, work, asvs, amp["unoise_minsize"], amp["unoise_alpha"], "z24_asv_")
    map_reads(pooled, asvs, work / "otutab.tsv", amp["map_id"])
    counts = otutab_to_long(work / "otutab.tsv")
    counts.to_csv(out_dir / "counts.csv", index=False)

    mapped = counts.groupby("run")["count"].sum().rename("reads_mapped")
    qc = qc.merge(mapped, left_on="run", right_index=True, how="left")
    qc.to_csv(out_dir / "read_qc.csv", index=False)
    pooled.unlink()
    shutil.rmtree(work)
    print(f"ASVs: {counts['asv_id'].nunique()}; runs: {counts['run'].nunique()}")
    print(qc[["pairs_in", "pairs_primer_trimmed", "reads_filtered", "reads_mapped"]].describe())


if __name__ == "__main__":
    main()
