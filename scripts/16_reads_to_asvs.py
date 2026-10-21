"""Phase 2, step 3: raw reads -> ASV table, one study at a time, one pipeline for all.

Per run (parallel, cached): download the run's FASTQ files from ENA, prepare the reads
according to the study's read_prep (D-024), filter by expected errors and the region's length
window, dereplicate. Raw FASTQs are deleted after each run.
  primers_in_reads   merge pairs, then keep the window between the primers on either strand
  primers_removed    merge pairs, then orient against the study's cultured reference
  pacbio             single reads, cut to the window of the study's reference primer pair (V4)
Per study (pooled over its runs): UNOISE3 + uchime3 de novo chimera removal, then every run's
reads are mapped back to the ASVs.

Usage: python scripts/16_reads_to_asvs.py STUDY_ID [STUDY_ID ...] [--per-run-only]
Outputs (data/interim/phase2/{study}/): derep/, qc/, asv_seqs.fasta, counts.csv, read_qc.csv
"""

import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from plate_truth.amplicon import (
    count_fastq_records,
    cut_window,
    denoise,
    fasta_size_total,
    fastq_urls,
    filter_derep,
    map_reads,
    merge_pairs,
    orient_reads,
    otutab_to_long,
)
from plate_truth.census import ena_table
from plate_truth.config import data_path, load_config
from plate_truth.http import download

WORKERS = 4


def prepare_reads(raw: list, mode: str, ref: dict, amp: dict, lo: int, hi: int, tmp) -> tuple:
    """Turn the downloaded files into window reads (FASTQ). Returns (path, reads_in)."""
    if mode == "pacbio":
        reads, n_in = raw[0], None
    else:
        if len(raw) != 2:
            raise ValueError(f"{mode} needs paired files, got {len(raw)}")
        reads = tmp / "merged.fq"
        merge_pairs(raw[0], raw[1], reads)
        n_in = count_fastq_records(reads)
    out = tmp / "window.fq"
    if mode in ("pacbio", "primers_in_reads"):
        rep = cut_window(reads, out, ref["primer_fwd"], ref["primer_rev"],
                         amp["primer_error_rate"], lo, hi, tmp / "cutadapt.json")
        n_in = n_in if n_in is not None else rep["read_counts"]["input"]
    elif mode == "primers_removed":
        orient_reads(reads, ref["fasta"], out)
    else:
        raise ValueError(f"unknown read_prep {mode!r}")
    return out, n_in


def process_run(row, ref: dict, amp: dict, lo: int, hi: int, out_dir, tmp_root) -> dict:
    """Download, prepare, filter and dereplicate one run; return read counts."""
    run, mode = row["run_accession"], row["read_prep"]
    derep = out_dir / "derep" / f"{run}.fasta"
    qc_file = out_dir / "qc" / f"{run}.json"
    if derep.exists() and qc_file.exists():
        return pd.read_json(qc_file, typ="series").to_dict()
    urls = fastq_urls(row["fastq_ftp"])
    if not urls:
        return {"run": run, "status": "no_fastq"}
    tmp = tmp_root / run
    tmp.mkdir(parents=True, exist_ok=True)
    try:
        raw = [download(u, tmp / f"raw_{i}.fastq.gz") for i, u in enumerate(urls, 1)]
        window, n_in = prepare_reads(raw, mode, ref, amp, lo, hi, tmp)
        n_window = count_fastq_records(window)
        filter_derep(window, tmp, run, derep, amp["max_ee"], lo, hi)
        qc = {"run": run, "status": "ok", "read_prep": mode, "reads_in": n_in,
              "reads_window": n_window, "reads_filtered": fasta_size_total(derep)}
        pd.Series(qc).to_json(qc_file)
        return qc
    except Exception as exc:                                   # noqa: BLE001
        return {"run": run, "status": f"error: {exc}"}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_sheet(study: dict) -> pd.DataFrame:
    """The study's runs from the sample sheet, with FASTQ URLs and the read_prep per run."""
    sheet = pd.read_csv(data_path("processed") / "sample_sheet.csv")
    sheet = sheet[sheet["study_id"] == study["id"]].copy()
    acc = study["accession"]
    urls = ena_table(acc, "read_run", "run_accession,fastq_ftp",
                     data_path("raw") / "census" / "ena" / f"{acc}_fastq.tsv")
    sheet = sheet.merge(urls, on="run_accession", how="left")
    sheet["read_prep"] = study["read_prep"]
    if study.get("read_prep_plate"):
        sheet.loc[sheet["kind"] == "plate", "read_prep"] = study["read_prep_plate"]
    return sheet


def process_study(study: dict, per_run_only: bool) -> None:
    """All runs of one study, then pooled denoising and mapping."""
    amp = load_config("thresholds")["amplicon"]
    refcfg = load_config("thresholds")["reference"]
    sets = pd.read_csv(data_path("processed") / "primer_sets.csv")
    ref = next(r for r in sets.to_dict("records") if study["id"] in r["studies"].split(";"))
    ref["fasta"] = (data_path("interim") / "reference" / "regions" /
                    f"cultured_{ref['ref_name']}.fasta")
    lo, hi = refcfg["region_length"][ref["region"]]
    sheet = run_sheet(study)
    out_dir = data_path("interim") / "phase2" / study["id"]
    tmp_root = out_dir / "tmp"
    for sub in ("derep", "qc"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)
    print(f"{study['id']}: {len(sheet)} runs, reference {ref['ref_name']}, window {lo}-{hi} bp",
          flush=True)

    results = []
    with ThreadPoolExecutor(WORKERS) as pool:
        futures = [pool.submit(process_run, r, ref, amp, lo, hi, out_dir, tmp_root)
                   for _, r in sheet.iterrows()]
        for i, fut in enumerate(as_completed(futures), 1):
            results.append(fut.result())
            if i % 25 == 0:
                print(f"  {i}/{len(futures)} runs", flush=True)
    qc = pd.DataFrame(results)
    bad = qc[qc["status"] != "ok"]
    if len(bad):
        print(bad.to_string(index=False))
    if per_run_only:
        return

    ok = qc.loc[qc["status"] == "ok", "run"]
    pooled = out_dir / "all_derep.fasta"
    with pooled.open("w") as handle:
        for run in ok:
            handle.write((out_dir / "derep" / f"{run}.fasta").read_text())
    work = out_dir / "work"
    work.mkdir(exist_ok=True)
    asvs = out_dir / "asv_seqs.fasta"
    denoise(pooled, work, asvs, amp["unoise_minsize"], amp["unoise_alpha"], f"{study['id']}_")
    map_reads(pooled, asvs, work / "otutab.tsv", amp["map_id"])
    counts = otutab_to_long(work / "otutab.tsv")
    counts.to_csv(out_dir / "counts.csv", index=False)
    mapped = counts.groupby("run")["count"].sum().rename("reads_mapped")
    qc = qc.merge(mapped, left_on="run", right_index=True, how="left")
    qc.to_csv(out_dir / "read_qc.csv", index=False)
    pooled.unlink()
    shutil.rmtree(work)
    print(f"{study['id']}: ASVs {counts['asv_id'].nunique()}; runs {counts['run'].nunique()}")
    print(qc[["reads_in", "reads_window", "reads_filtered", "reads_mapped"]].describe()
          .round(0).to_string())


def main() -> None:
    ids = [a for a in sys.argv[1:] if not a.startswith("--")]
    studies = {s["id"]: s for s in load_config("census_studies")["studies"]}
    for sid in ids:
        if studies[sid].get("read_prep") in (None, "unresolved"):
            print(f"{sid}: read_prep unresolved, skipped")
            continue
        process_study(studies[sid], "--per-run-only" in sys.argv)


if __name__ == "__main__":
    main()
