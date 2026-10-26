"""Phase 2, step 3: raw reads -> ASV table, one study at a time, one pipeline for all.

Per run (parallel, cached): download the run's FASTQ files from ENA, prepare the reads
according to the study's read_prep (D-024), filter by expected errors and the region's length
window, dereplicate. Raw FASTQs are deleted after each run.
  primers_in_reads   merge pairs, then keep the window between the primers on either strand
  primers_removed    merge pairs, then orient against the study's cultured reference
  pacbio             single reads, cut to the window of the study's reference primer pair (V4)
  merged_primers_removed  deposited reads already merged and primer-trimmed (one file per run):
                     orient against the study's cultured reference (D-028)
Per study (pooled over its runs): UNOISE3 + uchime3 de novo chimera removal, then every run's
reads are mapped back to the ASVs.

Usage: python scripts/16_reads_to_asvs.py STUDY_ID [STUDY_ID ...] [--per-run-only]
       [--tag=NAME --max-ee=X]   (a reprocessing with another expected-error filter, written to
                                  data/interim/phase2/{study}__{NAME}/ for a sensitivity analysis)
Outputs (data/interim/phase2/{study}/): derep/, qc/, asv_seqs.fasta, counts.csv, read_qc.csv
"""

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from plate_truth.amplicon import denoise_and_map, process_run
from plate_truth.census import ena_table
from plate_truth.config import cli_option, data_path, load_config, study_dir

WORKERS = 4


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


def process_study(study: dict, per_run_only: bool, tag: str | None = None,
                  max_ee: float | None = None) -> None:
    """All runs of one study, then pooled denoising and mapping. tag and max_ee: a
    reprocessing with another expected-error filter, in its own folder."""
    amp = load_config("thresholds")["amplicon"]
    if max_ee is not None:
        amp = {**amp, "max_ee": max_ee}
    refcfg = load_config("thresholds")["reference"]
    sets = pd.read_csv(data_path("processed") / "primer_sets.csv")
    ref = next(r for r in sets.to_dict("records") if study["id"] in r["studies"].split(";"))
    ref["fasta"] = (data_path("interim") / "reference" / "regions" /
                    f"cultured_{ref['ref_name']}.fasta")
    lo, hi = refcfg["region_length"][ref["region"]]
    sheet = run_sheet(study)
    out_dir = study_dir(study["id"], tag)
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
    counts = denoise_and_map(out_dir, list(ok), amp, f"{study['id']}_")
    mapped = counts.groupby("run")["count"].sum().rename("reads_mapped")
    qc = qc.merge(mapped, left_on="run", right_index=True, how="left")
    qc.to_csv(out_dir / "read_qc.csv", index=False)
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
        max_ee = cli_option(sys.argv, "max-ee")
        process_study(studies[sid], "--per-run-only" in sys.argv, cli_option(sys.argv, "tag"),
                      float(max_ee) if max_ee else None)


if __name__ == "__main__":
    main()
