"""Phase 5, step 2: amplicon surveys, raw reads to ASVs, identity and taxonomy
(ANALYSIS_PLAN.md section 21).

The training pipeline unchanged (plate_truth.amplicon.process_run and denoise_and_map, as in
script 16): per run, merge pairs, keep the window between the study's primer pair on either
strand, expected-error and length filter, dereplicate; per study, pooled UNOISE3 with de novo
chimera removal and mapping at 97%. Identity as in script 17 (reference 'all': id1 primary,
id2 and gap-blind mid for sensitivity analyses); taxonomy as in script 18 (SINTAX, SILVA 138.2,
cutoff 0.8).

Input: data/processed/survey_runs_amplicon.csv (script 39).
Outputs per study in data/interim/phase5/amplicon/<study>/: counts.csv, read_qc.csv,
asv_seqs.fasta, nearest_relative.csv, taxonomy.csv. Resumable (each step is skipped when its
output exists). The exhaustive identity searches are saved in chunks of CHUNK ASVs, so a
stopped run loses at most one chunk. Several copies can run at once (one with --reverse works
from the end of the list), so downloads overlap with searches: each study is claimed with a
claim.pid file, and a study claimed by another live copy is skipped. Each copy writes its
process id to data/interim/phase5/run_40_<pid>.pid while it runs (stop: taskkill /F /T /PID
<pid>; restart with the same command).

Usage: python scripts/40_survey_amplicon.py [--reverse] [STUDY ...]   (default: every study)
"""

import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import psutil

from plate_truth.amplicon import denoise_and_map, process_run
from plate_truth.config import data_path, load_config
from plate_truth.identity import nearest_relative_table, read_hits, search, search_chunked
from plate_truth.reference import read_fasta
from plate_truth.taxonomy import parse_sintax, run_sintax_chunks

WORKERS = 4
MIN_ID = 0.75
MAXACCEPTS, MAXREJECTS = 64, 256
SINTAX_CUTOFF = 0.8
N_CHUNKS = 8
CHUNK = 500                      # ASVs per saved piece of an exhaustive search


def study_folder(study: str):
    return data_path("interim") / "phase5" / "amplicon" / study


def reads_to_asvs(study: str, runs: pd.DataFrame, ref: dict, amp: dict, lo: int, hi: int):
    out = study_folder(study)
    if (out / "read_qc.csv").exists():
        return
    for sub in ("derep", "qc"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    runs = runs.assign(read_prep="primers_in_reads")
    results = []
    with ThreadPoolExecutor(WORKERS) as pool:
        futs = [pool.submit(process_run, r, ref, amp, lo, hi, out, out / "tmp")
                for _, r in runs.iterrows()]
        for fut in as_completed(futs):
            results.append(fut.result())
    qc = pd.DataFrame(results)
    ok = qc.loc[qc["status"] == "ok", "run"].tolist()
    ok = [r for r in ok if (out / "derep" / f"{r}.fasta").stat().st_size > 0]
    if not ok:
        qc.to_csv(out / "read_qc.csv", index=False)
        print(f"{study}: no run produced window reads", flush=True)
        return
    counts = denoise_and_map(out, ok, amp, f"{study}_")
    mapped = counts.groupby("run")["count"].sum().rename("reads_mapped")
    qc = qc.merge(mapped, left_on="run", right_index=True, how="left")
    qc.to_csv(out / "read_qc.csv", index=False)
    print(f"{study}: {counts['asv_id'].nunique()} ASVs from {len(ok)} runs", flush=True)


def identities(study: str, ref_fasta) -> None:
    out = study_folder(study)
    asvs = out / "asv_seqs.fasta"
    if (out / "nearest_relative.csv").exists() or not asvs.exists():
        return
    work = out / "nr"
    work.mkdir(exist_ok=True)
    parts = []
    for iddef in (1, 2):
        path = work / f"hits_exhaustive_id{iddef}.tsv"
        if not path.exists():
            search_chunked(asvs, ref_fasta, path, iddef, MIN_ID, 0, 0, CHUNK,
                           top_hits_only=True)
        parts.append(read_hits(path))
    path = work / "hits_bounded.tsv"
    if not path.exists():
        search(asvs, ref_fasta, path, 1, MIN_ID, MAXACCEPTS, MAXREJECTS)
    parts.append(read_hits(path))
    ids = [h.split(";")[0] for h, _ in read_fasta(asvs)]
    t = nearest_relative_table(ids, pd.concat(parts).drop_duplicates(["query", "target"]))
    t.to_csv(out / "nearest_relative.csv")
    print(f"{study}: median id1 {t['identity_id1'].median():.2f}; "
          f"no hit {t['identity_id1'].isna().mean():.1%}", flush=True)


def taxonomy(study: str, seed: int) -> None:
    out = study_folder(study)
    asvs = out / "asv_seqs.fasta"
    if (out / "taxonomy.csv").exists() or not asvs.exists():
        return
    tsv = out / "sintax.tsv"
    if not tsv.exists():
        run_sintax_chunks(asvs, data_path("interim") / "taxonomy" / "silva138.2_sintax.fasta",
                          tsv, SINTAX_CUTOFF, seed, N_CHUNKS, out)
    parse_sintax(tsv).to_csv(out / "taxonomy.csv")


def claim(study: str) -> bool:
    """Claim a study for this process; False if another live process holds it."""
    lock = study_folder(study) / "claim.pid"
    if lock.exists():
        other = int(lock.read_text() or 0)
        if other != os.getpid() and psutil.pid_exists(other):
            return False
    lock.parent.mkdir(parents=True, exist_ok=True)
    lock.write_text(str(os.getpid()))
    return True


def main(args: list[str]) -> None:
    pid = data_path("interim") / "phase5" / f"run_40_{os.getpid()}.pid"
    pid.write_text(str(os.getpid()))
    try:
        run([a for a in args if a != "--reverse"], "--reverse" in args)
    finally:
        pid.unlink(missing_ok=True)


def run(studies: list[str], reverse: bool = False) -> None:
    th = load_config("thresholds")
    amp, refcfg = th["amplicon"], th["reference"]
    seed = th["resampling"]["seed"]
    runs = pd.read_csv(data_path("processed") / "survey_runs_amplicon.csv", dtype=str)
    sets = pd.read_csv(data_path("processed") / "primer_sets.csv").set_index("ref_name")
    regions = data_path("interim") / "reference" / "regions"
    groups = list(runs.groupby("study", sort=False))
    for study, r in (groups[::-1] if reverse else groups):
        if studies and study not in studies:
            continue
        if (study_folder(study) / "taxonomy.csv").exists() or not claim(study):
            continue
        name = r["reference"].iloc[0]
        ref = {"ref_name": name, "primer_fwd": sets.loc[name, "primer_fwd"],
               "primer_rev": sets.loc[name, "primer_rev"],
               "fasta": regions / f"cultured_{name}.fasta"}
        lo, hi = refcfg["region_length"][sets.loc[name, "region"]]
        print(f"{study} ({r['habitat'].iloc[0]}): {len(r)} runs, {name}, window {lo}-{hi} bp",
              flush=True)
        try:
            reads_to_asvs(study, r, ref, amp, lo, hi)
            identities(study, ref["fasta"])
            taxonomy(study, seed)
        except Exception as exc:                                   # noqa: BLE001
            print(f"{study}: FAILED {exc}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
