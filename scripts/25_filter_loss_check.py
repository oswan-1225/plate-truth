"""Phase 2 QC: is the expected-error filter loss taxon-biased? (D-031)

verhoeven2023_sludge loses more than half of its merged reads at the maxEE filter (long V1-V3
amplicons). For two environmental and two plate runs: merge and orient as script 16 does, keep
reads inside the region's length window, split them by the filter (maxEE <= threshold or not),
assign the first 150 bases of each read to ASV start-groups, and compare passing with failing
reads within each run, against a random-loss baseline (as in script 22).

Usage: python scripts/25_filter_loss_check.py [STUDY_ID ...]   (default: verhoeven2023_sludge)
Output: results/phase2/filter_loss_check.csv (per run)
"""

import shutil
import sys

import numpy as np
import pandas as pd

from plate_truth.amplicon import asv_end_groups, merge_pairs, orient_reads, paired_fastq_urls
from plate_truth.config import data_path, load_config, project_root
from plate_truth.http import download
from plate_truth.qc import abundance_distortion
from plate_truth.reference import read_fasta, write_fasta
from plate_truth.tools import vsearch

PART = 150
MIN_READS = 50


def expected_errors(qual: str) -> float:
    """Sum of per-base error probabilities from a Phred+33 quality string."""
    return sum(10 ** (-(ord(c) - 33) / 10) for c in qual)


def pick_runs(sheet: pd.DataFrame, n: int = 2) -> pd.DataFrame:
    """n environment and n plate runs closest to the study's median ENA read count."""
    return pd.concat([sheet[sheet["kind"] == k].assign(
        d=lambda x: (x["read_count"] - x["read_count"].median()).abs()).nsmallest(n, "d")
        for k in ("environment", "plate")])


def assign(parts: list[str], groups_fa, tmp, tag: str) -> pd.Series:
    """Number of parts per ASV start-group (>= 97% identity)."""
    fa, hits = tmp / f"{tag}.fa", tmp / f"{tag}.tsv"
    write_fasta(((f"{tag}{i}", s) for i, s in enumerate(parts)), fa)
    vsearch("--usearch_global", fa, "--db", groups_fa, "--id", 0.97, "--strand", "plus",
            "--maxaccepts", 1, "--maxrejects", 32, "--userfields", "query+target",
            "--userout", hits)
    return pd.read_csv(hits, sep="\t", names=["query", "group_id"])["group_id"].value_counts()


def within(d: pd.DataFrame) -> float:
    """Share of reads in groups whose distortion is within 0.8-1.25."""
    w = d["total"] / d["total"].sum()
    return float(w[(d["distortion"] >= 0.8) & (d["distortion"] <= 1.25)].sum())


def main() -> None:
    th = load_config("thresholds")
    rng = np.random.default_rng(th["resampling"]["seed"])
    max_ee = th["amplicon"]["max_ee"]
    studies = {s["id"]: s for s in load_config("census_studies")["studies"]}
    sets = pd.read_csv(data_path("processed") / "primer_sets.csv")
    sheet = pd.read_csv(data_path("processed") / "sample_sheet.csv")
    rows = []
    for sid in sys.argv[1:] or ["verhoeven2023_sludge"]:
        st = studies[sid]
        lo, hi = th["reference"]["region_length"][st["region"]]
        ref = next(r.ref_name for r in sets.itertuples() if sid in r.studies.split(";"))
        ref_fa = data_path("interim") / "reference" / "regions" / f"cultured_{ref}.fasta"
        work = data_path("interim") / "phase2" / "filter_loss"
        work.mkdir(parents=True, exist_ok=True)
        asvs = [(h.split(";")[0], s) for h, s in read_fasta(
            data_path("interim") / "phase2" / sid / "asv_seqs.fasta")]
        g = asv_end_groups(asvs, PART, "start").drop_duplicates("group_id")
        groups_fa = work / f"{sid}_groups.fasta"
        write_fasta(g[["group_id", "part"]].itertuples(index=False, name=None), groups_fa)
        urls = pd.read_csv(data_path("raw") / "census" / "ena" / f"{st['accession']}_fastq.tsv",
                           sep="\t")
        runs = pick_runs(sheet[sheet["study_id"] == sid]).merge(urls, on="run_accession")
        for row in runs.to_dict("records"):
            tmp = work / row["run_accession"]
            tmp.mkdir(exist_ok=True)
            try:
                r1, r2 = [download(u, tmp / u.rsplit("/", 1)[1])
                          for u in paired_fastq_urls(row["fastq_ftp"])]
                merge_pairs(r1, r2, tmp / "merged.fq")
                orient_reads(tmp / "merged.fq", ref_fa, tmp / "oriented.fq")
                parts = {"pass": [], "fail": []}
                lines = (tmp / "oriented.fq").read_text().split("\n")
                for i in range(0, len(lines) - 3, 4):
                    seq, qual = lines[i + 1], lines[i + 3]
                    if lo <= len(seq) <= hi:
                        ok = expected_errors(qual) <= max_ee
                        parts["pass" if ok else "fail"].append(seq[:PART])
                c = pd.DataFrame({k: assign(v, groups_fa, tmp, k) for k, v in parts.items()})
                c = c.fillna(0)
                d = abundance_distortion(c["pass"], c["fail"], MIN_READS)
                total = c.sum(axis=1).astype(int)
                k = rng.binomial(total, c["pass"].sum() / total.sum())
                dn = abundance_distortion(pd.Series(k, c.index), total - k, MIN_READS)
                w = d["total"] / d["total"].sum()
                rows.append({"study_id": sid, "run": row["run_accession"], "kind": row["kind"],
                             "reads_pass": len(parts["pass"]), "reads_fail": len(parts["fail"]),
                             "within_1.25": within(d), "random_loss_within_1.25": within(dn),
                             "beyond_2": float(w[(d["distortion"] < 0.5) |
                                                 (d["distortion"] > 2)].sum())})
                print(pd.Series(rows[-1]).to_string(), flush=True)
            finally:
                shutil.rmtree(tmp, ignore_errors=True)
    pd.DataFrame(rows).round(4).to_csv(project_root() / "results" / "phase2" /
                                       "filter_loss_check.csv", index=False)


if __name__ == "__main__":
    main()
