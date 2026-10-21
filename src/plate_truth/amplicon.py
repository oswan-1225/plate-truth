"""Amplicon reads to ASVs: primer trimming (cutadapt) and VSEARCH merge/filter/denoise/map.

Per run: trim primers -> merge pairs -> quality filter -> dereplicate (cached).
Pooled over runs: dereplicate -> UNOISE3 denoise -> de novo chimera removal -> map reads back.
"""

import json
import re
from pathlib import Path

import pandas as pd

from plate_truth.tools import cutadapt, vsearch

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence, IUPAC codes included."""
    return seq.upper().translate(COMPLEMENT)[::-1]


def paired_fastq_urls(fastq_ftp: str) -> list[str]:
    """The two paired-end FASTQ URLs (_1, _2) from an ENA fastq_ftp field, as https URLs.

    ENA sometimes also lists an unpaired file without the _1/_2 suffix; it is ignored.
    Returns [] if the pair is not present.
    """
    urls = sorted(u for u in str(fastq_ftp).split(";") if re.search(r"_[12]\.fastq\.gz$", u))
    if len(urls) != 2:
        return []
    return ["https://" + u if not u.startswith("http") else u for u in urls]


def trim_primers(r1: Path, r2: Path, out1: Path, out2: Path, fwd: str, rev: str,
                 error_rate: float, report: Path) -> dict:
    """Remove the forward primer from R1 and the reverse primer from R2 with cutadapt.

    Read pairs in which either primer is not found are discarded, because a read without
    its primer is off-target or truncated. Returns the cutadapt JSON report as a dict.
    """
    cutadapt("-g", fwd, "-G", rev, "-e", error_rate, "--overlap", min(len(fwd), len(rev)) - 2,
             "--discard-untrimmed", "--pair-filter=any", "-j", 1,
             "-o", out1, "-p", out2, "--json", report, r1, r2)
    return json.loads(Path(report).read_text())


def merge_filter_derep(r1: Path, r2: Path, workdir: Path, run: str, out: Path,
                       max_ee: float, min_len: int) -> None:
    """Merge pairs, filter by expected errors, and dereplicate one run with VSEARCH.

    Output is a FASTA of unique sequences with ;size= annotations, labelled by run.
    """
    merged = workdir / "merged.fq"
    filt = workdir / "filt.fa"
    vsearch("--fastq_mergepairs", r1, "--reverse", r2, "--fastq_maxdiffs", 10,
            "--fastqout", merged)
    vsearch("--fastq_filter", merged, "--fastq_maxee", max_ee, "--fastq_minlen", min_len,
            "--fastq_maxns", 0, "--fastaout", filt, "--relabel", f"{run}.")
    part = out.with_name(out.name + ".part")
    vsearch("--derep_fulllength", filt, "--sizeout", "--relabel", f"{run}.u", "--output", part)
    part.replace(out)


def count_fastq_records(path: Path) -> int:
    """Number of records in an uncompressed FASTQ file."""
    with open(path, "rb") as handle:
        return sum(1 for _ in handle) // 4


def fasta_size_total(path: Path) -> int:
    """Sum of ;size=N annotations in a dereplicated FASTA (= reads represented)."""
    total = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                total += int(re.search(r";size=(\d+)", line).group(1))
    return total


def denoise(pooled_derep: Path, workdir: Path, asv_out: Path, minsize: int, alpha: float,
            prefix: str) -> None:
    """Pool-dereplicate, UNOISE3-denoise and remove de novo chimeras; write labelled ASVs."""
    pooled = workdir / "pooled.fa"
    denoised = workdir / "denoised.fa"
    nochim = workdir / "nochim.fa"
    vsearch("--derep_fulllength", pooled_derep, "--sizein", "--sizeout", "--output", pooled)
    vsearch("--cluster_unoise", pooled, "--minsize", minsize, "--unoise_alpha", alpha,
            "--centroids", denoised)
    vsearch("--uchime3_denovo", denoised, "--nonchimeras", nochim)
    vsearch("--fastx_filter", nochim, "--relabel", prefix, "--fastaout", asv_out)


def map_reads(derep: Path, asvs: Path, out_tsv: Path, map_id: float, threads: int = 0) -> None:
    """Map dereplicated reads (with sizes) to ASVs; write a VSEARCH OTU table."""
    vsearch("--usearch_global", derep, "--db", asvs, "--id", map_id, "--sizein",
            "--strand", "plus", "--threads", threads, "--otutabout", out_tsv)


def otutab_to_long(path: Path) -> pd.DataFrame:
    """Read a VSEARCH OTU table into long format: run, asv_id, count (zeros dropped)."""
    tab = pd.read_csv(path, sep="\t", index_col=0)
    tab.index.name = "asv_id"
    long = tab.reset_index().melt(id_vars="asv_id", var_name="run", value_name="count")
    return long[long["count"] > 0].reset_index(drop=True)


# --- Phase 2: one read-preparation path for every study ---------------------------------------
# Studies deposit reads differently (D-024): primers still in the reads in either orientation,
# primers already removed, or single PacBio reads. Every path ends in the same object: reads
# covering exactly the inter-primer window, in the forward orientation, then quality-filtered
# and dereplicated per run.

def fastq_urls(fastq_ftp: str) -> list[str]:
    """The run's read files as https URLs: the _1/_2 pair if present, else the single file."""
    pair = paired_fastq_urls(fastq_ftp)
    if pair:
        return pair
    files = [u for u in str(fastq_ftp).split(";") if u.endswith(".fastq.gz")]
    return ["https://" + u if not u.startswith("http") else u for u in files[:1]]


def merge_pairs(r1: Path, r2: Path, out: Path) -> None:
    """Merge read pairs; staggered pairs keep only the overlap (adapter read-through removed)."""
    vsearch("--fastq_mergepairs", r1, "--reverse", r2, "--fastq_maxdiffs", 10,
            "--fastq_allowmergestagger", "--fastqout", out)


def cut_window(reads: Path, out: Path, fwd: str, rev: str, error_rate: float, min_len: int,
               max_len: int, report: Path) -> dict:
    """Keep the part of each read between the two primers, in forward orientation.

    cutadapt linked adapter FWD...revcomp(REV), searched on both strands (--revcomp), so reads
    in either orientation, with barcodes or spacers before the primer, give the same window.
    Reads without both primer sites are discarded. Returns the cutadapt JSON report.
    """
    cutadapt("-g", f"{fwd}...{reverse_complement(rev)}", "-e", error_rate, "--revcomp",
             "--discard-untrimmed", "-m", min_len, "-M", max_len, "-j", 1,
             "--json", report, "-o", out, reads)
    return json.loads(Path(report).read_text())


def orient_reads(reads: Path, db: Path, out: Path) -> None:
    """Put primer-free reads in forward orientation by k-mer comparison with a reference;
    reads whose orientation cannot be decided are dropped."""
    vsearch("--orient", reads, "--db", db, "--fastqout", out)


def filter_derep(reads: Path, workdir: Path, run: str, out: Path, max_ee: float, min_len: int,
                 max_len: int) -> None:
    """Expected-error and length filter, then dereplicate one run (FASTA with ;size=)."""
    filt = workdir / "filt.fa"
    vsearch("--fastq_filter", reads, "--fastq_maxee", max_ee, "--fastq_minlen", min_len,
            "--fastq_maxlen", max_len, "--fastq_maxns", 0, "--fastaout", filt,
            "--relabel", f"{run}.")
    part = out.with_name(out.name + ".part")
    vsearch("--derep_fulllength", filt, "--sizeout", "--relabel", f"{run}.u", "--output", part)
    part.replace(out)
