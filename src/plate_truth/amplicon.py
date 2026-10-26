"""Amplicon reads to ASVs: primer trimming (cutadapt) and VSEARCH merge/filter/denoise/map.

Per run: trim primers -> merge pairs -> quality filter -> dereplicate (cached).
Pooled over runs: dereplicate -> UNOISE3 denoise -> de novo chimera removal -> map reads back.
"""

import gzip
import json
import re
import shutil
from pathlib import Path

import pandas as pd

from plate_truth.http import download
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


# --- Forward reads only (sensitivity for runs whose mates are missing, D-030) -----------------

def asv_end_groups(asvs: list[tuple[str, str]], length: int, end: str) -> pd.DataFrame:
    """Group ASVs that share their first (end='start') or last (end='end') `length` bases: a
    single read from that end of the amplicon cannot tell them apart.

    asvs: (asv_id, sequence). Returns asv_id, group_id, part; group ids {s|e}1, {s|e}2, ... in
    order of first appearance.
    """
    if end not in ("start", "end"):
        raise ValueError(f"end must be 'start' or 'end', not {end!r}")
    ids: dict[str, str] = {}
    rows = []
    for asv_id, seq in asvs:
        part = seq[:length] if end == "start" else seq[-length:]
        ids.setdefault(part, f"{end[0]}{len(ids) + 1}")
        rows.append({"asv_id": asv_id, "group_id": ids[part], "part": part})
    return pd.DataFrame(rows)


def reads_by_orientation(reads: Path, fwd: str, rev: str, error_rate: float, out_dir: Path
                         ) -> tuple[Path, Path]:
    """Split single reads by the primer they start with (spacers before it allowed) and remove
    it: (reads read from the forward primer, reads read from the reverse primer). Reads with
    neither primer are dropped."""
    cutadapt("-g", f"fwd={fwd}", "-g", f"rev={rev}", "-e", error_rate, "--discard-untrimmed",
             "-j", 1, "-o", out_dir / "start-{name}.fq", reads)
    return out_dir / "start-fwd.fq", out_dir / "start-rev.fq"


def window_end_reads(reads: Path, length: int, max_ee: float, out: Path, workdir: Path,
                     revcomp: bool) -> None:
    """Truncate primer-free reads to `length`, expected-error filter (no Ns), optionally reverse-
    complement (reads from the reverse primer then show the window's last bases in forward
    orientation), and dereplicate with sizes."""
    filt = workdir / f"{out.stem}_filt.fa"
    vsearch("--fastq_filter", reads, "--fastq_trunclen", length, "--fastq_maxee", max_ee,
            "--fastq_maxns", 0, "--fastaout", filt)
    if revcomp:
        rc = workdir / f"{out.stem}_rc.fa"
        vsearch("--fastx_revcomp", filt, "--fastaout", rc)
        filt = rc
    vsearch("--derep_fulllength", filt, "--sizeout", "--output", out)


def split_group_counts(group_counts: pd.DataFrame, groups: pd.DataFrame,
                       weights: pd.DataFrame) -> pd.DataFrame:
    """Share each (run, group) count among the group's ASVs in proportion to their counts in
    another table for the same run (weights: run, asv_id, count); equally if those are all
    zero. Integer counts by largest remainder, so each group's total is kept exactly.

    group_counts: run, group_id, count. groups: asv_id, group_id. Returns run, asv_id, count.
    """
    g = group_counts.merge(groups[["asv_id", "group_id"]], on="group_id")
    g = g.merge(weights.rename(columns={"count": "w"}), on=["run", "asv_id"], how="left")
    g["w"] = g["w"].fillna(0).astype(float)
    key = ["run", "group_id"]
    wsum = g.groupby(key)["w"].transform("sum")
    size = g.groupby(key)["asv_id"].transform("size")
    share = (g["w"] / wsum).where(wsum > 0, 1.0 / size)
    exact = g["count"] * share
    g["out"] = exact.astype(int)
    g["rem"] = exact - g["out"]
    short = g["count"] - g.groupby(key)["out"].transform("sum")
    g["rank"] = g.groupby(key)["rem"].rank(method="first", ascending=False)
    g["out"] += (g["rank"] <= short).astype(int)
    out = g[g["out"] > 0][["run", "asv_id", "out"]].rename(columns={"out": "count"})
    return out.reset_index(drop=True)


# --- Per run and per study (shared by the training studies and the surveys) ---------------------

def prepare_reads(raw: list, mode: str, ref: dict, amp: dict, lo: int, hi: int, tmp) -> tuple:
    """Turn the downloaded files into window reads (FASTQ). Returns (path, reads_in)."""
    if mode == "pacbio":
        reads, n_in = raw[0], None
    elif mode == "merged_primers_removed":
        if len(raw) != 1:
            raise ValueError(f"{mode} needs one merged file, got {len(raw)}")
        reads = tmp / "merged.fq"
        with gzip.open(raw[0], "rb") as src, open(reads, "wb") as dst:
            shutil.copyfileobj(src, dst)
        n_in = count_fastq_records(reads)
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
    elif mode in ("primers_removed", "merged_primers_removed"):
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


def denoise_and_map(out_dir: Path, runs: list[str], amp: dict, prefix: str) -> pd.DataFrame:
    """Pool the runs' dereplicated reads, denoise (UNOISE3, chimeras removed) into
    out_dir/asv_seqs.fasta and map every run's reads back; returns long counts and writes
    out_dir/counts.csv."""
    pooled = out_dir / "all_derep.fasta"
    with pooled.open("w") as handle:
        for run in runs:
            handle.write((out_dir / "derep" / f"{run}.fasta").read_text())
    work = out_dir / "work"
    work.mkdir(exist_ok=True)
    asvs = out_dir / "asv_seqs.fasta"
    denoise(pooled, work, asvs, amp["unoise_minsize"], amp["unoise_alpha"], prefix)
    map_reads(pooled, asvs, work / "otutab.tsv", amp["map_id"])
    counts = otutab_to_long(work / "otutab.tsv")
    counts.to_csv(out_dir / "counts.csv", index=False)
    pooled.unlink()
    shutil.rmtree(work)
    return counts
