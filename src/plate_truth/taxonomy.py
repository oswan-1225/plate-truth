"""ASV taxonomy (VSEARCH SINTAX against SILVA) and 16S copy-number lookup (rrnDB)."""

import re
from pathlib import Path

import pandas as pd

from plate_truth.reference import read_fasta
from plate_truth.tools import vsearch

RANKS = ["domain", "phylum", "class", "order", "family", "genus"]
SINTAX_CODES = dict(zip(RANKS, "dpcofg", strict=False))


def dada2_to_sintax_header(i: int, header: str) -> str:
    """'Bacteria;Bacillota;Bacilli;;' -> 'r1;tax=d:Bacteria,p:Bacillota,c:Bacilli;'.

    DADA2 training-set headers list ranks from domain down, ';'-separated; empty or missing
    trailing ranks are left out. Commas inside names would break SINTAX and are replaced."""
    names = [n.strip().replace(",", "_") for n in header.split(";")]
    parts = [f"{SINTAX_CODES[r]}:{n}" for r, n in zip(RANKS, names, strict=False) if n]
    return f"r{i};tax={','.join(parts)};"


def write_sintax_db(trainset: Path, out: Path) -> int:
    """Convert a DADA2-formatted SILVA training FASTA to a SINTAX database."""
    n = 0
    with open(out, "w") as handle:
        for i, (header, seq) in enumerate(read_fasta(trainset), 1):
            handle.write(f">{dada2_to_sintax_header(i, header)}\n{seq}\n")
            n += 1
    return n


def run_sintax(query: Path, db: Path, out: Path, cutoff: float, seed: int) -> None:
    """Classify queries with SINTAX (bootstrap cutoff applies to the reported ranks)."""
    vsearch("--sintax", query, "--db", db, "--tabbedout", out, "--sintax_cutoff", cutoff,
            "--strand", "plus", "--randseed", seed, "--threads", 1)


def parse_sintax(path: Path) -> pd.DataFrame:
    """SINTAX tabbed output -> one row per ASV with a column per rank, keeping only ranks
    that pass the cutoff (column 4). Unclassified ranks are NaN."""
    rows = []
    for line in Path(path).read_text().splitlines():
        cols = line.split("\t")
        row = {"asv_id": cols[0]}
        passed = cols[3] if len(cols) > 3 else ""
        for item in filter(None, passed.split(",")):
            code, name = item.split(":", 1)
            rank = next(r for r, c in SINTAX_CODES.items() if c == code)
            row[rank] = name
        rows.append(row)
    return pd.DataFrame(rows, columns=["asv_id", *RANKS]).set_index("asv_id")


def copy_number_lookup(pantaxa: pd.DataFrame) -> dict[tuple[str, str], float]:
    """(rank, name) -> mean 16S copy number from an rrnDB pan-taxa table."""
    t = pantaxa.dropna(subset=["name", "mean"])
    rows = zip(t["rank"], t["name"], t["mean"], strict=True)
    return {(r.lower(), n): float(m) for r, n, m in rows}


def clean_taxon(name: str) -> str:
    """Drop SILVA-style suffixes so names match NCBI: 'Bacillus_A' -> 'Bacillus',
    'Burkholderia-Caballeronia-Paraburkholderia' -> 'Burkholderia'."""
    return re.split(r"[_\-\s]", str(name))[0]


def asv_copy_number(tax: pd.DataFrame, lookup: dict, default: float) -> pd.DataFrame:
    """Copy number per ASV from its lowest classified rank that rrnDB knows.

    Returns columns copy_number and copy_number_rank ('default' when nothing matched).
    """
    cn, src = [], []
    for _, row in tax.iterrows():
        value, rank_used = default, "default"
        for rank in reversed(RANKS):
            name = row.get(rank)
            if isinstance(name, str) and (rank, clean_taxon(name)) in lookup:
                value, rank_used = lookup[(rank, clean_taxon(name))], rank
                break
        cn.append(value)
        src.append(rank_used)
    return pd.DataFrame({"copy_number": cn, "copy_number_rank": src}, index=tax.index)


def run_sintax_chunks(query: Path, db: Path, out: Path, cutoff: float, seed: int,
                      n_chunks: int, workdir: Path) -> None:
    """SINTAX in parallel: split the query into chunks, classify each with one thread and the
    same seed (so the result is reproducible), then concatenate the outputs."""
    from concurrent.futures import ThreadPoolExecutor

    from plate_truth.reference import write_fasta

    records = list(read_fasta(query))
    paths = []
    for i in range(n_chunks):
        chunk = workdir / f"sintax_chunk_{i}.fasta"
        write_fasta(records[i::n_chunks], chunk)
        paths.append(chunk)
    with ThreadPoolExecutor(n_chunks) as pool:
        list(pool.map(lambda p: run_sintax(p, db, p.with_suffix(".tsv"), cutoff, seed), paths))
    with open(out, "w") as handle:
        for p in paths:
            handle.write(p.with_suffix(".tsv").read_text())
            p.unlink()
            p.with_suffix(".tsv").unlink()


def base_phylum(name) -> str | None:
    """Phylum name without GTDB prefixes/suffixes: 'p__Bacillota_A' -> 'Bacillota'."""
    if not isinstance(name, str) or not name:
        return None
    return re.sub(r"_[A-Z]$", "", name.removeprefix("p__"))


def phylum_consistent(genome_phylum, ssu_phylum) -> bool:
    """False only when both phyla are known and differ; an unclassified 16S cannot be judged
    and is treated as consistent."""
    g, s = base_phylum(genome_phylum), base_phylum(ssu_phylum)
    return g is None or s is None or g == s
