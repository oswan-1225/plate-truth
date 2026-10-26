"""Nearest cultured relative: global alignment identity between ASVs and the reference.

VSEARCH identity definitions (manual, vsearch-userfields(7)):
  id1  matching columns / alignment length, ALL gaps counted (internal and terminal).
       Primary: the reference is cut to the same amplified region as the ASV first, so
       terminal gaps only appear where one sequence is genuinely longer.
  id2  as id1 but terminal gaps excluded (VSEARCH default).
  mid  matching columns / columns without gaps: ignores every gap (the gap-blind identity
       criticised by Steen et al. 2019).
"""

import shutil
from pathlib import Path

import pandas as pd

from plate_truth.reference import read_fasta, write_fasta
from plate_truth.tools import vsearch

FIELDS = ["query", "target", "id1", "id2", "mid", "alnlen", "ids", "mism", "gaps", "ql", "tl",
          "qcov"]


def search(query: Path, db: Path, out: Path, iddef: int, min_id: float, maxaccepts: int,
           maxrejects: int, threads: int = 0, top_hits_only: bool = False) -> None:
    """Global-align every query against the database; write tab-separated hits (FIELDS).

    maxaccepts/maxrejects = 0 makes the search exhaustive; otherwise candidate targets are
    tried in decreasing k-mer similarity and the search stops early. top_hits_only keeps only
    the hits tied at the best identity under iddef (needed to keep exhaustive output small).
    """
    extra = ["--top_hits_only"] if top_hits_only else []
    part = Path(out).with_name(Path(out).name + ".part")   # an interrupted search never looks done
    vsearch("--usearch_global", query, "--db", db, "--iddef", iddef, "--id", min_id,
            "--maxaccepts", maxaccepts, "--maxrejects", maxrejects, "--strand", "plus",
            "--threads", threads, "--userout", part, "--userfields", "+".join(FIELDS), *extra)
    part.replace(out)


def chunk_records(records, size: int) -> list[list]:
    """Split (id, sequence) pairs into consecutive lists of at most size records."""
    records = list(records)
    return [records[i:i + size] for i in range(0, len(records), size)]


def search_chunked(query: Path, db: Path, out: Path, iddef: int, min_id: float,
                   maxaccepts: int, maxrejects: int, chunk_size: int, threads: int = 0,
                   top_hits_only: bool = False) -> None:
    """search() run on chunks of chunk_size queries, each saved when it finishes, so a stopped
    search resumes at the first unfinished chunk. Each query is aligned independently, so the
    joined output equals one search over all queries."""
    out = Path(out)
    work = out.with_name(f"{out.stem}_chunks{chunk_size}")
    work.mkdir(exist_ok=True)
    done = []
    for i, recs in enumerate(chunk_records(read_fasta(query), chunk_size)):
        hits = work / f"hits_{i:04d}.tsv"
        if not hits.exists():
            q = work / f"query_{i:04d}.fasta"
            write_fasta(recs, q)
            search(q, db, hits, iddef, min_id, maxaccepts, maxrejects, threads, top_hits_only)
            q.unlink()
        done.append(hits)
    part = out.with_name(out.name + ".part")
    with open(part, "w") as handle:
        for hits in done:
            handle.write(hits.read_text())
    part.replace(out)
    shutil.rmtree(work)


def read_hits(path: Path) -> pd.DataFrame:
    """Read a hits file written by search()."""
    return pd.read_csv(path, sep="\t", names=FIELDS, header=None)


def best_hits(hits: pd.DataFrame, id_col: str) -> pd.DataFrame:
    """Per query: the best identity under id_col, the number of targets tied at that
    identity, and their ids (';'-joined). Ties are kept, not broken arbitrarily."""
    top = hits[hits[id_col] == hits.groupby("query")[id_col].transform("max")]
    return top.groupby("query").agg(identity=(id_col, "first"),
                                    n_tied=("target", "nunique"),
                                    targets=("target", lambda t: ";".join(sorted(set(t)))))


def nearest_relative_table(queries: list[str], hits: pd.DataFrame,
                           id_cols=("id1", "id2", "mid")) -> pd.DataFrame:
    """One row per query with the best identity under each definition (NaN = no hit)."""
    out = pd.DataFrame(index=pd.Index(queries, name="asv_id"))
    for col in id_cols:
        b = best_hits(hits, col)
        out[f"identity_{col}"] = b["identity"]
        out[f"n_tied_{col}"] = b["n_tied"]
        out[f"targets_{col}"] = b["targets"]
    return out
