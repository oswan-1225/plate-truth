"""Cultured 16S reference: collect sequences from cultured sources, cut out the amplified
region, dereplicate, and keep per-sequence provenance.

Sources (all from cultured organisms; MAG/SAG-derived sequences are excluded):
  bacdive_strain  16S accession listed in a strain's own BacDive record (fetched from ENA)
  refseq_type     NCBI RefSeq Targeted Loci 16S (curated, from type material)
  gtdb_isolate    16S found in a GTDB genome that NCBI does not flag as metagenome- or
                  single-cell-derived
"""

import gzip
import re
from collections.abc import Iterator
from pathlib import Path

import pandas as pd

from plate_truth.amplicon import reverse_complement
from plate_truth.tools import cutadapt

AMBIGUOUS = re.compile(r"[^ACGT]")
MAG_LIKE_STRAIN = re.compile(r"MAG|(^|[^A-Za-z])bin[._\s-]?\d", re.IGNORECASE)


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header without '>', sequence) from a FASTA file, gzipped or not."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line.upper())
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(records, path: Path) -> int:
    """Write (id, sequence) pairs to FASTA; return the number written."""
    n = 0
    with open(path, "w") as handle:
        for rid, seq in records:
            handle.write(f">{rid}\n{seq}\n")
            n += 1
    return n


def ambiguous_fraction(seq: str) -> float:
    """Fraction of positions that are not A, C, G or T."""
    return len(AMBIGUOUS.findall(seq.upper())) / len(seq) if seq else 1.0


def sequence_qc(seq: str | None, min_len: int, max_ambig: float) -> str:
    """'ok', 'missing', 'too_short' or 'too_ambiguous' for a full-length 16S sequence."""
    if not seq:
        return "missing"
    if len(seq) < min_len:
        return "too_short"
    if ambiguous_fraction(seq) > max_ambig:
        return "too_ambiguous"
    return "ok"


def is_isolate_genome(category: str, strain_ids: str) -> bool:
    """True for a GTDB genome from an isolate: NCBI category 'none' (not derived from a
    metagenome, single cell or environmental sample) and no MAG/bin-like strain name."""
    if category != "none":
        return False
    return not MAG_LIKE_STRAIN.search(str(strain_ids or ""))


def accession_base(acc: str) -> str:
    """Assembly accession without database prefix, GCA/GCF type or version:
    'RS_GCF_000331735.1' and 'GCA_000331735' both give '000331735'."""
    hit = re.search(r"GC[AF]_(\d{9})", str(acc))
    return hit.group(1) if hit else ""


def in_silico_pcr(fasta_in: Path, fasta_out: Path, fwd: str, rev: str, error_rate: float,
                  min_len: int, max_len: int, report: Path) -> None:
    """Cut each reference down to the region between the primers (primers removed).

    Uses a cutadapt linked adapter FWD...revcomp(REV); both primers must be found, and
    reverse-complemented references are also searched. Sequences without both primer sites
    (partial sequences, or primer mismatches beyond error_rate) are discarded.
    """
    linked = f"{fwd}...{reverse_complement(rev)}"
    cutadapt("-g", linked, "-e", error_rate, "--revcomp", "--discard-untrimmed",
             "-m", min_len, "-M", max_len, "-j", 0, "--json", report,
             "-o", fasta_out, fasta_in)


def dereplicate(records) -> tuple[list[tuple[str, str]], pd.DataFrame]:
    """Collapse identical sequences. Returns (unique records with ids ref_000001...,
    membership table: ref_id, member_id). Order of first appearance is kept."""
    seq_to_ref: dict[str, str] = {}
    members = []
    uniques = []
    for rid, seq in records:
        ref = seq_to_ref.get(seq)
        if ref is None:
            ref = f"ref_{len(seq_to_ref) + 1:06d}"
            seq_to_ref[seq] = ref
            uniques.append((ref, seq))
        members.append((ref, rid))
    return uniques, pd.DataFrame(members, columns=["ref_id", "member_id"])


def as_list(x) -> list:
    """BacDive fields are a dict when single and a list when repeated; always return a list."""
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def parse_bacdive_record(bid: int, rec: dict) -> dict:
    """Pull the fields this project needs from one BacDive v2 record.

    Returns {'strain': {...}, 's16': [...], 'genomes': [...]}. A strain counts as having a
    growth record if any culture-medium or culture-temperature entry does not say 'no'.
    """
    tax = rec.get("Name and taxonomic classification") or {}
    lpsn = tax.get("LPSN") or {}
    if isinstance(lpsn, list):
        lpsn = lpsn[0]
    cg = rec.get("Culture and growth conditions") or {}
    growth_entries = as_list(cg.get("culture medium")) + as_list(cg.get("culture temp"))
    has_growth = any(str(e.get("growth", "")).lower() != "no" for e in growth_entries)
    strain = {"bacdive_id": bid, "type_strain": tax.get("type strain"),
              "has_growth_record": has_growth}
    for rank in ["domain", "phylum", "genus", "species"]:
        strain[rank] = lpsn.get(rank) or tax.get(rank)
    seq = rec.get("Sequence information") or {}
    s16 = [{"bacdive_id": bid, "accession": s.get("accession"), "length": s.get("length")}
           for s in as_list(seq.get("16S sequences")) if s.get("accession")]
    genomes = [{"bacdive_id": bid, "assembly": g.get("INSDC accession")}
               for g in as_list(seq.get("Genome sequences")) if g.get("INSDC accession")]
    return {"strain": strain, "s16": s16, "genomes": genomes}
