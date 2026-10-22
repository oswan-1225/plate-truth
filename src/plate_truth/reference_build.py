"""Building the cultured reference: load the cultured sources, group studies by primer pair,
and flag genome assemblies whose 16S copies disagree with the genome's own phylum.

The loaders were moved here from scripts/03_build_reference.py unchanged, so that the Phase 0
build (script 03) and the per-region Phase 2 build (script 15) share one implementation.
"""

import json
import re
from pathlib import Path

import pandas as pd

from plate_truth.reference import is_isolate_genome, parse_bacdive_record, read_fasta, sequence_qc
from plate_truth.taxonomy import phylum_consistent


def load_bacdive(raw: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Parse every cached BacDive batch into strain, 16S-accession and genome tables."""
    strains, s16, genomes = [], [], []
    for path in sorted((raw / "bacdive").rglob("batch_*.json")):
        results = json.loads(path.read_text(encoding="utf-8")).get("results") or {}
        for bid, rec in results.items():
            parsed = parse_bacdive_record(int(bid), rec)
            strains.append(parsed["strain"])
            s16 += parsed["s16"]
            genomes += parsed["genomes"]
    st = pd.DataFrame(strains).drop_duplicates("bacdive_id")
    return st, pd.DataFrame(s16).drop_duplicates(), pd.DataFrame(genomes).drop_duplicates()


def load_ena_cache(raw: Path) -> dict[str, str]:
    """Accession (no version) -> sequence, from the cached ENA FASTA batches."""
    seqs = {}
    for path in sorted((raw / "ena_16s").glob("batch_*.fasta")):
        for header, seq in read_fasta(path):
            acc = header.split("|")[1] if "|" in header else header.split()[0]
            seqs[acc] = seq
    return seqs


def choose_strain_16s(s16: pd.DataFrame, ena: dict, cfg: dict) -> pd.DataFrame:
    """One 16S per strain: the longest listed accession (ties: first listed), then QC."""
    s = s16.copy()
    s["length"] = pd.to_numeric(s["length"], errors="coerce").fillna(0)
    best = s.sort_values(["bacdive_id", "length"], ascending=[True, False],
                         kind="stable").drop_duplicates("bacdive_id")
    best["seq"] = best["accession"].map(ena)
    best["status"] = [sequence_qc(x if isinstance(x, str) else None, cfg["full_min_len"],
                                  cfg["full_max_ambig"]) for x in best["seq"]]
    return best


def gtdb_isolates(raw: Path) -> pd.DataFrame:
    """GTDB genomes (bacteria + archaea) with isolate status and type-material flag."""
    cols = ["accession", "ncbi_genome_category", "ncbi_strain_identifiers",
            "ncbi_type_material_designation", "ncbi_genbank_assembly_accession", "gtdb_taxonomy"]
    meta = pd.concat([pd.read_csv(raw / "reference" / f"gtdb_r232_{d}_metadata.tsv.gz", sep="\t",
                                  usecols=cols, low_memory=False) for d in ("bac120", "ar53")])
    meta["is_isolate"] = [is_isolate_genome(c, s) for c, s in
                          zip(meta["ncbi_genome_category"], meta["ncbi_strain_identifiers"],
                              strict=True)]
    return meta


def full_length_records(best16, meta, raw: Path):
    """Yield (source|id, sequence) for every full-length cultured sequence."""
    for bid, seq in best16.loc[best16["status"] == "ok", ["bacdive_id", "seq"]].itertuples(
            index=False):
        yield f"bacdive_strain|{bid}", seq
    for dom in ("bacteria", "archaea"):
        for header, seq in read_fasta(raw / "reference" / f"{dom}.16SrRNA.fna.gz"):
            yield f"refseq_type|{header.split()[0]}", seq
    isolates = set(meta.loc[meta["is_isolate"], "accession"])
    for i, (header, seq) in enumerate(read_fasta(raw / "reference" / "gtdb_r232_ssu_all.fna.gz")):
        genome = header.split("~")[0]
        if genome in isolates:
            yield f"gtdb_isolate|{genome}|{i}", seq


def primer_sets(studies: list[dict], full_length_set: dict) -> pd.DataFrame:
    """One row per distinct (forward, reverse) primer pair among the given studies.

    Full-length studies (region 'full-length') are cut in silico to `full_length_set` (a dict
    with region, primers, primer_fwd, primer_rev) for the pooled analysis (D-017), so they share
    that set. A study whose deposited reads are already cut to the window of another study's
    primer pair names that study in `window_of` and shares its set (D-028). Returns columns:
    ref_name, region, primers, primer_fwd, primer_rev, studies (';'-joined ids). ref_name is
    '{region}_{primer names}', with a numeric suffix if two different sequence pairs share a name.
    """
    by_id = {s["id"]: s for s in studies}
    rows = {}
    for s in studies:
        if s["region"] == "full-length":
            src = full_length_set
        elif s.get("window_of"):
            src = by_id[s["window_of"]]
        else:
            src = s
        key = (src["primer_fwd"], src["primer_rev"])
        if key not in rows:
            rows[key] = {"region": src["region"], "primers": src["primers"],
                         "primer_fwd": key[0], "primer_rev": key[1], "studies": []}
        rows[key]["studies"].append(s["id"])
    out = pd.DataFrame(rows.values())
    base = (out["region"] + "_" + out["primers"]).map(lambda x: re.sub(r"[^A-Za-z0-9-]+", "_", x))
    out["ref_name"] = base + base.groupby(base).cumcount().map(lambda i: f"_{i + 1}" if i else "")
    out["studies"] = out["studies"].map(";".join)
    return out[["ref_name", "region", "primers", "primer_fwd", "primer_rev", "studies"]]


def flag_contaminated_genomes(members: pd.DataFrame, ssu_phylum: pd.Series,
                              genome_phylum: pd.Series) -> pd.DataFrame:
    """Per GTDB genome: how many of its 16S regions disagree with its own phylum.

    members: ref_id, member_id ('gtdb_isolate|{genome}|{i}' for genome-derived sequences),
    source. ssu_phylum: ref_id -> SINTAX phylum. genome_phylum: genome -> GTDB phylum.
    A genome is dropped if any of its regions is classified to a different phylum
    (phylum_consistent). Returns columns genome, n_ssu, n_inconsistent, dropped.
    """
    g = members[members["source"] == "gtdb_isolate"].copy()
    g["genome"] = g["member_id"].str.split("|").str[1]
    consistent = [phylum_consistent(gp, sp) for gp, sp in
                  zip(g["genome"].map(genome_phylum), g["ref_id"].map(ssu_phylum), strict=True)]
    g["inconsistent"] = [not c for c in consistent]
    out = g.groupby("genome").agg(n_ssu=("ref_id", "size"),
                                  n_inconsistent=("inconsistent", "sum")).reset_index()
    out["dropped"] = out["n_inconsistent"] > 0
    return out
