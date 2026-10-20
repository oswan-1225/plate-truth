"""Build the cultured 16S reference, cut to the study's amplified region (V4, 515F/806R).

Inputs (data/raw/):
  bacdive/**/batch_*.json      BacDive v2 records (full ID crawl, cached by the soil pilot)
  ena_16s/batch_*.fasta        ENA sequences for the 16S accessions BacDive lists
  reference/{bacteria,archaea}.16SrRNA.fna.gz      NCBI RefSeq Targeted Loci (type material)
  reference/gtdb_r232_ssu_all.fna.gz + *_metadata.tsv.gz   GTDB release 232
Outputs (data/interim/reference/):
  bacdive_strains.csv, bacdive_16s_choice.csv, bacdive_genomes.csv
  cultured_v4.fasta            unique amplified-region sequences (ref_000001, ...)
  cultured_v4_members.csv      ref_id -> every source sequence with that exact region
  pcr_report_*.json            cutadapt in-silico PCR reports
"""

import json
import re

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.reference import (
    dereplicate,
    in_silico_pcr,
    is_isolate_genome,
    parse_bacdive_record,
    read_fasta,
    sequence_qc,
    write_fasta,
)

STUDY = "zheng2024_soil"


def load_bacdive(raw) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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


def load_ena_cache(raw) -> dict[str, str]:
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


def gtdb_isolates(raw) -> pd.DataFrame:
    """GTDB genomes (bacteria + archaea) with isolate status and type-material flag."""
    cols = ["accession", "ncbi_genome_category", "ncbi_strain_identifiers",
            "ncbi_type_material_designation", "ncbi_genbank_assembly_accession", "gtdb_taxonomy"]
    meta = pd.concat([pd.read_csv(raw / "reference" / f"gtdb_r232_{d}_metadata.tsv.gz", sep="\t",
                                  usecols=cols, low_memory=False) for d in ("bac120", "ar53")])
    meta["is_isolate"] = [is_isolate_genome(c, s) for c, s in
                          zip(meta["ncbi_genome_category"], meta["ncbi_strain_identifiers"],
                              strict=True)]
    return meta


def full_length_records(best16, meta, raw):
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


def main() -> None:
    raw = data_path("raw")
    out = data_path("interim") / "reference"
    out.mkdir(parents=True, exist_ok=True)
    cfg = load_config("thresholds")["reference"]
    primers = next(s for s in load_config("studies")["studies"] if s["id"] == STUDY)["primers"]

    strains, s16, genomes = load_bacdive(raw)
    best16 = choose_strain_16s(s16, load_ena_cache(raw), cfg)
    strains.to_csv(out / "bacdive_strains.csv", index=False)
    best16.drop(columns="seq").to_csv(out / "bacdive_16s_choice.csv", index=False)
    genomes.to_csv(out / "bacdive_genomes.csv", index=False)
    print(f"BacDive strains: {len(strains)}; with growth record: "
          f"{strains['has_growth_record'].sum()}; 16S status:\n"
          f"{best16['status'].value_counts().to_string()}")

    meta = gtdb_isolates(raw)
    meta.drop(columns="gtdb_taxonomy").to_csv(out / "gtdb_genomes.csv", index=False)
    print(f"GTDB genomes: {len(meta)}; isolates: {meta['is_isolate'].sum()}")

    full = out / "cultured_full.fasta"
    n_full = write_fasta(full_length_records(best16, meta, raw), full)
    print(f"full-length cultured sequences: {n_full}")

    region = out / "cultured_region_raw.fasta"
    in_silico_pcr(full, region, primers["forward"], primers["reverse"], cfg["pcr_error_rate"],
                  cfg["region_min_len"], cfg["region_max_len"], out / "pcr_report_v4.json")
    clean = [(rid, seq) for rid, seq in read_fasta(region) if not re.search(r"[^ACGT]", seq)]
    uniques, members = dereplicate((rid.split()[0], seq) for rid, seq in clean)
    write_fasta(uniques, out / "cultured_v4.fasta")
    members["source"] = members["member_id"].str.split("|").str[0]
    members.to_csv(out / "cultured_v4_members.csv", index=False)

    full_src = pd.Series([h.split("|")[0] for h, _ in read_fasta(full)]).value_counts()
    summary = pd.DataFrame({"full_length": full_src,
                            "region_extracted": members["source"].value_counts()})
    summary["region_unique_seqs"] = members.drop_duplicates(["ref_id", "source"])[
        "source"].value_counts()
    summary["extraction_rate"] = summary["region_extracted"] / summary["full_length"]
    summary.to_csv(out / "reference_summary.csv")
    print(summary)
    print(f"unique region sequences: {len(uniques)} (from {len(clean)} extracted, "
          f"ambiguous-base regions dropped)")
    full.unlink()
    region.unlink()


if __name__ == "__main__":
    main()
