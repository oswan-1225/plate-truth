"""Build the cultured 16S reference, cut to the study's amplified region (V4, 515F/806R).

Inputs (data/raw/):
  bacdive/**/batch_*.json      BacDive v2 records (full ID crawl, cached by the soil pilot)
  ena_16s/batch_*.fasta        ENA sequences for the 16S accessions BacDive lists
  reference/{bacteria,archaea}.16SrRNA.fna.gz      NCBI RefSeq Targeted Loci (type material)
  reference/gtdb_r232_ssu_all.fna.gz + *_metadata.tsv.gz   GTDB release 232
Outputs (data/interim/reference/):
  bacdive_strains.csv, bacdive_16s_choice.csv, bacdive_genomes.csv
  cultured_v4_unfiltered.fasta unique amplified-region sequences (ref_000001, ...)
  cultured_v4_members_unfiltered.csv  ref_id -> every source sequence with that region
  (scripts/03b_filter_reference.py then removes contaminated genome assemblies)
  pcr_report_*.json            cutadapt in-silico PCR reports
"""

import re

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.reference import dereplicate, in_silico_pcr, read_fasta, write_fasta
from plate_truth.reference_build import (
    choose_strain_16s,
    full_length_records,
    gtdb_isolates,
    load_bacdive,
    load_ena_cache,
)

STUDY = "zheng2024_soil"


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
    write_fasta(uniques, out / "cultured_v4_unfiltered.fasta")
    members["source"] = members["member_id"].str.split("|").str[0]
    members.to_csv(out / "cultured_v4_members_unfiltered.csv", index=False)

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
