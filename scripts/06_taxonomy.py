"""Taxonomy and 16S copy number for the Zheng et al. soil ASVs.

Taxonomy: VSEARCH SINTAX (cutoff 0.8) against SILVA 138.2 (DADA2-formatted NR99 training set,
Zenodo record 14169026), converted to SINTAX format.
Copy number: rrnDB 5.10 pan-taxa means (NCBI taxonomy), taken from the lowest classified rank
whose (cleaned) name rrnDB knows; ASVs matching nothing get the rrnDB mean for Bacteria.

Output: data/processed/zheng2024_taxonomy.csv
"""

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.taxonomy import (
    asv_copy_number,
    copy_number_lookup,
    parse_sintax,
    run_sintax,
    write_sintax_db,
)

SINTAX_CUTOFF = 0.8


def main() -> None:
    ref = data_path("raw") / "reference"
    work = data_path("interim") / "taxonomy"
    work.mkdir(parents=True, exist_ok=True)
    db = work / "silva138.2_sintax.fasta"
    if not db.exists():
        n = write_sintax_db(ref / "silva_nr99_v138.2_toGenus_trainset.fa.gz", db)
        print(f"SINTAX database: {n} sequences")
    asvs = data_path("interim") / "zheng2024" / "asv_seqs.fasta"
    run_sintax(asvs, db, work / "sintax.tsv", SINTAX_CUTOFF,
               load_config("thresholds")["resampling"]["seed"])
    tax = parse_sintax(work / "sintax.tsv")

    pan = pd.read_csv(ref / "rrndb" / "rrnDB-5.10_pantaxa_stats_NCBI.tsv", sep="\t")
    lookup = copy_number_lookup(pan)
    default = lookup[("domain", "Bacteria")]
    tax = tax.join(asv_copy_number(tax, lookup, default))
    tax.to_csv(data_path("processed") / "zheng2024_taxonomy.csv")
    print(f"ASVs: {len(tax)}; classified to phylum {tax['phylum'].notna().mean():.1%}, "
          f"genus {tax['genus'].notna().mean():.1%}")
    print(tax["copy_number_rank"].value_counts(normalize=True).round(3).to_string())
    print(tax["phylum"].value_counts().head(12).to_string())


if __name__ == "__main__":
    main()
