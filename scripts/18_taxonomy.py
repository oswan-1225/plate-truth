"""Phase 2, step 5: taxonomy and 16S copy number for every ASV of each study.

Same method as script 06 (Phase 0): VSEARCH SINTAX (cutoff 0.8) against SILVA 138.2, and
rrnDB 5.10 pan-taxa mean copy numbers from the lowest classified rank rrnDB knows (Bacteria
mean if none). SINTAX runs in parallel chunks with a fixed seed, so results are reproducible.

Usage: python scripts/18_taxonomy.py STUDY_ID [STUDY_ID ...] [--tag=NAME]
Output: data/interim/phase2/{study}/taxonomy.csv
"""

import sys

import pandas as pd

from plate_truth.config import cli_option, data_path, load_config
from plate_truth.config import study_dir as folder
from plate_truth.taxonomy import asv_copy_number, copy_number_lookup, parse_sintax, \
    run_sintax_chunks

SINTAX_CUTOFF = 0.8
N_CHUNKS = 8


def main() -> None:
    db = data_path("interim") / "taxonomy" / "silva138.2_sintax.fasta"
    seed = load_config("thresholds")["resampling"]["seed"]
    pan = pd.read_csv(data_path("raw") / "reference" / "rrndb" /
                      "rrnDB-5.10_pantaxa_stats_NCBI.tsv", sep="\t")
    lookup = copy_number_lookup(pan)
    included = {s["id"] for s in load_config("census_studies")["studies"]
                if s["decision"] == "include"}
    for sid in [a for a in sys.argv[1:] if a in included]:
        study_dir = folder(sid, cli_option(sys.argv, "tag"))
        tsv = study_dir / "sintax.tsv"
        if not tsv.exists():
            run_sintax_chunks(study_dir / "asv_seqs.fasta", db, tsv, SINTAX_CUTOFF, seed,
                              N_CHUNKS, study_dir)
        tax = parse_sintax(tsv)
        tax = tax.join(asv_copy_number(tax, lookup, lookup[("domain", "Bacteria")]))
        tax.insert(0, "study_id", sid)
        tax.to_csv(study_dir / "taxonomy.csv")
        print(f"{sid}: ASVs {len(tax)}; phylum {tax['phylum'].notna().mean():.1%}, "
              f"genus {tax['genus'].notna().mean():.1%}; "
              f"chloroplast/mitochondria {tax['order'].eq('Chloroplast').sum()}/"
              f"{tax['family'].eq('Mitochondria').sum()}")


if __name__ == "__main__":
    main()
