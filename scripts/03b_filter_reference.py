"""Remove contaminated genome assemblies from the cultured reference.

Some GTDB genomes from isolates carry 16S sequences from unrelated organisms (assembly
contamination); GTDB's FILE_DESCRIPTIONS.txt warns about this for ssu_all. Such a sequence
would give an environmental ASV a 'cultured relative' that was never cultured.

Rule: every unique reference region is classified with SINTAX against SILVA 138.2. A GTDB
genome is dropped entirely if any of its 16S regions is classified (cutoff 0.8) to a phylum
different from the genome's own GTDB phylum. BacDive and RefSeq sequences are kept as they are.

Inputs: cultured_v4_unfiltered.fasta, cultured_v4_members_unfiltered.csv (from script 03)
Outputs (data/interim/reference/): cultured_v4.fasta, cultured_v4_members.csv (retained only),
cultured_v4_members_flagged.csv (all members with the check), gtdb_genome_check.csv;
results/phase0/reference_contamination.csv
"""

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.reference import read_fasta, write_fasta
from plate_truth.taxonomy import parse_sintax, phylum_consistent, run_sintax_chunks

N_CHUNKS = 8


def main() -> None:
    ref = data_path("interim") / "reference"
    work = ref / "sintax"
    work.mkdir(exist_ok=True)
    seed = load_config("thresholds")["resampling"]["seed"]
    db = data_path("interim") / "taxonomy" / "silva138.2_sintax.fasta"
    tsv = work / "reference_sintax.tsv"
    if not tsv.exists():
        run_sintax_chunks(ref / "cultured_v4_unfiltered.fasta", db, tsv, 0.8, seed, N_CHUNKS,
                          work)
    ssu_phylum = parse_sintax(tsv)["phylum"]

    members = pd.read_csv(ref / "cultured_v4_members_unfiltered.csv")
    members["ssu_phylum"] = members["ref_id"].map(ssu_phylum)
    is_gtdb = members["source"] == "gtdb_isolate"
    members["genome"] = members["member_id"].str.split("|").str[1].where(is_gtdb)
    raw = data_path("raw") / "reference"
    tax = pd.concat([pd.read_csv(raw / f"gtdb_r232_{d}_metadata.tsv.gz", sep="\t",
                                 usecols=["accession", "gtdb_taxonomy"])
                     for d in ("bac120", "ar53")])
    genome_phylum = tax.set_index("accession")["gtdb_taxonomy"].str.split(";").str[1]
    members["genome_phylum"] = members["genome"].map(genome_phylum)
    members["consistent"] = [not g or phylum_consistent(gp, sp) for g, gp, sp in
                             zip(members["genome"].notna(), members["genome_phylum"],
                                 members["ssu_phylum"], strict=True)]

    g = members[is_gtdb].groupby("genome").agg(n_ssu=("ref_id", "size"),
                                               n_inconsistent=("consistent", lambda c: (~c).sum()))
    g["dropped"] = g["n_inconsistent"] > 0
    g.to_csv(ref / "gtdb_genome_check.csv")
    dropped = set(g.index[g["dropped"]])
    members["retained"] = ~members["genome"].isin(dropped)
    members.to_csv(ref / "cultured_v4_members_flagged.csv", index=False)
    kept = members[members["retained"]]
    kept[["ref_id", "member_id", "source"]].to_csv(ref / "cultured_v4_members.csv", index=False)
    keep_refs = set(kept["ref_id"])
    write_fasta(((r, s) for r, s in read_fasta(ref / "cultured_v4_unfiltered.fasta")
                 if r in keep_refs), ref / "cultured_v4.fasta")

    summary = pd.DataFrame({
        "gtdb_genomes_with_16s": [len(g)],
        "gtdb_genomes_dropped": [len(dropped)],
        "gtdb_16s_regions_dropped": [int((~members["retained"]).sum())],
        "unique_regions_before": [members["ref_id"].nunique()],
        "unique_regions_after": [len(keep_refs)],
        "max_16s_per_dropped_genome": [int(g.loc[g["dropped"], "n_ssu"].max())],
    })
    out = data_path("results") / "phase0"
    summary.to_csv(out / "reference_contamination.csv", index=False)
    print(summary.T.to_string())
    print(g[g["dropped"]].sort_values("n_ssu", ascending=False).head(10).to_string())


if __name__ == "__main__":
    main()
