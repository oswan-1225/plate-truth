"""How many cultured BacDive strains have a 16S sequence, before and after expanding the
reference (pilot: 52% of cultured strains lacked their own 16S).

Two strain sets: any culture-medium or culture-temperature entry not saying 'no'; and,
closer to the pilot, both a medium and a temperature entry. A strain is covered by:
  own_16s        its own BacDive 16S accession passed QC (the pilot's only route)
  own_genome     a genome BacDive lists for it is a GTDB isolate genome with a 16S gene
  type_species   it is a type strain and RefSeq type-material 16S exists for its species name
Output: results/phase0/strain_16s_coverage.csv (overall) and *_by_phylum.csv
"""

import pandas as pd

from plate_truth.config import data_path
from plate_truth.reference import accession_base, read_fasta


def refseq_species(raw) -> set[str]:
    """Species names in the RefSeq type-material 16S headers ('NR_1.1 Genus species 16S ...')."""
    names = set()
    for dom in ("bacteria", "archaea"):
        for header, _ in read_fasta(raw / "reference" / f"{dom}.16SrRNA.fna.gz"):
            words = header.split(" 16S")[0].split()[1:]
            names.add(" ".join(words[:2]))
    return names


def main() -> None:
    raw, interim = data_path("raw"), data_path("interim") / "reference"
    strains = pd.read_csv(interim / "bacdive_strains.csv")
    choice = pd.read_csv(interim / "bacdive_16s_choice.csv")
    genomes = pd.read_csv(interim / "bacdive_genomes.csv")
    members = pd.read_csv(interim / "cultured_v4_members.csv")

    cultured = strains[strains["has_growth_record"]].copy()
    ok16 = set(choice.loc[choice["status"] == "ok", "bacdive_id"])
    cultured["own_16s"] = cultured["bacdive_id"].isin(ok16)

    gtdb_ids = members.loc[members["source"] == "gtdb_isolate", "member_id"].str.split("|").str[1]
    with_ssu = {accession_base(a) for a in gtdb_ids.unique()}
    genomes["base"] = genomes["assembly"].map(accession_base)
    strain_genome = set(genomes.loc[genomes["base"].isin(with_ssu), "bacdive_id"])
    cultured["own_genome"] = cultured["bacdive_id"].isin(strain_genome)

    species = refseq_species(raw)
    cultured["type_species"] = (cultured["type_strain"].eq("yes")
                                & cultured["species"].isin(species))
    cultured["any_strain_level"] = cultured["own_16s"] | cultured["own_genome"]
    cultured["any"] = cultured["any_strain_level"] | cultured["type_species"]

    routes = ["own_16s", "own_genome", "type_species", "any_strain_level", "any"]
    tables = []
    for set_name, mask in [("any_growth_record", cultured["has_growth_record"]),
                           ("medium_and_temperature", cultured["has_medium_and_temp"])]:
        sub = cultured[mask]
        t = pd.DataFrame({"strain_set": set_name, "route": routes,
                          "n_cultured_strains": len(sub),
                          "n_covered": sub[routes].sum().values,
                          "share_covered": sub[routes].mean().values})
        tables.append(t)
    overall = pd.concat(tables, ignore_index=True)
    overall["share_lacking"] = 1 - overall["share_covered"]
    out = data_path("results") / "phase0"
    out.mkdir(parents=True, exist_ok=True)
    overall.to_csv(out / "strain_16s_coverage.csv", index=False)
    by_phylum = cultured.groupby("phylum")[["own_16s", "any"]].agg(["size", "mean"])
    by_phylum.columns = ["n", "share_own_16s", "_", "share_any"]
    by_phylum.drop(columns="_").sort_values("n", ascending=False).to_csv(
        out / "strain_16s_coverage_by_phylum.csv")
    print(overall.to_string())


if __name__ == "__main__":
    main()
