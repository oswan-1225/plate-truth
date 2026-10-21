"""Phase 2, step 2: one cultured reference per primer pair used by the included studies.

Same sources and rules as scripts 03 and 03b (Phase 0), generalized from V4 to every primer
pair in config/census_studies.yaml: collect full-length cultured 16S, cut each primer pair's
region in silico (both primer sites required, so every reference covers the whole window;
identity counts terminal gaps), drop regions with ambiguous bases, dereplicate, and remove GTDB
genomes whose 16S regions disagree with the genome's own phylum. A genome flagged in any region
is removed from every region (D-023). Full-length studies use the V4 reference (D-017).

Outputs:
  data/interim/reference/regions/cultured_{ref_name}.fasta, ..._members.csv.gz
  data/interim/reference/regions/gtdb_genome_check.csv     (all regions, union of drops)
  data/processed/primer_sets.csv                           (study -> reference)
  results/phase2/reference_regions.csv                     (coverage per region and source)
"""

import re

import pandas as pd

from plate_truth.config import data_path, load_config
from plate_truth.reference import dereplicate, in_silico_pcr, read_fasta, write_fasta
from plate_truth.reference_build import (
    choose_strain_16s,
    flag_contaminated_genomes,
    full_length_records,
    gtdb_isolates,
    load_bacdive,
    load_ena_cache,
    primer_sets,
)
from plate_truth.taxonomy import parse_sintax, run_sintax_chunks

N_CHUNKS = 8


def build_full_length(raw, out, cfg) -> tuple:
    """Write every full-length cultured sequence to one FASTA (reused if present)."""
    full = out / "cultured_full.fasta"
    meta = gtdb_isolates(raw)
    if not full.exists():
        strains, s16, _ = load_bacdive(raw)
        best16 = choose_strain_16s(s16, load_ena_cache(raw), cfg)
        n = write_fasta(full_length_records(best16, meta, raw), full)
        print(f"full-length cultured sequences: {n}")
    return full, meta


def cut_region(full, row, cfg, out) -> pd.DataFrame:
    """In-silico PCR, ambiguity filter and dereplication for one primer set."""
    lo, hi = cfg["region_length"][row.region]
    raw_region = out / f"tmp_{row.ref_name}.fasta"
    in_silico_pcr(full, raw_region, row.primer_fwd, row.primer_rev, cfg["pcr_error_rate"],
                  lo, hi, out / f"pcr_report_{row.ref_name}.json")
    clean = [(rid.split()[0], seq) for rid, seq in read_fasta(raw_region)
             if not re.search(r"[^ACGT]", seq)]
    raw_region.unlink()
    uniques, members = dereplicate(clean)
    write_fasta(uniques, out / f"cultured_{row.ref_name}_unfiltered.fasta")
    members["source"] = members["member_id"].str.split("|").str[0]
    return members


def main() -> None:
    raw = data_path("raw")
    out = data_path("interim") / "reference" / "regions"
    out.mkdir(parents=True, exist_ok=True)
    cfg = load_config("thresholds")["reference"]
    seed = load_config("thresholds")["resampling"]["seed"]
    studies = [s for s in load_config("census_studies")["studies"] if s["decision"] == "include"]
    fl = next(s for s in studies if s["id"] == cfg["full_length_cut_to"])
    sets = primer_sets(studies, fl)
    sets.to_csv(data_path("processed") / "primer_sets.csv", index=False)
    print(sets.to_string(index=False))

    full, meta = build_full_length(raw, out, cfg)
    full_src = pd.Series([h.split("|")[0] for h, _ in read_fasta(full)]).value_counts()
    genome_phylum = meta.set_index("accession")["gtdb_taxonomy"].str.split(";").str[1]
    db = data_path("interim") / "taxonomy" / "silva138.2_sintax.fasta"

    members, checks = {}, []
    for row in sets.itertuples(index=False):
        m = cut_region(full, row, cfg, out)
        tsv = out / f"sintax_{row.ref_name}.tsv"
        if not tsv.exists():
            run_sintax_chunks(out / f"cultured_{row.ref_name}_unfiltered.fasta", db, tsv, 0.8,
                              seed, N_CHUNKS, out)
        checks.append(flag_contaminated_genomes(m, parse_sintax(tsv)["phylum"], genome_phylum)
                      .assign(ref_name=row.ref_name))
        members[row.ref_name] = m
        print(f"{row.ref_name}: {m['ref_id'].nunique()} unique regions from {len(m)} sequences")

    check = pd.concat(checks)
    check.to_csv(out / "gtdb_genome_check.csv", index=False)
    dropped = set(check.loc[check["dropped"], "genome"])
    print(f"GTDB genomes dropped (union over regions): {len(dropped)}")

    rows = []
    for row in sets.itertuples(index=False):
        m = members[row.ref_name]
        genome = m["member_id"].str.split("|").str[1].where(m["source"] == "gtdb_isolate")
        kept = m[~genome.isin(dropped)]
        kept.to_csv(out / f"cultured_{row.ref_name}_members.csv.gz", index=False)
        keep = set(kept["ref_id"])
        seqs = [(r, s) for r, s in read_fasta(out / f"cultured_{row.ref_name}_unfiltered.fasta")
                if r in keep]
        write_fasta(seqs, out / f"cultured_{row.ref_name}.fasta")
        (out / f"cultured_{row.ref_name}_unfiltered.fasta").unlink()
        lengths = pd.Series([len(s) for _, s in seqs])
        for source, n_full in full_src.items():
            src = kept[kept["source"] == source]
            rows.append({"ref_name": row.ref_name, "region": row.region, "studies": row.studies,
                         "source": source, "full_length": n_full,
                         "region_extracted": len(src),
                         "extraction_rate": round(len(src) / n_full, 4),
                         "region_unique_seqs": src["ref_id"].nunique(),
                         "unique_regions_total": len(keep),
                         "len_p01": int(lengths.quantile(0.01)),
                         "len_median": int(lengths.median()),
                         "len_p99": int(lengths.quantile(0.99))})
    summary = pd.DataFrame(rows)
    res = data_path("results") / "phase2"
    res.mkdir(parents=True, exist_ok=True)
    summary.to_csv(res / "reference_regions.csv", index=False)
    print(summary.to_string(index=False))
    full.unlink()


if __name__ == "__main__":
    main()
