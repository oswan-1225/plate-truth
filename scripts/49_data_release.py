"""Phase 6: the harmonized paired dataset for a Zenodo release (data/release/plate-truth-data/).

Exports, from data/plate_truth.sqlite (script 19), data/processed/census.csv and the Phase 3
analysis rows (script 26): studies, sequencing runs, plating conditions, media and their
classes, ASVs with sequences, read counts, taxonomy, nearest-cultured-relative identities and
the reference members of each best hit, the study census, and the (ASV, plating condition)
table the model was fitted on. Writes data_dictionary.csv (every column described; the script
stops if one is not), README.md and MD5SUMS.txt. Nothing is uploaded.
"""

import hashlib
import shutil
import sqlite3
import sys

import pandas as pd

from plate_truth.config import data_path

OUT = data_path("raw").parent / "release" / "plate-truth-data"

TABLES = {  # file -> (description, source: sql table or callable name)
    "studies.csv": "One row per study dataset analysed (the seawater study is held out).",
    "runs.csv": "One row per sequencing run: environmental, inoculum, plate or control.",
    "conditions.csv": "One row per plating condition (replicate plates pooled).",
    "media.csv": "Medium label per study, its harmonized name and medium class.",
    "medium_classes.csv": "Definition of each medium class.",
    "asvs.csv.gz": "Amplicon sequence variants, denoised per study, with sequences.",
    "counts.csv.gz": "Read counts per run and ASV (zeros not stored); primary counts.",
    "counts_every_read.csv.gz": "Zheng et al. 2024 datasets only: counts from every read, "
                                "including reads deposited without their mate (the primary "
                                "counts for those datasets in the model).",
    "taxonomy.csv.gz": "SINTAX taxonomy against SILVA 138.2 and rrnDB copy numbers per ASV.",
    "nearest_relative.csv.gz": "Identity of each ASV to its nearest cultured relative, per "
                               "reference set and identity definition.",
    "reference_members.csv.gz": "Cultured strains and genomes behind each reference sequence "
                                "that is some ASV's best hit.",
    "census.csv": "Every candidate study screened, with the inclusion decision and reason.",
    "analysis_rows.csv.gz": "The unit of analysis: one row per (ASV, usable plating "
                            "condition), with the outcomes and covariates of the model.",
}

DICTIONARY = {
    "studies.csv": {
        "study_id": "Study dataset identifier (first author, year, habitat).",
        "study_group": "Independent study; datasets of one group are held out together.",
        "accession": "BioProject or SRA study accession of the raw reads.",
        "doi": "DOI of the paper describing the study.",
        "citation": "Short citation.",
        "habitat": "Habitat of the source samples.",
        "host_associated": "True for host-associated habitats (sponge, plant, food, rumen).",
        "tier": "Census data-quality grade (A best).",
        "linkage": "How plates link to source samples (per sample, per compartment, ...).",
        "region": "Amplified 16S region.",
        "primers": "Primer pair names.",
        "primer_fwd": "Forward primer sequence (IUPAC).",
        "primer_rev": "Reverse primer sequence (IUPAC).",
        "platform": "Sequencing platform.",
        "read_prep": "How reads entered the pipeline (merged pairs, single reads, long reads).",
        "reference": "Cultured reference (primer pair) the study's ASVs were compared with.",
        "incubation_days": "Incubation time in days, when constant across the study.",
    },
    "runs.csv": {
        "run_accession": "ENA/SRA run accession.",
        "study_id": "Study dataset.",
        "sample_accession": "BioSample accession.",
        "library_name": "Library name as deposited.",
        "kind": "environment, inoculum, plate or control.",
        "source_sample": "Source sample the run belongs to or traces back to (the enrichment "
                         "link between plates and their environmental sample).",
        "replicate": "Replicate label, where given.",
        "condition_id": "Plating condition (plate runs only).",
        "read_count": "Read or spot count reported by ENA.",
        "reads_in": "Reads entering primer cutting (merged pairs or single reads).",
        "reads_window": "Reads with the amplified window found.",
        "reads_filtered": "Reads passing the expected-error and length filters.",
        "reads_mapped": "Reads mapped to an ASV at 97% identity.",
    },
    "conditions.csv": {
        "condition_id": "Plating condition: study|source|medium|oxygen|temperature|detail.",
        "study_id": "Study dataset.",
        "source_sample": "Source sample plated.",
        "medium": "Medium label as used in the study.",
        "oxygen": "aerobic, anaerobic or microaerophilic (empty if not stated).",
        "temperature_c": "Incubation temperature in degrees C (empty if not resolved).",
        "incubation_days": "Incubation time in days.",
        "condition_detail": "Any other study-specific factor (antibiotic, dilution, ...).",
        "n_plates": "Number of plate runs pooled in the condition.",
    },
    "media.csv": {
        "study_id": "Study dataset.",
        "medium": "Medium label as used in the study.",
        "medium_name": "Harmonized medium name.",
        "medium_class": "Medium class by the ordered recipe rule (see medium_classes.csv).",
        "class_rule": "Recipe facts that decided the class.",
    },
    "medium_classes.csv": {
        "medium_class": "Medium class.",
        "definition": "How the class is assigned from recipe facts.",
    },
    "asvs.csv.gz": {
        "asv_id": "ASV identifier (study-prefixed; unique across studies).",
        "study_id": "Study dataset in which the ASV was denoised.",
        "sequence": "ASV sequence over the study's amplified window.",
        "length": "Sequence length (bases).",
    },
    "counts.csv.gz": {
        "run_accession": "Run.", "asv_id": "ASV.", "count": "Reads.",
    },
    "counts_every_read.csv.gz": {
        "run_accession": "Run.", "asv_id": "ASV.",
        "count": "Reads, counting every read once (shared among ASVs a read cannot separate, "
                 "in proportion to their paired counts).",
    },
    "taxonomy.csv.gz": {
        "asv_id": "ASV.",
        **{r: f"SILVA 138.2 {r} (SINTAX, cutoff 0.8; empty if not assigned)."
           for r in ("domain", "phylum", "class", "order", "family", "genus")},
        "copy_number": "16S copy number from rrnDB 5.10 at the lowest matched rank.",
        "copy_number_rank": "Rank at which the copy number was matched.",
        "is_organelle": "1 for chloroplast, mitochondrial or non-prokaryotic ASVs (excluded "
                        "from all analyses).",
    },
    "nearest_relative.csv.gz": {
        "asv_id": "ASV.",
        "reference_set": "all (primary: every cultured source), bacdive_strain or refseq_type.",
        "reference": "Cultured reference (primer pair).",
        "identity_id1": "Percent identity to the nearest cultured relative, all gaps counted "
                        "(VSEARCH --iddef 1; primary). Empty: no hit at 75% or above.",
        "n_tied_id1": "Number of reference sequences tied at that identity.",
        "targets_id1": "Tied reference sequences (';'-joined ref_id).",
        "identity_id2": "Identity with terminal gaps ignored (--iddef 2).",
        "n_tied_id2": "Ties under id2.", "targets_id2": "Tied reference sequences under id2.",
        "identity_mid": "Gap-blind identity (approximate; sensitivity analyses only).",
        "n_tied_mid": "Ties under the gap-blind identity.",
        "targets_mid": "Tied reference sequences under the gap-blind identity.",
    },
    "reference_members.csv.gz": {
        "reference": "Cultured reference (primer pair).",
        "ref_id": "Reference sequence (a unique region sequence).",
        "member_id": "A strain or genome carrying it: bacdive_strain|{BacDive ID}, "
                     "refseq_type|{accession} or gtdb_isolate|{genome}|{copy}.",
        "source": "bacdive_strain, refseq_type or gtdb_isolate.",
    },
    "census.csv": {
        "study_id": "Study dataset identifier.", "study_group": "Independent study.",
        "accession": "Raw-read accession.", "doi": "DOI.", "citation": "Short citation.",
        "habitat": "Habitat.", "host_associated": "Host-associated habitat.",
        "tier": "Data-quality grade.", "linkage": "Plate-to-source link.",
        "tier_note": "Reason for the grade.", "media": "Media used (';'-separated).",
        "temperature_c": "Incubation temperature(s).", "oxygen": "Oxygen regime(s).",
        "incubation_days": "Incubation time.", "region": "Amplified region.",
        "primers": "Primer pair.", "decision": "include or exclude.",
        "reason": "Reason for the decision.", "n_runs": "Runs in the accession.",
        "n_amplicon_runs": "Amplicon runs.", "n_environment_runs": "Environmental runs.",
        "n_plate_runs": "Plate runs.", "n_other_runs": "Other runs.",
        "platform_observed": "Platform recorded by ENA.", "note": "Further notes.",
    },
    "analysis_rows.csv.gz": {
        "study_id": "Study dataset.", "study_group": "Independent study (validation fold).",
        "source_sample": "Source sample.", "condition_id": "Plating condition.",
        "asv_id": "ASV (present in the source sample with at least 2 reads).",
        "env_reads": "Reads of the ASV in the source sample's pooled environmental runs.",
        "env_rel": "Its relative abundance among the source sample's prokaryotic reads.",
        "plate_reads": "Reads of the ASV in the condition's pooled plate runs.",
        "plate_rel": "Its relative abundance among the condition's prokaryotic plate reads.",
        "fold": "plate_rel / env_rel.",
        "enriched_on_plates": "Primary outcome: at least 5 plate reads and fold at least 2.",
        "detected_on_plates": "Secondary outcome: at least 2 plate reads.",
        "plate_depth": "Prokaryotic reads in the condition's pooled plates.",
        "identity": "Identity to the nearest cultured relative (identity_id1, reference set "
                    "'all'); empty = no hit at 75% or above.",
        "phylum": "SILVA phylum.",
        "medium_class": "Medium class of the plate (as used in the model; unknown if not "
                        "classed when the model was registered).",
        "temperature": "Temperature band: cold (<=20 C), moderate, warm (>=33 C) or unknown.",
        "oxygen": "Oxygen regime.", "incubation_days": "Incubation days.",
        "habitat": "Habitat.", "host_associated": "Host-associated habitat.",
    },
}

ROW_COLUMNS = ["study_id", "study_group", "source_sample", "condition_id", "asv_id",
               "env_reads", "env_rel", "plate_reads", "plate_rel", "fold", "enriched_on_plates",
               "detected_on_plates", "plate_depth", "identity", "phylum", "medium",
               "temperature", "oxygen", "incubation_days", "habitat", "host_associated"]


def frames(con) -> dict[str, pd.DataFrame]:
    q = lambda sql: pd.read_sql_query(sql, con)  # noqa: E731
    rows = pd.read_parquet(data_path("interim") / "phase3" / "rows" / "primary.parquet",
                           columns=ROW_COLUMNS).rename(columns={"medium": "medium_class"})
    rows["enriched_on_plates"] = rows["enriched_on_plates"].astype(int)
    rows["detected_on_plates"] = rows["detected_on_plates"].astype(int)
    return {
        "studies.csv": q("select * from studies"),
        "runs.csv": q("select * from samples"),
        "conditions.csv": q("select * from conditions"),
        "media.csv": q("select * from media"),
        "medium_classes.csv": q("select * from medium_class"),
        "asvs.csv.gz": q("select * from asvs"),
        "counts.csv.gz": q("select * from abundances"),
        "counts_every_read.csv.gz": q("select * from abundances_every_read"),
        "taxonomy.csv.gz": q("select * from taxonomy"),
        "nearest_relative.csv.gz": q("select * from nearest_relative"),
        "reference_members.csv.gz": q("select * from reference_strains"),
        "census.csv": pd.read_csv(data_path("processed") / "census.csv"),
        "analysis_rows.csv.gz": rows[[c if c != "medium" else "medium_class"
                                      for c in ROW_COLUMNS]],
    }


def md5(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


README = """# plate-truth: paired environmental and plate-grown 16S communities

Harmonized data from {n_studies} public studies (10 independent study groups) that sequenced
both an environmental sample and the communities grown from it on agar plates: soil,
activated sludge, a marine sponge, oak, broiler carcasses and rumen. All raw reads were
reprocessed with one pipeline (cutadapt, VSEARCH, UNOISE3 denoising, SINTAX against SILVA
138.2), and every ASV's identity to its nearest cultured relative was computed over the study's
own amplified region with all gaps counted.

- {n_runs} sequencing runs, {n_conditions} plating conditions, {n_asvs} ASVs.
- {n_rows} (ASV, usable plating condition) rows, the unit of the calibration analysis.

Raw reads remain at their archive accessions (studies.csv). The code that produced these files,
the pre-registered analysis plan and the results are at
https://github.com/oswan-1225/plate-truth.

## Files

{files}

Every column is described in data_dictionary.csv. MD5SUMS.txt lists checksums.

## Notes

- ASVs are denoised per study, so one sequence found in two studies has two identifiers.
- "Enriched on plates" (at least 5 plate reads and at least twice the source relative
  abundance) is an observation from sequencing; nothing in these files is a model prediction.
- Medium classes follow one ordered rule from recipe facts: selective, habitat extract,
  defined, dilute complex (at most 4 g/L organic matter), rich complex.
- The cultured reference combines BacDive strain 16S sequences, NCBI RefSeq type-material
  16S sequences and 16S genes of GTDB release 232 isolate genomes; reference sequences
  themselves are not redistributed here, only their identifiers.
"""


def main() -> None:
    con = sqlite3.connect(data_path("database"))
    data = frames(con)
    missing = [(f, c) for f, df in data.items() for c in df.columns
               if c not in DICTIONARY.get(f, {})]
    if missing:
        sys.exit(f"columns without a dictionary entry: {missing}")
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    for f, df in data.items():
        df.to_csv(OUT / f, index=False)
        print(f"{f}: {len(df):,} rows")
    pd.DataFrame([{"file": f, "column": c, "description": DICTIONARY[f][c]}
                  for f, df in data.items() for c in df.columns]).to_csv(
        OUT / "data_dictionary.csv", index=False)
    runs = data["runs.csv"]
    text = README.format(
        n_studies=len(data["studies.csv"]), n_runs=f"{len(runs):,}",
        n_conditions=f"{len(data['conditions.csv']):,}", n_asvs=f"{len(data['asvs.csv.gz']):,}",
        n_rows=f"{len(data['analysis_rows.csv.gz']):,}",
        files="\n".join(f"- `{f}`: {d}" for f, d in TABLES.items()))
    (OUT / "README.md").write_text(text, encoding="utf-8")
    sums = [f"{md5(p)}  {p.name}" for p in sorted(OUT.iterdir()) if p.name != "MD5SUMS.txt"]
    (OUT / "MD5SUMS.txt").write_text("\n".join(sums) + "\n", encoding="utf-8")
    size = sum(p.stat().st_size for p in OUT.iterdir()) / 1e6
    print(f"wrote {OUT} ({size:.0f} MB)")


if __name__ == "__main__":
    main()
