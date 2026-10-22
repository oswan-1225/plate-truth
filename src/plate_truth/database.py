"""The harmonized SQLite database: table definitions and the functions that build its rows.

The database is rebuilt from the per-study pipeline outputs by scripts/19_build_database.py
and is not committed (D-003). Table and column meanings are in notes/schema.md; downstream
queries live in plate_truth.queries.
"""

import sqlite3
from pathlib import Path

import pandas as pd

SCHEMA = """
CREATE TABLE studies (
    study_id TEXT PRIMARY KEY, study_group TEXT, accession TEXT, doi TEXT, citation TEXT,
    habitat TEXT, host_associated INTEGER, tier TEXT, linkage TEXT, region TEXT, primers TEXT,
    primer_fwd TEXT, primer_rev TEXT, platform TEXT, read_prep TEXT, reference TEXT,
    incubation_days REAL);
CREATE TABLE samples (
    run_accession TEXT PRIMARY KEY, study_id TEXT REFERENCES studies, sample_accession TEXT,
    library_name TEXT, kind TEXT, source_sample TEXT, replicate REAL, condition_id TEXT,
    read_count INTEGER, reads_in INTEGER, reads_window INTEGER, reads_filtered INTEGER,
    reads_mapped INTEGER);
CREATE TABLE conditions (
    condition_id TEXT PRIMARY KEY, study_id TEXT REFERENCES studies, source_sample TEXT,
    medium TEXT, oxygen TEXT, temperature_c REAL, incubation_days REAL, condition_detail TEXT,
    n_plates INTEGER);
CREATE TABLE asvs (
    asv_id TEXT PRIMARY KEY, study_id TEXT REFERENCES studies, sequence TEXT, length INTEGER);
CREATE TABLE abundances (
    run_accession TEXT REFERENCES samples, asv_id TEXT REFERENCES asvs, count INTEGER);
CREATE TABLE abundances_every_read (
    run_accession TEXT REFERENCES samples, asv_id TEXT REFERENCES asvs, count INTEGER);
CREATE TABLE taxonomy (
    asv_id TEXT PRIMARY KEY REFERENCES asvs, domain TEXT, phylum TEXT, class TEXT,
    "order" TEXT, family TEXT, genus TEXT, copy_number REAL, copy_number_rank TEXT,
    is_organelle INTEGER);
CREATE TABLE nearest_relative (
    asv_id TEXT REFERENCES asvs, reference_set TEXT, reference TEXT,
    identity_id1 REAL, n_tied_id1 REAL, targets_id1 TEXT,
    identity_id2 REAL, n_tied_id2 REAL, targets_id2 TEXT,
    identity_mid REAL, n_tied_mid REAL, targets_mid TEXT,
    PRIMARY KEY (asv_id, reference_set));
CREATE TABLE reference_strains (
    reference TEXT, ref_id TEXT, member_id TEXT, source TEXT);
CREATE TABLE media (
    study_id TEXT REFERENCES studies, medium TEXT, medium_name TEXT, medium_class TEXT,
    class_rule TEXT, PRIMARY KEY (study_id, medium));
CREATE TABLE medium_class (
    medium_class TEXT PRIMARY KEY, definition TEXT);
CREATE INDEX idx_abund_run ON abundances (run_accession);
CREATE INDEX idx_abund_asv ON abundances (asv_id);
CREATE INDEX idx_abund_er_run ON abundances_every_read (run_accession);
CREATE INDEX idx_samples_study ON samples (study_id, kind);
CREATE INDEX idx_refstrains ON reference_strains (reference, ref_id);
"""


def condition_ids(sheet: pd.DataFrame) -> pd.Series:
    """A readable plating-condition id for plate rows (NaN otherwise):
    study|source|medium|oxygen|temperature|detail, missing parts written as 'NA'.

    Plates sharing all six parts are replicate plates of one condition and are pooled.
    """
    parts = ["study_id", "source_sample", "medium", "oxygen", "temperature_c",
             "condition_detail"]
    def fmt(v):
        if pd.isna(v):
            return "NA"
        return f"{v:g}" if isinstance(v, float) else str(v)
    cid = sheet[parts].apply(lambda r: "|".join(fmt(v) for v in r), axis=1)
    return cid.where(sheet["kind"] == "plate")


def conditions_table(sheet: pd.DataFrame) -> pd.DataFrame:
    """One row per plating condition with its fields and number of plate runs."""
    plates = sheet[sheet["kind"] == "plate"]
    cols = ["condition_id", "study_id", "source_sample", "medium", "oxygen", "temperature_c",
            "incubation_days", "condition_detail"]
    out = plates[cols].drop_duplicates("condition_id")
    n = plates.groupby("condition_id").size().rename("n_plates")
    return out.merge(n, left_on="condition_id", right_index=True).reset_index(drop=True)


def is_organelle(tax: pd.DataFrame) -> pd.Series:
    """True for chloroplast or mitochondrial 16S, and for ASVs not classified as Bacteria or
    Archaea at domain level (host or other non-prokaryotic sequences)."""
    organelle = tax["order"].eq("Chloroplast") | tax["family"].eq("Mitochondria")
    prokaryote = tax["domain"].isin(["Bacteria", "Archaea"])
    return organelle | ~prokaryote


def write_tables(db: Path, tables: dict[str, pd.DataFrame]) -> None:
    """Create a fresh database with SCHEMA and append each DataFrame to its table."""
    if db.exists():
        db.unlink()
    with sqlite3.connect(db) as con:
        con.executescript(SCHEMA)
        for name, df in tables.items():
            cols = [r[1] for r in con.execute(f"PRAGMA table_info({name})")]
            df[[c for c in cols if c in df.columns]].to_sql(name, con, if_exists="append",
                                                           index=False)
