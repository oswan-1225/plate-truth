"""Queries on data/plate_truth.sqlite used by the analysis phases.

Every function takes an open sqlite3 connection and returns a DataFrame. Read counts are raw;
relative abundances and the grew/detected definitions are computed in plate_truth.growth so
that the same code serves Phase 0 and the pooled analysis.
"""

import sqlite3
from pathlib import Path

import pandas as pd

from plate_truth.config import project_root

COUNT_TABLES = ("abundances", "abundances_every_read")


def connect(db: Path | None = None) -> sqlite3.Connection:
    """Open the harmonized database (default data/plate_truth.sqlite)."""
    return sqlite3.connect(db or project_root() / "data" / "plate_truth.sqlite")


def _count_table(table: str) -> str:
    """The read-count table to use: paired-read counts (default) or, for the studies that have
    it, counts from every read (D-030)."""
    if table not in COUNT_TABLES:
        raise ValueError(f"table must be one of {COUNT_TABLES}, not {table!r}")
    return table


def environment_counts(con: sqlite3.Connection, prokaryotes_only: bool = True,
                       table: str = "abundances") -> pd.DataFrame:
    """Environmental reads pooled over each source sample's runs (PCR or sample replicates).

    Columns: study_id, source_sample, asv_id, count. Organelle and non-prokaryote ASVs are
    dropped unless prokaryotes_only is False. table: 'abundances' or 'abundances_every_read'
    (only the studies whose reads were partly deposited without mates have the latter).
    """
    org = "AND t.is_organelle = 0" if prokaryotes_only else ""
    return pd.read_sql(f"""
        SELECT s.study_id, s.source_sample, a.asv_id, SUM(a.count) AS count
        FROM {_count_table(table)} a JOIN samples s USING (run_accession)
        JOIN taxonomy t USING (asv_id)
        WHERE s.kind = 'environment' {org}
        GROUP BY s.study_id, s.source_sample, a.asv_id""", con)


def plate_counts(con: sqlite3.Connection, prokaryotes_only: bool = True,
                 table: str = "abundances") -> pd.DataFrame:
    """Plate reads pooled over the replicate plates of each plating condition.

    Columns: study_id, condition_id, asv_id, count. table as in environment_counts.
    """
    org = "AND t.is_organelle = 0" if prokaryotes_only else ""
    return pd.read_sql(f"""
        SELECT s.study_id, s.condition_id, a.asv_id, SUM(a.count) AS count
        FROM {_count_table(table)} a JOIN samples s USING (run_accession)
        JOIN taxonomy t USING (asv_id)
        WHERE s.kind = 'plate' {org}
        GROUP BY s.study_id, s.condition_id, a.asv_id""", con)


def plating_conditions(con: sqlite3.Connection, linked_only: bool = True) -> pd.DataFrame:
    """Plating conditions with study fields and medium class; by default only conditions whose
    source sample has at least one environmental run."""
    link = """AND EXISTS (SELECT 1 FROM samples e WHERE e.kind = 'environment'
              AND e.study_id = c.study_id AND e.source_sample = c.source_sample)""" \
        if linked_only else ""
    return pd.read_sql(f"""
        SELECT c.*, st.study_group, st.habitat, st.host_associated, st.tier, st.region,
               m.medium_name, m.medium_class
        FROM conditions c JOIN studies st USING (study_id)
        LEFT JOIN media m ON m.study_id = c.study_id AND m.medium = c.medium
        WHERE 1 = 1 {link}""", con)


def nearest_relative(con: sqlite3.Connection, reference_set: str = "all") -> pd.DataFrame:
    """Nearest cultured relative per ASV for one reference set (identities in percent, 0-100;
    NaN = no hit at or above 75%)."""
    return pd.read_sql("""SELECT n.*, a.study_id FROM nearest_relative n JOIN asvs a USING (asv_id)
                          WHERE n.reference_set = ?""", con, params=(reference_set,))


def asv_taxonomy(con: sqlite3.Connection) -> pd.DataFrame:
    """Taxonomy and 16S copy number per ASV."""
    return pd.read_sql("SELECT t.*, a.study_id FROM taxonomy t JOIN asvs a USING (asv_id)", con)
