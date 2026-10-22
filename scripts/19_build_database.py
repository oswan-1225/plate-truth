"""Phase 2, step 6: build data/plate_truth.sqlite from the per-study pipeline outputs.

Reads config/census_studies.yaml, data/processed/sample_sheet.csv and primer_sets.csv,
config/media.yaml, and for every study that has been processed (scripts 16-18):
data/interim/phase2/{study}/asv_seqs.fasta, counts.csv, read_qc.csv, taxonomy.csv,
nearest_relative.csv, and counts_r1all.csv where it exists (every-read counts, D-030). Studies
not yet processed are skipped and listed. The reference_strains
table keeps only reference sequences that are some ASV's nearest relative.

Output: data/plate_truth.sqlite (not committed; rebuilt from the pipeline, D-003)
With --tag=NAME, studies that have a tagged reprocessing (data/interim/phase2/{study}__{NAME}/)
use it, and the database is written to data/interim/plate_truth__{NAME}.sqlite instead.
"""

import sys

import pandas as pd

from plate_truth.config import cli_option, data_path, load_config, project_root, study_dir
from plate_truth.database import condition_ids, conditions_table, is_organelle, write_tables
from plate_truth.media import media_table
from plate_truth.reference import read_fasta

STUDY_COLS = ["study_id", "study_group", "accession", "doi", "citation", "habitat",
              "host_associated", "tier", "linkage", "region", "primers", "primer_fwd",
              "primer_rev", "platform", "read_prep", "incubation_days"]


def study_outputs(sid: str, tag: str | None = None) -> dict | None:
    """The pipeline outputs of one study (its tagged reprocessing if one exists), or None if it
    has not been fully processed."""
    d = study_dir(sid, tag) if tag and study_dir(sid, tag).exists() else study_dir(sid)
    needed = ["asv_seqs.fasta", "counts.csv", "read_qc.csv", "taxonomy.csv",
              "nearest_relative.csv"]
    if not all((d / f).exists() for f in needed):
        return None
    asvs = pd.DataFrame([(h.split(";")[0], s) for h, s in read_fasta(d / "asv_seqs.fasta")],
                        columns=["asv_id", "sequence"])
    asvs["study_id"], asvs["length"] = sid, asvs["sequence"].str.len()
    counts = pd.read_csv(d / "counts.csv").rename(columns={"run": "run_accession"})
    qc = pd.read_csv(d / "read_qc.csv").rename(columns={"run": "run_accession"})
    tax = pd.read_csv(d / "taxonomy.csv")
    tax["is_organelle"] = is_organelle(tax).astype(int)
    every = (pd.read_csv(d / "counts_r1all.csv").rename(columns={"run": "run_accession"})
             if (d / "counts_r1all.csv").exists() else counts.iloc[0:0])
    return {"asvs": asvs, "abundances": counts, "abundances_every_read": every, "qc": qc,
            "taxonomy": tax, "nearest_relative": pd.read_csv(d / "nearest_relative.csv")}


def reference_members(nr: pd.DataFrame) -> pd.DataFrame:
    """Members (source sequences) of every reference sequence that is a best hit."""
    regions = data_path("interim") / "reference" / "regions"
    rows = []
    for ref, g in nr.groupby("reference"):
        targets = set()
        for col in ("targets_id1", "targets_id2"):
            for t in g[col].dropna():
                targets.update(t.split(";"))
        m = pd.read_csv(regions / f"cultured_{ref}_members.csv.gz")
        rows.append(m[m["ref_id"].isin(targets)].assign(reference=ref))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def media_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """The media mapping and class definitions from config/media.yaml."""
    cfg = load_config("media")
    classes = pd.DataFrame([{"medium_class": k, "definition": v}
                            for k, v in cfg["classes"].items()])
    return media_table(cfg.get("media") or []), classes


def main() -> None:
    tag = cli_option(sys.argv, "tag")
    studies = [s for s in load_config("census_studies")["studies"] if s["decision"] == "include"]
    sets = pd.read_csv(data_path("processed") / "primer_sets.csv")
    ref_of = {sid: r.ref_name for r in sets.itertuples() for sid in r.studies.split(";")}
    sheet = pd.read_csv(data_path("processed") / "sample_sheet.csv")
    sheet["condition_id"] = condition_ids(sheet)

    done, skipped, parts = [], [], {k: [] for k in
                                    ("asvs", "abundances", "abundances_every_read", "qc",
                                     "taxonomy", "nearest_relative")}
    for s in studies:
        out = study_outputs(s["id"], tag)
        if out is None:
            skipped.append(s["id"])
            continue
        done.append(s["id"])
        for k, v in out.items():
            parts[k].append(v)
    if not done:
        raise SystemExit("no processed studies yet")
    t = {k: pd.concat(v, ignore_index=True) for k, v in parts.items()}

    st = pd.DataFrame([{**{k: s.get(k) for k in STUDY_COLS if k != "study_id"},
                        "study_id": s["id"], "reference": ref_of.get(s["id"])}
                       for s in studies if s["id"] in done])
    st["host_associated"] = st["host_associated"].astype(int)
    samples = sheet[sheet["study_id"].isin(done)].merge(
        t["qc"][["run_accession", "reads_in", "reads_window", "reads_filtered", "reads_mapped"]],
        on="run_accession", how="left")
    media, classes = media_tables()
    tables = {
        "studies": st,
        "samples": samples,
        "conditions": conditions_table(samples),
        "asvs": t["asvs"],
        "abundances": t["abundances"][t["abundances"]["run_accession"].isin(samples["run_accession"])],
        "abundances_every_read": t["abundances_every_read"][
            t["abundances_every_read"]["run_accession"].isin(samples["run_accession"])],
        "taxonomy": t["taxonomy"],
        "nearest_relative": t["nearest_relative"],
        "reference_strains": reference_members(t["nearest_relative"]),
        "media": media[media["study_id"].isin(done)],
        "medium_class": classes,
    }
    db = (data_path("interim") / f"plate_truth__{tag}.sqlite" if tag
          else project_root() / "data" / "plate_truth.sqlite")
    write_tables(db, tables)
    print(f"built {db.name}: studies {len(done)}; skipped (not processed): {skipped}")
    for name, df in tables.items():
        print(f"  {name:18s} {len(df):>10,d} rows")


if __name__ == "__main__":
    main()
