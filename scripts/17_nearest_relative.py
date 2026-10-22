"""Phase 2, step 4: nearest cultured relative for every ASV of each study.

Same method as script 05 (Phase 0), per study and against the reference cut with that
study's own primer pair (script 15): exhaustive searches ranked by id1 (primary, gaps counted)
and by id2 (terminal gaps ignored), top hits only, plus a bounded id1 search whose hits give
the approximate gap-blind identity (mid) used only in sensitivity analyses.

Reference sets: 'all' (every cultured source after the contaminated-genome filter; primary),
'bacdive_strain' and 'refseq_type' (single sources, for the D-016 question of whether the
relative's source matters).

Usage: python scripts/17_nearest_relative.py STUDY_ID [STUDY_ID ...] [--tag=NAME]
Output: data/interim/phase2/{study}/nearest_relative.csv (asv_id x reference set); hit files
are cached in data/interim/phase2/{study}/nr/ (delete it after changing ASVs or reference).
"""

import sys

import pandas as pd

from plate_truth.config import cli_option, data_path, load_config, study_dir
from plate_truth.identity import nearest_relative_table, read_hits, search
from plate_truth.reference import read_fasta, write_fasta

MIN_ID = 0.75          # hits below 75% identity are reported as 'no hit' (as in Phase 0)
MAXACCEPTS, MAXREJECTS = 64, 256
SETS = [("all", None), ("bacdive_strain", "bacdive_strain"), ("refseq_type", "refseq_type")]


def reference_for(study_id: str) -> str:
    """The reference name (primer set) script 15 assigned to a study."""
    sets = pd.read_csv(data_path("processed") / "primer_sets.csv")
    return next(r.ref_name for r in sets.itertuples() if study_id in r.studies.split(";"))


def nearest_relatives(study_id: str, tag: str | None = None) -> pd.DataFrame:
    """Search one study's ASVs against each reference set; return the combined table."""
    ref_name = reference_for(study_id)
    regions = data_path("interim") / "reference" / "regions"
    uniques = dict(read_fasta(regions / f"cultured_{ref_name}.fasta"))
    members = pd.read_csv(regions / f"cultured_{ref_name}_members.csv.gz")
    folder = study_dir(study_id, tag)
    asvs = folder / "asv_seqs.fasta"
    asv_ids = [h.split(";")[0] for h, _ in read_fasta(asvs)]
    work = folder / "nr"
    work.mkdir(exist_ok=True)

    tables = []
    for name, source in SETS:
        ids = uniques.keys() if source is None else \
            set(members.loc[members["source"] == source, "ref_id"])
        db = work / f"ref_{name}.fasta"
        write_fasta(((r, s) for r, s in uniques.items() if r in ids), db)
        parts = []
        for iddef in (1, 2):
            path = work / f"hits_{name}_exhaustive_id{iddef}.tsv"
            if not path.exists():
                search(asvs, db, path, iddef, MIN_ID, 0, 0, top_hits_only=True)
            parts.append(read_hits(path))
        path = work / f"hits_{name}_bounded.tsv"
        if not path.exists():
            search(asvs, db, path, 1, MIN_ID, MAXACCEPTS, MAXREJECTS)
        parts.append(read_hits(path))
        t = nearest_relative_table(asv_ids, pd.concat(parts).drop_duplicates(["query", "target"]))
        t["reference_set"] = name
        tables.append(t.reset_index())
        print(f"{study_id} {name}: median id1 {t['identity_id1'].median():.2f}; "
              f"no hit {t['identity_id1'].isna().mean():.3%}", flush=True)
    out = pd.concat(tables, ignore_index=True)
    out.insert(0, "study_id", study_id)
    out["reference"] = ref_name
    return out


def main() -> None:
    included = {s["id"] for s in load_config("census_studies")["studies"]
                if s["decision"] == "include"}
    for sid in [a for a in sys.argv[1:] if a in included]:
        tag = cli_option(sys.argv, "tag")
        nr = nearest_relatives(sid, tag)
        nr.to_csv(study_dir(sid, tag) / "nearest_relative.csv", index=False)


if __name__ == "__main__":
    main()
