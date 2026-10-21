"""Nearest cultured relative for every Zheng et al. soil ASV.

Reference sets: all sources after removing contaminated genome assemblies (primary; script
03b), all sources before that filter ('all_unfiltered', sensitivity), BacDive strains only (the
pilot's reference) and RefSeq type material only.

For each set: an exhaustive search ranked by id1 and another ranked by id2 (top hits only),
so the best id1 and id2 identities are exact. The gap-blind identity (mid) is the best over
the union of those hits plus a bounded search's 64 accepted hits; it is an approximation and
is used only in the sensitivity analysis.

Hit files in data/interim/nearest_relative/ are reused if present; delete that folder after
changing the ASVs or the reference.

Output: data/processed/zheng2024_nearest_relative.csv (asv_id x reference set)
"""

import pandas as pd

from plate_truth.config import data_path
from plate_truth.identity import nearest_relative_table, read_hits, search
from plate_truth.reference import read_fasta, write_fasta

MIN_ID = 0.75          # hits below 75% identity are reported as 'no hit'
MAXACCEPTS, MAXREJECTS = 64, 256


def subset_fasta(uniques: dict, members: pd.DataFrame, source: str | None, path) -> None:
    """Write the unique reference sequences that have at least one member from source."""
    ids = uniques.keys() if source is None else \
        set(members.loc[members["source"] == source, "ref_id"])
    write_fasta(((r, uniques[r]) for r in uniques if r in ids), path)


def main() -> None:
    ref_dir = data_path("interim") / "reference"
    work = data_path("interim") / "nearest_relative"
    work.mkdir(parents=True, exist_ok=True)
    asvs = data_path("interim") / "zheng2024" / "asv_seqs.fasta"
    asv_ids = [h.split(";")[0] for h, _ in read_fasta(asvs)]
    uniques = dict(read_fasta(ref_dir / "cultured_v4.fasta"))
    members = pd.read_csv(ref_dir / "cultured_v4_members.csv")
    unfiltered = dict(read_fasta(ref_dir / "cultured_v4_unfiltered.fasta"))

    tables = []
    for name, source in [("all", None), ("all_unfiltered", None),
                         ("bacdive_strain", "bacdive_strain"), ("refseq_type", "refseq_type")]:
        db = work / f"ref_{name}.fasta"
        subset_fasta(unfiltered if name == "all_unfiltered" else uniques, members, source, db)
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
        hits = pd.concat(parts).drop_duplicates(["query", "target"])
        t = nearest_relative_table(asv_ids, hits)
        t["reference_set"] = name
        tables.append(t.reset_index())
        print(f"{name}: median id1 {t['identity_id1'].median():.2f}; "
              f"no hit {t['identity_id1'].isna().mean():.3%}", flush=True)
    nr = pd.concat(tables, ignore_index=True)
    nr.to_csv(data_path("processed") / "zheng2024_nearest_relative.csv", index=False)


if __name__ == "__main__":
    main()
