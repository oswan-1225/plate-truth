"""Phase 5, step 3: shotgun surveys, SSU fragments placed on the V4 window, identity and
taxonomy (ANALYSIS_PLAN.md section 21).

Per study: download each selected run's MGnify v5.0 SSU fragment FASTA; place every fragment
(both strands) against SILVA 138.2 cut in silico to the 515F/806R window (VSEARCH, iddef 2,
>= 70%); keep fragments covering >= 100 window bases, oriented and trimmed to the covered part
(from the alignment's CIGAR string); identity of each trimmed fragment to the cultured V4
reference with terminal gaps excluded (iddef 2; internal gaps counted; best hit, no hit below
75%); SINTAX taxonomy (SILVA 138.2, cutoff 0.8).

The SILVA window database is built once (in-silico PCR with the V4 reference's primers,
dereplicated) in data/interim/phase5/silva_window/.

Input: data/processed/survey_runs_shotgun.csv (script 39).
Output per study in data/interim/phase5/shotgun/<study>/: fragments.csv (one row per window
fragment: run, overlap, identity, taxonomy) and placement_qc.csv (fragments per run).
Usage: python scripts/41_survey_shotgun.py [STUDY ...]
"""

import gzip
import sys

import pandas as pd

from plate_truth.amplicon import reverse_complement
from plate_truth.config import data_path, load_config
from plate_truth.http import download
from plate_truth.identity import nearest_relative_table, read_hits, search
from plate_truth.phase5 import window_span
from plate_truth.reference import dereplicate, in_silico_pcr, read_fasta, write_fasta
from plate_truth.taxonomy import parse_sintax, run_sintax_chunks
from plate_truth.tools import vsearch

MIN_ID = 0.75
SINTAX_CUTOFF = 0.8
N_CHUNKS = 8


def silva_window(cfg: dict) -> "Path":
    """SILVA 138.2 cut to the V4 window and dereplicated (built once)."""
    d = data_path("interim") / "phase5" / "silva_window"
    out = d / "silva_v4.fasta"
    if out.exists():
        return out
    d.mkdir(parents=True, exist_ok=True)
    sets = pd.read_csv(data_path("processed") / "primer_sets.csv").set_index("ref_name")
    ref = sets.loc[cfg["shotgun_window"]["reference"]]
    th = load_config("thresholds")["reference"]
    raw = d / "silva_v4_raw.fasta"
    if not raw.exists():
        in_silico_pcr(data_path("interim") / "taxonomy" / "silva138.2_sintax.fasta", raw,
                      ref["primer_fwd"], ref["primer_rev"], th["pcr_error_rate"],
                      *th["region_length"][ref["region"]], d / "pcr_report.json")
    uniques, _ = dereplicate((h.split(";")[0], s.upper()) for h, s in read_fasta(raw))
    write_fasta(uniques, out)
    return out


def fetch_fragments(study: str, runs: pd.DataFrame, out) -> int:
    """Download each run's SSU FASTA and write all fragments with ids <run>.<i>."""
    path = out / "fragments_raw.fasta"
    if path.exists():
        return 0
    recs = []
    for r in runs.itertuples():
        f = download(r.ssu_url, data_path("raw") / "mgnify" / "ssu" / study /
                     f"{r.run_accession}_SSU.fasta.gz")
        with gzip.open(f, "rt") as handle:
            lines = handle.read().splitlines()
        tmp = out / "tmp_run.fasta"
        tmp.write_text("\n".join(lines) + "\n")
        recs += [(f"{r.run_accession}.{i}", s.upper()) for i, (_, s) in
                 enumerate(read_fasta(tmp), 1)]
        tmp.unlink()
    return write_fasta(recs, path)


def place(out, silva, cfg: dict) -> pd.DataFrame:
    """Place fragments on the window; write the oriented, trimmed window fragments."""
    w = cfg["shotgun_window"]
    hits = out / "placement.tsv"
    if not hits.exists():
        vsearch("--usearch_global", out / "fragments_raw.fasta", "--db", silva, "--id",
                w["locator_min_id"], "--iddef", 2, "--strand", "both", "--maxaccepts", 8,
                "--top_hits_only", "--threads", 0, "--userout", hits, "--userfields",
                "query+target+qstrand+caln+id")
    h = pd.read_csv(hits, sep="\t", names=["query", "target", "qstrand", "caln", "id"])
    h = h.drop_duplicates("query")
    seqs = dict(read_fasta(out / "fragments_raw.fasta"))
    rows, recs = [], []
    for r in h.itertuples():
        s = seqs[r.query] if r.qstrand == "+" else reverse_complement(seqs[r.query])
        a, b, covered = window_span(r.caln, len(s))
        rows.append({"fragment": r.query, "placement_id": r.id, "overlap": covered,
                     "length": len(s)})
        if covered >= w["min_overlap"]:
            recs.append((r.query, s[a:b]))
    write_fasta(recs, out / "window.fasta")
    p = pd.DataFrame(rows)
    p["run"] = p["fragment"].str.rsplit(".", n=1).str[0]
    return p


def identities(out, ref_fasta) -> pd.DataFrame:
    path = out / "hits_id2.tsv"
    if not path.exists():
        search(out / "window.fasta", ref_fasta, path, 2, MIN_ID, 0, 0, top_hits_only=True)
    ids = [h for h, _ in read_fasta(out / "window.fasta")]
    return nearest_relative_table(ids, read_hits(path), id_cols=("id2",))


def main(studies: list[str]) -> None:
    cfg = load_config("phase5")
    seed = load_config("thresholds")["resampling"]["seed"]
    silva = silva_window(cfg)
    ref_fasta = (data_path("interim") / "reference" / "regions" /
                 f"cultured_{cfg['shotgun_window']['reference']}.fasta")
    sintax_db = data_path("interim") / "taxonomy" / "silva138.2_sintax.fasta"
    runs = pd.read_csv(data_path("processed") / "survey_runs_shotgun.csv", dtype=str)
    for study, r in runs.groupby("study", sort=False):
        if studies and study not in studies:
            continue
        out = data_path("interim") / "phase5" / "shotgun" / study
        if (out / "fragments.csv").exists():
            continue
        out.mkdir(parents=True, exist_ok=True)
        try:
            fetch_fragments(study, r, out)
            p = place(out, silva, cfg)
            nr = identities(out, ref_fasta)
            tsv = out / "sintax.tsv"
            if not tsv.exists():
                run_sintax_chunks(out / "window.fasta", sintax_db, tsv, SINTAX_CUTOFF, seed,
                                  N_CHUNKS, out)
            tax = parse_sintax(tsv)
            frag = p[p["overlap"] >= cfg["shotgun_window"]["min_overlap"]].set_index("fragment")
            frag = frag.join(nr.rename(columns={"identity_id2": "identity"})[["identity"]])
            frag = frag.join(tax)
            frag.index.name = "fragment"
            frag.to_csv(out / "fragments.csv")
            n_raw = pd.Series([h.rsplit(".", 1)[0] for h, _ in
                               read_fasta(out / "fragments_raw.fasta")]).value_counts()
            qc = pd.DataFrame({"fragments": n_raw,
                               "placed": p["run"].value_counts(),
                               "window": frag["run"].value_counts()}).fillna(0).astype(int)
            qc.index.name = "run"
            qc.to_csv(out / "placement_qc.csv")
            print(f"{study} ({r['habitat'].iloc[0]}): window fragments per run "
                  f"{qc['window'].to_dict()}; no cultured hit {frag['identity'].isna().mean():.1%}",
                  flush=True)
        except Exception as exc:                                   # noqa: BLE001
            print(f"{study}: FAILED {exc}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
