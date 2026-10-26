"""Phase 5, step 1: survey census and selection (ANALYSIS_PLAN.md section 20).

For each habitat and arm (amplicon, shotgun), the MGnify studies of the habitat's biome lineages
are screened in a seeded random order until the target number of studies is included or the
screening cap is reached. Every screened study is recorded with its decision and reason.

Amplicon: the study has a MGnify amplicon analysis; ENA runs that are Illumina, paired,
AMPLICON, >= 10,000 read pairs and not controls; one run per biosample, up to 8 chosen at
random; the first two chosen runs must carry a training primer pair in >= 50% of their first
2,000 read pairs (streamed), with reads long enough to merge.
Shotgun: the study has MGnify v5.0 metagenomic analyses; one analysis per biosample, taken in
a seeded random order, kept if MGnify reports >= 10,000 predicted SSU sequences and publishes
the SSU FASTA; up to 5 per study.

Outputs: data/processed/survey_census_<arm>.csv (every screened study), data/processed/
survey_runs_<arm>.csv (selected runs). API and ENA responses are cached in data/raw/mgnify/.

Usage: python scripts/39_survey_census.py [amplicon|shotgun ...]
"""

import sys
from statistics import median

import numpy as np
import pandas as pd

from plate_truth.amplicon import paired_fastq_urls
from plate_truth.census import ena_table
from plate_truth.config import data_path, load_config
from plate_truth.mgnify import (
    V2,
    cached_json,
    analysis_rows,
    analysis_summary,
    ssu_fasta_url,
    stream_fastq_head,
    studies_by_lineage,
    study_analyses,
)
from plate_truth.phase5 import (
    best_reference,
    choose_runs,
    habitat_of,
    has_word,
    merge_possible,
    primer_pair_shares,
    screening_order,
)

ENA_FIELDS = ("run_accession,sample_accession,study_accession,library_strategy,library_layout,"
              "instrument_platform,read_count,fastq_ftp,sample_alias,sample_title,library_name")
FORMAT_CHECKED = {"MGYS00003941"}          # its SSU file was opened to check the format


def raw(*parts):
    return data_path("raw").joinpath("mgnify", *parts)


def training_projects() -> set[str]:
    """ENA project accessions of every census study (plating studies are not surveys)."""
    return {s["accession"] for s in load_config("census_studies")["studies"] if s.get("accession")}


def candidates(habitat: str, h: dict, habitats: dict) -> pd.DataFrame:
    """Studies of a habitat's lineages whose own lineage maps to this habitat."""
    rows = {}
    for lineage in h["lineages"]:
        for s in studies_by_lineage(lineage, raw("api")):
            lin = (s.get("biome") or {}).get("lineage")
            if habitat_of(lin, habitats) == habitat:
                rows[s["accession"]] = {"study": s["accession"], "title": s["title"],
                                        "biome": lin,
                                        "ena": ";".join(s.get("ena_accessions") or [])}
    return pd.DataFrame(rows.values())


def ena_project(ena: str) -> str | None:
    """The PRJ accession among a study's ENA accessions (else the first)."""
    accs = [a for a in ena.split(";") if a]
    prj = [a for a in accs if a.startswith("PRJ")]
    return (prj or accs or [None])[0]


def ena_runs(project: str) -> pd.DataFrame:
    t = ena_table(project, "read_run", ENA_FIELDS, raw("ena", f"{project}.tsv"))
    t["read_count"] = pd.to_numeric(t["read_count"], errors="coerce").fillna(0)
    return t


def not_control(t: pd.DataFrame, words: list[str]) -> pd.Series:
    text = t["sample_alias"] + " " + t["sample_title"] + " " + t["library_name"]
    return ~text.map(lambda s: has_word(s, words))


def common_screen(row, sel: dict, training: set) -> str | None:
    """Study-level exclusions shared by both arms; returns a reason or None."""
    if row.study in FORMAT_CHECKED:
        return "format of its SSU file was checked before registration"
    if has_word(row.title, sel["exclude_study_words"]):
        return "title indicates cultures, isolates, enrichments, mock communities or plates"
    if set(row.ena.split(";")) & training:
        return "a training (plating) study"
    return None


def is_candidate(study: str, arm: str) -> bool:
    """A study is a candidate for the amplicon arm if MGnify has an amplicon analysis of it,
    and for the shotgun arm if it has a v5.0 metagenomic (raw reads) analysis. Studies of the
    other kinds are recorded but do not count towards the screening cap."""
    def wanted(r):
        if arm == "amplicon":
            return r["experiment_type"] == "Amplicon"
        return r["experiment_type"] == "Metagenomic" and r["pipeline_version"] == "V5"
    seen, page = 0, 1
    while page <= 30:
        d = cached_json(f"{V2}/studies/{study}/analyses/", {"page": page, "page_size": 100},
                        raw("api"))
        if any(wanted(r) for r in analysis_rows(d["items"])):
            return True
        seen += len(d["items"])
        if seen >= d["count"] or not d["items"]:
            return False
        page += 1
    return False


def screen_amplicon(row, sel: dict, seed: int) -> tuple[str, str, pd.DataFrame | None, dict]:
    """(decision, reason, chosen runs, details) for one amplicon candidate."""
    a = sel["amplicon"]
    project = ena_project(row.ena)
    t = ena_runs(project)
    ok = ((t["library_strategy"] == "AMPLICON") & (t["instrument_platform"] == "ILLUMINA")
          & (t["library_layout"] == "PAIRED") & (t["read_count"] >= a["min_ena_read_count"])
          & (t["fastq_ftp"].map(lambda f: len(paired_fastq_urls(f)) == 2))
          & not_control(t, sel["exclude_run_words"]))
    elig = t[ok]
    n_samples = elig["sample_accession"].nunique()
    if n_samples < a["min_runs_per_study"]:
        return "exclude", f"{n_samples} eligible runs (Illumina paired amplicon, >= " \
            f"{a['min_ena_read_count']} pairs, not controls); need {a['min_runs_per_study']}", \
            None, {"project": project}
    chosen = choose_runs(elig, a["runs_per_study"], seed)
    refs = a["references"]
    best, details = [], {"project": project, "eligible_runs": n_samples}
    for i, r in chosen.head(a["screen_runs"]).iterrows():
        u1, u2 = paired_fastq_urls(r["fastq_ftp"])
        try:
            s1 = stream_fastq_head(u1, a["screen_reads"])
            s2 = stream_fastq_head(u2, a["screen_reads"])
        except Exception as exc:                                   # noqa: BLE001
            return "exclude", f"reads could not be streamed ({exc})", None, details
        pairs = list(zip(s1, s2))
        shares = primer_pair_shares(pairs, refs, a["primer_max_offset"],
                                    a["primer_max_mismatches"])
        ref = best_reference(shares, a["min_primer_pair_share"])
        details[f"run{i + 1}"] = r["run_accession"]
        details[f"run{i + 1}_shares"] = ";".join(f"{k}={v:.2f}" for k, v in shares.items())
        len1, len2 = median(map(len, s1)), median(map(len, s2))
        details[f"run{i + 1}_lengths"] = f"{len1:g}+{len2:g}"
        if ref is None:
            return "exclude", "no training primer pair at the read starts", None, details
        if not merge_possible(len1, len2, a["amplicon_length"][ref], a["min_merge_overlap"]):
            return "exclude", f"reads ({len1:g}+{len2:g} bp) too short to merge {ref}", None, \
                details
        best.append(ref)
    if len(set(best)) != 1:
        return "exclude", f"screened runs disagree on the primer pair ({best})", None, details
    details["reference"] = best[0]
    return "include", f"primer pair {best[0]}", chosen.assign(reference=best[0]), details


def screen_shotgun(row, sel: dict, seed: int) -> tuple[str, str, pd.DataFrame | None, dict]:
    """(decision, reason, chosen runs, details) for one shotgun candidate."""
    s = sel["shotgun"]
    an = pd.DataFrame(analysis_rows(study_analyses(row.study, raw("api"))))
    an = an[(an["experiment_type"] == "Metagenomic") & (an["pipeline_version"] == "V5")]
    project = ena_project(row.ena)
    t = ena_runs(project)
    t = t[(t["instrument_platform"] == "ILLUMINA") & not_control(t, sel["exclude_run_words"])]
    an = an.merge(t[["run_accession", "sample_accession"]].rename(
        columns={"sample_accession": "ena_sample"}), on="run_accession")
    an = an.sort_values("run_accession").drop_duplicates("ena_sample")
    order = an.iloc[np.random.default_rng(seed).permutation(len(an))]
    kept, tried = [], 0
    for _, r in order.iterrows():
        if len(kept) >= s["runs_per_study"] or tried >= s["max_queried_runs"]:
            break
        tried += 1
        ssu = int(analysis_summary(r["analysis"], raw("api")).get("Predicted SSU sequences", 0))
        if ssu < s["min_predicted_ssu"]:
            continue
        url = ssu_fasta_url(r["analysis"], raw("api"))
        if url:
            kept.append({**r.to_dict(), "predicted_ssu": ssu, "ssu_url": url})
    details = {"project": project, "eligible_runs": len(an), "queried": tried}
    if len(kept) < s["min_runs_per_study"]:
        return "exclude", f"{len(kept)} runs with >= {s['min_predicted_ssu']} predicted SSU " \
            f"sequences and an SSU FASTA (of {tried} queried); need {s['min_runs_per_study']}", \
            None, details
    chosen = pd.DataFrame(kept).sort_values("run_accession")
    return "include", f"{len(chosen)} runs", chosen, details


def census_arm(arm: str, cfg: dict) -> tuple[list[dict], list[pd.DataFrame]]:
    sel, habitats = cfg["selection"], cfg["habitats"]
    target = sel[arm]["studies_per_habitat"]
    training = training_projects()
    rows, runs = [], []
    for habitat, h in habitats.items():
        cand = candidates(habitat, h, habitats)
        order = screening_order(cand["study"].tolist(), sel["seed"]) if len(cand) else []
        cand = cand.set_index("study")
        included = screened = 0
        for acc in order:
            if included >= target or screened >= sel[arm]["max_screened_per_habitat"]:
                break
            row = cand.loc[acc]
            row = row.copy()
            row.name = acc
            row["study"] = acc
            try:
                candidate = is_candidate(acc, arm)
            except Exception as exc:                               # noqa: BLE001
                print(f"{arm} {habitat} {acc}: analyses not retrieved ({exc})", flush=True)
                candidate = False
            if not candidate:
                rows.append({"arm": arm, "habitat": habitat, "screen_rank": None, "study": acc,
                             "title": row["title"], "biome": row["biome"], "ena": row["ena"],
                             "decision": "not_candidate",
                             "reason": f"no MGnify {arm} analysis of the required kind"})
                continue
            screened += 1
            reason = common_screen(row, sel, training)
            chosen, details = None, {}
            if reason:
                decision = "exclude"
            else:
                try:
                    fn = screen_amplicon if arm == "amplicon" else screen_shotgun
                    decision, reason, chosen, details = fn(row, sel, sel["seed"])
                except Exception as exc:                           # noqa: BLE001
                    decision, reason = "exclude", f"screen failed: {exc}"
            included += decision == "include"
            rows.append({"arm": arm, "habitat": habitat, "screen_rank": screened, "study": acc,
                         "title": row["title"], "biome": row["biome"], "ena": row["ena"],
                         "decision": decision, "reason": reason, **details})
            if chosen is not None:
                runs.append(chosen.assign(arm=arm, habitat=habitat, study=acc))
            print(f"{arm} {habitat} #{screened} {acc}: {decision} ({reason})", flush=True)
        print(f"{arm} {habitat}: {included} included of {screened} screened "
              f"({len(cand)} candidates)", flush=True)
    return rows, runs


def main(arms: list[str]) -> None:
    cfg = load_config("phase5")
    out = data_path("processed")
    for arm in arms:
        rows, runs = census_arm(arm, cfg)
        pd.DataFrame(rows).to_csv(out / f"survey_census_{arm}.csv", index=False)
        if runs:
            pd.concat(runs, ignore_index=True).to_csv(out / f"survey_runs_{arm}.csv",
                                                      index=False)


if __name__ == "__main__":
    main(sys.argv[1:] or ["amplicon", "shotgun"])
