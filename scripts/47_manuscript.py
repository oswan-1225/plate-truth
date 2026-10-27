"""Phase 6: build the preprint from manuscript/preprint_template.md.

Every {{placeholder}} is filled from results/ (and the census and survey tables in
data/processed/), so no result is typed into the text. Every [@key] citation becomes an
author-year label, and the reference list is built from each DOI's Crossref record (cached in
data/raw/manuscript/crossref/). DOIs: config/manuscript.yaml and, for the census studies,
config/census_studies.yaml (checked by script 21, which must pass first).

Outputs: manuscript/preprint.md, manuscript/preprint.docx (pandoc, if available) and
results/phase6/manuscript_references.csv (key, DOI, label, Crossref title) for auditing.
Run after scripts 21 and 46.
"""

import datetime
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from plate_truth.config import data_path, load_config, project_root
from plate_truth.http import get
from plate_truth.manuscript import (
    citation_labels,
    cited_keys,
    clean,
    fill_placeholders,
    format_reference,
    record_year,
    replace_citations,
    surnames,
)
from plate_truth.report import (
    fmt_count,
    fmt_interval,
    fmt_pct,
    fmt_share,
    study_region_labels,
)
from plate_truth.tools import cutadapt_exe, vsearch_exe

RES = data_path("results")
HABITAT_NAMES = {"soil": "Soil", "activated_sludge": "Activated sludge", "plant": "Plant",
                 "marine_water": "Marine water", "marine_sediment": "Marine sediment",
                 "freshwater": "Freshwater", "human_gut": "Human gut"}


def csv(*parts) -> pd.DataFrame:
    return pd.read_csv(RES.joinpath(*parts))


def est(x, lo, hi, digits=2) -> str:
    """'0.50 (0.35-0.59)'."""
    return f"{fmt_share(x, digits)} {fmt_interval(lo, hi, digits)}"


def est_pct(x, lo, hi, digits=1) -> str:
    """'1.1% (0.2%-5.5%)'."""
    return f"{fmt_pct(x, digits)} {fmt_interval(lo, hi, digits, pct=True)}"


def signed(x, lo, hi) -> str:
    """'+0.76 (0.24 to 1.29)' for a log-odds effect."""
    return f"{x:+.2f} ({lo:.2f} to {hi:.2f})"


def parse_interval(text: str) -> str:
    """'mean AUC gain 0.073 (0.029-0.123)' -> '0.073 (0.029-0.123)' (last number + interval)."""
    m = re.search(r"(-?\d+\.\d+) \((-?\d+\.\d+)-(-?\d+\.\d+)\)\s*$", text)
    if not m:
        raise ValueError(f"no estimate with interval in {text!r}")
    return f"{m.group(1)} ({m.group(2)}-{m.group(3)})"


# ---------------------------------------------------------------- values, section by section
def census_values() -> dict:
    interim = data_path("interim") / "census"
    census = pd.read_csv(data_path("processed") / "census.csv")
    links = csv("phase2", "qc_asvs_links.csv")
    rows = csv("phase3", "rows_summary.csv")
    rows = rows[rows["variant"] == "primary"]
    support = csv("phase3", "primary", "level_support.csv")
    medium = support[(support["factor"] == "medium") & (support["level"] != "unknown")]
    runs = links[["runs_control", "runs_environment", "runs_inoculum", "runs_plate"]].sum()
    return {
        "n_articles": fmt_count(len(pd.read_csv(interim / "candidates.csv"))),
        "n_projects": fmt_count(len(pd.read_csv(interim / "accession_evidence.csv"))),
        "n_census_groups": str(census.loc[census["decision"] == "include",
                                          "study_group"].nunique()),
        "n_datasets": str(rows["study_id"].nunique()),
        "n_groups": str(rows["study_group"].nunique()),
        "n_runs": fmt_count(runs.sum()),
        "n_env_runs": fmt_count(runs["runs_environment"]),
        "n_plate_runs": fmt_count(runs["runs_plate"]),
        "n_asvs": fmt_count(links["asvs"].sum()),
        "n_conditions": fmt_count(links["conditions"].sum()),
        "n_usable_conditions": fmt_count(rows["conditions"].sum()),
        "n_rows": fmt_count(rows["rows"].sum()),
        "n_enriched": fmt_count(rows["events"].sum()),
        "min_conditions": str(rows["conditions"].min()),
        "max_conditions": str(rows["conditions"].max()),
        "min_sources": str(links["sources_linked"].min()),
        "max_sources": str(links["sources_linked"].max()),
        "n_medium_classes": str(len(medium)),
        "warm_studies": str(len(support.query("factor == 'temperature' and level == 'warm'")
                                ["studies"].iloc[0].split(";"))),
    }


def curve_at(curve: pd.DataFrame, identity: float) -> pd.Series:
    c = curve[(curve["no_hit"] == 0) & np.isclose(curve["identity"], identity)]
    if len(c) != 1:
        raise ValueError(f"curve has {len(c)} rows at identity {identity}")
    return c.iloc[0]


def phase3_values() -> dict:
    res = ("phase3", "primary")
    summ = csv(*res, "heldout_summary.csv").set_index("model")
    by = csv(*res, "heldout_by_group.csv")
    by = by[by["model"] == "primary"].set_index("study_group")
    hyp = csv(*res, "hypotheses.csv").set_index("hypothesis")
    rec = csv(*res, "recovery_by_group.csv").set_index("study_group")
    curve = csv(*res, "curve.csv")
    host = csv(*res, "curve_by_host.csv")
    calib = csv(*res, "calibration_by_identity_bin.csv")
    nohit = calib[(calib["model"] == "primary") & (calib["identity_bin"] == "none")].iloc[0]
    sens = csv("phase3", "sensitivity_criteria.csv")
    sens = sens[sens["h1_model"] == "primary"].set_index("variant")
    others = sens.drop(index="primary")
    crit = load_config("phase3")["criteria"]
    oe_lo, oe_hi = crit["oe_range"]
    recal = csv("phase3", "exploratory", "local_recalibration_summary.csv")
    recal = recal.set_index(["pilot_fraction", "method"])
    v = {}
    for key, ident in (("curve_90", 90), ("curve_945", 94.5), ("curve_97", 97), ("curve_99", 99),
                       ("curve_100", 100)):
        r = curve_at(curve, ident)
        v[key] = est_pct(r["p"], r["p_lo"], r["p_hi"])
    r99 = curve_at(curve, 99)
    v["curve_99_p10"], v["curve_99_p90"] = fmt_pct(r99["study_p10"]), fmt_pct(r99["study_p90"])
    for name, h in (("env", -1.0), ("host", 1.0)):
        for ident in (90, 99):
            r = curve_at(host[host["host"] == h], ident)
            v[f"{name}_{ident}"] = est_pct(r["p"], r["p_lo"], r["p_hi"])
    io = summ.loc["identity_only"]
    v.update({
        "nohit_k": fmt_count(nohit["k"]), "nohit_n": fmt_count(nohit["n"]),
        "h1_groups": str(int(io["groups_auc_above_0.5"])),
        "h1_auc": est(io["mean_auc"], io["mean_auc_lo"], io["mean_auc_hi"], 3),
        "h4_gain": parse_interval(hyp.loc["H4 identity adds to covariates", "statistic"]),
        "primary_auc": est(*summ.loc["primary", ["mean_auc", "mean_auc_lo", "mean_auc_hi"]], 3),
        "lgbm_auc": est(*summ.loc["lightgbm", ["mean_auc", "mean_auc_lo", "mean_auc_hi"]], 3),
        "auc_yu2020": f"{by.loc['yu2020', 'auc']:.3f}",
        "auc_li2023": f"{by.loc['li2023', 'auc']:.3f}",
        "oe_min": f"{by['oe'].min():.2f}", "oe_max": f"{by['oe'].max():.2f}",
        "oe_low_n": str(int((by["oe"] < oe_lo).sum())),
        "oe_high_n": str(int((by["oe"] > oe_hi).sum())),
        "h2_groups": str(int(summ.loc["primary", "groups_calibrated"])),
        "h3_groups": str(int((rec["proxy_over_observed"] > 1).sum())),
        "h3_ratio": parse_interval(hyp.loc["H3 shortcut overstates recovery", "statistic"]),
        "ratio_min": f"{rec['proxy_over_observed'].min():.1f}",
        "ratio_max": f"{rec['proxy_over_observed'].max():.0f}",
        "n_sens": str(len(others)),
        "n_sens_same": str(int((~others["changes_conclusion"]).sum())),
        "detected_ratio": f"{sens.loc['detected', 'median_ratio']:.2f}",
        "detected_ratio_lo": f"{sens.loc['detected', 'median_ratio_lo']:.2f}",
        "notierc_groups": str(int(sens.loc["no_tier_c", "groups"])),
        "notierc_auc_groups": str(int(sens.loc["no_tier_c", "groups_auc_above_0.5"])),
        "slope_min": f"{by['slope'].min():.2f}", "slope_max": f"{by['slope'].max():.2f}",
        "recal_level_oe": str(int(recal.loc[(0.1, "intercept"), "groups_oe_in_range"])),
        "recal_level_slope_out": str(int(recal.loc[(0.1, "intercept"), "groups"]
                                         - recal.loc[(0.1, "intercept"),
                                                     "groups_slope_in_range"])),
        "recal_logistic_25": str(int(recal.loc[(0.25, "logistic"), "groups_calibrated"])),
        "recal_logistic_50": str(int(recal.loc[(0.5, "logistic"), "groups_calibrated"])),
    })
    eff = csv("phase6", "condition_effects.csv").set_index(["factor", "level"])
    for key, idx in (("or_dilute", ("medium", "dilute_oligotrophic")),
                     ("or_extract", ("medium", "environment_extract")),
                     ("or_selective", ("medium", "selective")),
                     ("or_warm", ("temperature", "warm")),
                     ("or_days", ("continuous", "incubation time (per doubling)")),
                     ("or_abundance", ("continuous", "source abundance (per 10-fold)"))):
        r = eff.loc[idx]
        v[key] = est(r["odds_ratio"], r["or_lo"], r["or_hi"])
    v["table_hypotheses"] = hypothesis_rows(v, summ, crit)
    return v


def hypothesis_rows(v: dict, summ: pd.DataFrame, crit: dict) -> str:
    """Markdown rows of Table 1 from the computed values."""
    n = summ.loc["primary", "groups"].astype(int)
    rows = [
        ("H1 identity ranks growth",
         f"identity-only held-out AUC > 0.5 in at least {crit['auc_supported_groups']} of {n} "
         "groups and the interval of the mean AUC above 0.5",
         f"{v['h1_groups']} of {n} groups; mean AUC {v['h1_auc']}", "supported"),
        ("H2 calibrated across studies",
         f"calibration slope {crit['calibration_slope_range'][0]}-"
         f"{crit['calibration_slope_range'][1]} and O/E {crit['oe_range'][0]}-"
         f"{crit['oe_range'][1]} in at least {crit['calibrated_groups']} of {n} groups",
         f"{v['h2_groups']} of {n} groups; O/E {v['oe_min']}-{v['oe_max']}", "not supported"),
        ("H3 the shortcut overstates recovery",
         f"shortcut / observed > 1 in at least {crit['proxy_over_groups']} of {n} groups and "
         "the interval of the median ratio above 1",
         f"{v['h3_groups']} of {n} groups; median {v['h3_ratio']}", "supported"),
        ("H4 identity adds to covariates",
         "mean held-out AUC higher with identity than without, interval of the difference "
         "above 0", f"AUC gain {v['h4_gain']}", "supported"),
    ]
    return "\n".join(f"| {a} | {b} | {c} | {d} |" for a, b, c, d in rows)


def phase4_values() -> dict:
    q5 = csv("phase4", "q5_tests.csv").set_index("taxon").loc["family"]
    fam = csv("phase4", "q5_family_medium.csv").set_index(["family", "medium"])
    q6 = csv("phase4", "q6_summary.csv").set_index("analysis").loc["primary"]
    ex = csv("phase4", "exploratory", "curation_by_study.csv").set_index("analysis")
    q7 = csv("phase4", "q7_family.csv").set_index("family")
    q7p = csv("phase4", "q7_phylum.csv").set_index("phylum_full")
    under = q7[q7["call"] == "under"]
    v = {
        "q5_groups": str(int(q5["groups"])), "q5_sd": f"{q5['sd_taxon_medium']:.2f}",
        "q5_p": "< 0.001" if q5["p_value"] < 0.001 else f"= {q5['p_value']:.3f}",
        "q5_cells_sig": str(int(fam["interval_excludes_0"].sum())),
        "q5_cells": str(len(fam)),
        "q5_improved": str(int(q5["groups_improved"])),
        "q5_dll": f"{q5['mean_logloss_diff']:.4f} "
                  f"({q5['mean_logloss_diff_lo']:.4f} to {q5['mean_logloss_diff_hi']:.4f})",
        "q6_or": est(q6["odds_ratio"], q6["or_lo"], q6["or_hi"]),
        "q6_improved": str(int(q6["groups_improved"])),
        "q6_bystudy_or": est(*ex.loc["bystudy", ["odds_ratio", "or_lo", "or_hi"]]),
        "q6_li": est(*ex.loc["li2023", ["odds_ratio", "or_lo", "or_hi"]]),
        "q6_zheng": est(*ex.loc["zheng2024", ["odds_ratio", "or_lo", "or_hi"]]),
        "q7_families": str(len(q7)), "q7_under": str(len(under)),
        "q7_over": str(int((q7["call"] == "over").sum())),
        "q7_zero": str(int((under["observed"] == 0).sum())),
        "q7_phyla": str(len(q7p)), "q7_phyla_under": str(int((q7p["call"] == "under").sum())),
    }
    for family, medium in (("Bacillaceae", "rich_complex"), ("Planococcaceae", "rich_complex")):
        r = fam.loc[(family, medium)]
        v[f"{family.lower()}_rich"] = signed(r["effect"], r["lo"], r["hi"])
    for name in ("Pseudomonadaceae", "Aeromonadaceae", "Rhizobiaceae"):
        r = q7.loc[name]
        v[f"oe_{name.lower()}"] = est(r["oe"], r["oe_lo"], r["oe_hi"])
    for name in ("Acidobacteriota", "Verrucomicrobiota", "Planctomycetota"):
        r = q7p.loc[name]
        v[f"oe_{name.lower()}"] = (f"{r['oe']:.2f}" if r["oe_hi"] == r["oe_lo"] else
                                   est(r["oe"], r["oe_lo"], r["oe_hi"]))
    return v


def habitat_row(r, training: bool) -> str:
    name = HABITAT_NAMES[r["habitat"]] + ("" if training else " *")
    return (f"| {name} | {int(r['studies'])} / {int(r['samples'])} | "
            f"{est(r['shortcut'], r['shortcut_lo'], r['shortcut_hi'])} | "
            f"{fmt_share(r['r1_p10'], 3)}-{fmt_share(r['r1_p90'], 3)} | "
            f"{est(r['d_mean_r1'], r['d_lo_r1'], r['d_hi_r1'])} |")


def phase5_values() -> dict:
    cfg = load_config("phase5")
    amp = csv("phase5", "amplicon", "habitats.csv")
    shot = csv("phase5", "shotgun", "habitats.csv")
    amp_studies = csv("phase5", "amplicon", "studies.csv")
    shot_studies = csv("phase5", "shotgun", "studies.csv")
    tc = csv("phase5", "training_check_criteria.csv")
    pub = csv("phase5", "published_estimates.csv")
    order = [h for h in cfg["habitats"] if h in set(amp["habitat"])]
    amp = amp.set_index("habitat").loc[order].reset_index()
    a = amp.set_index("habitat")
    lo_hab, hi_hab = amp.loc[amp["d_mean_r1"].idxmin()], amp.loc[amp["d_mean_r1"].idxmax()]
    both = pd.concat([amp["shortcut"], shot["shortcut"]])
    shot_text = []
    for r in shot.set_index("habitat").loc[[h for h in cfg["habitats"]
                                             if h in set(shot["habitat"])]].reset_index()\
            .itertuples():
        shot_text.append(f"{HABITAT_NAMES[r.habitat].lower()} ({r.studies} / {r.samples}) "
                         f"shortcut {fmt_share(r.shortcut)}, difference "
                         f"{est(r.d_mean_r1, r.d_lo_r1, r.d_hi_r1)}")
    tcs = tc.set_index(tc["check"].str.split(" \\(").str[0])
    lloyd = pub[pub["source"].str.startswith("Lloyd") & (pub["habitat"] == "soil")]
    martiny = pub[pub["source"].str.startswith("Martiny") & pub["measure"].str.contains(">97%")]
    marine = [f"{fmt_share(a.loc[h, 'family_not_evaluated'])} in {HABITAT_NAMES[h].lower()}"
              for h in ("marine_water", "marine_sediment")]
    return {
        "n_survey_studies": str(len(amp_studies) + len(shot_studies)),
        "n_amp_studies": str(len(amp_studies)), "n_shot_studies": str(len(shot_studies)),
        "n_amp_samples": str(int(amp["samples"].sum())),
        "n_shot_samples": str(int(shot["samples"].sum())),
        "shortcut_min": fmt_pct(both.min(), 0), "shortcut_max": fmt_pct(both.max(), 0),
        "tc_upper": str(int(tcs.loc["upper end", "groups"])),
        "tc_range": str(int(tcs.loc["range describes the spread", "groups"])),
        "tc_range_needed": str(int(tcs.loc["range describes the spread", "needed"])),
        "d_min_habitat": HABITAT_NAMES[lo_hab["habitat"]].lower(),
        "d_min": est(lo_hab["d_mean_r1"], lo_hab["d_lo_r1"], lo_hab["d_hi_r1"]),
        "d_max_habitat": HABITAT_NAMES[hi_hab["habitat"]].lower(),
        "d_max": est(hi_hab["d_mean_r1"], hi_hab["d_lo_r1"], hi_hab["d_hi_r1"]),
        "d8_gut": est(*a.loc["human_gut", ["d_mean_r8", "d_lo_r8", "d_hi_r8"]]),
        "table_habitats_amplicon": "\n".join(
            habitat_row(r, cfg["habitats"][r["habitat"]]["training"])
            for _, r in amp.iterrows()),
        "shotgun_sentence": "; ".join(shot_text) + " (studies / samples; the soil row is one "
                                                   "soil sequenced several times).",
        "comp_soil": fmt_share(a.loc["soil", "family_under"]),
        "comp_shortcut_soil": fmt_share(a.loc["soil", "shortcut_in_under_family"]),
        "comp_gut": fmt_share(a.loc["human_gut", "family_under"]),
        "comp_shortcut_gut": fmt_share(a.loc["human_gut", "shortcut_in_under_family"]),
        "comp_ne_marine": " and ".join(marine),
        "lloyd_soil": f"{lloyd['percent'].iloc[0]:.0f}%",
        "shortcut_soil": fmt_pct(a.loc["soil", "shortcut"], 0),
        "martiny_97": f"{martiny['percent'].iloc[0]:.0f}%",
    }


def reference_values() -> dict:
    ref = pd.read_csv(RES / "phase0" / "reference_contamination.csv").iloc[0]
    return {"gtdb_dropped": fmt_count(ref["gtdb_genomes_dropped"]),
            "gtdb_genomes": fmt_count(ref["gtdb_genomes_with_16s"]),
            "gtdb_dropped_pct": fmt_pct(ref["gtdb_genomes_dropped"]
                                        / ref["gtdb_genomes_with_16s"])}


def study_table() -> str:
    census = pd.read_csv(data_path("processed") / "census.csv").set_index("study_id")
    rows = csv("phase3", "rows_summary.csv")
    rows = rows[rows["variant"] == "primary"].set_index("study_id")
    links = csv("phase2", "qc_asvs_links.csv").set_index("study_id")
    regions = study_region_labels(pd.read_csv(data_path("processed") / "primer_sets.csv"),
                                  census.reset_index())
    out = []
    for s in rows.index:
        c = census.loc[s]
        out.append(f"| [@{c['study_group']}] | {c['habitat']} | "
                   f"{c['accession']} | {regions[s]} | {int(links.loc[s, 'sources_linked'])} | "
                   f"{int(rows.loc[s, 'conditions'])} |")
    return "\n".join(out)


def tool_version(cmd: list[str], pattern: str) -> str:
    res = subprocess.run(cmd, capture_output=True, text=True)
    m = re.search(pattern, res.stdout + res.stderr)
    if not m:
        raise RuntimeError(f"no version in output of {cmd}")
    return m.group(1)


def versions() -> dict:
    import lightgbm
    import matplotlib
    import scipy

    from plate_truth.rfit import r_env, rscript_exe
    r = subprocess.run([rscript_exe(), "-e",
                        "cat(paste(R.version$major, R.version$minor, sep='.'), "
                        "as.character(packageVersion('glmmTMB')))"],
                       capture_output=True, text=True, check=True, env=r_env()).stdout.split()
    return {"ver_python": sys.version.split()[0], "ver_pandas": pd.__version__,
            "ver_numpy": np.__version__, "ver_scipy": scipy.__version__,
            "ver_lightgbm": lightgbm.__version__, "ver_matplotlib": matplotlib.__version__,
            "ver_r": r[0], "ver_glmmtmb": r[1],
            "ver_vsearch": tool_version([vsearch_exe(), "--version"], r"vsearch v(\S+?)[_ ,]"),
            "ver_cutadapt": tool_version([cutadapt_exe(), "--version"], r"(\d+\.\d+\S*)")}


# ---------------------------------------------------------------- references
def crossref(doi: str) -> dict:
    cache = data_path("raw") / "manuscript" / "crossref" / (re.sub(r"[^\w.-]", "_", doi) + ".json")
    if not cache.exists():
        resp = get("https://api.crossref.org/works/" + doi)
        resp.raise_for_status()
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(resp.json()["message"]), encoding="utf-8")
    return json.loads(cache.read_text(encoding="utf-8"))


def reference_dois() -> dict[str, str]:
    dois = dict(load_config("manuscript")["references"])
    for s in load_config("census_studies")["studies"]:
        dois.setdefault(s["study_group"], s["doi"])
    return dois


def build_references(text: str) -> tuple[str, pd.DataFrame]:
    """Replace citations in text; return the text and the reference table (cited works only)."""
    check = pd.read_csv(RES / "reference_check.csv")
    if not check["ok"].all():
        sys.exit("scripts/21_verify_references.py reports failures; fix them first")
    dois = reference_dois()
    keys = cited_keys(text)
    missing = [k for k in keys if k not in dois]
    if missing:
        sys.exit(f"no DOI for citation keys: {', '.join(missing)}")
    records = {k: crossref(dois[k]) for k in keys}
    labels = citation_labels(records)
    text = replace_citations(text, labels)
    table = pd.DataFrame([{"key": k, "doi": dois[k], "label": labels[k],
                           "first_author": surnames(records[k])[0],
                           "year": record_year(records[k]),
                           "title": clean((records[k].get("title") or [""])[0]),
                           "entry": format_reference(records[k], labels[k].split()[-1])}
                          for k in keys])
    table = table.sort_values(["first_author", "label"], key=lambda s: s.str.lower())
    return text, table


def main() -> None:
    root = project_root()
    template = (root / "manuscript" / "preprint_template.md").read_text(encoding="utf-8")
    template = re.sub(r"<!--.*?-->\n*", "", template, flags=re.S)
    values = {"date": datetime.date.today().isoformat()}
    for part in (census_values, phase3_values, phase4_values, phase5_values,
                 reference_values, versions):
        values.update(part())
    values["table_studies"] = study_table()
    values["references"] = "{{references}}"
    text = fill_placeholders(template, values)
    text, table = build_references(text)
    text = text.replace("{{references}}", "\n\n".join(table["entry"]))
    out = root / "manuscript" / "preprint.md"
    out.write_text(text, encoding="utf-8")
    (RES / "phase6").mkdir(exist_ok=True)
    table.drop(columns="entry").to_csv(RES / "phase6" / "manuscript_references.csv",
                                       index=False)
    print(f"wrote {out.relative_to(root)} ({len(text.split())} words, {len(table)} references)")
    pandoc = shutil.which("pandoc") or next(
        iter(sorted(Path(sys.executable).parent.glob("Library/bin/pandoc.exe"))), None)
    if pandoc:
        subprocess.run([pandoc, str(out), "-o", str(out.with_suffix(".docx"))], check=True,
                       cwd=root / "manuscript")
        print("wrote manuscript/preprint.docx")


if __name__ == "__main__":
    main()
