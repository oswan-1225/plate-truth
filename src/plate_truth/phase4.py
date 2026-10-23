"""Phase 4: medium classes by taxon, the curation hypothesis, under-recovered taxa.

ANALYSIS_PLAN.md sections 11-17. Three parts:
- Collection media: the medium class of every culture-collection record (BacDive strain ->
  MediaDive medium), from the medium's MediaDive recipe by the plate-media rule.
- Links from each ASV to the culture-collection strains behind its nearest cultured relative,
  and the share of those strains recorded on each medium class.
- Helpers for the Q5-Q7 analyses (taxon labels, likelihood-ratio test on a variance, held-out
  log loss, the matched-pair statistic, per-taxon observed/expected).
Mixed models are fitted in R (plate_truth.rfit); everything that transforms data or computes a
metric is here, so it is tested in Python.
"""

import re

import numpy as np
import pandas as pd
from scipy.stats import chi2

from plate_truth.media import DILUTE_MAX_G_PER_L, classify_medium
from plate_truth.reference import as_list

GROUP = "study_group"
NO_RECIPE = "no_recipe"


# --- collection media ----------------------------------------------------------------------------

def culture_media_records(bid: int, rec: dict) -> list[dict]:
    """The culture-medium entries of one BacDive record: medium name, growth as recorded, and the
    MediaDive medium id when the entry links to MediaDive (None otherwise)."""
    cg = rec.get("Culture and growth conditions") or {}
    out = []
    for m in as_list(cg.get("culture medium")):
        hit = re.search(r"mediadive\.dsmz\.de/medium/([^/?#\s]+)", m.get("link") or "")
        out.append({"bacdive_id": bid, "medium_name": m.get("name"), "growth": m.get("growth"),
                    "mediadive_id": hit.group(1) if hit else None})
    return out


def has_organic_carbon(formula: str | None) -> bool:
    """True if a sum formula contains carbon that is not inorganic carbon.

    Carbon is an upper-case C not followed by a lower-case letter (so Cl, Ca, Co, Cu are not
    carbon). Carbonate and bicarbonate (CO3), CO and CO2 on their own, and cyanide or
    thiocyanate (CN) count as inorganic. A missing formula is False."""
    if not formula or not isinstance(formula, str):
        return False
    if not re.search(r"C(?![a-z])", formula):
        return False
    if re.search(r"CO3|CN", formula) or re.fullmatch(r"CO2?", formula.strip()):
        return False
    return True


def ingredient_roles(ingredients: pd.DataFrame, rules: dict) -> pd.DataFrame:
    """Per MediaDive ingredient (columns id, name, complex_compound, formula): the facts the
    medium rule needs, from the name lists in config/phase4.yaml (exact, case-insensitive).

    antibacterial, habitat_material: on the respective list.
    complex_source: flagged complex by MediaDive and not on the not_complex list.
    organic: flagged complex by MediaDive (so starch or glycerol count although they are not
    complex sources), or organic carbon by formula, and not on organic_excluded.
    composite_agar: a complex source whose name contains 'agar' (a complete dehydrated agar
    medium, whose weight includes its agar)."""
    def on(names):
        lower = {n.lower() for n in names}
        return ingredients["name"].str.lower().isin(lower)

    out = ingredients[["id", "name"]].copy()
    out["antibacterial"] = on(rules["antibacterial"])
    out["habitat_material"] = on(rules["habitat_material"])
    flagged = ingredients["complex_compound"].fillna(0).astype(int) == 1
    out["complex_source"] = flagged & ~on(rules["not_complex"])
    carbon = ingredients["formula"].map(has_organic_carbon)
    out["organic"] = (flagged | carbon) & ~on(rules["organic_excluded"])
    out["composite_agar"] = out["complex_source"] & out["name"].str.contains(
        "agar", case=False)
    return out


def medium_facts(composition: pd.DataFrame, roles: pd.DataFrame, agar_g_per_l: float
                 ) -> pd.DataFrame:
    """Recipe facts per medium from its flattened composition.

    composition: medium_id, ingredient_id, g_l, optional (one row per ingredient).
    Optional ingredients are ignored. organic_known_g: summed g/L of organic ingredients with a
    stated amount (a composite agar medium minus agar_g_per_l, not below 0); organic_unknown:
    some organic ingredient has no g/L (e.g. given in mL)."""
    c = composition[composition["optional"].fillna(0).astype(int) == 0]
    c = c.merge(roles.rename(columns={"id": "ingredient_id"}), on="ingredient_id", how="left")
    for col in ["antibacterial", "habitat_material", "complex_source", "organic",
                "composite_agar"]:
        c[col] = c[col].fillna(False).astype(bool)
    g = c["g_l"].astype(float)
    g = g.where(~c["composite_agar"], (g - agar_g_per_l).clip(lower=0))
    c["organic_g"] = np.where(c["organic"], g, 0.0)
    c["organic_missing"] = c["organic"] & g.isna()
    agg = c.groupby("medium_id").agg(
        antibacterial=("antibacterial", "any"), habitat_material=("habitat_material", "any"),
        complex_source=("complex_source", "any"),
        organic_known_g=("organic_g", lambda s: float(np.nansum(s))),
        organic_unknown=("organic_missing", "any"), n_ingredients=("ingredient_id", "size"))
    return agg.reset_index()


def classify_collection_media(facts: pd.DataFrame) -> pd.DataFrame:
    """medium_class and rule per medium by the plate-media rule (classify_medium).

    When some organic amount is unknown, the known part decides only if it already exceeds the
    dilute cut-off (the medium is rich whatever the rest is); otherwise the class is
    'unresolved'."""
    out = facts.copy()
    res = []
    for r in facts.itertuples(index=False):
        if r.organic_unknown and r.organic_known_g <= DILUTE_MAX_G_PER_L:
            organic = None
        else:
            organic = r.organic_known_g
        res.append(classify_medium(r.antibacterial, r.habitat_material, r.complex_source,
                                   organic))
    out["medium_class"] = [c for c, _ in res]
    out["class_rule"] = [r for _, r in res]
    return out


def strain_classes(records: pd.DataFrame, media_classes: pd.DataFrame, classes: list[str]
                   ) -> pd.DataFrame:
    """Distinct (bacdive_id, medium_class) for growth-positive records whose MediaDive medium
    falls in one of `classes`. Growth counts as positive unless recorded as 'no' (the same
    rule as the reference build)."""
    pos = records[records["growth"].astype(str).str.lower() != "no"].dropna(
        subset=["mediadive_id"])
    m = pos.merge(media_classes[["medium_id", "medium_class"]].astype({"medium_id": str}),
                  left_on="mediadive_id", right_on="medium_id")
    m = m[m["medium_class"].isin(classes)]
    return m[["bacdive_id", "medium_class"]].drop_duplicates().reset_index(drop=True)


# --- ASV -> culture-collection strains -----------------------------------------------------------

def genome_base(accession: str) -> str:
    """GTDB or NCBI assembly accession -> 'GCA_<digits>' without version or RS_/GB_ prefix, so
    GenBank and RefSeq copies of one assembly match."""
    a = re.sub(r"^(RS_|GB_)", "", accession)
    a = re.sub(r"^GCF_", "GCA_", a)
    return re.sub(r"\.\d+$", "", a)


def member_strains(members: pd.DataFrame, genomes: pd.DataFrame) -> pd.DataFrame:
    """(reference, ref_id, bacdive_id) for reference members that are BacDive strains.

    members: reference, ref_id, member_id ('bacdive_strain|<id>', 'gtdb_isolate|<acc>|<i>',
    'refseq_type|<acc>'). genomes: bacdive_id, assembly (from BacDive records). A BacDive strain
    16S links to that strain; a GTDB isolate genome links to every BacDive strain listing that
    assembly; RefSeq type-material 16S records carry no strain identifier and are not linked."""
    src = members["member_id"].str.split("|", n=1).str[0]
    direct = members[src == "bacdive_strain"].assign(
        bacdive_id=lambda d: d["member_id"].str.split("|").str[1].astype(int))
    gtdb = members[src == "gtdb_isolate"].assign(
        base=lambda d: d["member_id"].str.split("|").str[1].map(genome_base))
    gmap = genomes.assign(base=genomes["assembly"].map(genome_base))[["base", "bacdive_id"]]
    via = gtdb.merge(gmap.drop_duplicates(), on="base")
    cols = ["reference", "ref_id", "bacdive_id"]
    return pd.concat([direct[cols], via[cols]]).drop_duplicates().reset_index(drop=True)


def asv_strains(nr: pd.DataFrame, links: pd.DataFrame) -> pd.DataFrame:
    """(asv_id, bacdive_id) through every tied best hit. nr: asv_id, reference, targets
    (';'-joined ref_ids; missing = no hit)."""
    t = nr.dropna(subset=["targets"]).assign(ref_id=lambda d: d["targets"].str.split(";"))
    t = t.explode("ref_id")[["asv_id", "reference", "ref_id"]]
    return t.merge(links, on=["reference", "ref_id"])[["asv_id", "bacdive_id"]] \
        .drop_duplicates().reset_index(drop=True)


def recorded_shares(asv_links: pd.DataFrame, classes_by_strain: pd.DataFrame,
                    classes: list[str]) -> pd.DataFrame:
    """Per ASV: n_strains (linked strains with at least one classed growth record), and for
    each class X, rec_X = share of those strains recorded on X and any_X = 1 if any is.
    ASVs with no such strain are absent."""
    s = asv_links.merge(classes_by_strain, on="bacdive_id")
    n = s.groupby("asv_id")["bacdive_id"].nunique().rename("n_strains")
    k = s.groupby(["asv_id", "medium_class"])["bacdive_id"].nunique().unstack(fill_value=0)
    k = k.reindex(columns=classes, fill_value=0)
    out = pd.DataFrame({"n_strains": n})
    for c in classes:
        out[f"rec_{c}"] = k[c].reindex(out.index, fill_value=0) / out["n_strains"]
        out[f"any_{c}"] = (out[f"rec_{c}"] > 0).astype(float)
    return out.rename_axis("asv_id").reset_index()


def match_term(rows: pd.DataFrame, classes: list[str], prefix: str = "rec") -> pd.Series:
    """The recorded share (or indicator) for the class of each row's own plating medium; 0 when
    the row's class has no recorded-class column."""
    out = pd.Series(0.0, index=rows.index)
    for c in classes:
        sel = rows["medium"] == c
        out[sel] = rows.loc[sel, f"{prefix}_{c}"]
    return out


# --- Q5 ------------------------------------------------------------------------------------------

def taxon_label(tax: pd.DataFrame, rank: str) -> pd.Series:
    """Taxon at `rank` (e.g. 'family'); unassigned ASVs get 'unassigned <rank> (<phylum>)', or
    'unassigned' when the phylum is unassigned too."""
    phylum = tax["phylum"].fillna("")
    fallback = np.where(phylum != "", "unassigned " + rank + " (" + phylum + ")", "unassigned")
    return tax[rank].where(tax[rank].notna() & (tax[rank] != ""), fallback)


def combine(*cols: pd.Series) -> pd.Series:
    """Key for an interaction random effect: the values joined with ':'."""
    out = cols[0].astype(str)
    for c in cols[1:]:
        out = out + ":" + c.astype(str)
    return out


def boundary_lrt_p(loglik_null: float, loglik_full: float) -> float:
    """p-value of a likelihood-ratio test for one variance component, which sits on the boundary
    (0) under the null: the statistic follows a 50:50 mixture of chi-square(0) and
    chi-square(1), so p = 0.5 * P(chi2_1 > stat), and 1 if stat <= 0."""
    stat = 2.0 * (loglik_full - loglik_null)
    return 1.0 if stat <= 0 else float(0.5 * chi2.sf(stat, 1))


def mean_log_loss(y: np.ndarray, p: np.ndarray) -> float:
    """Mean negative Bernoulli log-likelihood (natural log); lower is better."""
    y = np.asarray(y, dtype=float)
    p = np.clip(np.asarray(p, dtype=float), 1e-12, 1 - 1e-12)
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log1p(-p)))


def reportable_taxa(rows: pd.DataFrame, taxon: str, min_rows: int, min_classes: int
                    ) -> list[str]:
    """Taxa with at least min_rows rows in total and rows in at least min_classes medium
    classes (row counts only; outcomes are not used)."""
    g = rows.groupby(taxon).agg(n=("y", "size"), k=("medium", "nunique"))
    return sorted(g[(g["n"] >= min_rows) & (g["k"] >= min_classes)].index)


# --- Q6 matched pairs ----------------------------------------------------------------------------

PAIR_KEYS = ["study_id", "source_key", "oxygen", "temperature", "incubation_days"]


def matched_condition_pairs(conds: pd.DataFrame) -> pd.DataFrame:
    """Pairs of plating conditions from the same study and source sample, with the same oxygen,
    temperature band and incubation time, but different medium classes. Oriented so that
    class_x < class_y alphabetically. conds: condition_id, medium and the PAIR_KEYS."""
    a = conds[["condition_id", "medium", *PAIR_KEYS]].drop_duplicates()
    p = a.merge(a, on=PAIR_KEYS, suffixes=("_x", "_y"))
    p = p[p["medium_x"] < p["medium_y"]]
    p = p.rename(columns={"medium_x": "class_x", "medium_y": "class_y"})
    return p[["study_id", "source_key", "condition_id_x", "condition_id_y", "class_x",
              "class_y"]].reset_index(drop=True)


def discordant_pairs(rows: pd.DataFrame, pairs: pd.DataFrame, recorded: pd.DataFrame
                     ) -> pd.DataFrame:
    """Per matched condition pair and ASV: keep ASVs enriched on exactly one of the two plates
    whose linked strains are recorded on exactly one of the two classes (any_<class>).
    Returns study_id, group, asv_id, stratum ('class_x|class_y') and the cell: 'a' recorded x,
    grew x; 'b' recorded x, grew y; 'c' recorded y, grew x; 'd' recorded y, grew y."""
    r = rows[["condition_id", "asv_id", GROUP, "y"]]
    m = pairs.merge(r.rename(columns={"condition_id": "condition_id_x", "y": "y_x"}),
                    on="condition_id_x")
    m = m.merge(r[["condition_id", "asv_id", "y"]].rename(
        columns={"condition_id": "condition_id_y", "y": "y_y"}),
        on=["condition_id_y", "asv_id"])
    m = m[m["y_x"] != m["y_y"]].merge(recorded, on="asv_id")
    rec_x = np.array([m.at[i, f"any_{c}"] for i, c in zip(m.index, m["class_x"])]) > 0
    rec_y = np.array([m.at[i, f"any_{c}"] for i, c in zip(m.index, m["class_y"])]) > 0
    keep = rec_x != rec_y
    m, rec_x = m[keep], rec_x[keep]
    grew_x = (m["y_x"] == 1).to_numpy()
    cell = np.where(rec_x, np.where(grew_x, "a", "b"), np.where(grew_x, "c", "d"))
    return pd.DataFrame({"study_id": m["study_id"].to_numpy(), GROUP: m[GROUP].to_numpy(),
                         "asv_id": m["asv_id"].to_numpy(),
                         "stratum": (m["class_x"] + "|" + m["class_y"]).to_numpy(),
                         "cell": cell})


def mantel_haenszel_or(counts: pd.DataFrame) -> float:
    """Mantel-Haenszel odds ratio over strata. counts: one row per stratum with columns a, b,
    c, d. OR = sum(a d / n) / sum(b c / n), n = a + b + c + d per stratum."""
    n = counts[["a", "b", "c", "d"]].sum(axis=1).replace(0, np.nan)
    num = (counts["a"] * counts["d"] / n).sum()
    den = (counts["b"] * counts["c"] / n).sum()
    return float(num / den) if den > 0 else np.nan


def pair_counts(disc: pd.DataFrame, weights: pd.Series | None = None) -> pd.DataFrame:
    """2 x 2 cell counts per stratum (study, class pair), optionally weighting each ASV (for
    the bootstrap)."""
    d = disc.assign(w=1.0 if weights is None else disc["asv_id"].map(weights).fillna(0.0))
    t = d.pivot_table(index=["study_id", "stratum"], columns="cell", values="w",
                      aggfunc="sum", fill_value=0.0)
    return t.reindex(columns=["a", "b", "c", "d"], fill_value=0.0).reset_index()


def bootstrap_pair_or(disc: pd.DataFrame, n_boot: int, seed: int) -> np.ndarray:
    """Bootstrap distribution of the Mantel-Haenszel OR, resampling ASVs with replacement
    within each study group (an ASV's weight is how often it was drawn)."""
    rng = np.random.default_rng(seed)
    asvs = disc[[GROUP, "asv_id"]].drop_duplicates()
    out = np.empty(n_boot)
    for b in range(n_boot):
        w = []
        for _, s in asvs.groupby(GROUP):
            draw = rng.integers(0, len(s), len(s))
            w.append(pd.Series(np.bincount(draw, minlength=len(s)), index=s["asv_id"].values))
        out[b] = mantel_haenszel_or(pair_counts(disc, pd.concat(w)))
    return out


# --- Q7 ------------------------------------------------------------------------------------------

def taxon_oe(df: pd.DataFrame, taxon: str, pred: str) -> pd.DataFrame:
    """Observed and expected events per (taxon, group): O = sum(y), E = sum(pred)."""
    return df.groupby([taxon, GROUP]).agg(observed=("y", "sum"), expected=(pred, "sum"),
                                          rows=("y", "size")).reset_index()


def classify_taxa(per_group: pd.DataFrame, taxon: str, min_expected: float, min_groups: int,
                  consistency: float, n_boot: int, seed: int) -> pd.DataFrame:
    """Pooled O/E per taxon with a bootstrap interval over study groups, and the pre-set call.

    A taxon is evaluated if its expected events reach min_expected in at least min_groups
    groups (only those groups count). 'under' if pooled O/E < 1, the interval's upper bound
    < 1 and O/E < 1 in at least `consistency` of its groups; 'over' is the mirror image;
    otherwise 'neither'."""
    rng = np.random.default_rng(seed)
    rows = []
    for t, s in per_group.groupby(taxon):
        s = s[s["expected"] >= min_expected]
        if len(s) < min_groups:
            continue
        o, e = s["observed"].to_numpy(float), s["expected"].to_numpy(float)
        idx = rng.integers(0, len(s), size=(n_boot, len(s)))
        boots = o[idx].sum(axis=1) / e[idx].sum(axis=1)
        lo, hi = np.quantile(boots, [0.025, 0.975])
        ratio = o / e
        est = o.sum() / e.sum()
        below, above = np.mean(ratio < 1), np.mean(ratio > 1)
        call = "neither"
        if est < 1 and hi < 1 and below >= consistency:
            call = "under"
        elif est > 1 and lo > 1 and above >= consistency:
            call = "over"
        rows.append({taxon: t, "groups": len(s), "observed": o.sum(), "expected": e.sum(),
                     "oe": est, "oe_lo": lo, "oe_hi": hi, "share_groups_below_1": below,
                     "share_groups_above_1": above, "call": call})
    return pd.DataFrame(rows)
