"""Phase 5: what the calibration model can say about public surveys (ANALYSIS_PLAN.md s19-25).

Three parts:
- Survey selection: habitat membership of a MGnify biome lineage, excluded words, the seeded
  screening order, run choice, and the primer screen of streamed read pairs.
- Applying the Phase 3 model to survey sequences: the plating-condition mix of the training
  studies, the predicted recovery share for one plate and for any of K plates at a given study
  effect, and the shortcut share.
- Aggregation: habitat replicates (studies resampled x fixed-effect draws), the shortcut verdict
  and composition shares.
Downloads, VSEARCH and R calls live in the numbered scripts; everything here is tested.
"""

import re
import subprocess

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm

from plate_truth.rfit import r_env, rscript_exe

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT", "S": "CG",
         "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
         "N": "ACGT"}
RANKS = ["domain", "phylum", "class", "order", "family", "genus"]


# --- survey selection --------------------------------------------------------------------------

def habitat_of(lineage: str, habitats: dict) -> str | None:
    """The habitat whose lineages prefix this MGnify biome lineage and whose exclude words it
    does not contain (case-insensitive); None if no habitat matches."""
    low = (lineage or "").lower()
    for name, h in habitats.items():
        if any(low == p.lower() or low.startswith(p.lower() + ":") for p in h["lineages"]) \
                and not any(w.lower() in low for w in h.get("exclude", [])):
            return name
    return None


def has_word(text: str, words: list[str]) -> bool:
    """True if any of the words occurs in text as a whole word (case-insensitive)."""
    t = str(text or "")
    return any(re.search(rf"(?<![A-Za-z]){re.escape(w)}(?![A-Za-z])", t, re.IGNORECASE)
               for w in words)


def screening_order(accessions: list[str], seed: int) -> list[str]:
    """Candidates sorted by accession, then shuffled with a fixed seed (reproducible order)."""
    acc = sorted(set(accessions))
    rng = np.random.default_rng(seed)
    return [acc[i] for i in rng.permutation(len(acc))]


def choose_runs(runs: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    """One run per biosample (the first by run accession), then up to n at random (fixed seed),
    returned sorted by run accession."""
    one = runs.sort_values("run_accession").drop_duplicates("sample_accession")
    if len(one) > n:
        one = one.sample(n, random_state=seed)
    return one.sort_values("run_accession").reset_index(drop=True)


def primer_mismatches(seq: str, primer: str) -> int:
    """Mismatches between the start of seq and an IUPAC primer (a base outside the primer's
    code is a mismatch; missing bases count as mismatches)."""
    mm = max(0, len(primer) - len(seq))
    for s, p in zip(seq, primer):
        if s.upper() not in IUPAC.get(p.upper(), p.upper()):
            mm += 1
    return mm


def primer_at_start(read: str, primers: list[str], max_offset: int, max_mm: int) -> bool:
    """True if any primer matches the read starting within the first max_offset bases."""
    return any(primer_mismatches(read[o:o + len(p)], p) <= max_mm
               for p in primers for o in range(max_offset + 1))


def pair_has_primers(r1: str, r2: str, fwd: list[str], rev: list[str], max_offset: int,
                     max_mm: int) -> bool:
    """A read pair carries the primer pair if one read starts with a forward primer and the
    other with a reverse primer, in either orientation."""
    def at(read, ps):
        return primer_at_start(read, ps, max_offset, max_mm)
    return (at(r1, fwd) and at(r2, rev)) or (at(r1, rev) and at(r2, fwd))


def primer_pair_shares(pairs: list[tuple[str, str]], references: dict, max_offset: int,
                       max_mm: int) -> dict:
    """Share of read pairs carrying each reference's primer pair."""
    out = {}
    for name, ref in references.items():
        hits = sum(pair_has_primers(a, b, ref["fwd"], ref["rev"], max_offset, max_mm)
                   for a, b in pairs)
        out[name] = hits / len(pairs) if pairs else 0.0
    return out


def best_reference(shares: dict, min_share: float) -> str | None:
    """The reference with the largest share if it reaches min_share (ties: first listed)."""
    if not shares:
        return None
    name = max(shares, key=lambda k: shares[k])
    return name if shares[name] >= min_share else None


def merge_possible(len1: int, len2: int, amplicon_len: int, min_overlap: int) -> bool:
    """Paired reads can span the amplicon (primers included) with at least min_overlap bases."""
    return len1 + len2 - amplicon_len >= min_overlap


# --- the plating-condition mix -----------------------------------------------------------------

def condition_weights(groups: pd.Series) -> pd.Series:
    """Weight of each condition so that every study group carries the same total (summing to
    1): 1 / (number of groups x conditions in the condition's group)."""
    n = groups.map(groups.value_counts())
    return 1.0 / (groups.nunique() * n)


def mixture_points(g: np.ndarray, w: np.ndarray, sigma: float, n_points: int) -> np.ndarray:
    """Equally weighted points representing the distribution sum_k w_k Normal(g_k, sigma^2):
    its quantiles at (j - 0.5) / n_points, j = 1..n_points (by bisection on the CDF)."""
    g, w = np.asarray(g, float), np.asarray(w, float) / np.sum(w)
    t = (np.arange(n_points) + 0.5) / n_points
    if sigma == 0:
        order = np.argsort(g)
        cum = np.cumsum(w[order])
        return g[order][np.minimum(np.searchsorted(cum, t), len(g) - 1)]
    lo = np.full(n_points, g.min() - 10 * sigma)
    hi = np.full(n_points, g.max() + 10 * sigma)
    for _ in range(80):
        mid = (lo + hi) / 2
        cdf = norm.cdf((mid[:, None] - g[None, :]) / sigma) @ w
        below = cdf < t
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    return (lo + hi) / 2


def gh_nodes(sigma: float, n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Hermite nodes and weights for the expectation over Normal(0, sigma^2)."""
    if sigma == 0:
        return np.zeros(1), np.ones(1)
    x, w = np.polynomial.hermite_e.hermegauss(n_nodes)
    return sigma * x, w / w.sum()


def taxon_probabilities(eta: np.ndarray, points: np.ndarray, u: float, sigma_source: float,
                        k: int, n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Per-sequence predicted probabilities at study effect u (ANALYSIS_PLAN.md s22): p1, of
    being enriched on one plate drawn from the condition mix, and pK, on any of k plates
    (independent given the shared study and source effects).

    eta: linear predictor of each sequence without condition terms; points: the condition mix
    plus condition effect (mixture_points)."""
    eta = np.asarray(eta, float)
    s, ws = gh_nodes(sigma_source, n_nodes)
    pbar = expit(eta[:, None, None] + u + s[None, :, None] + points[None, None, :]).mean(axis=2)
    return pbar @ ws, 1.0 - ((1.0 - pbar) ** k) @ ws


def recovery_shares(eta: np.ndarray, a: np.ndarray, points: np.ndarray, u: float,
                    sigma_source: float, k: int, n_nodes: int) -> tuple[float, float]:
    """Predicted recovery shares of one sample at study effect u: (R1, RK), the a-weighted
    share (a = relative abundance) of taxa enriched on one plate and on any of k plates
    (taxon_probabilities)."""
    a = np.asarray(a, float)
    p1, pk = taxon_probabilities(eta, points, u, sigma_source, k, n_nodes)
    return float(a @ p1 / a.sum()), float(a @ pk / a.sum())


def study_effect(q: float, sd_study: float) -> float:
    """The study random effect at quantile q of Normal(0, sd_study^2)."""
    return float(norm.ppf(q) * sd_study)


def shortcut_share(identity: np.ndarray, a: np.ndarray, threshold: float) -> float:
    """a-weighted share of sequences with a cultured relative at or above threshold (no hit =
    NaN counts as below)."""
    identity, a = np.asarray(identity, float), np.asarray(a, float)
    return float(a @ np.nan_to_num(identity >= threshold) / a.sum())


# --- identity spline basis (R, so it is exactly the basis of the fit) --------------------------

def spline_basis(values: np.ndarray, cfg: dict) -> np.ndarray:
    """R splines::ns basis of identity with the knots of config/phase3.yaml (5 columns), the
    basis the model was fitted with, evaluated at values."""
    s = cfg["identity_spline"]
    vals = "\n".join(f"{v:.6f}" for v in np.asarray(values, float))
    knots = ",".join(f"{k:g}" for k in s["interior_knots"])
    lo, hi = s["boundary_knots"]
    code = (f"library(splines); x <- scan(file('stdin'), quiet = TRUE); "
            f"b <- ns(x, knots = c({knots}), Boundary.knots = c({lo:g}, {hi:g})); "
            f"write.table(format(b, digits = 17), stdout(), sep = ',', row.names = FALSE, "
            f"col.names = FALSE, quote = FALSE)")
    res = subprocess.run([rscript_exe(), "-e", code], env=r_env(), capture_output=True,
                         text=True, check=True, input=vals + "\n")
    out = np.loadtxt(res.stdout.strip().splitlines(), delimiter=",", ndmin=2)
    return out


# --- shotgun fragments ---------------------------------------------------------------------------

def pseudo_taxon(tax: pd.DataFrame, rank: str) -> pd.Series:
    """Lineage key down to rank, stopping at the first unassigned rank (SINTAX columns)."""
    keep = RANKS[:RANKS.index(rank) + 1]

    def key(row):
        parts = []
        for r in keep:
            v = row.get(r)
            if not isinstance(v, str) or not v:
                break
            parts.append(f"{r[0]}:{v}")
        return ";".join(parts) or "unassigned"
    return tax.apply(key, axis=1)


def group_share(keys: pd.Series) -> pd.Series:
    """For each element, the share of all elements with the same key."""
    return keys.map(keys.value_counts()) / len(keys)


# --- aggregation and verdicts --------------------------------------------------------------------

def replicate_means(values: np.ndarray, n_rep: int, seed: int) -> np.ndarray:
    """Habitat replicates: values is (studies, draws). Replicate r resamples the studies with
    replacement and takes draw r mod draws; returns the mean over the resampled studies. With
    one study only the draws vary."""
    values = np.atleast_2d(np.asarray(values, float))
    n_st, n_dr = values.shape
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n_st, size=(n_rep, n_st)) if n_st > 1 else np.zeros((n_rep, 1), int)
    d = np.arange(n_rep) % n_dr
    return values[idx, d[:, None]].mean(axis=1)


def required_share(n: int, num: int, den: int) -> int:
    """Smallest count reaching num/den of n (integer arithmetic, rounded up)."""
    return -(-n * num // den)


def shortcut_verdict(d_by_study: np.ndarray, d_replicates: np.ndarray, num: int, den: int
                     ) -> dict:
    """The shortcut test of one habitat: the mean difference (shortcut share minus the upper
    end of the range) over studies, its 95% interval from the replicates, the number of
    studies with a positive difference, and the verdict (lower bound > 0 and enough studies)."""
    d = np.asarray(d_by_study, float)
    lo, hi = np.quantile(d_replicates, [0.025, 0.975])
    need = required_share(len(d), num, den)
    positive = int((d > 0).sum())
    return {"studies": len(d), "d_mean": float(d.mean()), "d_lo": float(lo), "d_hi": float(hi),
            "studies_positive": positive, "studies_needed": need,
            "shortcut_exceeds": bool(lo > 0 and positive >= need)}


def composition_shares(a: np.ndarray, labels: pd.Series) -> pd.Series:
    """a-weighted share of a sample in each label (e.g. Q7 call of the sequence's family)."""
    a = pd.Series(np.asarray(a, float), index=labels.index)
    return a.groupby(labels).sum() / a.sum()


# --- the fitted model's parts ------------------------------------------------------------------

def spline_terms(names) -> list[str]:
    """Coefficient names of the identity spline basis, in order."""
    return [n for n in names if str(n).startswith("ns(identity")]


def condition_columns(columns) -> list[str]:
    """Design columns that describe a plating condition (depth, days, medium, temperature,
    oxygen), as opposed to the sequence (identity, no-hit, abundance, host)."""
    return [c for c in columns if c not in ("identity", "no_hit", "log_env", "host")]


def sequence_eta(basis: np.ndarray, no_hit: np.ndarray, log_env: np.ndarray,
                 phylum_effect: np.ndarray, beta: pd.Series, env_center: float) -> np.ndarray:
    """Linear predictor of each sequence without the condition terms: intercept + identity
    spline + no-hit + centred log source abundance + phylum effect. beta may be a Series of
    coefficients or a DataFrame of draws (one column per draw) indexed by term."""
    b = beta if isinstance(beta, pd.DataFrame) else beta.to_frame()
    sp = spline_terms(b.index)
    eta = (np.asarray(basis) @ b.loc[sp].to_numpy()
           + np.outer(no_hit, b.loc["no_hit"].to_numpy())
           + np.outer(np.asarray(log_env) - env_center, b.loc["log_env"].to_numpy())
           + b.loc["(Intercept)"].to_numpy()[None, :]
           + np.asarray(phylum_effect)[:, None])
    return eta[:, 0] if not isinstance(beta, pd.DataFrame) else eta


def condition_effects(x_cond: pd.DataFrame, beta: pd.Series | pd.DataFrame) -> np.ndarray:
    """Contribution of each plating condition (rows of x_cond) to the linear predictor; one
    column per draw when beta is a DataFrame of draws."""
    b = beta if isinstance(beta, pd.DataFrame) else beta.to_frame()
    cols = list(x_cond.columns)
    g = x_cond.to_numpy(float) @ b.loc[cols].to_numpy()
    return g[:, 0] if not isinstance(beta, pd.DataFrame) else g


def window_span(caln: str, query_len: int) -> tuple[int, int, int]:
    """Part of an (oriented) query that lies in the target window, from a VSEARCH caln string
    (M: both, D: query base absent from the target, I: target base absent from the query).
    Returns (start, end) of the query slice [start:end] and the number of target (window)
    bases spanned. Terminal D/I runs are overhangs and are excluded."""
    ops = [(int(n) if n else 1, op) for n, op in re.findall(r"(\d*)([MDI])", caln)]
    if not any(op == "M" for _, op in ops):
        return 0, 0, 0
    first = next(i for i, (_, op) in enumerate(ops) if op == "M")
    last = max(i for i, (_, op) in enumerate(ops) if op == "M")
    start = sum(n for n, op in ops[:first] if op == "D")
    trail = sum(n for n, op in ops[last + 1:] if op == "D")
    covered = sum(n for n, op in ops[first:last + 1] if op in "MI")
    return start, query_len - trail, covered


# --- survey samples -------------------------------------------------------------------------------

def is_prokaryote(tax: pd.DataFrame) -> pd.Series:
    """Bacteria or Archaea, not chloroplast or mitochondria (the training rule)."""
    organelle = tax["order"].eq("Chloroplast") | tax["family"].eq("Mitochondria")
    return tax["domain"].isin(["Bacteria", "Archaea"]) & ~organelle


def amplicon_samples(counts: pd.DataFrame, identity: pd.Series, tax: pd.DataFrame,
                     present_min: int, min_reads: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per run (sample): the present prokaryotic ASVs (>= present_min reads) with relative
    abundance among ALL the run's prokaryotic reads, identity and taxonomy; and a depth table.
    Runs with fewer than min_reads prokaryotic reads are dropped (listed in the depth table)."""
    prok = is_prokaryote(tax.reindex(counts["asv_id"].unique()).fillna(""))
    c = counts[counts["asv_id"].map(prok).fillna(False).astype(bool)]
    depth = c.groupby("run")["count"].sum().rename("prokaryotic_reads").to_frame()
    depth["used"] = depth["prokaryotic_reads"] >= min_reads
    c = c[c["run"].map(depth["used"])].copy()
    c["a"] = c["count"] / c["run"].map(depth["prokaryotic_reads"])
    c = c[c["count"] >= present_min]
    c["identity"] = c["asv_id"].map(identity)
    return c.join(tax, on="asv_id").reset_index(drop=True), depth.reset_index()


def shotgun_samples(frags: pd.DataFrame, rank: str, min_fragments: int
                    ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per run: prokaryotic window fragments, each with weight a = 1 / fragments in the run
    and a source abundance = share of the run's fragments in the same lineage down to rank."""
    f = frags[is_prokaryote(frags.fillna(""))].copy()
    depth = f.groupby("run").size().rename("prokaryotic_reads").to_frame()
    depth["used"] = depth["prokaryotic_reads"] >= min_fragments
    f = f[f["run"].map(depth["used"])].copy()
    f["taxon_key"] = pseudo_taxon(f, rank)
    f["a"] = 1.0 / f["run"].map(depth["prokaryotic_reads"])
    f["abundance"] = f.groupby("run")["taxon_key"].transform(lambda k: group_share(k))
    return f.reset_index(drop=True), depth.reset_index()


def run_window_qc(qc: pd.DataFrame, min_share: float,
                  read_pairs: pd.Series | None = None) -> pd.DataFrame:
    """Per run: the share of merged reads that lie between the study's primers
    (reads_window / reads_in), flagged when below min_share (the census primer screen's
    threshold). A flagged run was probably amplified with another primer pair; the census
    screens only the first runs of a study, so this catches studies that mix regions.
    With read_pairs (ENA read pairs per run, indexed by run), also the share of pairs that
    merged (reads_in / read_pairs): low values mean short overlaps or poor read ends."""
    out = qc[["run", "status", "reads_in", "reads_window", "reads_mapped"]].copy()
    if read_pairs is not None:
        out["merged_share"] = out["reads_in"] / out["run"].map(read_pairs)
    out["window_share"] = out["reads_window"] / out["reads_in"]
    out["low_window_share"] = out["window_share"] < min_share
    return out
