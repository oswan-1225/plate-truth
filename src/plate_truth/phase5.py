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


def recovery_shares(eta: np.ndarray, a: np.ndarray, points: np.ndarray, u: float,
                    sigma_source: float, k: int, n_nodes: int) -> tuple[float, float]:
    """Predicted recovery shares of one sample at study effect u (ANALYSIS_PLAN.md s22).

    eta: linear predictor of each sequence without condition terms; a: its weight (relative
    abundance); points: the condition mix plus condition effect (mixture_points). Returns
    (R1, RK): the a-weighted share of taxa enriched on one plate, and on any of k plates drawn
    from the mix (independent given the shared study and source effects)."""
    eta, a = np.asarray(eta, float), np.asarray(a, float)
    s, ws = gh_nodes(sigma_source, n_nodes)
    pbar = expit(eta[:, None, None] + u + s[None, :, None] + points[None, None, :]).mean(axis=2)
    p1 = pbar @ ws
    pk = 1.0 - ((1.0 - pbar) ** k) @ ws
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
    vals = ",".join(f"{v:.6f}" for v in np.asarray(values, float))
    knots = ",".join(f"{k:g}" for k in s["interior_knots"])
    lo, hi = s["boundary_knots"]
    code = (f"library(splines); x <- c({vals}); "
            f"b <- ns(x, knots = c({knots}), Boundary.knots = c({lo:g}, {hi:g})); "
            f"write.table(b, stdout(), sep = ',', row.names = FALSE, col.names = FALSE)")
    res = subprocess.run([rscript_exe(), "-e", code], env=r_env(), capture_output=True,
                         text=True, check=True)
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
