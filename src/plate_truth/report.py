"""Shared helpers for the Phase 6 outputs: study labels for figures and text, the habitat of a
study group, presence of each condition level per study, and number formatting for the
manuscript (every number in the text is formatted here from a value read from results/).
"""

import re

import numpy as np
import pandas as pd


def study_label(group: str) -> str:
    """'indraningrat2019' -> 'Indraningrat 2019' (author and year from a study-group id)."""
    m = re.fullmatch(r"([a-z]+)(\d{4})", group)
    if not m:
        raise ValueError(f"not an author+year study id: {group!r}")
    return f"{m.group(1).capitalize()} {m.group(2)}"


def group_habitats(census: pd.DataFrame) -> dict[str, str]:
    """Study group -> its habitat(s), in census order, joined with ' + ' when a group spans
    several (one paper sampling two habitats is one study group)."""
    out: dict[str, list[str]] = {}
    for g, h in zip(census["study_group"], census["habitat"], strict=True):
        out.setdefault(g, [])
        if h not in out[g]:
            out[g].append(h)
    return {g: " + ".join(hs) for g, hs in out.items()}


def region_label(region: str, primers: str, ref_name: str, full_length: bool) -> str:
    """How a study's compared window is described: 'V4, 515F/806R'. A primer set whose
    reference name ends in '_2' is a modified version of a same-named pair (more degenerate
    bases or a shorter primer) and is marked; full-length reads were cut in silico."""
    label = f"{region}, {primers}" + (" (modified)" if ref_name.endswith("_2") else "")
    return f"full length, cut to {label}" if full_length else label


def study_region_labels(primer_sets: pd.DataFrame, census: pd.DataFrame) -> dict[str, str]:
    """Study -> region_label of the primer set it was compared with (primer_sets.csv lists
    each set's studies, ';'-separated; census.csv says which studies are full length)."""
    full = set(census.loc[census["region"] == "full-length", "study_id"])
    return {s: region_label(r.region, r.primers, r.ref_name, s in full)
            for r in primer_sets.itertuples() for s in r.studies.split(";")}


def level_presence(level_support: pd.DataFrame, studies: list[str]) -> pd.DataFrame:
    """Boolean table (study x factor:level): does the study have analysed rows at that level?

    level_support rows hold a factor, a level and the ';'-separated studies with rows there
    (results/phase3/primary/level_support.csv). Studies not listed for a level get False.
    """
    cols = {}
    for r in level_support.itertuples():
        listed = set(str(r.studies).split(";"))
        cols[f"{r.factor}:{r.level}"] = [s in listed for s in studies]
    return pd.DataFrame(cols, index=studies)


def omitted_level(estimates: pd.Series, vcov: pd.DataFrame, cols: list[str]) -> tuple[float,
                                                                                       float]:
    """Estimate and standard error of the level left out of a sum-to-zero coded factor.

    With deviation coding the omitted level's effect is minus the sum of the coded levels'
    effects; its variance is the sum of every entry of their covariance block (1' V 1).
    """
    est = -float(estimates[cols].sum())
    var = float(vcov.loc[cols, cols].to_numpy().sum())
    return est, float(np.sqrt(var))


def fmt_share(x: float, digits: int = 2) -> str:
    """A share (0-1) as a fixed-decimal string: 0.5012 -> '0.50'."""
    return f"{x:.{digits}f}"


def fmt_pct(x: float, digits: int = 1) -> str:
    """A share (0-1) as a percentage: 0.0331 -> '3.3%'."""
    return f"{100 * x:.{digits}f}%"


def fmt_interval(lo: float, hi: float, digits: int = 2, pct: bool = False) -> str:
    """'(lo-hi)' with the same formatting as the estimate; a negative bound keeps its sign and
    the separator becomes ' to ' so '-0.07-0.36' cannot be misread."""
    f = (lambda v: fmt_pct(v, digits)) if pct else (lambda v: fmt_share(v, digits))
    sep = " to " if lo < 0 or hi < 0 else "-"
    return f"({f(lo)}{sep}{f(hi)})"


def fmt_count(n: float) -> str:
    """An integer count with thousands separators: 665023 -> '665,023'."""
    if float(n) != int(n):
        raise ValueError(f"not a whole number: {n}")
    return f"{int(n):,}"
