"""Phase 6: the six manuscript figures, from results/ (scripts 20, 26, 28, 36, 43).

Each panel title states what the panel shows in plain words; axes say which direction is
better or worse; jargon (AUC, O/E, calibration slope, log-odds) is explained on the figure.

fig1_census: the ten study groups (eleven studies). (a) plating conditions and source samples
  per study; (b) medium classes, temperature bands and oxygen regimes covered, and the amplified
  region with its primer pair (the reference each study was compared with).
fig2_calibration: the headline. (a) P(enriched) against identity for a new study, with the
  held-out observed share per identity bin; (b) held-out AUC per study group; (c) held-out
  observed / predicted enrichments per study group.
fig3_recovery: (a) per study group, the shortcut share against the share enriched in any
  condition, in the median condition, and predicted for the median condition; (b) shortcut /
  observed (any condition) with 95% intervals.
fig4_shifts: (a) the curve for environmental and host-associated habitats (habitat model);
  (b) held-out calibration slope per study group; (c) condition effects as odds ratios relative
  to the average level of each factor, every level shown (the level omitted by sum-to-zero
  coding is recovered from the covariance matrix; table results/phase6/condition_effects.csv).
fig5_medium_taxa: Q5 medium-by-taxon effects: every reported phylum, and the families with an
  effect of at least FAMILY_MIN_EFFECT log-odds on some class (display rule; the full family map
  is fig_phase4_medium_taxa, Supplementary Fig. S5).
fig6_application: Phase 5 surveys: predicted one-plate and any-of-8 ranges against the
  shortcut, per habitat with data, both survey types; the shortcut / one-plate upper end ratio
  is printed beside each habitat.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import norm

from plate_truth.config import data_path, load_config
from plate_truth.report import (
    group_habitats,
    level_presence,
    omitted_level,
    study_label,
    study_region_labels,
)

# Reference palette (validated categorical slots in fixed order), as in the phase figures.
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
BLUE, ORANGE, AQUA = SLOTS[:3]
LIGHT, RED, MIDGREY = "#a9c9ef", "#e34948", "#f0efec"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
DIVERGING = LinearSegmentedColormap.from_list("red_grey_blue", [RED, MIDGREY, BLUE])
HABITATS = ["soil", "activated sludge", "marine sponge", "plant (oak)", "food (broiler carcass)",
            "rumen"]
HAB_COLOR = dict(zip(HABITATS, SLOTS, strict=True))
LEVEL_LABELS = {
    "medium": {"rich_complex": "rich", "dilute_oligotrophic": "dilute",
               "defined_minimal": "defined", "environment_extract": "habitat extract",
               "selective": "selective", "unknown": "not classed"},
    "temperature": {"cold": "cold (≤20 °C)", "moderate": "moderate", "warm": "warm (≥33 °C)",
                    "unknown": "not resolved"},
    "oxygen": {"aerobic": "aerobic", "anaerobic": "anaerobic",
               "microaerophilic": "microaerophilic", "unknown": "not stated"},
}
FACTOR_TITLES = {"medium": "Medium class", "temperature": "Temperature", "oxygen": "Oxygen"}
CLASS_SHORT = {"rich_complex": "Rich", "dilute_oligotrophic": "Dilute",
               "defined_minimal": "Defined", "environment_extract": "Habitat\nextract",
               "selective": "Selective"}
FIT_DIR = ("interim", "phase3/fits_rev1/primary/primary/all")
FAMILY_MIN_EFFECT = 0.4
HABITAT_NAMES = {"soil": "Soil", "activated_sludge": "Activated sludge", "plant": "Plant",
                 "marine_water": "Marine water", "marine_sediment": "Marine sediment",
                 "freshwater": "Freshwater", "human_gut": "Human gut"}


def style(ax, grid_axis="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    if grid_axis:
        getattr(ax, f"{grid_axis}axis").grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def title(ax, letter, text, sub=None):
    """Panel letter and a one-line message; an optional muted line explains a term."""
    ax.set_title(f"{letter}  {text}", loc="left", color=INK, fontsize=11, fontweight="bold",
                 pad=30 if sub else 8)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=8.5, color=MUTED, va="bottom")


def save(fig, name):
    out = data_path("figures")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"{name}.{ext}", dpi=300)
    plt.close(fig)
    print(f"wrote figures/{name}.png/.pdf")


def census_table() -> pd.DataFrame:
    return pd.read_csv(data_path("processed") / "census.csv")


def group_labels(groups) -> list[str]:
    hab = group_habitats(census_table())
    return [f"{study_label(g)} · {hab[g]}" for g in groups]


def region_labels() -> dict[str, str]:
    return study_region_labels(pd.read_csv(data_path("processed") / "primer_sets.csv"),
                               census_table())


def log_axis(ax, ticks, lim):
    ax.set_xscale("log")
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlim(*lim)


def habitat_label(habitat: str, cfg: dict) -> str:
    """Survey habitat name; * marks a habitat with no training study (extrapolation)."""
    name = HABITAT_NAMES[habitat]
    return name if cfg["habitats"][habitat]["training"] else f"{name}*"


# ---------------------------------------------------------------- Figure 1: census overview
def fig1() -> None:
    census = census_table().set_index("study_id")
    rows = pd.read_csv(data_path("results") / "phase3" / "rows_summary.csv")
    rows = rows[rows["variant"] == "primary"].set_index("study_id")
    links = pd.read_csv(data_path("results") / "phase2" / "qc_asvs_links.csv").set_index(
        "study_id")
    regions = region_labels()
    studies = sorted(rows.index, key=lambda s: (HABITATS.index(census.loc[s, "habitat"]), s))
    pres = level_presence(pd.read_csv(data_path("results") / "phase3" / "primary" /
                                      "level_support.csv"), studies)
    factors = list(LEVEL_LABELS)
    cols = sorted([c for c in pres.columns if c.split(":")[0] in LEVEL_LABELS
                   and c.split(":")[1] in LEVEL_LABELS[c.split(":")[0]]],
                  key=lambda c: (factors.index(c.split(":")[0]),
                                 list(LEVEL_LABELS[c.split(":")[0]]).index(c.split(":")[1])))
    fig, (a, b) = plt.subplots(1, 2, figsize=(15, 6.2), sharey=True,
                               gridspec_kw={"width_ratios": [0.75, 1.7], "wspace": 0.04})
    style(a, "x")
    y = np.arange(len(studies))
    a.barh(y, rows.loc[studies, "conditions"],
           color=[HAB_COLOR[census.loc[s, "habitat"]] for s in studies], height=0.6)
    for yi, s in zip(y, studies, strict=True):
        n = int(links.loc[s, "sources_linked"])
        a.text(rows.loc[s, "conditions"] * 1.12, yi, f"{n} source sample{'s' * (n != 1)}",
               va="center", fontsize=8, color=MUTED)
    a.set_xscale("log")
    a.set_xlim(1, 5000)
    a.set_xlabel("Plating conditions analysed (log scale)", color=INK)
    shared = census.loc[studies, "study_group"].duplicated(keep=False)
    a.set_yticks(y, [study_label(census.loc[s, "study_group"])
                     + (f" ({census.loc[s, 'habitat'].split(' ')[-1]})" if shared[s] else "")
                     for s in studies])
    a.invert_yaxis()
    fig.text(0.01, 0.955, "a  Eleven datasets from ten studies", color=INK, fontsize=11,
             fontweight="bold", ha="left")
    fig.text(0.01, 0.92, "bar colour = habitat; text = environmental samples the plates "
                         "came from", color=MUTED, fontsize=8.5, ha="left")

    style(b, None)
    b.tick_params(left=False, bottom=False, labelbottom=False, top=True, labeltop=True)
    b.spines["left"].set_visible(False)
    b.spines["bottom"].set_visible(False)
    xs, x, prev = [], 0.0, None
    for c in cols:
        f = c.split(":")[0]
        if prev is not None and f != prev:
            x += 0.9
        xs.append(x)
        x += 1
        prev = f
    for xi, c in zip(xs, cols, strict=True):
        known = c.split(":")[1] != "unknown"
        for yi, s in zip(y, studies, strict=True):
            if pres.loc[s, c]:
                b.scatter(xi, yi, s=80, color=INK if known else "white", edgecolors=INK,
                          linewidths=1.3, zorder=3)
            else:
                b.scatter(xi, yi, s=8, color=GRID, zorder=2)
    b.set_xticks(xs, [LEVEL_LABELS[c.split(":")[0]][c.split(":")[1]] for c in cols],
                 rotation=40, ha="left", fontsize=8.5)
    for f in factors:
        fx = [xi for xi, c in zip(xs, cols, strict=True) if c.startswith(f + ":")]
        b.annotate(FACTOR_TITLES[f], xy=(np.mean(fx), len(studies) - 0.2),
                   ha="center", va="top", fontsize=9, color=INK, fontweight="bold",
                   annotation_clip=False)
    reg_x = xs[-1] + 1.4
    b.text(reg_x, len(studies) - 0.2, "16S region and primers", ha="left", va="top",
           fontsize=9, color=INK, fontweight="bold")
    for yi, s in zip(y, studies, strict=True):
        b.text(reg_x, yi, regions[s], va="center", fontsize=8, color=INK)
    b.set_xlim(-0.6, reg_x + 6.2)
    b.set_ylim(len(studies) + 0.5, -0.6)
    fig.text(b.get_position().x0, 0.955, "b  What the analysed plates covered", color=INK,
             fontsize=11, fontweight="bold", ha="left")
    handles = [plt.Rectangle((0, 0), 1, 1, color=HAB_COLOR[h], label=h) for h in HABITATS]
    handles += [plt.Line2D([], [], marker="o", color=INK, lw=0, ms=8,
                           label="condition present"),
                plt.Line2D([], [], marker="o", mfc="white", mec=INK, lw=0, ms=8,
                           label="present but class not known")]
    fig.legend(handles=handles, loc="lower center", ncol=8, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.125, right=0.99, top=0.78, bottom=0.12)
    save(fig, "fig1_census")


# ---------------------------------------------------------------- Figure 2: headline
def bin_x(label: str) -> float:
    """x position of an identity-bin label ('none' plotted at 77, '<90' at 85)."""
    if label == "none":
        return 77.0
    if label.startswith("<"):
        return 85.0
    if label.startswith(">="):
        return 99.5
    lo, hi = map(float, label.split("-"))
    return (lo + hi) / 2


def curve_panel(ax, curve, colour=BLUE, label="Model prediction for a new study",
                spread=True, band_label="95% interval"):
    c = curve[curve["no_hit"] == 0]
    ax.fill_between(c["identity"], 100 * c["p_lo"], 100 * c["p_hi"], color=colour, alpha=0.13,
                    lw=0, label=band_label)
    if spread:
        for col in ("study_p10", "study_p90"):
            ax.plot(c["identity"], 100 * c[col], color=colour, lw=1.1, ls="--",
                    label="low- and high-recovery study (10th, 90th percentile)"
                    if col == "study_p10" else None)
    ax.plot(c["identity"], 100 * c["p"], color=colour, lw=2.2, label=label)


def fig2() -> None:
    res = data_path("results") / "phase3" / "primary"
    by = pd.read_csv(res / "heldout_by_group.csv")
    by = by[by["model"] == "primary"].sort_values("oe").reset_index(drop=True)
    calib = pd.read_csv(res / "calibration_by_identity_bin.csv")
    calib = calib[calib["model"] == "primary"]
    fig = plt.figure(figsize=(15.5, 5.6))
    axes = [fig.add_axes(r) for r in ([0.05, 0.11, 0.38, 0.72], [0.64, 0.11, 0.15, 0.72],
                                      [0.83, 0.11, 0.155, 0.72])]
    a = axes[0]
    style(a)
    curve_panel(a, pd.read_csv(res / "curve.csv"))
    x = calib["identity_bin"].map(bin_x)
    a.errorbar(x, 100 * calib["observed"],
               yerr=[100 * (calib["observed"] - calib["observed_lo"]),
                     100 * (calib["observed_hi"] - calib["observed"])],
               fmt="o", ms=8, color=ORANGE, elinewidth=2, capsize=0, mec="white", mew=1,
               label="Observed share enriched (each study predicted by the others)", zorder=4)
    a.set_xlim(75, 100.5)
    a.set_ylim(0, 26)
    a.set_xlabel("16S identity to the nearest cultured relative (%)", color=INK)
    a.set_ylabel("ASVs enriched on plates (%)", color=INK)
    a.legend(frameon=False, fontsize=8, loc="upper left")
    i99 = calib.loc[calib["identity_bin"].str.startswith(">="), "observed"].iloc[0]
    a.annotate(f"even at ≥99% identity,\n{100 * (1 - i99):.0f}% of ASVs were not enriched",
               xy=(99.5, 100 * i99), xytext=(86.5, 14), fontsize=8.5, color=INK,
               arrowprops={"arrowstyle": "->", "color": MUTED, "lw": 1})
    title(a, "a", "Closer cultured relatives mean more growth, but most still do not grow")

    y = np.arange(len(by))
    b, c = axes[1], axes[2]
    for ax in (b, c):
        style(ax, "x")
    b.axvline(0.5, color=MUTED, lw=1)
    b.plot(by["auc"], y, "o", ms=8, color=BLUE, mec="white", mew=1)
    b.set_yticks(y, group_labels(by["study_group"]))
    b.set_xlim(0.45, 1.0)
    b.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    b.set_xlabel("Held-out AUC", color=INK)
    b.text(0.505, -0.9, "chance", fontsize=8, color=MUTED, va="bottom")
    b.text(0.99, -0.9, "perfect →", fontsize=8, color=MUTED, ha="right", va="bottom")
    title(b, "b", "Ranking works in a new study",
          "AUC: chance that a grower ranks\nabove a non-grower")
    c.axvspan(0.67, 1.5, color=GRID, alpha=0.7, lw=0)
    c.axvline(1, color=MUTED, lw=1)
    c.plot(by["oe"], y, "o", ms=8, color=BLUE, mec="white", mew=1)
    c.set_yticks(y, [])
    log_axis(c, [0.2, 0.5, 1, 2, 3], (0.15, 4))
    c.set_xlabel("Observed ÷ predicted (log)", color=INK)
    c.text(0.16, -0.9, "← model too high", fontsize=8, color=MUTED, va="bottom")
    c.text(3.9, -0.9, "too low →", fontsize=8, color=MUTED, ha="right", va="bottom")
    c.text(1.0, len(by) - 0.45, "acceptable\nrange", fontsize=7.5, color=MUTED, ha="center",
           va="bottom")
    title(c, "c", "The overall level does not",
          "enriched ASVs observed ÷\npredicted, per study")
    for ax in (b, c):
        ax.set_ylim(len(by) - 0.4, -1.1)
    save(fig, "fig2_calibration")


# ---------------------------------------------------------------- Figure 3: recovery per study
def fig3() -> None:
    rec = pd.read_csv(data_path("results") / "phase3" / "primary" / "recovery_by_group.csv")
    rec = rec.sort_values("proxy_over_observed").reset_index(drop=True)
    y = np.arange(len(rec))
    fig, (a, b) = plt.subplots(1, 2, figsize=(13.5, 5.6), sharey=True,
                               gridspec_kw={"width_ratios": [1.35, 0.8], "wspace": 0.06})
    style(a, "x")
    for yi, r in zip(y, rec.itertuples(), strict=True):
        a.plot([100 * r.observed_any_condition, 100 * r.proxy_has_relative], [yi, yi],
               color=GRID, lw=3, zorder=0, solid_capstyle="butt")
    a.plot(100 * rec["proxy_has_relative"], y, "o", ms=9, color=ORANGE, mec="white", mew=1,
           label="Shortcut: reads with a cultured relative (≥94.5%)")
    a.plot(100 * rec["observed_any_condition"], y, "o", ms=9, color=BLUE, mec="white", mew=1,
           label="Observed: reads in ASVs enriched on any plate")
    a.plot(100 * rec["observed_median"], y, "o", ms=6, mfc="white", mec=BLUE, mew=1.5, lw=0,
           label="Observed on a typical plate (median condition)")
    a.plot(100 * rec["predicted_median"], y, "x", ms=6, color=AQUA, mew=2, lw=0,
           label="Predicted for a typical plate (study held out)")
    a.set_yticks(y, group_labels(rec["study_group"]))
    a.set_xlim(0, 100)
    a.set_xlabel("Share of the environmental community's reads (%)", color=INK)
    a.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.45, -0.12),
             ncol=2)
    title(a, "a", "The shortcut counts far more of the community than grew")
    style(b, "x")
    b.axvline(1, color=MUTED, lw=1)
    one = rec["ratio_lo"].isna()
    b.errorbar(rec.loc[~one, "proxy_over_observed"], y[~one],
               xerr=[rec.loc[~one, "proxy_over_observed"] - rec.loc[~one, "ratio_lo"],
                     rec.loc[~one, "ratio_hi"] - rec.loc[~one, "proxy_over_observed"]],
               fmt="o", ms=8, color=ORANGE, mec="white", mew=1, elinewidth=1.5, capsize=0,
               label="95% interval (source samples resampled)")
    b.plot(rec.loc[one, "proxy_over_observed"], y[one], "o", ms=8, mfc="white", mec=ORANGE,
           mew=1.5, lw=0, label="one source sample: no interval")
    for yi, r in zip(y, rec.itertuples(), strict=True):
        right = r.proxy_over_observed if pd.isna(r.ratio_hi) else r.ratio_hi
        b.text(right * 1.2, yi, f"{r.proxy_over_observed:.0f}×" if r.proxy_over_observed >= 10
               else f"{r.proxy_over_observed:.1f}×", va="center", fontsize=8, color=INK)
    log_axis(b, [1, 3, 10, 30, 100, 300], (0.8, 900))
    b.set_xlabel("Shortcut ÷ observed (log)", color=INK)
    b.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    lo, hi = rec["proxy_over_observed"].min(), rec["proxy_over_observed"].max()
    title(b, "b", f"By {lo:.0f} to {hi:.0f} times")
    b.invert_yaxis()
    fig.subplots_adjust(left=0.23, right=0.98, top=0.92, bottom=0.25)
    save(fig, "fig3_recovery")


# ---------------------------------------------------------------- Figure 4: shifts
def condition_effects() -> pd.DataFrame:
    """Odds ratio of every medium, temperature and oxygen level relative to the factor's
    average level, plus the continuous covariates, from the full-data primary fit."""
    fit = data_path(FIT_DIR[0]) / FIT_DIR[1]
    coef = pd.read_csv(fit / "coef.csv").set_index("term")
    vcov = pd.read_csv(fit / "vcov.csv", index_col=0)
    est = coef["estimate"]
    support = pd.read_csv(data_path("results") / "phase3" / "primary" / "level_support.csv")
    out = []
    for f in LEVEL_LABELS:
        coded = [t for t in est.index if t.startswith(f + "_")]
        levels = [lev for lev in support.loc[support["factor"] == f, "level"]
                  if lev != "unknown"]
        for lev in levels:
            term = f"{f}_{lev}"
            if term in coded:
                e, se = float(est[term]), float(coef.loc[term, "se"])
                omitted = False
            else:
                e, se = omitted_level(est, vcov, coded)
                omitted = True
            out.append({"factor": f, "level": lev, "log_or": e, "se": se, "omitted": omitted})
    per = {"log_env": ("source abundance", "per 10-fold"),
           "log_depth": ("plate sequencing depth", "per 10-fold"),
           "log2_days": ("incubation time", "per doubling")}
    for term, (lev, unit) in per.items():
        out.append({"factor": "continuous", "level": f"{lev} ({unit})",
                    "log_or": float(est[term]), "se": float(coef.loc[term, "se"]),
                    "omitted": False})
    t = pd.DataFrame(out)
    z = norm.ppf(0.975)
    t["odds_ratio"] = np.exp(t["log_or"])
    t["or_lo"] = np.exp(t["log_or"] - z * t["se"])
    t["or_hi"] = np.exp(t["log_or"] + z * t["se"])
    return t


def bands_overlap(host: pd.DataFrame) -> bool:
    """Do the environmental and host-associated 95% bands overlap at every identity?"""
    h = host[host["no_hit"] == 0]
    e = h[h["host"] == -1].set_index("identity")
    o = h[h["host"] == 1].set_index("identity")
    return bool(((e["p_lo"] <= o["p_hi"]) & (o["p_lo"] <= e["p_hi"])).all())


def fig4() -> None:
    res = data_path("results") / "phase3" / "primary"
    eff = condition_effects()
    out = data_path("results") / "phase6"
    out.mkdir(parents=True, exist_ok=True)
    eff.to_csv(out / "condition_effects.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.8),
                             gridspec_kw={"width_ratios": [1.05, 0.75, 0.85], "wspace": 0.75})
    a = axes[0]
    style(a)
    host = pd.read_csv(res / "curve_by_host.csv")
    for h, colour, lab in ((-1.0, BLUE, "Environmental (soil, sludge)"),
                           (1.0, ORANGE, "Host-associated (sponge, oak, broiler, rumen)")):
        curve_panel(a, host[host["host"] == h], colour, lab, spread=False, band_label=None)
    a.set_xlim(75, 100.5)
    a.set_xlabel("16S identity to the nearest cultured relative (%)", color=INK)
    a.set_ylabel("ASVs enriched on plates, new study (%)", color=INK)
    a.legend(frameon=False, fontsize=8, loc="upper left")
    overlap = bands_overlap(host)
    a.text(76, 24, "shaded: 95% intervals", fontsize=8, color=MUTED, va="top")
    title(a, "a", "Habitat type: intervals overlap throughout" if overlap else
          "Habitat type shifts the curve")

    b = axes[1]
    style(b, "x")
    by = pd.read_csv(res / "heldout_by_group.csv")
    by = by[by["model"] == "primary"].sort_values("slope").reset_index(drop=True)
    hosted = census_table().groupby("study_group")["host_associated"].first()
    y = np.arange(len(by))
    b.axvline(1, color=MUTED, lw=1)
    for yi, r in zip(y, by.itertuples(), strict=True):
        colour = ORANGE if hosted[r.study_group] else BLUE
        b.errorbar(r.slope, yi, xerr=[[r.slope - r.slope_lo], [r.slope_hi - r.slope]],
                   fmt="o", ms=7, color=colour, mec="white", mew=1, elinewidth=1.5)
    b.set_yticks(y, group_labels(by["study_group"]), fontsize=8)
    b.set_ylim(len(by) - 0.4, -1.1)
    b.text(0.95, -0.9, "← weaker", fontsize=8, color=MUTED, ha="right", va="bottom")
    b.text(1.05, -0.9, "stronger →", fontsize=8, color=MUTED, va="bottom")
    b.set_xlabel("Held-out calibration slope", color=INK)
    title(b, "b", "Identity matters more in some studies",
          "1 = identity separates growers as\nstrongly as the model expects")

    c = axes[2]
    style(c, "x")
    order = [(f, lev) for f in ("medium", "temperature", "oxygen")
             for lev in LEVEL_LABELS[f] if lev in set(eff.loc[eff["factor"] == f, "level"])]
    order += [("continuous", lev) for lev in eff.loc[eff["factor"] == "continuous", "level"]]
    labels, yy, pos, prev = [], [], 0.0, None
    for f, lev in order:
        if prev is not None and f != prev:
            pos += 0.8
        r = eff[(eff["factor"] == f) & (eff["level"] == lev)].iloc[0]
        c.errorbar(r["odds_ratio"], pos, xerr=[[r["odds_ratio"] - r["or_lo"]],
                                                [r["or_hi"] - r["odds_ratio"]]],
                   fmt="o", ms=6, color=INK, elinewidth=1.5)
        name = LEVEL_LABELS[f][lev] if f in LEVEL_LABELS else lev
        labels.append(f"{FACTOR_TITLES[f]}: {name}" if f != "continuous"
                      else name[0].upper() + name[1:])
        yy.append(pos)
        pos += 1
        prev = f
    c.axvline(1, color=MUTED, lw=1)
    log_axis(c, [0.3, 0.5, 0.7, 1, 1.4, 2], (0.25, 2))
    c.set_yticks(yy, labels, fontsize=8)
    c.set_ylim(pos - 0.4, -1.1)
    c.text(0.95, -0.9, "← less growth", fontsize=8, color=MUTED, ha="right", va="bottom")
    c.text(1.05, -0.9, "more →", fontsize=8, color=MUTED, va="bottom")
    c.set_xlabel("Odds ratio (log scale)", color=INK)
    title(c, "c", "Conditions shift growth modestly",
          "each level vs the average of its factor")
    fig.subplots_adjust(left=0.05, right=0.99, top=0.84, bottom=0.11)
    save(fig, "fig4_shifts")


# ---------------------------------------------------------------- Figure 5: medium x taxon
def heat_panel(ax, t: pd.DataFrame, taxon: str, norm_, order: list[str]):
    classes = list(CLASS_SHORT)
    grid = t.pivot(index=taxon, columns="medium", values="effect").reindex(index=order,
                                                                           columns=classes)
    sig = t.pivot(index=taxon, columns="medium", values="interval_excludes_0").reindex(
        index=order, columns=classes)
    ax.imshow(np.ma.masked_invalid(grid.to_numpy(float)), cmap=DIVERGING, norm=norm_,
              aspect="auto")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v, s = grid.iat[i, j], sig.iat[i, j] is True
            if pd.isna(v):
                ax.text(j, i, "–", ha="center", va="center", fontsize=8, color=MUTED)
            elif abs(v) >= 0.25 or s:
                ax.text(j, i, f"{v:+.1f}" + ("*" if s else ""), ha="center", va="center",
                        fontsize=8, color="white" if abs(v) > 0.9 else INK,
                        fontweight="bold" if s else "normal")
    ax.set_xticks(range(len(classes)), [CLASS_SHORT[c] for c in classes], fontsize=8.5)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(order)), order, fontsize=8.5)
    ax.tick_params(length=0, colors=MUTED, labelcolor=INK)
    for side in ax.spines.values():
        side.set_visible(False)


def fig5() -> None:
    res = data_path("results") / "phase4"
    fam = pd.read_csv(res / "q5_family_medium.csv")
    phy = pd.read_csv(res / "q5_phylum_medium.csv")
    keep = fam.groupby("family")["effect"].apply(lambda e: e.abs().max() >= FAMILY_MIN_EFFECT)
    fam = fam[fam["family"].isin(keep[keep].index)]
    fam_order = fam.groupby("family")["effect"].max().sort_values(ascending=False).index.tolist()
    phy_order = (phy.groupby("phylum_full")["effect"].apply(lambda e: e.abs().max())
                 .sort_values(ascending=False).index.tolist())
    norm_ = TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 6.8),
                               gridspec_kw={"width_ratios": [1, 1], "wspace": 0.75})
    heat_panel(a, phy, "phylum_full", norm_, phy_order)
    heat_panel(b, fam, "family", norm_, fam_order)
    a.set_title("a  Phyla: Bacillota prefer rich media", loc="left", color=INK, fontsize=11,
                fontweight="bold", pad=34)
    b.set_title(f"b  Families with an effect of at least {FAMILY_MIN_EFFECT:g}",
                loc="left", color=INK, fontsize=11, fontweight="bold", pad=34)
    sm = plt.cm.ScalarMappable(norm=norm_, cmap=DIVERGING)
    cb = fig.colorbar(sm, ax=[a, b], orientation="horizontal", fraction=0.04, pad=0.04,
                      aspect=45)
    cb.set_label("grew less often on this class  ←  effect (log-odds)  →  more often",
                 color=INK, fontsize=9)
    cb.outline.set_visible(False)
    fig.text(0.01, 0.015, "* 95% interval excludes 0 (bold).  – no rows on that class.  "
                          "Numbers shown where |effect| ≥ 0.25.  Each effect is relative to the "
                          "taxon's own average and the class average, within studies.",
             fontsize=8, color=MUTED)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.88, bottom=0.2)
    save(fig, "fig5_medium_taxa")


# ---------------------------------------------------------------- Figure 6: surveys
def fig6() -> None:
    cfg = load_config("phase5")
    arms = {arm: pd.read_csv(data_path("results") / "phase5" / arm / "habitats.csv")
            .set_index("habitat") for arm in ("amplicon", "shotgun")}
    order = [h for h in cfg["habitats"] if h in arms["amplicon"].index]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.2), sharey=True,
                             gridspec_kw={"wspace": 0.12})
    for ax, (arm, h) in zip(axes, arms.items(), strict=True):
        style(ax, "x")
        for y, hab in enumerate(order):
            if hab not in h.index:
                ax.text(0.01, y, "no survey study", va="center", fontsize=8, color=MUTED)
                continue
            r = h.loc[hab]
            ax.plot([r["r8_p10"], r["r8_p90"]], [y, y], color=LIGHT, lw=11,
                    solid_capstyle="butt")
            ax.plot([r["r1_p10"], r["r1_p90"]], [y, y], color=BLUE, lw=11,
                    solid_capstyle="butt")
            ax.plot([r["r1_p90_lo"], r["r1_p90_hi"]], [y - 0.21, y - 0.21], color=INK, lw=1,
                    marker="|", ms=5, markevery=[0, 1])
            ax.plot([r["shortcut_lo"], r["shortcut_hi"]], [y, y], color=RED, lw=1.2)
            ax.scatter(r["shortcut"], y, marker="D", s=55, zorder=3, color=RED,
                       edgecolors="white", linewidths=0.8)
            ax.text(r["shortcut"], y - 0.17, f"{r['shortcut'] / r['r1_p90']:.0f}×",
                    ha="center", va="bottom", fontsize=8.5, color=INK, fontweight="bold")
            n = int(r["studies"])
            ax.text(0.0, y + 0.36, f"{n} stud{'ies' if n > 1 else 'y'}, "
                                   f"{int(r['samples'])} samples", fontsize=7, color=MUTED,
                    va="center")
        ax.set_xlim(0, 1.08)
        ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0", "20", "40", "60", "80", "100"])
        ax.set_xlabel("Share of the sample's reads (%)", color=INK)
        ax.set_title({"amplicon": "16S amplicon surveys",
                      "shotgun": "Shotgun metagenomes (16S fragments, V4 window)"}[arm],
                     loc="left", color=INK, fontsize=10.5)
    axes[0].set_yticks(range(len(order)), [habitat_label(h, cfg) for h in order])
    axes[0].set_ylim(len(order) - 0.4, -0.7)
    fig.suptitle("In every habitat, the shortcut counts several times more of the community "
                 "as culturable than one plating is predicted to recover",
                 x=0.01, ha="left", fontsize=11.5, fontweight="bold", color=INK)
    handles = [plt.Line2D([], [], color=BLUE, lw=9,
                          label="predicted recovery, one plate: low- to high-recovery study"),
               plt.Line2D([], [], color=LIGHT, lw=9,
                          label="the same on any of 8 plates (generous)"),
               plt.Line2D([], [], color=INK, lw=1, marker="|", ms=6,
                          label="95% interval of the one-plate high end"),
               plt.Line2D([], [], marker="D", color=RED, lw=1.2,
                          label="shortcut: reads with a cultured relative ≥94.5% (95% interval)"),
               plt.Line2D([], [], lw=0,
                          label="N× = shortcut ÷ high end of the one-plate range;  "
                                "* no training study")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8.5)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.88, bottom=0.19)
    save(fig, "fig6_application")


def main() -> None:
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
    fig6()


if __name__ == "__main__":
    main()
