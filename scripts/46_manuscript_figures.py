"""Phase 6: the six manuscript figures, from results/ (scripts 20, 26, 28, 36, 43).

fig1_census.(png|pdf): the ten study groups (eleven studies) analysed. (a) plating conditions
  and source samples per study; (b) which medium classes, temperature bands and oxygen regimes
  each study's analysed rows cover, with the amplified 16S region.
fig2_calibration.(png|pdf): the headline. (a) P(enriched) against identity for a new study,
  with the held-out observed share per identity bin; (b) held-out AUC per study group;
  (c) held-out observed / expected per study group.
fig3_recovery.(png|pdf): (a) per study group, the share of source reads with a cultured
  relative (>= 94.5%), enriched in any condition, enriched in the median condition, and
  predicted for the median condition; (b) shortcut / observed (any condition) with 95% intervals.
fig4_shifts.(png|pdf): (a) the curve for environmental and host-associated habitats (habitat
  model); (b) held-out calibration slope per study group (how sharply identity separates growers
  in that study, relative to the model); (c) condition effects as odds ratios relative to the
  average level of each factor, every level shown (the level omitted by sum-to-zero coding is
  recovered from the covariance matrix). The table is written to
  results/phase6/condition_effects.csv.
fig5_medium_taxa / fig6_application: copies of fig_phase4_medium_taxa (script 37) and
  fig_phase5_application (script 45), so the manuscript numbering has one file per figure.
"""

import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

from plate_truth.config import data_path
from plate_truth.report import group_habitats, level_presence, omitted_level, study_label

# Reference palette (validated categorical slots in fixed order), as in the phase figures.
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
BLUE, ORANGE, AQUA = SLOTS[:3]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
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
FIT_DIR = ("interim", "phase3/fits_rev1/primary/primary/all")


def style(ax, grid_axis="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    if grid_axis:
        getattr(ax, f"{grid_axis}axis").grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def title(ax, text):
    ax.set_title(text, loc="left", color=INK, fontsize=11)


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


# ---------------------------------------------------------------- Figure 1: census overview
def fig1() -> None:
    census = census_table().set_index("study_id")
    rows = pd.read_csv(data_path("results") / "phase3" / "rows_summary.csv")
    rows = rows[rows["variant"] == "primary"].set_index("study_id")
    links = pd.read_csv(data_path("results") / "phase2" / "qc_asvs_links.csv").set_index(
        "study_id")
    studies = sorted(rows.index, key=lambda s: (HABITATS.index(census.loc[s, "habitat"]), s))
    pres = level_presence(pd.read_csv(data_path("results") / "phase3" / "primary" /
                                      "level_support.csv"), studies)
    factors = [f for f in LEVEL_LABELS]
    cols = [c for f in factors for c in pres.columns if c.split(":")[0] == f
            and c.split(":")[1] in LEVEL_LABELS[f]]
    cols = sorted(cols, key=lambda c: (factors.index(c.split(":")[0]),
                                       list(LEVEL_LABELS[c.split(":")[0]]).index(
                                           c.split(":")[1])))
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True,
                               gridspec_kw={"width_ratios": [0.8, 1.6], "wspace": 0.05})
    style(a, "x")
    y = np.arange(len(studies))
    colors = [HAB_COLOR[census.loc[s, "habitat"]] for s in studies]
    a.barh(y, rows.loc[studies, "conditions"], color=colors, height=0.6)
    for yi, s in zip(y, studies, strict=True):
        a.text(rows.loc[s, "conditions"] * 1.12, yi,
               f"{int(links.loc[s, 'sources_linked'])} source"
               f"{'s' if links.loc[s, 'sources_linked'] != 1 else ''}",
               va="center", fontsize=8, color=MUTED)
    a.set_xscale("log")
    a.set_xlim(1, 3000)
    a.set_xlabel("Plating conditions analysed (log scale)", color=INK)
    a.set_yticks(y, [f"{study_label(census.loc[s, 'study_group'])} · {census.loc[s, 'habitat']}"
                     for s in studies])
    a.invert_yaxis()
    title(a, "a  Studies, conditions and source samples")

    style(b, None)
    b.tick_params(left=False)
    b.spines["left"].set_visible(False)
    xs, x, prev = [], 0.0, None
    for c in cols:
        f = c.split(":")[0]
        if prev is not None and f != prev:
            x += 0.8
        xs.append(x)
        x += 1
        prev = f
    for xi, c in zip(xs, cols, strict=True):
        f, lev = c.split(":")
        known = lev != "unknown"
        for yi, s in zip(y, studies, strict=True):
            if pres.loc[s, c]:
                b.scatter(xi, yi, s=70, color=INK if known else "white", edgecolors=INK,
                          linewidths=1.2, zorder=3)
            else:
                b.scatter(xi, yi, s=8, color=GRID, zorder=2)
    for f in factors:
        fx = [xi for xi, c in zip(xs, cols, strict=True) if c.startswith(f + ":")]
        b.text(np.mean(fx), -1.25, FACTOR_TITLES[f], ha="center", fontsize=9, color=INK)
    b.set_xticks(xs, [LEVEL_LABELS[c.split(":")[0]][c.split(":")[1]] for c in cols],
                 rotation=45, ha="right", fontsize=8)
    reg_x = xs[-1] + 1.6
    b.text(reg_x, -1.25, "16S region", ha="left", fontsize=9, color=INK)
    for yi, s in zip(y, studies, strict=True):
        b.text(reg_x, yi, census.loc[s, "region"], va="center", fontsize=8, color=INK)
    b.set_xlim(-0.6, reg_x + 2.2)
    b.set_ylim(len(studies) - 0.5, -1.8)
    title(b, "b  Conditions covered (filled) and not classed (open)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=HAB_COLOR[h], label=h) for h in HABITATS]
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.9, bottom=0.2)
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


def curve_panel(ax, curve, colour=BLUE, label="Model, new study", spread=True):
    c = curve[curve["no_hit"] == 0]
    ax.fill_between(c["identity"], 100 * c["p_lo"], 100 * c["p_hi"], color=colour, alpha=0.15,
                    lw=0, label="95% interval" if spread else None)
    if spread:
        for col in ("study_p10", "study_p90"):
            ax.plot(c["identity"], 100 * c[col], color=colour, lw=1.1, ls="--",
                    label="10th / 90th percentile study" if col == "study_p10" else None)
    ax.plot(c["identity"], 100 * c["p"], color=colour, lw=2, label=label)


def fig2() -> None:
    res = data_path("results") / "phase3" / "primary"
    by = pd.read_csv(res / "heldout_by_group.csv")
    by = by[by["model"] == "primary"].sort_values("oe").reset_index(drop=True)
    calib = pd.read_csv(res / "calibration_by_identity_bin.csv")
    calib = calib[calib["model"] == "primary"]
    fig = plt.figure(figsize=(15, 5))
    axes = [fig.add_axes(r) for r in ([0.05, 0.12, 0.40, 0.8], [0.66, 0.12, 0.15, 0.8],
                                      [0.83, 0.12, 0.15, 0.8])]
    a = axes[0]
    style(a)
    curve_panel(a, pd.read_csv(res / "curve.csv"))
    x = calib["identity_bin"].map(bin_x)
    a.errorbar(x, 100 * calib["observed"],
               yerr=[100 * (calib["observed"] - calib["observed_lo"]),
                     100 * (calib["observed_hi"] - calib["observed"])],
               fmt="o", ms=8, color=ORANGE, elinewidth=2, capsize=0, mec="white", mew=1,
               label="Observed, each study held out", zorder=4)
    a.plot(x, 100 * calib["predicted"], "_", ms=16, mew=2.5, color=INK, zorder=5,
           label="Mean held-out prediction")
    a.set_xlim(75, 100.5)
    a.set_xlabel("16S identity to nearest cultured relative (%)", color=INK)
    a.set_ylabel("ASVs enriched on plates (%)", color=INK)
    a.legend(frameon=False, fontsize=8, loc="upper left")
    title(a, "a  Growth on plates by identity")

    y = np.arange(len(by))
    labels = group_labels(by["study_group"])
    for ax, col, ref, xlab, ttl in ((axes[1], "auc", 0.5, "Held-out AUC", "b  Ranking"),
                                    (axes[2], "oe", 1.0, "Observed / expected (log)",
                                     "c  Level")):
        style(ax, "x")
        ax.axvline(ref, color=MUTED, lw=1)
        ax.plot(by[col], y, "o", ms=8, color=BLUE, mec="white", mew=1)
        ax.set_yticks(y, labels if ax is axes[1] else [])
        ax.set_xlabel(xlab, color=INK)
        title(ax, ttl)
    axes[1].set_xlim(0.45, 1.0)
    axes[2].axvspan(0.67, 1.5, color=GRID, alpha=0.6, lw=0)
    axes[2].set_xscale("log")
    ticks = [0.2, 0.5, 1, 2, 3]
    axes[2].set_xticks(ticks, [f"{t:g}" for t in ticks])
    axes[2].xaxis.set_minor_formatter(plt.NullFormatter())
    axes[2].set_xlim(0.15, 4)
    save(fig, "fig2_calibration")


# ---------------------------------------------------------------- Figure 3: recovery per study
def fig3() -> None:
    rec = pd.read_csv(data_path("results") / "phase3" / "primary" / "recovery_by_group.csv")
    rec = rec.sort_values("proxy_over_observed").reset_index(drop=True)
    y = np.arange(len(rec))
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True,
                               gridspec_kw={"width_ratios": [1.3, 0.8], "wspace": 0.08})
    style(a, "x")
    series = [("proxy_has_relative", ORANGE, "o", "Has a cultured relative (≥94.5%): shortcut"),
              ("observed_any_condition", BLUE, "o", "Enriched, any condition (observed)"),
              ("observed_median", BLUE, "D", "Enriched, median condition (observed)"),
              ("predicted_median", AQUA, "D", "Predicted, median condition (held out)")]
    cols = [s[0] for s in series]
    for yi, (_, r) in zip(y, rec.iterrows(), strict=True):
        a.plot([100 * r[cols].min(), 100 * r[cols].max()], [yi, yi], color=GRID, lw=2, zorder=0)
    for col, colour, marker, lab in series:
        a.plot(100 * rec[col], y, marker, ms=7, color=colour, mec="white", mew=1, label=lab)
    a.set_yticks(y, group_labels(rec["study_group"]))
    a.set_xlim(0, 100)
    a.set_xlabel("Share of the source community's reads (%)", color=INK)
    a.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.45, -0.13),
             ncol=2)
    title(a, "a  Shortcut, observed and predicted recovery")
    style(b, "x")
    b.axvline(1, color=MUTED, lw=1)
    one = rec["ratio_lo"].isna()
    b.errorbar(rec.loc[~one, "proxy_over_observed"], y[~one],
               xerr=[rec.loc[~one, "proxy_over_observed"] - rec.loc[~one, "ratio_lo"],
                     rec.loc[~one, "ratio_hi"] - rec.loc[~one, "proxy_over_observed"]],
               fmt="o", ms=7, color=ORANGE, mec="white", mew=1, elinewidth=1.5, capsize=0,
               label="95% interval, sources resampled")
    b.plot(rec.loc[one, "proxy_over_observed"], y[one], "o", ms=7, mfc="white", mec=ORANGE,
           mew=1.5, lw=0, label="one source sample: no interval")
    b.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    b.set_xscale("log")
    ticks = [1, 3, 10, 30, 100, 300]
    b.set_xticks(ticks, [f"{t:g}" for t in ticks])
    b.xaxis.set_minor_formatter(plt.NullFormatter())
    b.set_xlim(0.8, 500)
    b.set_xlabel("Shortcut / enriched in any condition (log)", color=INK)
    title(b, "b  Overstatement factor")
    b.invert_yaxis()
    fig.subplots_adjust(left=0.22, right=0.98, top=0.92, bottom=0.24)
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


def fig4() -> None:
    res = data_path("results") / "phase3" / "primary"
    eff = condition_effects()
    out = data_path("results") / "phase6"
    out.mkdir(parents=True, exist_ok=True)
    eff.to_csv(out / "condition_effects.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2),
                             gridspec_kw={"width_ratios": [1.1, 0.75, 0.8], "wspace": 0.75})
    a = axes[0]
    style(a)
    host = pd.read_csv(res / "curve_by_host.csv")
    for h, colour, lab in ((-1.0, BLUE, "Environmental habitats"),
                           (1.0, ORANGE, "Host-associated habitats")):
        curve_panel(a, host[host["host"] == h], colour, lab, spread=False)
    a.set_xlim(75, 100.5)
    a.set_xlabel("16S identity to nearest cultured relative (%)", color=INK)
    a.set_ylabel("ASVs enriched on plates, new study (%)", color=INK)
    handles, _ = a.get_legend_handles_labels()
    handles.append(plt.Rectangle((0, 0), 1, 1, color=MUTED, alpha=0.2, label="95% interval"))
    a.legend(handles=handles, frameon=False, fontsize=8, loc="upper left")
    title(a, "a  Environmental vs host-associated")

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
    b.legend(handles=[plt.Line2D([], [], marker="o", color=BLUE, lw=0, label="environmental"),
                      plt.Line2D([], [], marker="o", color=ORANGE, lw=0,
                                 label="host-associated")],
             frameon=False, fontsize=8, loc="lower right")
    b.set_xlabel("Held-out calibration slope", color=INK)
    title(b, "b  Strength of identity per study")

    c = axes[2]
    style(c, "x")
    order = [(f, lev) for f in ("medium", "temperature", "oxygen")
             for lev in LEVEL_LABELS[f] if lev in set(eff.loc[eff["factor"] == f, "level"])]
    order += [("continuous", lev) for lev in eff.loc[eff["factor"] == "continuous", "level"]]
    labels, yy, pos = [], [], 0.0
    prev = None
    for f, lev in order:
        if prev is not None and f != prev:
            pos += 0.7
        r = eff[(eff["factor"] == f) & (eff["level"] == lev)].iloc[0]
        c.errorbar(r["odds_ratio"], pos, xerr=[[r["odds_ratio"] - r["or_lo"]],
                                                [r["or_hi"] - r["odds_ratio"]]],
                   fmt="o", ms=6, color=INK, elinewidth=1.5)
        name = LEVEL_LABELS[f][lev] if f in LEVEL_LABELS else lev
        labels.append(f"{FACTOR_TITLES.get(f, 'Per unit')}: {name}" if f != "continuous"
                      else name[0].upper() + name[1:])
        yy.append(pos)
        pos += 1
        prev = f
    c.axvline(1, color=MUTED, lw=1)
    c.set_xscale("log")
    ticks = [0.3, 0.5, 0.7, 1, 1.4, 2]
    c.set_xticks(ticks, [f"{t:g}" for t in ticks])
    c.set_xlim(0.25, 2)
    c.xaxis.set_minor_formatter(plt.NullFormatter())
    c.set_yticks(yy, labels, fontsize=8)
    c.invert_yaxis()
    c.set_xlabel("Odds ratio (log scale)", color=INK)
    title(c, "c  Conditions and covariates")
    fig.subplots_adjust(left=0.05, right=0.99, top=0.92, bottom=0.12)
    save(fig, "fig4_shifts")


def copy_figures() -> None:
    out = data_path("figures")
    for src, dst in (("fig_phase4_medium_taxa", "fig5_medium_taxa"),
                     ("fig_phase5_application", "fig6_application")):
        for ext in ("png", "pdf"):
            shutil.copyfile(out / f"{src}.{ext}", out / f"{dst}.{ext}")
        print(f"copied figures/{src} -> {dst}")


def main() -> None:
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig1()
    fig2()
    fig3()
    fig4()
    copy_figures()


if __name__ == "__main__":
    main()
