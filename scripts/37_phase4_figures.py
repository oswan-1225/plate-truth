"""Phase 4 figures (from results/phase4/, script 36).

fig_phase4_medium_taxa.(png|pdf):
  (a) family x medium-class effects (log-odds deviation of a family on a class from its own
      average and from the class average; conditional modes of the Q5 model) for the families
      reported in Q5. Blue: grew more often on that class; red: less often. A dot marks an
      effect whose 95% interval excludes 0. Blank: the family had no rows on that class.
  (b) the same for phyla.
fig_phase4_curation.(png|pdf):
  Q6: odds ratio for "relatives kept on the plate's own medium class" (match term; share
      from 0 to 1), primary analysis and sensitivity analyses, and the matched-plate
      Mantel-Haenszel odds ratio, with 95% intervals. 1 = no medium-specific information.
      Below the line (grey, exploratory): recorded-class effects estimated per study, and
      the registered model fitted within each of the two largest study groups.
fig_phase4_under_recovered.(png|pdf):
  Q7: observed / expected enriched rows per family against the relatedness model (held
      out, each study's level shifted to its own total), with 95% group-bootstrap intervals.
      Red: called systematically under-recovered; blue: over-recovered; grey: neither. Families
      with no enriched rows at all (O/E = 0) are drawn as triangles at 0.003.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from plate_truth.config import data_path, load_config

# Reference palette: diverging blue <-> red with a neutral grey midpoint; ink and grid tokens
# as in the Phase 0 and Phase 3 figures.
BLUE, RED, MID = "#2a78d6", "#e34948", "#f0efec"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
DIVERGING = LinearSegmentedColormap.from_list("blue_grey_red", [RED, MID, BLUE])
CLASS_LABELS = {"rich_complex": "Rich", "dilute_oligotrophic": "Dilute",
                "defined_minimal": "Defined", "environment_extract": "Habitat\nextract",
                "selective": "Selective"}


def style(ax, grid_axis="y"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    if grid_axis:
        getattr(ax, f"{grid_axis}axis").grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def heatmap(ax, t: pd.DataFrame, taxon: str, classes: list[str], vmax: float):
    """Taxon x class effects; NaN where the taxon has no rows on a class."""
    order = t.groupby(taxon)["effect"].apply(lambda s: s.abs().max()).sort_values().index
    cols = [c for c in classes if c in set(t["medium"])]
    eff = t.pivot(index=taxon, columns="medium", values="effect").reindex(order)[cols]
    sig = t.pivot(index=taxon, columns="medium", values="interval_excludes_0").reindex(
        order)[cols].fillna(False).astype(bool)
    im = ax.imshow(eff.to_numpy(float), cmap=DIVERGING, aspect="auto",
                   norm=TwoSlopeNorm(0.0, -vmax, vmax), interpolation="nearest")
    yy, xx = np.nonzero(sig.to_numpy())
    ax.scatter(xx, yy, s=10, color=INK, linewidths=0)
    ax.set_xticks(range(len(cols)), [CLASS_LABELS[c] for c in cols], fontsize=8)
    ax.set_yticks(range(len(order)), order, fontsize=7)
    ax.tick_params(length=0, colors=MUTED, labelcolor=INK)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def forest(ax, labels, est, lo, hi, colors, xlabel):
    y = np.arange(len(labels))[::-1]
    ax.hlines(y, lo, hi, color=colors, lw=2)
    ax.scatter(est, y, s=36, color=colors, zorder=3, edgecolor="white", linewidth=1)
    ax.axvline(1.0, color=MUTED, lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(y, labels, fontsize=8)
    ax.set_xlabel(xlabel)
    style(ax, "x")


def main() -> None:
    res = data_path("results") / "phase4"
    out = data_path("figures")
    classes = load_config("phase4")["medium_classes"]

    fam = pd.read_csv(res / "q5_family_medium.csv")
    phy = pd.read_csv(res / "q5_phylum_medium.csv")
    vmax = float(np.ceil(max(fam["effect"].abs().max(), phy["effect"].abs().max()) * 2) / 2)
    fig, axes = plt.subplots(1, 2, figsize=(11, max(6, 0.17 * fam["family"].nunique() + 1.5)),
                             gridspec_kw={"width_ratios": [1, 1], "wspace": 0.75})
    im = heatmap(axes[0], fam, "family", classes, vmax)
    axes[0].set_title("(a) Families", loc="left", fontsize=10, color=INK)
    heatmap(axes[1], phy, "phylum_full", classes, vmax)
    axes[1].set_title("(b) Phyla", loc="left", fontsize=10, color=INK)
    cb = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.02)
    cb.set_label("Taxon x medium-class effect (log-odds)\nblue: grew more often on the class",
                 fontsize=8, color=INK)
    cb.outline.set_visible(False)
    fig.savefig(out / "fig_phase4_medium_taxa.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / "fig_phase4_medium_taxa.pdf", bbox_inches="tight")
    plt.close(fig)

    q6 = pd.read_csv(res / "q6_summary.csv")
    q6 = q6[q6["odds_ratio"].notna()]
    names = {"primary": "Match term, primary", "identity_ge_97": "Identity at least 97%",
             "indicators": "Any strain recorded (indicator)",
             "strain_16s_links": "Strain 16S links only"}
    labels = [names[a] for a in q6["analysis"]]
    est, lo, hi = (list(q6[c]) for c in ("odds_ratio", "or_lo", "or_hi"))
    p = q6[q6["analysis"] == "primary"].iloc[0]
    labels.append("Matched plates (Mantel-Haenszel)")
    est.append(p["pairs_or"]), lo.append(p["pairs_or_lo"]), hi.append(p["pairs_or_hi"])
    n_registered = len(labels)
    ex = pd.read_csv(res / "exploratory" / "curation_by_study.csv").set_index("analysis")
    for key, lab in (("bystudy_no_sponges", "Recorded effects per study*"),
                     ("li2023", "Li et al. 2023 alone*"), ("zheng2024", "Zheng et al. 2024 alone*")):
        labels.append(lab)
        est.append(ex.at[key, "odds_ratio"]), lo.append(ex.at[key, "or_lo"])
        hi.append(ex.at[key, "or_hi"])
    colors = [BLUE] * (n_registered - 1) + [INK] + [MUTED] * (len(labels) - n_registered)

    fig, ax = plt.subplots(figsize=(6.5, 0.42 * len(labels) + 1.2))
    forest(ax, labels, est, lo, hi, colors, "Odds ratio (relatives kept on the plate's class)")
    ax.axhline(len(labels) - n_registered - 0.5, color=GRID, lw=1)
    ax.set_xticks([0.5, 0.75, 1, 1.5, 2, 3, 5], ["0.5", "0.75", "1", "1.5", "2", "3", "5"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_title("Curation hypothesis: medium-specific information in collection records",
                 loc="left", fontsize=10, color=INK)
    ax.text(0.0, -0.18, "* exploratory, not registered", transform=ax.transAxes, fontsize=7,
            color=MUTED)
    fig.savefig(out / "fig_phase4_curation.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / "fig_phase4_curation.pdf", bbox_inches="tight")
    plt.close(fig)

    q7 = pd.read_csv(res / "q7_family.csv").sort_values(["oe", "expected"],
                                                        ascending=[True, False])
    call_color = {"under": RED, "over": BLUE, "neither": MUTED}
    floor = 0.003                      # log axis: O/E of 0 is drawn at this floor
    zero = q7["oe"] == 0
    fig, ax = plt.subplots(figsize=(6.5, 0.15 * len(q7) + 1.5))
    forest(ax, list(q7["family"]), q7["oe"].clip(lower=floor), q7["oe_lo"].clip(lower=floor),
           q7["oe_hi"].clip(lower=floor), [call_color[c] for c in q7["call"]],
           "Observed / expected enriched rows")
    y = np.arange(len(q7))[::-1]
    ax.scatter(np.full(zero.sum(), floor), y[zero.to_numpy()], marker="<", s=40,
               color=RED, zorder=4, edgecolor="white", linewidth=0.8)
    ax.set_xlim(floor * 0.7, 20)
    ax.set_xticks([0.01, 0.1, 1, 10], ["0.01", "0.1", "1", "10"])
    ax.tick_params(axis="y", labelsize=6)
    ax.set_title("Families against the relatedness model", loc="left", fontsize=10, color=INK)
    handles = [plt.Line2D([], [], color=c, marker="o", lw=2, label=lab) for lab, c in
               (("Under-recovered", RED), ("Over-recovered", BLUE), ("Neither", MUTED))]
    handles.append(plt.Line2D([], [], color=RED, marker="<", lw=0,
                              label="None observed (drawn at 0.003)"))
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="upper left",
              bbox_to_anchor=(1.01, 1.0))
    fig.savefig(out / "fig_phase4_under_recovered.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / "fig_phase4_under_recovered.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
