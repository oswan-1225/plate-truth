"""Phase 5 figures (from results/phase5/, scripts 42-44).

fig_phase5_application.(png|pdf) (Figure 6):
  One row per habitat (training habitats first; extrapolated habitats marked "extrapolation";
  human gut last), one panel per arm (amplicon, shotgun).
  - Dark bar: predicted share of the sample's reads in taxa enriched on ONE plate, from a
    low-recovery study (10th percentile of the study effect) to a high-recovery one (90th),
    averaged over the training plating conditions; whiskers at each end: 95% intervals.
  - Light bar: the same for ANY of 8 plates (generous; plates assumed independent).
  - Diamond: the database shortcut, share of reads with a cultured relative at >= 94.5%
    identity, with its 95% interval over studies.
  - Filled diamond: the shortcut exceeds the upper end of the one-plate range (registered
    test); a ring around it: it also exceeds the any-of-8 upper end.
  Numbers of studies and samples are printed beside each habitat, with a note when a study
  is flagged in config/phase5.yaml study_notes (one sample sequenced several times; a
  laboratory-reactor study).
fig_phase5_composition.(png|pdf):
  Share of reads in families called under-recovered, over-recovered or neither in Phase 4
  (and families not evaluated), per habitat and arm (mean over studies).
fig_phase5_training_check.(png|pdf):
  The procedure applied to the ten training studies with held-out models: predicted one-plate
  range (bar), observed per-plate recovery (dot), observed any-plate recovery (open circle)
  against the any-of-own-K range (light bar), and the shortcut share (diamond).
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plate_truth.config import data_path, load_config

BLUE, LIGHT, RED, GREY = "#2a78d6", "#a9c9ef", "#e34948", "#9c9a94"
MID = "#6fa3e0"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CALL_COLORS = {"under": RED, "neither": "#c9c7c1", "over": BLUE, "not_evaluated": "#efeeea"}
ARMS = ["amplicon", "shotgun"]
FLAG_TEXT = {"one_sample_replicates": "all runs from one sample",
             "lab_reactor": "includes a lab-reactor study"}


def style(ax, grid_axis="x"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    if grid_axis:
        getattr(ax, f"{grid_axis}axis").grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def habitat_order(cfg: dict) -> list[str]:
    h = cfg["habitats"]
    train = [k for k in h if h[k]["training"]]
    extra = [k for k in h if not h[k]["training"] and h[k]["report"] != "human_gut"]
    gut = [k for k in h if h[k]["report"] == "human_gut"]
    return train + extra + gut


def label(habitat: str, cfg: dict) -> str:
    name = habitat.replace("_", " ")
    return name if cfg["habitats"][habitat]["training"] else f"{name} (extrapolation)"


def study_flags(arm: str, cfg: dict) -> dict[str, str]:
    """Habitat -> note on the arm's flagged studies (config study_notes), for the figure."""
    f = data_path("results") / "phase5" / arm / "studies.csv"
    if not f.exists():
        return {}
    notes = cfg.get("study_notes", {})
    out = {}
    for hab, st in pd.read_csv(f).groupby("habitat"):
        text = [FLAG_TEXT[notes[s]["flag"]] for s in st["study"]
                if notes.get(s, {}).get("flag") in FLAG_TEXT]
        if text:
            out[hab] = "; ".join(text)
    return out


def load(arm: str) -> pd.DataFrame | None:
    f = data_path("results") / "phase5" / arm / "habitats.csv"
    return pd.read_csv(f) if f.exists() else None


def application(cfg: dict, out) -> None:
    order = habitat_order(cfg)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 0.6 * len(order) + 1.8), sharey=True,
                             gridspec_kw={"wspace": 0.42})
    for ax, arm in zip(axes, ARMS):
        h = load(arm)
        style(ax)
        ax.set_title(f"{arm.capitalize()} surveys", color=INK, loc="left", fontsize=11)
        ax.set_xlim(0, 1)
        ax.set_xlabel("Share of reads", color=INK)
        ax.text(1.03, -0.75, "studies / samples", fontsize=7, color=MUTED, va="center",
                transform=ax.get_yaxis_transform(), clip_on=False)
        if h is None:
            continue
        h = h.set_index("habitat")
        flags = study_flags(arm, cfg)
        for y, hab in enumerate(order):
            if hab not in h.index:
                ax.text(0.01, y, "no study included", va="center", fontsize=8, color=MUTED)
                continue
            r = h.loc[hab]
            ax.plot([r["r8_p10"], r["r8_p90"]], [y, y], color=LIGHT, lw=9, solid_capstyle="butt")
            ax.plot([r["r1_p10"], r["r1_p90"]], [y, y], color=BLUE, lw=9, solid_capstyle="butt")
            for (lo, hi), dy, col in (((r["r1_p10_lo"], r["r1_p10_hi"]), 0.24, BLUE),
                                      ((r["r1_p90_lo"], r["r1_p90_hi"]), 0.24, BLUE),
                                      ((r["r8_p90_lo"], r["r8_p90_hi"]), 0.36, MID)):
                ax.plot([lo, hi], [y + dy, y + dy], color=col, lw=1, marker="|", ms=4,
                        markevery=[0, 1])
            ax.plot([r["shortcut_lo"], r["shortcut_hi"]], [y, y], color=INK, lw=1)
            filled = bool(r["shortcut_exceeds_r1"])
            ax.scatter(r["shortcut"], y, marker="D", s=46, zorder=3, color=RED if filled else
                       "white", edgecolors=RED, linewidths=1.4)
            if bool(r["shortcut_exceeds_r8"]):
                ax.scatter(r["shortcut"], y, marker="o", s=150, zorder=2, facecolors="none",
                           edgecolors=RED, linewidths=1)
            tr = ax.get_yaxis_transform()
            ax.text(1.03, y, f"{int(r['studies'])} / {int(r['samples'])}", va="center",
                    ha="left", fontsize=7.5, color=MUTED, transform=tr, clip_on=False)
            if hab in flags:
                ax.text(1.03, y + 0.33, flags[hab].replace("; ", "\n"), va="top", ha="left",
                        fontsize=6.3, color=MUTED, style="italic", transform=tr,
                        clip_on=False, linespacing=1.0)
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels([label(h, cfg) for h in order])
    axes[0].invert_yaxis()
    handles = [plt.Line2D([], [], color=BLUE, lw=8, label="predicted, one plate (10th-90th "
                                                         "percentile study)"),
               plt.Line2D([], [], color=LIGHT, lw=8, label="predicted, any of 8 plates"),
               plt.Line2D([], [], color=BLUE, lw=1, marker="|", ms=4,
                          label="95% interval, one-plate ends"),
               plt.Line2D([], [], color=MID, lw=1, marker="|", ms=4,
                          label="95% interval, any-of-8 upper end"),
               plt.Line2D([], [], marker="D", color=RED, lw=0, label="shortcut (>= 94.5% "
                                                                     "identity); filled = "
                                                                     "exceeds one-plate range"),
               plt.Line2D([], [], marker="o", markerfacecolor="none", color=RED, lw=0, ms=10,
                          label="also exceeds any-of-8 range")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.2, right=0.9, top=0.93, bottom=0.17)
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig_phase5_application.{ext}", dpi=300)
    plt.close(fig)


def composition(cfg: dict, out) -> None:
    order = habitat_order(cfg)
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.45 * len(order) + 1.4), sharey=True)
    calls = ["under", "neither", "over", "not_evaluated"]
    for ax, arm in zip(axes, ARMS):
        h = load(arm)
        style(ax)
        ax.set_title(f"{arm.capitalize()} surveys", color=INK, loc="left", fontsize=11)
        ax.set_xlim(0, 1)
        ax.set_xlabel("Share of reads by Phase 4 family call", color=INK)
        if h is None:
            continue
        h = h.set_index("habitat")
        for y, hab in enumerate(order):
            if hab not in h.index:
                ax.text(0.01, y, "no study included", va="center", fontsize=8, color=MUTED)
                continue
            left = 0.0
            for c in calls:
                v = h.loc[hab, f"family_{c}"]
                ax.barh(y, v, left=left, color=CALL_COLORS[c], height=0.7,
                        edgecolor="white", linewidth=0.5)
                left += v
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels([label(h, cfg) for h in order])
    axes[0].invert_yaxis()
    handles = [plt.Rectangle((0, 0), 1, 1, color=CALL_COLORS[c], label=c.replace("_", " "))
               for c in calls]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig_phase5_composition.{ext}", dpi=300)
    plt.close(fig)


def training_check(out) -> None:
    t = pd.read_csv(data_path("results") / "phase5" / "training_check.csv")
    fig, ax = plt.subplots(figsize=(8, 5.2))
    style(ax)
    for y, r in enumerate(t.itertuples()):
        ax.plot([r.rown_p10, r.rown_p90], [y, y], color=LIGHT, lw=9, solid_capstyle="butt")
        ax.plot([r.r1_p10, r.r1_p90], [y, y], color=BLUE, lw=9, solid_capstyle="butt")
        ax.scatter(r.observed_per_condition, y, color=INK, s=22, zorder=3)
        ax.scatter(r.observed_any_condition, y, facecolors="white", edgecolors=INK, s=26,
                   zorder=3)
        ax.scatter(r.shortcut, y, marker="D", color=RED, s=30, zorder=3)
    ax.set_yticks(range(len(t)))
    ax.set_yticklabels([f"{g} ({h}; K={k})" for g, h, k in
                        zip(t["study_group"], t["habitat"], t["k_own"])], fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("symlog", linthresh=0.01)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share of reads (symmetric log scale below 0.01)", color=INK)
    handles = [plt.Line2D([], [], color=BLUE, lw=8, label="predicted one plate, 10th-90th"),
               plt.Line2D([], [], color=LIGHT, lw=8, label="predicted any of K plates"),
               plt.Line2D([], [], marker="o", color=INK, lw=0, label="observed per plate"),
               plt.Line2D([], [], marker="o", markerfacecolor="white", color=INK, lw=0,
                          label="observed any plate"),
               plt.Line2D([], [], marker="D", color=RED, lw=0, label="shortcut")]
    fig.legend(handles=handles, fontsize=7.5, frameon=False, loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig_phase5_training_check.{ext}", dpi=300)
    plt.close(fig)


def main() -> None:
    cfg = load_config("phase5")
    out = data_path("figures")
    application(cfg, out)
    composition(cfg, out)
    training_check(out)


if __name__ == "__main__":
    main()
