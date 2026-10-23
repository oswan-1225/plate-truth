"""Phase 3 figures (from results/phase3/primary/, script 28).

fig_phase3_calibration.(png|pdf):
  (a) P(enriched on plates) against identity to the nearest cultured relative for a study the
      model has not seen, at average conditions: line = estimate, dark band = 95% interval of the
      fixed effects, dashed lines = 10th and 90th percentile across studies. Points: observed share
      enriched per identity bin over all held-out groups (95% Wilson interval) and the mean
      held-out prediction for that bin (black tick).
  (b) held-out calibration per study group: observed / expected (O/E). Shaded: the registered
      "calibrated" range 0.67-1.5.
  (c) per study group, read-weighted share of the source community that has a cultured relative
      (>= 94.5%), that was enriched in at least one condition, and the per-condition medians of
      observed and held-out predicted shares.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plate_truth.config import data_path, load_config

# Reference palette, categorical slots 1-3 (validated, light mode), as in the Phase 0 figure.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def bin_centres(labels: list[str]) -> dict[str, float]:
    """x position for each identity bin label ('none' plotted at 77)."""
    out = {}
    for lab in labels:
        if lab == "none":
            out[lab] = 77.0
        elif lab.startswith("<"):
            out[lab] = 85.0
        elif lab.startswith(">="):
            out[lab] = 99.5
        else:
            lo, hi = map(float, lab.split("-"))
            out[lab] = (lo + hi) / 2
    return out


def panel_curve(ax, curve: pd.DataFrame, calib: pd.DataFrame):
    c = curve[curve["no_hit"] == 0]
    ax.fill_between(c["identity"], 100 * c["p_lo"], 100 * c["p_hi"], color=BLUE, alpha=0.18,
                    lw=0, label="95% interval")
    for col, lab in (("study_p10", "Study 10th / 90th percentile"), ("study_p90", None)):
        ax.plot(c["identity"], 100 * c[col], color=BLUE, lw=1.2, ls="--", label=lab)
    ax.plot(c["identity"], 100 * c["p"], color=BLUE, lw=2, label="Model, new study")
    d = calib[calib["model"] == "primary"]
    x = d["identity_bin"].map(bin_centres(list(d["identity_bin"])))
    ax.errorbar(x, 100 * d["observed"], yerr=[100 * (d["observed"] - d["observed_lo"]),
                                              100 * (d["observed_hi"] - d["observed"])],
                fmt="o", ms=8, color=ORANGE, ecolor=ORANGE, elinewidth=2, capsize=0,
                mec="white", mew=1, label="Observed, held-out groups", zorder=4)
    ax.plot(x, 100 * d["predicted"], "_", ms=16, mew=2.5, color=INK, zorder=5,
            label="Mean held-out prediction")
    ax.set_xlim(75, 100.5)
    ax.set_xlabel("16S identity to nearest cultured relative (%)", color=INK)
    ax.set_ylabel("Enriched on plates (%)", color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_title("a  Growth on plates by identity", loc="left", color=INK, fontsize=11)


def panel_oe(ax, by_group: pd.DataFrame, oe_range):
    d = by_group[by_group["model"] == "primary"].sort_values("oe")
    y = np.arange(len(d))
    ax.axvspan(*oe_range, color=GRID, alpha=0.6, lw=0)
    ax.axvline(1, color=MUTED, lw=1)
    ax.plot(d["oe"], y, "o", ms=8, color=BLUE, mec="white", mew=1)
    ax.set_xscale("log")
    ticks = [0.2, 0.5, 1, 2, 3]
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlim(0.15, 4)
    ax.set_yticks(y, d["study_group"])
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_xlabel("Observed / expected (held out)", color=INK)
    ax.set_title("b  Calibration per study", loc="left", color=INK, fontsize=11)


def panel_recovery(ax, rec: pd.DataFrame):
    d = rec.set_index("study_group").sort_values("observed_any_condition")
    y = np.arange(len(d))
    series = [("proxy_has_relative", ORANGE, "o", "Has cultured relative (≥94.5%)"),
              ("observed_any_condition", BLUE, "o", "Enriched, any condition"),
              ("observed_median", BLUE, "D", "Enriched, median condition"),
              ("predicted_median", AQUA, "D", "Predicted, median condition")]
    for col, colour, marker, lab in series:
        ax.plot(100 * d[col], y, marker, ms=7, color=colour, mec="white", mew=1, label=lab)
    cols = [s[0] for s in series]
    for yi, (_, r) in zip(y, d.iterrows(), strict=True):
        ax.plot([100 * r[cols].min(), 100 * r[cols].max()], [yi, yi], color=GRID, lw=2,
                zorder=0)
    ax.set_yticks(y, d.index)
    ax.set_xlim(0, 100)
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_xlabel("Share of the source community's reads (%)", color=INK)
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.45, -0.14),
              ncol=2)
    ax.set_title("c  Recoverable share per study", loc="left", color=INK, fontsize=11)


def main() -> None:
    res = data_path("results") / "phase3" / "primary"
    oe_range = load_config("phase3")["criteria"]["oe_range"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8),
                             gridspec_kw={"width_ratios": [1.4, 0.8, 1.1]})
    for ax in axes:
        style(ax)
    panel_curve(axes[0], pd.read_csv(res / "curve.csv"),
                pd.read_csv(res / "calibration_by_identity_bin.csv"))
    panel_oe(axes[1], pd.read_csv(res / "heldout_by_group.csv"), oe_range)
    panel_recovery(axes[2], pd.read_csv(res / "recovery_by_group.csv"))
    fig.tight_layout()
    out = data_path("figures")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig_phase3_calibration.{ext}", dpi=300)
    print("wrote figures/fig_phase3_calibration.png/.pdf")


if __name__ == "__main__":
    main()
