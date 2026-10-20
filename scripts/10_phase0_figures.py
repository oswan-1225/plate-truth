"""Phase 0 figures.

fig_phase0_calibration.(png|pdf):
  (a) observed share of soil ASVs that grew, per nearest-relative identity bin, for the
      enriched (primary) and detected definitions, with soil-bootstrap 95% intervals, and the
      mean held-out prediction for that bin (open marker).
  (b) per soil: read-weighted share of the soil community that has a cultured relative
      (>= 94.5%, the database proxy), that the held-out calibration predicts will grow, and
      that was observed enriched on plates.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plate_truth.config import data_path

# Reference palette, categorical slots 1-3 (validated all-pairs, light mode).
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
LABELS = {"enriched_on_plates": "Enriched on plates (primary)",
          "detected_on_plates": "Detected on plates (≥2 reads)"}


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def panel_bins(ax, calib):
    sub = calib[calib["unit"] == "asv_soil"]
    bins = list(dict.fromkeys(sub["identity_bin"]))
    x = np.arange(len(bins))
    for (outcome, colour), shift in zip([("enriched_on_plates", BLUE),
                                         ("detected_on_plates", ORANGE)], (-0.12, 0.12),
                                        strict=True):
        d = sub[sub["outcome"] == outcome].set_index("identity_bin").reindex(bins)
        err = np.vstack([d["observed"] - d["boot_lo"], d["boot_hi"] - d["observed"]])
        ax.errorbar(x + shift, d["observed"], yerr=np.clip(err, 0, None), fmt="o", ms=7,
                    color=colour, ecolor=colour, elinewidth=2, capsize=0,
                    label=LABELS[outcome] + ", observed")
        ax.plot(x + shift, d["predicted"], "o", ms=8, mfc="none", mec=colour, mew=1.5,
                label=LABELS[outcome].split(" (")[0] + ", held-out prediction")
    ax.set_xticks(x, [b if b != "none" else "no hit" for b in bins])
    ax.set_xlabel("16S identity to nearest cultured relative (%)", color=INK)
    ax.set_ylabel("Share of soil ASVs", color=INK)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_title("a  Growth on plates by identity", loc="left", color=INK, fontsize=11)


def panel_recovery(ax, rec):
    d = rec[(rec["unit"] == "asv_soil") & (rec["outcome"] == "enriched_on_plates")]
    d = d.set_index("source_sample").sort_values("observed_grew")
    y = np.arange(len(d))
    for col, colour, lab in [("proxy_has_relative", ORANGE, "Has cultured relative (≥94.5%)"),
                             ("predicted_calibrated", AQUA, "Predicted (held-out calibration)"),
                             ("observed_grew", BLUE, "Observed enriched on plates")]:
        ax.plot(d[col], y, "o", ms=8, color=colour, label=lab, mec="white", mew=1)
    for yi, (_, r) in zip(y, d.iterrows(), strict=True):
        lo = r[["proxy_has_relative", "predicted_calibrated", "observed_grew"]].min()
        hi = r[["proxy_has_relative", "predicted_calibrated", "observed_grew"]].max()
        ax.plot([lo, hi], [yi, yi], color=GRID, lw=2, zorder=0)
    ax.set_yticks(y, d.index)
    ax.set_xlabel("Read-weighted share of the soil community", color=INK)
    ax.set_xlim(0, 1)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.yaxis.grid(False)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("b  Recoverable share per soil", loc="left", color=INK, fontsize=11)


def main() -> None:
    res = data_path("results") / "phase0"
    calib = pd.read_csv(res / "calibration_bins.csv")
    rec = pd.read_csv(res / "recovery_by_soil.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.3, 1]})
    for ax in axes:
        style(ax)
    panel_bins(axes[0], calib)
    panel_recovery(axes[1], rec)
    fig.tight_layout()
    out = data_path("figures")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"fig_phase0_calibration.{ext}", dpi=300)
    print("wrote figures/fig_phase0_calibration.png/.pdf")


if __name__ == "__main__":
    main()
