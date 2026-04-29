import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import numpy as np

# ── load data ──────────────────────────────────────────────────────────────────
file_path = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/tt_vs_hv_fits.csv"
plot_folder = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/plots/mdom_transit/"

# ── load data ──────────────────────────────────────────────────────────────────
df = pd.read_csv(file_path)

modules = list(df["module"].unique())
colors  = {"mDOM_D009_v1": "#4a7fc1",
           "mDOM_D192_v1": "#5aab8a",
           "mDOM_D230_v1": "#d4735a",
           "mDOM_D232_v1": "#a06db8"}
labels  = {m: m.replace("mDOM_", "").replace("_v1", "") for m in modules}

# ── binning helper ─────────────────────────────────────────────────────────────
def fine_bins(series, n_bins=40):
    lo, hi = series.min(), series.max()
    pad = (hi - lo) * 0.05
    return np.linspace(lo - pad, hi + pad, n_bins + 1)

# ── shared bin edges (same scale across all modules per parameter) ─────────────
bins_slope     = fine_bins(df["slope"],     n_bins=40)
bins_intercept = fine_bins(df["intercept"], n_bins=40)

# ── tick-inside helper ─────────────────────────────────────────────────────────
def ticks_inside(ax):
    ax.tick_params(axis="both", direction="in",
                   top=True, right=True,
                   which="both", length=4)
    ax.tick_params(axis="both", which="minor", length=2)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

# ── generic 2×2 figure builder ─────────────────────────────────────────────────
def make_figure(col, bins, xlabel, fname, suptitle):
    global_mean = df[col].mean()

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=False)
    fig.suptitle(suptitle, fontsize=13)

    for ax, mod in zip(axes.flat, modules):
        vals = df.loc[df["module"] == mod, col]

        ax.hist(vals, bins=bins,
                color=colors[mod], alpha=0.80,
                edgecolor="white", linewidth=0.5)

        # per-module mean
        mod_mean = vals.mean()
        ax.axvline(mod_mean, color=colors[mod], linewidth=1.4,
                   linestyle="--", label=f"mean = {mod_mean:.4f}")

        # global mean for reference
        ax.axvline(global_mean, color="black", linewidth=1.0,
                   linestyle=":", alpha=0.6, label=f"global = {global_mean:.4f}")

        ax.set_title(labels[mod], fontsize=11, color=colors[mod], fontweight="bold")
        ax.set_ylabel("channels", fontsize=10)
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", linewidth=0.4, alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)
        ticks_inside(ax)

    # shared x-label on bottom row only
    for ax in axes[1]:
        ax.set_xlabel(xlabel, fontsize=10)

    plt.tight_layout()
    plt.savefig(fname, dpi=150, bbox_inches="tight")
    print(f"Saved {fname}")
    # plt.show()
    plt.close()

# ── produce the two figures ────────────────────────────────────────────────────
make_figure("slope",     bins_slope,
            "slope  a  (ns/V)",
            plot_folder + "hist_slope.png",
            "HV fit — slope per module")

make_figure("intercept", bins_intercept,
            "intercept  b  (ns)",
            plot_folder+"hist_intercept.png",
            "HV fit — intercept per module")
