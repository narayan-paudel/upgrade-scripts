"""
fit_gaussian_rise_exp_decay.py
──────────────────────────────
Fits a Gaussian-CDF rise × Exponential decay model to a histogram.

Model
─────
    f(t) = A · Φ((t − μ) / σ) · exp(−(t − μ) / τ)   for t ≥ μ
         = A · Φ((t − μ) / σ)                          for t < μ

where Φ is the standard normal CDF (= 0.5·(1 + erf(·/√2))).

Usage
─────
    # With your own data (1-D array of event times / values):
    #   set USE_SYNTHETIC = False and point DATA_FILE at your file.
    #
    # Out-of-the-box the script generates synthetic noisy data and fits it.

Dependencies: numpy, scipy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import erf
from scipy.optimize import curve_fit
from scipy.stats import chi2

# ══════════════════════════════════════════════════════════════════════════════
#  CONFIG
# ══════════════════════════════════════════════════════════════════════════════
USE_SYNTHETIC = True          # False → load DATA_FILE
DATA_FILE     = "data.txt"    # one value per line (used when USE_SYNTHETIC=False)
N_BINS        = 60            # histogram bins
SEED          = 42

# True parameters for the synthetic dataset (used only for labelling)
TRUE = dict(A=4000, mu=3.0, sigma=0.5, tau=1.5)

# Dark theme palette
C = dict(
    bg      = "#0f1117",
    panel   = "#161b27",
    grid    = "#1e2533",
    text    = "#e2e8f0",
    muted   = "#64748b",
    hist    = "#334155",
    fit     = "#a78bfa",   # violet  – total fit
    rise    = "#38bdf8",   # sky     – Gaussian CDF rise envelope
    decay   = "#fb923c",   # amber   – exponential decay envelope
    resid   = "#34d399",   # emerald – residuals
    err     = "#f87171",   # red     – ±1σ band
)


# ══════════════════════════════════════════════════════════════════════════════
#  MODEL
# ══════════════════════════════════════════════════════════════════════════════
def gaussian_cdf(t, mu, sigma):
    """Smooth 0→1 rise (error-function CDF)."""
    return 0.5 * (1.0 + erf((t - mu) / (np.sqrt(2.0) * sigma)))


def exp_decay(t, mu, tau):
    """Exponential decay anchored at t = mu."""
    d = np.ones_like(t, dtype=float)
    mask = t >= mu
    d[mask] = np.exp(-(t[mask] - mu) / tau)
    return d


def model(t, A, mu, sigma, tau):
    """
    Full waveform: amplitude × Gaussian-CDF rise × Exponential decay.
    Parameters
    ----------
    A     : peak amplitude scale
    mu    : onset / peak position
    sigma : rise width  (Gaussian σ)
    tau   : decay time-constant
    """
    return A * gaussian_cdf(t, mu, sigma) * exp_decay(t, mu, tau)


# ══════════════════════════════════════════════════════════════════════════════
#  DATA  (synthetic or loaded)
# ══════════════════════════════════════════════════════════════════════════════
rng = np.random.default_rng(SEED)

if USE_SYNTHETIC:
    # Sample from the model by rejection sampling
    t_fine = np.linspace(0, 12, 5000)
    pdf    = model(t_fine, **TRUE)
    pdf   /= pdf.max()                       # normalise for rejection

    samples = []
    while len(samples) < 8000:
        ts = rng.uniform(t_fine[0], t_fine[-1], 20_000)
        us = rng.uniform(0, 1, 20_000)
        pdf_at_ts = model(ts, **TRUE)
        pdf_at_ts /= pdf_at_ts.max()
        samples.extend(ts[us < pdf_at_ts].tolist())

    data = np.array(samples[:8000])
    t_range = (0, 12)
else:
    data    = np.loadtxt(DATA_FILE)
    t_range = (data.min(), data.max())


# ══════════════════════════════════════════════════════════════════════════════
#  HISTOGRAM
# ══════════════════════════════════════════════════════════════════════════════
counts, edges = np.histogram(data, bins=N_BINS, range=t_range)
centres       = 0.5 * (edges[:-1] + edges[1:])
bin_width     = edges[1] - edges[0]
errors        = np.sqrt(np.maximum(counts, 1))   # Poisson √N errors


# ══════════════════════════════════════════════════════════════════════════════
#  FIT  (weighted least-squares via curve_fit)
# ══════════════════════════════════════════════════════════════════════════════
# Initial guesses — robust even for real data
peak_bin  = centres[np.argmax(counts)]
p0 = [counts.max() * 2, peak_bin * 0.8,
      (centres[-1] - centres[0]) * 0.05,
      (centres[-1] - peak_bin)   * 0.5]

bounds = (
    [0,   t_range[0], 1e-6, 1e-6],          # lower
    [np.inf, t_range[1], np.inf, np.inf],   # upper
)

popt, pcov = curve_fit(
    model, centres, counts,
    p0=p0, sigma=errors, absolute_sigma=True,
    bounds=bounds, maxfev=20_000
)
perr = np.sqrt(np.diag(pcov))
A_fit, mu_fit, sigma_fit, tau_fit = popt

# ── Goodness-of-fit ──────────────────────────────────────────────────────────
residuals = counts - model(centres, *popt)
chi2_val  = np.sum((residuals / errors) ** 2)
ndof      = len(counts) - len(popt)
chi2_red  = chi2_val / ndof
p_value   = 1.0 - chi2.cdf(chi2_val, ndof)


# ══════════════════════════════════════════════════════════════════════════════
#  PRINT RESULTS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "═"*52)
print("  FIT RESULTS")
print("═"*52)
labels = ["A (amplitude)", "μ (onset)", "σ (rise width)", "τ (decay const)"]
for lbl, val, err in zip(labels, popt, perr):
    print(f"  {lbl:<20s}  {val:>10.4f}  ±  {err:.4f}")
print(f"\n  χ² / ndof  =  {chi2_val:.1f} / {ndof}  =  {chi2_red:.3f}")
print(f"  p-value    =  {p_value:.4f}")
print("═"*52 + "\n")


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT
# ══════════════════════════════════════════════════════════════════════════════
t_fine = np.linspace(t_range[0], t_range[1], 2000)

fit_fine   = model(t_fine, *popt)
rise_env   = A_fit * gaussian_cdf(t_fine, mu_fit, sigma_fit)   # rise envelope
decay_env  = A_fit * exp_decay(t_fine, mu_fit, tau_fit)        # decay envelope

# ±1σ uncertainty band via parameter covariance
n_mc   = 500
mc_curves = []
for sample in rng.multivariate_normal(popt, pcov, n_mc):
    try:
        mc_curves.append(model(t_fine, *sample))
    except Exception:
        pass
mc_curves  = np.array(mc_curves)
band_lo    = np.percentile(mc_curves, 16, axis=0)
band_hi    = np.percentile(mc_curves, 84, axis=0)

fig = plt.figure(figsize=(12, 10), facecolor=C["bg"])
gs  = gridspec.GridSpec(
    2, 1, height_ratios=[3, 1], hspace=0.08,
    top=0.91, bottom=0.08, left=0.10, right=0.96
)


def style_ax(ax, ylabel, xlim):
    ax.set_facecolor(C["panel"])
    ax.tick_params(colors=C["text"], labelsize=9)
    for sp in ax.spines.values():
        sp.set_edgecolor("#2d3748")
    ax.set_ylabel(ylabel, color=C["text"], fontsize=10)
    ax.grid(color=C["grid"], linewidth=0.6, linestyle="--", alpha=0.8)
    ax.set_xlim(*xlim)


# ── Main panel ────────────────────────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0])
style_ax(ax1, "Counts / bin", t_range)
ax1.tick_params(labelbottom=False)

# Histogram bars
ax1.bar(centres, counts, width=bin_width * 0.88,
        color=C["hist"], alpha=0.75, label="Data histogram", zorder=2)
ax1.errorbar(centres, counts, yerr=errors,
             fmt="none", ecolor=C["muted"], elinewidth=0.8, capsize=2, zorder=3)

# ±1σ uncertainty band
ax1.fill_between(t_fine, band_lo, band_hi,
                 color=C["fit"], alpha=0.15, label="±1σ fit uncertainty")

# Component envelopes
ax1.plot(t_fine, rise_env,  color=C["rise"],  lw=1.4, ls="--",
         alpha=0.85, label=f"Rise envelope  (Gaussian CDF, σ={sigma_fit:.3f})")
ax1.plot(t_fine, decay_env, color=C["decay"], lw=1.4, ls="--",
         alpha=0.85, label=f"Decay envelope  (exp, τ={tau_fit:.3f})")

# Total fit
ax1.plot(t_fine, fit_fine, color=C["fit"], lw=2.5,
         label="Total fit  f(t) = A · Rise · Decay")

# Peak marker
peak_t = t_fine[np.argmax(fit_fine)]
peak_y = fit_fine.max()
ax1.plot(peak_t, peak_y, "o", color="white", ms=7, zorder=6)
ax1.annotate(f"peak  t={peak_t:.2f}",
             xy=(peak_t, peak_y),
             xytext=(peak_t + 0.4, peak_y * 0.92),
             arrowprops=dict(arrowstyle="->", color="white", lw=0.9),
             color="white", fontsize=8)

# Parameter box
box_txt = (
    f"A  = {A_fit:.1f} ± {perr[0]:.1f}\n"
    f"μ  = {mu_fit:.3f} ± {perr[1]:.3f}\n"
    f"σ  = {sigma_fit:.3f} ± {perr[2]:.3f}\n"
    f"τ  = {tau_fit:.3f} ± {perr[3]:.3f}\n"
    f"\nχ²/ndof = {chi2_red:.3f}\np-value  = {p_value:.4f}"
)
ax1.text(0.975, 0.97, box_txt,
         transform=ax1.transAxes, va="top", ha="right",
         fontsize=8.5, color=C["text"], family="monospace",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="#0f1117",
                   edgecolor="#2d3748", alpha=0.85))

ax1.legend(fontsize=8, facecolor=C["bg"], labelcolor=C["text"],
           framealpha=0.7, loc="upper right",
           bbox_to_anchor=(0.975, 0.60))

# ── Residual panel ────────────────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[1], sharex=ax1)
style_ax(ax2, "Pull  (σ)", t_range)

pulls = residuals / errors
ax2.bar(centres, pulls, width=bin_width * 0.88,
        color=[C["resid"] if p >= 0 else C["err"] for p in pulls],
        alpha=0.75)
ax2.axhline(0,  color=C["muted"], lw=1.0)
ax2.axhline(+2, color=C["muted"], lw=0.7, ls=":")
ax2.axhline(-2, color=C["muted"], lw=0.7, ls=":")
ax2.set_xlabel("t  (a.u.)", color=C["text"], fontsize=10)
ax2.set_ylim(-4.5, 4.5)

fig.suptitle(
    "Gaussian-CDF Rise  ×  Exponential Decay  —  Histogram Fit",
    color=C["text"], fontsize=13, fontweight="bold"
)

out_png = "/mnt/user-data/outputs/fit_gaussian_rise_exp_decay.png"
plt.savefig(out_png, dpi=150, bbox_inches="tight",
            facecolor=fig.get_facecolor())
plt.show()
print(f"Plot saved → {out_png}")