import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf,erfc
from scipy.stats import chi2


from pathlib import Path
home = str(Path.home())


import glob
import json



from transit_time_dependence_plots import (extract_json, extract_json_tts_filter, extract_json_min_tts_filter, extract_json_tts_chi2_filter,plot_single_transit_time_histogram,
                                            prod_id_to_icm_id,extract_channel,extract_fit_params,gauss,get_pmt_uid,extract_run_number)


from refit_bad_tt_fits import merge_bins_reduceat
# -------------------------------------------------------------------
# Model: Gaussian-convolved rise * Exponential decay
# Physically: instrument response (Gaussian) convolved with
# an exponential decay — common in fluorescence lifetime, kinetics
# -------------------------------------------------------------------


plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})


# colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#1f78b4','#33a02c','#e31a1c','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#cab2d6','#ffff99','#6a3d9a','#b15928','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#1f78b4']

#Dark theme palattee
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




def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))


def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
    return A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) + A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2))


def triple_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3):
    return A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) + A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2)) + A3 * np.exp(-(x - mu3)**2 / (2 * sigma3**2))


def gaussian_rise_exp_decay(t, A, t0, sigma, tau):
    """
    Gaussian rise with exponential decay (analytically convolved).

    G(t) = A * exp(-((t - t0)^2) / (2*sigma^2))  [rise]
    E(t) = exp(-t / tau)                           [decay]
    Combined = convolution of both, evaluated analytically.

    Parameters:
        t     : time array
        A     : amplitude
        t0    : center/onset of Gaussian rise
        sigma : width of Gaussian (rise time)
        tau   : exponential decay time constant
    """
    # Analytical convolution of Gaussian IRF with exponential decay
    # Result: A/2 * exp(-t/tau + sigma^2/(2*tau^2)) * erfc(-(t-t0-sigma^2/tau) / (sqrt(2)*sigma))
    factor = np.exp(-((t - t0) / tau) + (sigma ** 2) / (2 * tau ** 2))
    erf_term = erfc(-((t - t0) - sigma ** 2 / tau) / (np.sqrt(2) * sigma))
    return (A / 2) * factor * erf_term


def gaussian_rise_exp_decay_simple(t, A, t0, sigma, tau, offset):
    """
    Simpler empirical version (not analytically convolved):
    Rise shaped by Gaussian CDF, decay by exponential.

    Useful when you want an intuitive parameterization.

    Parameters:
        t      : time array
        A      : amplitude
        t0     : onset time
        sigma  : rise width
        tau    : decay time constant
        offset : baseline offset
    """
    from scipy.special import erf
    rise = 0.5 * (1 + erf((t - t0) / (np.sqrt(2) * sigma)))  # Gaussian CDF
    decay = np.exp(-np.maximum(t - t0, 0) / tau)              # causal exponential
    return A * rise * decay + offset


def gaussian_rise(t,A,t0,sigma):
    """
    Simpler empirical version (not analytically convolved):
    Rise shaped by Gaussian CDF, decay by exponential.

    Useful when you want an intuitive parameterization.

    Parameters:
        t      : time array
        A      : amplitude
        t0     : onset time
        sigma  : rise width
        tau    : decay time constant
        offset : baseline offset
    """
    from scipy.special import erf 
    return 0.5 * (1 + erf((t - t0) / (np.sqrt(2) * sigma)))

def exp_decay_simple(t, A, t0, tau):
    """
    Simpler empirical version (not analytically convolved):
    Rise shaped by Gaussian CDF, decay by exponential.

    Useful when you want an intuitive parameterization.

    Parameters:
        t      : time array
        A      : amplitude
        t0     : onset time
        sigma  : rise width
        tau    : decay time constant
        offset : baseline offset
    """
    return A*np.exp(-np.maximum(t - t0, 0) / tau)


def gaussian_cdf(t, mu, sigma):
    """Smooth 0→1 rise (error-function CDF)."""
    return 0.5 * (1.0 + erf((t - mu) / (np.sqrt(2.0) * sigma)))
 
 
def exp_decay(t, mu, tau):
    """Exponential decay anchored at t = mu."""
    d = np.ones_like(t, dtype=float)
    mask = t >= mu
    d[mask] = np.exp(-(t[mask] - mu) / tau)
    return d
 
 
def gaussian_rise_exponential_decay(t, A, mu, sigma, tau):
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


def plot_gaussian_rise_exponential_decay_fit(x_values, data,plotFolder):
    x = np.array(x_values)
    data = np.array(data)
    peak_bin = x[np.argmax(data)]

    p0 = [data.max()*1.2, 55,3, 5]  # A, mu, sigma, tau

    # Bounds: all positive, t0 within data range
    bounds = (
        [data.max()*0.3, 45, 1,- np.inf],   # lower
        [data.max()*1.2,  65, 5, np.inf],  # upper
    )
    poisson_errors = np.sqrt(np.maximum(np.asarray(data),1))

    popt, pcov = curve_fit(
        gaussian_rise_exponential_decay,
        x, data,
        p0=p0,
        bounds=bounds,
        maxfev=10000,
        sigma=poisson_errors,
        absolute_sigma=True
    )

    # popt, pcov = curve_fit(
    #     gaussian_rise_exp_decay,
    #     p0=p0,
    #     maxfev=10000,
    # )

    perr = np.sqrt(np.diag(pcov))  # 1-sigma uncertainties

    A_fit, t0_fit, sigma_fit, tau_fit = popt

    ############################
    perr = np.sqrt(np.diag(pcov))
    perr_dof = perr / np.sqrt(len(x_values) - len(popt))
    chi2_val = np.sum(((data - gaussian_rise_exponential_decay(x_values, *popt))/poisson_errors)**2)
    ndof = len(x_values) - len(popt)
    reduced_chi2 = chi2_val / ndof
    p_value   = 1.0 - chi2.cdf(chi2_val, ndof)


    # -------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------
    y_fit = gaussian_rise_exponential_decay(x_values, *popt)
    residuals = data - y_fit
    RSS = np.sum(residuals**2)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Main plot
    ax1.errorbar(x_values, data, yerr=poisson_errors, fmt='o', color=colorsCustom[0], label=f"transit time")
    fit_x_values = np.linspace(x_values.min(), x_values.max(), 400)
    ax1.plot(fit_x_values, gaussian_rise_exponential_decay(fit_x_values, *popt),  '-',c=colorsCustom[1],  linewidth=2,   label=f'gauss rise, expon decay \nRSS = {RSS:.0f}, χ² = {reduced_chi2:.1f}')
    ax1.plot(fit_x_values, A_fit * gaussian_cdf(fit_x_values, t0_fit, sigma_fit),  '-',c=colorsCustom[2],  linewidth=2,   label=rf"A:{A_fit:.0f}, $\mu$:{t0_fit:.0f}, $\sigma$:{sigma_fit:.0f}")
    ax1.plot(fit_x_values, A_fit * exp_decay(fit_x_values, t0_fit, tau_fit),  '-',c=colorsCustom[3],  linewidth=2,   label=rf"A:{A_fit:.0f}, $\mu$:{t0_fit:.0f}, $\tau$:{tau_fit:.0f}")
    ax1.set_ylabel('Count', fontsize=22)
    ax1.set_title('Gaussian Rise Exponential Decay Fit', fontsize=22)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax1.set_ylabel('Count', fontsize=22)
    ax1.set_title('Gaussian Rise Exponential Decay Fit', fontsize=22)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # Annotate fitted params
    # param_text = (f"A = {A_fit:.3f} ± {perr[0]:.3f}\n"
    #             f"t₀ = {t0_fit:.3f} ± {perr[1]:.3f}\n"
    #             f"σ = {sigma_fit:.3f} ± {perr[2]:.3f}\n"
    #             f"τ = {tau_fit:.3f} ± {perr[3]:.3f}")
    # ax1.text(0.97, 0.95, param_text, transform=ax1.transAxes,
    #         fontsize=10, va='top', ha='right',
    #         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Residuals
    ax2.errorbar(x_values, residuals, yerr=poisson_errors, fmt='o', color=colorsCustom[0], alpha=1)
    ax2.axhline(0, color=colorsCustom[1], linewidth=2)
    ax2.set_xlabel('Transit times [ns]', fontsize=22)
    ax2.set_ylabel('Residuals', fontsize=22)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both',which='both', direction='in', labelsize=22)
    plt.tight_layout()
    plt.savefig(plotFolder + '/gaussian_rise_exponential_decay_residuals.pdf', dpi=150)
    plt.close()






def get_transit_time_data(mDOM_prod_id, channel,mdom_tt_dir,run_picks_json,need_refits_json,empty_meas_json,site_str="_M",filter_non_zero=False,check_outliers=None):
    refit_tt = []
    refit_tt2 = []
    original_tt = []
    chi2_values = []
    chi2_refit = []
    tt_diff = []
    chi2_diff = []
    with open(run_picks_json, 'r') as f:
        run_picks_data = json.load(f)
    with open(need_refits_json, 'r') as f:
        need_refits_data = json.load(f)
    with open(empty_meas_json, 'r') as f:
        empties_data = json.load(f)
    run_picks_data_mdoms = [ielt["mDOM"] for ielt in run_picks_data]
    need_refits_data_mdoms = [ielt["mDOM"] for ielt in need_refits_data]
    need_refits_data_mdoms_pmt = [ielt["mDOM"]+"_"+ielt["channel"] for ielt in need_refits_data]
    empties_data_mdoms = [ielt["mDOM"] for ielt in empties_data]
    empties_data_mdoms_pmt = [ielt["mDOM"]+"_"+ielt["channel"] for ielt in empties_data]

    meas_files = glob.glob(f"{mdom_tt_dir}/{mDOM_prod_id}*/*.json")
    available_channels = [extract_channel(ifile) for ifile in meas_files]
    # print(f"available channels for {mDOM_prod_id}: {list(set(available_channels))}")
    meas_files = [ifile for ifile in meas_files if extract_channel(ifile) == channel]
    # print(f"meas files for {mDOM_prod_id} {meas_files}")
    remove_pedestal = None
    merge_bins = None
    if mDOM_prod_id in ["mDOM_D092"] and channel == 1:        
        merge_bins = 2
        print(f"{merge_bins} consecutive bins merged for {mDOM_prod_id} channel {channel}")
    if mDOM_prod_id in ["mDOM_D029"] and channel == 1 or mDOM_prod_id in ["mDOM_D034"] and channel == 21:
        remove_pedestal = 0
        print(f"pedestal {remove_pedestal} removed for {mDOM_prod_id} channel {channel}")
    if mDOM_prod_id in run_picks_data_mdoms:
        select_mdom_run_data = [ielt for ielt in run_picks_data if ielt["mDOM"] == mDOM_prod_id][0]
        select_channel_run_data = [ielt for ielt in select_mdom_run_data["select_run"] if int(ielt["channel"]) == int(channel)][0]
        select_run = select_channel_run_data["run"]
        meas_files = [ifile for ifile in meas_files if int(extract_run_number(ifile)) == int(select_run)]
    for i,ifile in enumerate(meas_files):
        with open(ifile, 'r') as f:
            data = json.load(f)
            transit_times = []
            y_values = data["meas_data"][0]['y_values']
            x_min = data["meas_data"][0]["x_min"]
            x_max = data["meas_data"][0]["x_max"]
            n_bins = data["meas_data"][0]["n_bins"]
            x_values = np.linspace(x_min,x_max,n_bins)
            x_label = data["meas_data"][0]["x_label"]
            y_label = data["meas_data"][0]["y_label"]
            if remove_pedestal is not None:
                x_values_wo_ped = []
                y_values_wo_ped = []
                for ix,iy in zip(x_values,y_values):
                    if abs(iy - remove_pedestal) > 0.001:
                        x_values_wo_ped.append(ix)
                        y_values_wo_ped.append(iy)

                x_values = x_values_wo_ped
                y_values = y_values_wo_ped
                # print(x_values)
            if merge_bins is not None:
                y_values, x_values = merge_bins_reduceat(y_values, x_values, n=merge_bins)
    return x_values, y_values



def plot_model_fit(x, data,plotFolder):

    # t = np.linspace(-2, 20, 500)
    # true_params = dict(A=3.0, t0=2.0, sigma=0.8, tau=4.0)

    # y_true = gaussian_rise_exp_decay(t, **true_params)
    # y_noise = y_true + np.random.normal(0, 0.05, size=len(t))


    # -------------------------------------------------------------------
    # Fit
    # -------------------------------------------------------------------
    # Initial guesses: [A, t0, sigma, tau]
    p0 = [300.0, 55.0, 1.0, 2.0]

    # Bounds: all positive, t0 within data range
    bounds = (
        [0,    45, 1, -500.0],   # lower
        [1e6,  65, 5, 500.0],  # upper
    )
    poisson_errors = np.sqrt(np.asarray(data))
    poisson_errors[poisson_errors == 0] = 1

    popt, pcov = curve_fit(
        gaussian_rise_exp_decay,
        x, data,
        p0=p0,
        bounds=bounds,
        maxfev=10000,
        sigma=poisson_errors,
        absolute_sigma=True
    )

    # popt, pcov = curve_fit(
    #     gaussian_rise_exp_decay,
    #     x, data,
    #     p0=p0,
    #     maxfev=10000,
    # )

    perr = np.sqrt(np.diag(pcov))  # 1-sigma uncertainties

    A_fit, t0_fit, sigma_fit, tau_fit = popt


    # print("=" * 45)
    # print(f"{'Parameter':<10} {'True':>8} {'Fitted':>10} {'±1σ':>10}")
    # print("-" * 45)
    # print(f"{'A':<10} {true_params['A']:>8.3f} {A_fit:>10.3f} {perr[0]:>10.3f}")
    # print(f"{'t0':<10} {true_params['t0']:>8.3f} {t0_fit:>10.3f} {perr[1]:>10.3f}")
    # print(f"{'sigma':<10} {true_params['sigma']:>8.3f} {sigma_fit:>10.3f} {perr[2]:>10.3f}")
    # print(f"{'tau':<10} {true_params['tau']:>8.3f} {tau_fit:>10.3f} {perr[3]:>10.3f}")
    # print("=" * 45)


    # -------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------
    y_fit = gaussian_rise_exp_decay(x, *popt)
    residuals = data - y_fit

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Main plot
    ax1.errorbar(x, data, yerr=poisson_errors, fmt='o', color='steelblue', alpha=0.5, label="transit time")
    ax1.plot(x, y_fit,  'r-',  linewidth=2,   label='Fit')
    ax1.plot(x, gaussian_rise(x, A_fit, t0_fit, sigma_fit),  'g-',  linewidth=2,   label='Gaussian Rise')
    ax1.plot(x, exp_decay_simple(x, A_fit, t0_fit, tau_fit),  'b-',  linewidth=2,   label='Exponential Decay')
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Gaussian Rise + Exponential Decay Fit', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0,400)

    # Annotate fitted params
    param_text = (f"A = {A_fit:.3f} ± {perr[0]:.3f}\n"
                f"t₀ = {t0_fit:.3f} ± {perr[1]:.3f}\n"
                f"σ = {sigma_fit:.3f} ± {perr[2]:.3f}\n"
                f"τ = {tau_fit:.3f} ± {perr[3]:.3f}")
    ax1.text(0.97, 0.95, param_text, transform=ax1.transAxes,
            fontsize=10, va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Residuals
    ax2.scatter(x, residuals, s=4, color='gray', alpha=0.5)
    ax2.axhline(0, color='red', linewidth=1)
    ax2.set_xlabel('Time', fontsize=12)
    ax2.set_ylabel('Residuals', fontsize=12)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(plotFolder + '/gaussian_rise_exp_decay_fit.pdf', dpi=150)
    plt.close()



def plot_triple_gaussian_fit_residual(x_values, data,plotFolder):

    initial_guess = [np.max(data), 55, 3, np.max(data)/1.5, 55, 3,np.max(data)/1.5, 55, 3]
    # bounds = ([0, 30, 0.1,0, 30, 0.1], #lower limits
    #            [np.max(y_values)+100, 80, 6,np.max(y_values)+100, 80, 6])      # upper limits
    bounds = ([0.5*np.max(data), 45, 0.1,0, 45, 0.1,0, 45, 0.1], #lower limits
                [np.max(data)+100, 65, 5,np.max(data)+100, 65, 5,np.max(data)+100, 65, 5])      # upper limits


    poisson_errors = np.sqrt(np.asarray(data))
    poisson_errors[poisson_errors == 0] = 1

    popt, pcov = curve_fit(
        triple_gaussian,
        x_values, data,
        p0=initial_guess,
        bounds=bounds,
        maxfev=10000,
        sigma=poisson_errors,
        absolute_sigma=True
    )

    perr = np.sqrt(np.diag(pcov))
    perr_dof = perr / np.sqrt(len(x_values) - len(popt))
    chi2 = np.sum(((data - triple_gaussian(x_values, *popt))/poisson_errors)**2)
    ndof = len(x_values) - len(popt)
    reduced_chi2 = chi2 / ndof
    A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3 = popt

    # -------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------
    y_fit = triple_gaussian(x_values, *popt)
    residuals = data - y_fit
    RSS = np.sum(residuals**2)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Main plot
    ax1.errorbar(x_values, data, yerr=poisson_errors, fmt='o', color=colorsCustom[0], label=f"transit time")
    fit_x_values = np.linspace(x_values.min(), x_values.max(), 400)
    ax1.plot(fit_x_values, triple_gaussian(fit_x_values, *popt),  '-',c=colorsCustom[1],  linewidth=2,   label=f'triple Gaussian Fit\nRSS = {RSS:.0f}, χ² = {reduced_chi2:.1f}')
    ax1.plot(fit_x_values, gaussian(fit_x_values, A1, mu1, sigma1),  '-',c=colorsCustom[2],  linewidth=2,   label=rf"A$_{{{1}}}$:{A1:.0f}, $\mu_{{{1}}}$:{mu1:.0f}, $\sigma_{{{1}}}$:{sigma1:.0f}")
    ax1.plot(fit_x_values, gaussian(fit_x_values, A2, mu2, sigma2),  '-',c=colorsCustom[3],  linewidth=2,   label=rf"A$_{{{2}}}$:{A2:.0f}, $\mu_{{{2}}}$:{mu2:.0f}, $\sigma_{{{2}}}$:{sigma2:.0f}")
    ax1.plot(fit_x_values, gaussian(fit_x_values, A3, mu3, sigma3),  '-',c=colorsCustom[4],  linewidth=2,   label=rf"A$_{{{3}}}$:{A3:.0f}, $\mu_{{{3}}}$:{mu3:.0f}, $\sigma_{{{3}}}$:{sigma3:.0f}")
    ax1.set_ylabel('Count', fontsize=22)
    ax1.set_title('triple gaussian Fit', fontsize=22)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # Annotate fitted params
    # param_text = (f"A = {A_fit:.3f} ± {perr[0]:.3f}\n"
    #             f"t₀ = {t0_fit:.3f} ± {perr[1]:.3f}\n"
    #             f"σ = {sigma_fit:.3f} ± {perr[2]:.3f}\n"
    #             f"τ = {tau_fit:.3f} ± {perr[3]:.3f}")
    # ax1.text(0.97, 0.95, param_text, transform=ax1.transAxes,
    #         fontsize=10, va='top', ha='right',
    #         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Residuals
    ax2.errorbar(x_values, residuals, yerr=poisson_errors, fmt='o', color=colorsCustom[0], alpha=1)
    ax2.axhline(0, color=colorsCustom[1], linewidth=2)
    ax2.set_xlabel('Transit times [ns]', fontsize=22)
    ax2.set_ylabel('Residuals', fontsize=22)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both',which='both', direction='in', labelsize=22)
    plt.tight_layout()
    plt.savefig(plotFolder + '/triple_gaussian_fit_residuals.pdf', dpi=150)
    plt.close()

def plot_double_gaussian_fit_residual(x_values, data,plotFolder):

    initial_guess = [np.max(data), 55, 3, np.max(data)/1.5, 55, 3]
    # bounds = ([0, 30, 0.1,0, 30, 0.1], #lower limits
    #            [np.max(y_values)+100, 80, 6,np.max(y_values)+100, 80, 6])      # upper limits
    bounds = ([0.5*np.max(data), 45, 0.1,0, 45, 0.1], #lower limits
                [np.max(data)+100, 65, 5,np.max(data)+100, 65, 5])      # upper limits


    poisson_errors = np.sqrt(np.asarray(data))
    poisson_errors[poisson_errors == 0] = 1

    popt, pcov = curve_fit(
        double_gaussian,
        x_values, data,
        p0=initial_guess,
        bounds=bounds,
        maxfev=10000,
        sigma=poisson_errors,
        absolute_sigma=True
    )

    perr = np.sqrt(np.diag(pcov))
    perr_dof = perr / np.sqrt(len(x_values) - len(popt))
    chi2 = np.sum(((data - double_gaussian(x_values, *popt))/poisson_errors)**2)
    ndof = len(x_values) - len(popt)
    reduced_chi2 = chi2 / ndof
    A1, mu1, sigma1, A2, mu2, sigma2 = popt

    # -------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------
    y_fit = double_gaussian(x_values, *popt)
    residuals = data - y_fit
    RSS = np.sum(residuals**2)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Main plot
    ax1.errorbar(x_values, data, yerr=poisson_errors, fmt='o', color=colorsCustom[0], label=f"transit time")
    fit_x_values = np.linspace(x_values.min(), x_values.max(), 400)
    ax1.plot(fit_x_values, double_gaussian(fit_x_values, *popt),  '-',c=colorsCustom[1],  linewidth=2,   label=f'double Gaussian Fit\nRSS = {RSS:.0f}, χ² = {reduced_chi2:.1f}')
    ax1.plot(fit_x_values, gaussian(fit_x_values, A1, mu1, sigma1),  '-',c=colorsCustom[2],  linewidth=2,   label=rf"A$_{{{1}}}$:{A1:.0f}, $\mu_{{{1}}}$:{mu1:.0f}, $\sigma_{{{1}}}$:{sigma1:.0f}")
    ax1.plot(fit_x_values, gaussian(fit_x_values, A2, mu2, sigma2),  '-',c=colorsCustom[3],  linewidth=2,   label=rf"A$_{{{2}}}$:{A2:.0f}, $\mu_{{{2}}}$:{mu2:.0f}, $\sigma_{{{2}}}$:{sigma2:.0f}")
    ax1.set_ylabel('Count', fontsize=22)
    ax1.set_title('double gaussian Fit', fontsize=22)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # Annotate fitted params
    # param_text = (f"A = {A_fit:.3f} ± {perr[0]:.3f}\n"
    #             f"t₀ = {t0_fit:.3f} ± {perr[1]:.3f}\n"
    #             f"σ = {sigma_fit:.3f} ± {perr[2]:.3f}\n"
    #             f"τ = {tau_fit:.3f} ± {perr[3]:.3f}")
    # ax1.text(0.97, 0.95, param_text, transform=ax1.transAxes,
    #         fontsize=10, va='top', ha='right',
    #         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Residuals
    ax2.errorbar(x_values, residuals, yerr=poisson_errors, fmt='o', color=colorsCustom[0], alpha=1)
    ax2.axhline(0, color=colorsCustom[1], linewidth=2)
    ax2.set_xlabel('Transit times [ns]', fontsize=22)
    ax2.set_ylabel('Residuals', fontsize=22)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both',which='both', direction='in', labelsize=22)
    plt.tight_layout()
    plt.savefig(plotFolder + '/double_gaussian_fit_residual.pdf', dpi=150)
    plt.close()


def plot_single_gaussian_fit_residual(x_values, data,plotFolder):

    initial_guess = [np.max(data), 55, 3]
    # bounds = ([0, 30, 0.1,0, 30, 0.1], #lower limits
    #            [np.max(y_values)+100, 80, 6,np.max(y_values)+100, 80, 6])      # upper limits
    bounds = ([0.5*np.max(data), 45, 0.1], #lower limits
                [np.max(data)+100, 65, 5])      # upper limits


    poisson_errors = np.sqrt(np.asarray(data))
    poisson_errors[poisson_errors == 0] = 1

    popt, pcov = curve_fit(
        gaussian,
        x_values, data,
        p0=initial_guess,
        bounds=bounds,
        maxfev=10000,
        sigma=poisson_errors,
        absolute_sigma=True
    )

    perr = np.sqrt(np.diag(pcov))
    perr_dof = perr / np.sqrt(len(x_values) - len(popt))
    chi2 = np.sum(((data - gaussian(x_values, *popt))/poisson_errors)**2)
    ndof = len(x_values) - len(popt)
    reduced_chi2 = chi2 / ndof
    A1, mu1, sigma1 = popt

    # -------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------
    y_fit = gaussian(x_values, *popt)
    residuals = data - y_fit
    RSS = np.sum(residuals**2)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # Main plot
    ax1.errorbar(x_values, data, yerr=poisson_errors, fmt='o', color=colorsCustom[0], label=f"transit time")
    fit_x_values = np.linspace(x_values.min(), x_values.max(), 400)
    ax1.plot(fit_x_values, gaussian(fit_x_values, *popt),  '-',c=colorsCustom[1],  linewidth=2,   label=rf'Gaussian Fit A$_{{{1}}}$:{A1:.0f}, $\mu_{{{1}}}$:{mu1:.0f}, $\sigma_{{{1}}}$:{sigma1:.0f}'+'\n'+f'RSS = {RSS:.0f}, χ² = {reduced_chi2:.1f}')
    ax1.set_ylabel('Count', fontsize=22)
    ax1.set_title('gaussian Fit', fontsize=22)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # Annotate fitted params
    # param_text = (f"A = {A_fit:.3f} ± {perr[0]:.3f}\n"
    #             f"t₀ = {t0_fit:.3f} ± {perr[1]:.3f}\n"
    #             f"σ = {sigma_fit:.3f} ± {perr[2]:.3f}\n"
    #             f"τ = {tau_fit:.3f} ± {perr[3]:.3f}")
    # ax1.text(0.97, 0.95, param_text, transform=ax1.transAxes,
    #         fontsize=10, va='top', ha='right',
    #         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Residuals
    ax2.errorbar(x_values, residuals, yerr=poisson_errors, fmt='o', color=colorsCustom[0], alpha=1)
    ax2.axhline(0, color=colorsCustom[1], linewidth=2)
    ax2.set_xlabel('Transit times [ns]', fontsize=22)
    ax2.set_ylabel('Residuals', fontsize=22)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both',which='both', direction='in', labelsize=22)
    plt.tight_layout()
    plt.savefig(plotFolder + '/gaussian_fit_residuals.pdf', dpi=150)
    plt.close()




def main():
    upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"
    geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))
    mdom_tt_dir = home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    ################################################################################
    #main file
    # transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    #selected files
    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time_select.json"
    ################################################################################
    ##########run picks file###########
    run_picks_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_run_picks.json"
    refit_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_tt_needing_refit.json"
    empty_meas_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_empty_meas.json"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    # -------------------------------------------------------------------
    # Generate synthetic data
    # -------------------------------------------------------------------
    np.random.seed(42)
    x_values, y_values = get_transit_time_data("mDOM_M055", 13,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_M",filter_non_zero=False,check_outliers=None)
    plot_model_fit(x_values, y_values,plotFolder)
    x_values, y_values = get_transit_time_data("mDOM_M134", 5,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_M",filter_non_zero=False,check_outliers=None)
    plot_model_fit(x_values, y_values,plotFolder)
    x_values, y_values = get_transit_time_data("mDOM_M083", 13,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_M",filter_non_zero=False,check_outliers=None)
    plot_model_fit(x_values, y_values,plotFolder)


if __name__ == "__main__":
    main()