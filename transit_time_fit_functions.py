import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erfc


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
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Gaussian Rise + Exponential Decay Fit', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

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
    plt.savefig(plotFolder + '/gaussian_rise_exp_decay_fit.png', dpi=150)
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


if __name__ == "__main__":
    main()