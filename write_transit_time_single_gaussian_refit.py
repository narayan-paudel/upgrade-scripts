# !/usr/bin/env python
import glob
from subprocess import run

import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from transit_time_dependence_plots import (extract_json, extract_json_tts_filter, extract_json_min_tts_filter, extract_json_tts_chi2_filter,plot_single_transit_time_histogram,
                                           prod_id_to_icm_id,extract_channel,extract_fit_params)

from refit_bad_tt_fits import merge_bins_reduceat

import json

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})


# colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#1f78b4','#33a02c','#e31a1c','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#cab2d6','#ffff99','#6a3d9a','#b15928','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#1f78b4']


def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def get_single_gaussian_fit(mDOM_prod_id, channel, mdom_tt_dir) -> None:
    '''Plots the transit time histogram for a single mDOM and all its channels'''
    meas_files = glob.glob(f"{mdom_tt_dir}/{mDOM_prod_id}*/*.json")
    meas_files = [ifile for ifile in meas_files if extract_channel(ifile) == channel]
    remove_pedestal = None
    merge_bins = None
    A1, mu1, sigma1, reduced_chi2, RSS = np.nan, np.nan, np.nan, np.nan, np.nan
    if mDOM_prod_id in ["mDOM_D092"] and channel == 1:        
        merge_bins = 2
        print(f"{merge_bins} consecutive bins merged for {mDOM_prod_id} channel {channel}")
    if mDOM_prod_id in ["mDOM_D029"] and channel == 1 or mDOM_prod_id in ["mDOM_D034"] and channel == 21:
        remove_pedestal = 0
        print(f"pedestal {remove_pedestal} removed for {mDOM_prod_id} channel {channel}")

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
            initial_guess = [np.max(y_values), 55, 5]
            bounds = ([0, 45, 0.1], #lower limits
                       [np.max(y_values)+100, 65, 5])      # upper limits
            poisson_errors = np.sqrt(np.asarray(y_values))
            poisson_errors[poisson_errors == 0] = 1

            try:
                popt, pcov = curve_fit(gaussian,x_values,y_values,p0=initial_guess,bounds=bounds,sigma=poisson_errors,absolute_sigma=True)
                A1, mu1, sigma1 = popt
                y_fit = gaussian(x_values, A1, mu1, sigma1)
                residuals = y_values - y_fit
                RSS = np.sum(residuals**2)
                chi2 = np.sum(((y_values - gaussian(x_values, *popt))/poisson_errors)**2)
                ndof = len(x_values) - len(popt)
                reduced_chi2 = chi2 / ndof
                run = data["run_number"]
                # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f" tt {mu:.1f}\u00B1{sigma:.1f} ns PMT HV {hv:.1f} V {temperature:.1f} \u00b0C Run {run}",alpha=1)
                tt, tt_spread, a, b, c, chi2, hv = extract_fit_params(ifile)
            except RuntimeError:
                print(f"Gaussian fit did not converge for {mDOM_prod_id} channel {channel}")
                A1, mu1, sigma1, reduced_chi2, RSS = np.nan, np.nan, np.nan, np.nan, np.nan
    # print(f"mdom {mDOM_prod_id} channel {channel}: tt {tt:.1f} ns, tts {tt_spread:.1f} ns, a {a:.1f}, b {b:.1f}, c {c:.1f}, chi2 {chi2:.1f}, hv {hv:.1f} V")
    return A1, mu1, sigma1, reduced_chi2, RSS


def rewrite_json_unique_gaussian_refit_tt(transit_time_file,mdom_tt_dir ) -> list:
    '''
    Reads the json file with the list of DOMs and transit time measurements and returns a list of transit times
    and other data values
    for chi2, object key : "chi2"
    for transit time, object key:"mu"

    '''
    mdom_list = []
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            mdom_list.append(mdom)

    print(f"MDOMs: {len(mdom_list)}")
    channels = [int(ichannel) for ichannel in np.linspace(0,23,24)]
    mdom_dict = {}
    for mdom in mdom_list:
        pmt_dict = {}
        icm_id = data[mdom]["icm_id"]
        for channel in channels:
            tt_data = data[f"{mdom}"]
            if len(tt_data["transit_times"][f"channel_{channel}"])> 1:
                print(f"refitting {mdom} channel {channel} but has multiple runs {tt_data["transit_times"][f"channel_{channel}"]}")
            tt_info_list = tt_data["transit_times"][f"channel_{channel}"]
            if len(tt_info_list) > 0:
                tt_info = tt_data["transit_times"][f"channel_{channel}"][0]
                A1, mu1, sigma1, reduced_chi2, RSS = get_single_gaussian_fit(mdom, int(channel), mdom_tt_dir)
                tt_info["mu"] = mu1 - 12
                tt_info["sigma"] = sigma1
                tt_info["a"] = A1
                tt_info["b"] = mu1
                tt_info["c"] = sigma1
                tt_info["chi2"] = reduced_chi2
                tt_info["RSS"] = RSS
                tt_info_elt = [tt_info]
                pmt_dict[f"channel_{channel}"] = tt_info_elt
        mdom_dict[mdom] = {"icm_id": icm_id, "transit_times": pmt_dict}
        channels_list = [f"channel_{i}" for i in range(0,24)]
        pmt_dict_ordered = {channel: pmt_dict.get(channel, []) for channel in channels_list}                
        mdom_dict[mdom] = {"icm_id": icm_id, "transit_times": pmt_dict_ordered}
    with open('/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time_gaussian_refit.json', 'w') as f:
        json.dump(mdom_dict, f, indent=4)      






def main() -> None:
    upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"
    geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))
    mdom_tt_dir = home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    run_picks_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_run_picks.json"
    refit_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_tt_needing_refit.json"
    empty_meas_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_empty_meas.json"

    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    transit_times = extract_json(transit_time_file,obj_key="mu",filter_non_zero=False)
    transit_times_single_per_pmt_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time_single_pick.json"
    # rewrite_json_unique_tt(transit_time_file, mdom_tt_dir, run_picks_json, refit_json, empty_meas_json, filter_non_zero=False, check_outliers=[40,60])
    rewrite_json_unique_gaussian_refit_tt(transit_times_single_per_pmt_file, mdom_tt_dir)

if __name__ == "__main__":
    main()