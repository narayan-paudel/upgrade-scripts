# !/usr/bin/env python
import glob
from subprocess import run

import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from transit_time_dependence_plots import extract_json, extract_json_tts_filter, extract_json_min_tts_filter, extract_json_tts_chi2_filter,plot_single_transit_time_histogram
from transit_time_dependence_plots import prod_id_to_icm_id,extract_channel,extract_fit_params,gauss,get_pmt_uid

import json

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})


# colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
colorsCustom = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#1f78b4','#33a02c','#e31a1c','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#cab2d6','#ffff99','#6a3d9a','#b15928','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#1f78b4']

mis_fit_list = [{"mDOM":"mDOM_D007","channel":18,"run":22},{"mDOM":"mDOM_D015","channel":0,"run":2}]

def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def plot_single_transit_time_histogram_refit(mDOM_prod_id, channel,run, mdom_tt_dir, plotFolder,fit_line=False,fit_xlim=[-30,100],exclude_runs=[]) -> None:
    '''Plots the transit time histogram for a single mDOM and all its channels'''
    meas_files = glob.glob(f"{mdom_tt_dir}/{mDOM_prod_id}*/*_{run}.json")
    meas_files = [ifile for ifile in meas_files if extract_channel(ifile) == channel]
    print(f"meas files for {mDOM_prod_id} channel {channel} run {run}: {meas_files}")
    # print(f"meas files for {mDOM_prod_id} {meas_files}")
    # --- Initial guesses ---

    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for i,ifile in enumerate(meas_files):
        with open(ifile, 'r') as f:
            data = json.load(f)
            transit_times = []
            y_values = data["meas_data"][0]['y_values']
            x_min = data["meas_data"][0]["x_min"]
            x_max = data["meas_data"][0]["x_max"]
            n_bins = data["meas_data"][0]["n_bins"]
            x_values = np.linspace(x_min,x_max,n_bins)
            x_values_short = [ix for ix in x_values if ix > 55 and ix < 65]
            y_values_short = [iy for ix,iy in zip(x_values,y_values) if ix > 55 and ix < 65]
            x_label = data["meas_data"][0]["x_label"]
            y_label = data["meas_data"][0]["y_label"]
            A_init = np.max(y_values)
            mu_init = np.mean(55)
            sigma_init = np.std(4)
            popt, pcov = curve_fit(gaussian, x_values, y_values, p0=[A_init, mu_init, sigma_init])
            # popt, pcov = curve_fit(gaussian, x_values_short, y_values_short, p0=[A_init, mu_init, sigma_init])
            A_fit, mu_fit, sigma_fit = popt

            run = data["run_number"]
            if run in exclude_runs:
                continue
            # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f" tt {mu:.1f}\u00B1{sigma:.1f} ns PMT HV {hv:.1f} V {temperature:.1f} \u00b0C Run {run}",alpha=1)
            tt, tt_spread, a, b, c, chi2, hv = extract_fit_params(ifile)
            print(f"mdom {mDOM_prod_id} channel {channel} Run {run}: tt {tt:.1f} ns, tts {tt_spread:.1f} ns, a {a:.1f}, b {b:.1f}, c {c:.1f}, chi2 {chi2:.1f}, hv {hv:.1f} V")
            # if c < 20:
            # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f"Run {run}: tt {b:.0f} ns ({tt:.0f} ns) tts {tt_spread:.0f} ns chi2 {chi2:.1f}",alpha=1)
            if fit_line:
                print(f"Fitting Gaussian to {mDOM_prod_id} channel {channel} run {run} with fit parameters A {A_fit:.1f}, mu {mu_fit:.1f}, sigma {sigma_fit:.1f}")
                print(f"Fitting Gaussian to {mDOM_prod_id} channel {channel} run {run} with previous fit A {a:.1f}, mu {b:.1f}, sigma {c:.1f}")
                x_values_fit = np.linspace(fit_xlim[0],fit_xlim[1],1000)                
                # ax.plot(x_values_fit,gauss(x_values_fit,a,b,c),ls='--',lw = 2.5,c=colorsCustom[i],alpha=1)
                ax.plot(x_values_fit,gaussian(x_values_fit,A_fit,mu_fit,sigma_fit),ls='--',lw = 2.5,c=colorsCustom[i],alpha=1)

    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel} {get_pmt_uid(mDOM_prod_id, int(channel), mdom_tt_dir)}", transform=ax.transAxes, ha='left', fontsize=10)
    ax.set_xlabel(f"{x_label}", fontsize=22)
    ax.set_ylabel(f"{y_label}", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=8,ncols=1)
    plot_name = f"{plotFolder}/{mDOM_prod_id}_channel_{channel}_refit"
    print(plot_name)
    # plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()

def main() -> None:
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
    # transit_times
    # plot_single_transit_time_histogram_refit("mDOM_D007", 18, 22, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[])#Run 151 very off
    # plot_single_transit_time_histogram_refit("mDOM_D015", 0, 2, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[])#Run 151 very off
    plot_single_transit_time_histogram_refit("mDOM_D084", 20, 143, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[])#Run 151 very off


if __name__ == "__main__":
    main()