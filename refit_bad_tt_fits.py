# !/usr/bin/env python
import glob
from subprocess import run

import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from transit_time_dependence_plots import (extract_json, extract_json_tts_filter, extract_json_min_tts_filter, extract_json_tts_chi2_filter,plot_single_transit_time_histogram,
                                            prod_id_to_icm_id,extract_channel,extract_fit_params,gauss,get_pmt_uid)

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


def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
    return A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) + A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2))

def merge_bins_reduceat(counts, bin_edges, n=2):
    """Fast merging using reduceat."""
    indices = np.arange(0, len(counts), n)
    new_counts = np.add.reduceat(counts, indices)
    bin_edges = np.asarray(bin_edges)
    new_edges = bin_edges[indices]
    return new_counts, new_edges



def plot_refit_pmt(mdom_channel, mdom_tt_dir, plotFolder,fit_line=False,fit_xlim=[-30,100],ylim=None,exclude_runs=[],remove_pedestal=None,merge_bins=None) -> None:
    # print(f"refitting {ipmt}")
    mDOM_prod_id, channel = mdom_channel.split("_Ch")
    channel = int(channel)
    meas_files = glob.glob(f"{mdom_tt_dir}/{mDOM_prod_id}*/*.json")
    meas_files = [ifile for ifile in meas_files if extract_channel(ifile) == channel]
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
                # print(x_values)
                # print(y_values)
            initial_guess = [np.max(y_values), 55, 5]
            bounds = ([0, 45, 0.1], #lower limits
                    [np.max(y_values)+100, 65, 5])      # upper limits
            poisson_errors = np.sqrt(np.asarray(y_values))
            poisson_errors[poisson_errors == 0] = 1

            popt, pcov = curve_fit(gaussian,x_values,y_values,p0=initial_guess,bounds=bounds,sigma=poisson_errors,absolute_sigma=True)
            A1, mu1, sigma1 = popt

            chi2 = np.sum(((y_values - gaussian(x_values, *popt))/poisson_errors)**2)
            ndof = len(x_values) - len(popt)
            reduced_chi2 = chi2 / ndof

            run = data["run_number"]
            # if run in exclude_runs:
            #     continue
            # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f" tt {mu:.1f}\u00B1{sigma:.1f} ns PMT HV {hv:.1f} V {temperature:.1f} \u00b0C Run {run}",alpha=1)
            tt, tt_spread, a, b, c, chi2, hv = extract_fit_params(ifile)
            print(f"mdom {mDOM_prod_id} channel {channel} Run {run}: tt {tt:.1f} ns, tts {tt_spread:.1f} ns, a {a:.1f}, b {b:.1f}, c {c:.1f}, chi2 {chi2:.1f}, hv {hv:.1f} V")
            # if c < 20:
            # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f"Run {run}: tt {b:.0f} ns ({tt:.0f} ns) tts {tt_spread:.0f} ns chi2 {chi2:.1f}",alpha=1)
            # ax.plot(x_values,y_values,'o',lw = 2.5,c=colorsCustom[i],label=f"Run {run}: tt {b:.0f} ns ({tt:.0f} ns) tts {tt_spread:.0f} ns chi2 {chi2:.1f}",alpha=1)
            ax.errorbar(x_values,y_values,yerr=poisson_errors,fmt='o',lw = 2.5,c=colorsCustom[i],label=f"Run:{run}",alpha=1)
            if fit_line:
                x_values_fit = np.linspace(fit_xlim[0],fit_xlim[1],1000)                
                ax.plot(x_values_fit,gaussian(x_values_fit,*popt),ls='--',lw = 2.5,c=colorsCustom[i+1],label=f"refit A:{A1:.0f}, \u03bc:{mu1:.0f}, \u03C3:{sigma1 :.0f}, \u03C7\u00B2:{reduced_chi2:.1f}",alpha=1)
                ax.plot(x_values_fit,gaussian(x_values_fit,a,b,c),ls='--',lw = 2.5,c=colorsCustom[i+2],label=f"fit A:{a:.0f}, \u03bc:{b:.0f}, \u03C3:{c :.0f}, \u03C7\u00B2:{chi2:.1f}",alpha=1)

    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel}", transform=ax.transAxes, ha='left', fontsize=10)
    # ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel} {get_pmt_uid(mDOM_prod_id, int(channel), mdom_tt_dir)}", transform=ax.transAxes, ha='left', fontsize=10)
    ax.set_xlabel(f"{x_label}", fontsize=22)
    ax.set_ylabel(f"{y_label}", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=12,ncols=1)
    if ylim:
        ax.set_ylim(ylim)
    # ax.set_yscale('log')
    plot_name = f"{plotFolder}/{mDOM_prod_id}_channel_{channel}_badfit_refit_single_gauss"
    # plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()

def refit_bad_tt_fits(need_refits_file, transit_time_file, mdom_tt_dir, plotFolder,fit_line=False,fit_xlim=[-30,100],ylim=None,exclude_runs=[],remove_pedestal=None,merge_bins=None):
    with open(need_refits_file, 'r') as f:
        need_refits_data = json.load(f)
    need_refits_data_mdoms = [ielt["mDOM"] for ielt in need_refits_data]
    need_refits_data_mdoms_pmt = list(set([ielt["mDOM"]+"_Ch"+ielt["channel"] for ielt in need_refits_data]))
    print(f"there are {len(need_refits_data_mdoms_pmt)} PMTs needing refit: {need_refits_data_mdoms_pmt}")
    for ipmt in need_refits_data_mdoms_pmt[:]:
        plot_refit_pmt(ipmt, mdom_tt_dir, plotFolder,fit_line=fit_line,fit_xlim=fit_xlim,ylim=ylim,exclude_runs=exclude_runs,remove_pedestal=remove_pedestal,merge_bins=merge_bins)





def main():
    upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"
    geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))
    mdom_tt_dir = home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    run_picks_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_run_picks.json"
    refit_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_tt_needing_refit.json"
    empty_meas_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_empty_meas.json"

    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    # transit_times = extract_json(transit_time_file,obj_key="mu",filter_non_zero=False)
    # # transit_times_single_per_pmt = extract_json_unique_tt(transit_time_file,run_picks_json,refit_json,empty_meas_json,obj_key="mu",filter_non_zero=True)
    # transit_times_single_per_pmt_msu, sigma_list, chi2_list = extract_json_unique_tt(transit_time_file,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_M",filter_non_zero=False)
    # transit_times_single_per_pmt_desy, sigma_list, chi2_list = extract_json_unique_tt(transit_time_file,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_D",filter_non_zero=False)
    #for refitting all pmts needing refit
    refit_all = True
    if refit_all:
        refit_bad_tt_fits(refit_json, transit_time_file, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[30,80],ylim=None,exclude_runs=[])
    # plot_refit_pmt("mDOM_D036_Ch18", mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[30,80],ylim=None,exclude_runs=[])
    # plot_refit_pmt("mDOM_D036_Ch4", mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[30,80],ylim=None,exclude_runs=[])
    plot_refit_pmt("mDOM_M162_Ch13", mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[50,60],ylim=[0,5],exclude_runs=[])
    plot_refit_pmt("mDOM_D029_Ch1", mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],ylim=None,exclude_runs=[],remove_pedestal=0)
    plot_refit_pmt("mDOM_D034_Ch21", mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],ylim=None,exclude_runs=[],remove_pedestal=0)
    plot_refit_pmt("mDOM_D092_Ch1", mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],ylim=None,exclude_runs=[],remove_pedestal=None,merge_bins=2)
    #slightly off
    plot_refit_pmt("mDOM_D089_Ch18", mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],ylim=None,exclude_runs=[],remove_pedestal=None,merge_bins=None)
    plot_refit_pmt("mDOM_M182_Ch13", mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],ylim=None,exclude_runs=[],remove_pedestal=None,merge_bins=None)

if __name__ == "__main__":
    main()