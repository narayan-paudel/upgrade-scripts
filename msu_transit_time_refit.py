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

def plot_single_transit_time_histogram_refit_double_gaussian(mDOM_prod_id, channel, mdom_tt_dir, plotFolder,fit_line=False,fit_xlim=[-30,100],exclude_runs=[]) -> None:
    '''Plots the transit time histogram for a single mDOM and all its channels'''
    meas_files = glob.glob(f"{mdom_tt_dir}/{mDOM_prod_id}*/*.json")
    available_channels = [extract_channel(ifile) for ifile in meas_files]
    print(f"available channels for {mDOM_prod_id}: {list(set(available_channels))}")
    meas_files = [ifile for ifile in meas_files if extract_channel(ifile) == channel]
    print(f"meas files for {mDOM_prod_id} {meas_files}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    x_label, y_label = None, None
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
            initial_guess = [np.max(y_values), 55, 3, np.max(y_values)/1.5, 55, 3]
            # bounds = ([0, 30, 0.1,0, 30, 0.1], #lower limits
            #            [np.max(y_values)+100, 80, 6,np.max(y_values)+100, 80, 6])      # upper limits
            bounds = ([0, 45, 0.1,0, 45, 0.1], #lower limits
                       [np.max(y_values)+100, 65, 5,np.max(y_values)+100, 65, 5])      # upper limits
            
            poisson_errors = np.sqrt(np.asarray(y_values))
            poisson_errors[poisson_errors == 0] = 1
            try:
                popt, pcov = curve_fit(double_gaussian,x_values,y_values,p0=initial_guess,bounds=bounds,sigma=poisson_errors,absolute_sigma=True)
                perr = np.sqrt(np.diag(pcov))
                perr_dof = perr / np.sqrt(len(x_values) - len(popt))
                chi2 = np.sum(((y_values - double_gaussian(x_values, *popt))/poisson_errors)**2)
                ndof = len(x_values) - len(popt)
                reduced_chi2 = chi2 / ndof
                A1, mu1, sigma1, A2, mu2, sigma2 = popt

                run = data["run_number"]
                if run in exclude_runs:
                    continue
                # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f" tt {mu:.1f}\u00B1{sigma:.1f} ns PMT HV {hv:.1f} V {temperature:.1f} \u00b0C Run {run}",alpha=1)
                tt, tt_spread, a, b, c, chi2, hv = extract_fit_params(ifile)
                print(f"mdom {mDOM_prod_id} channel {channel} Run {run}: tt {tt:.1f} ns, tts {tt_spread:.1f} ns, a {a:.1f}, b {b:.1f}, c {c:.1f}, chi2 {chi2:.1f}, hv {hv:.1f} V")
                # if c < 20:
                # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f"Run {run}: tt {b:.0f} ns ({tt:.0f} ns) tts {tt_spread:.0f} ns chi2 {chi2:.1f}",alpha=1)
                # ax.plot(x_values,y_values,'o',lw = 2.5,c=colorsCustom[i],label=f"Run {run}: tt {b:.0f} ns ({tt:.0f} ns) tts {tt_spread:.0f} ns chi2 {chi2:.1f}",alpha=1)
                ax.errorbar(x_values,y_values,yerr=poisson_errors,fmt='o',lw = 2.5,c=colorsCustom[i],label=f"Run:{run}\ntt {b:.0f} ns, tts {tt_spread:.0f} ns, \u03C7\u00B2 {chi2:.1f}",alpha=1)
                if fit_line:
                    x_values_fit = np.linspace(fit_xlim[0],fit_xlim[1],1000)                
                    ax.plot(x_values_fit,double_gaussian(x_values_fit,*popt),ls='--',lw = 2.5,c=colorsCustom[i+1],label=f"double gaussian, \u03C7\u00B2:{reduced_chi2:.1f}",alpha=1)
                    ax.plot(x_values_fit,gaussian(x_values_fit,A1,mu1,sigma1),ls='--',lw = 2.5,c=colorsCustom[i+2],label=rf"A$_{{{1}}}$:{A1:.0f}, $\mu_{{{1}}}$:{mu1:.0f}, $\sigma_{{{1}}}$:{sigma1:.0f}",alpha=1)
                    ax.plot(x_values_fit,gaussian(x_values_fit,A2,mu2,sigma2),ls='--',lw = 2.5,c=colorsCustom[i+3],label=rf"A$_{{{2}}}$:{A2:.0f}, $\mu_{{{2}}}$:{mu2:.0f}, $\sigma_{{{2}}}$:{sigma2:.0f}",alpha=1)

            except RuntimeError:
                print(f"Failed to fit double gaussian for {mDOM_prod_id} channel {channel}")
                continue

    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel}", transform=ax.transAxes, ha='left', fontsize=10)
    # ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel} {get_pmt_uid(mDOM_prod_id, int(channel), mdom_tt_dir)}", transform=ax.transAxes, ha='left', fontsize=10)
    ax.set_xlabel(f"{x_label}", fontsize=22)
    ax.set_ylabel(f"{y_label}", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=12)
    # ax.set_yscale('log')
    plot_name = f"{plotFolder}/{mDOM_prod_id}_channel_{channel}_refit"
    # plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_single_transit_time_histogram_refit_single_gaussian(mDOM_prod_id, channel, mdom_tt_dir, plotFolder,fit_line=False,fit_xlim=[-30,100],exclude_runs=[]) -> None:
    '''Plots the transit time histogram for a single mDOM and all its channels'''
    meas_files = glob.glob(f"{mdom_tt_dir}/{mDOM_prod_id}*/*.json")
    available_channels = [extract_channel(ifile) for ifile in meas_files]
    print(f"available channels for {mDOM_prod_id}: {list(set(available_channels))}")
    meas_files = [ifile for ifile in meas_files if extract_channel(ifile) == channel]
    print(f"meas files for {mDOM_prod_id} {meas_files}")
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
            if run in exclude_runs:
                continue
            # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f" tt {mu:.1f}\u00B1{sigma:.1f} ns PMT HV {hv:.1f} V {temperature:.1f} \u00b0C Run {run}",alpha=1)
            tt, tt_spread, a, b, c, chi2, hv = extract_fit_params(ifile)
            print(f"mdom {mDOM_prod_id} channel {channel} Run {run}: tt {tt:.1f} ns, tts {tt_spread:.1f} ns, a {a:.1f}, b {b:.1f}, c {c:.1f}, chi2 {chi2:.1f}, hv {hv:.1f} V")
            # if c < 20:
            # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f"Run {run}: tt {b:.0f} ns ({tt:.0f} ns) tts {tt_spread:.0f} ns chi2 {chi2:.1f}",alpha=1)
            # ax.plot(x_values,y_values,'o',lw = 2.5,c=colorsCustom[i],label=f"Run {run}: tt {b:.0f} ns ({tt:.0f} ns) tts {tt_spread:.0f} ns chi2 {chi2:.1f}",alpha=1)
            ax.errorbar(x_values,y_values,yerr=poisson_errors,fmt='o',lw = 2.5,c=colorsCustom[i],label=f"Run:{run}\ntt:{b:.0f} ns, tts:{tt_spread:.0f} ns, \u03C7\u00B2:{chi2:.1f}",alpha=1)
            if fit_line:
                x_values_fit = np.linspace(fit_xlim[0],fit_xlim[1],1000)                
                ax.plot(x_values_fit,gaussian(x_values_fit,*popt),ls='--',lw = 2.5,c=colorsCustom[i+1],label=f"A:{A1:.0f}, \u03bc:{mu1:.0f}, \u03C3:{sigma1 :.0f}, \u03C7\u00B2:{reduced_chi2:.1f}",alpha=1)

    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel}", transform=ax.transAxes, ha='left', fontsize=10)
    # ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel} {get_pmt_uid(mDOM_prod_id, int(channel), mdom_tt_dir)}", transform=ax.transAxes, ha='left', fontsize=10)
    ax.set_xlabel(f"{x_label}", fontsize=22)
    ax.set_ylabel(f"{y_label}", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=12,ncols=1)
    # ax.set_yscale('log')
    plot_name = f"{plotFolder}/{mDOM_prod_id}_channel_{channel}_refit_single_gauss"
    # plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()

def extract_json_mdom_pmt_list(transit_time_file, site) -> list:
    '''
    Reads the json file with the list of DOMs and transit time measurements and returns a list of (mdom, pmts)

    '''
    pmt_all = []
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            # print(tt_data.keys())
            for ichannel in tt_data["transit_times"]:
                pmt_all.append((mdom,ichannel))
    pmt_all = list(set(pmt_all))
    if site == None:
        pmt_select = pmt_all
    else:
        if site == "MSU":
            pmt_select = [pmt for pmt in pmt_all if "_M" in pmt[0]]
        elif site == "Desy":
            pmt_select = [pmt for pmt in pmt_all if "_D" in pmt[0]]
        else:
            raise ValueError(f"Invalid site {site}. Valid options are 'MSU', 'Desy', or None.")
    return pmt_select

def plot_transit_time_hist_after_refit(refit_tt,refit_tt2, original_tt, plotFolder) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(original_tt,histtype='step', bins=1000, alpha=0.5, label='Original TT')
    ax.hist(refit_tt, histtype='step', bins=200, alpha=0.5, label='Refit TT')
    ax.hist(refit_tt2, histtype='step', bins=200, alpha=0.5, label='Refit TT2')

    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    # ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel} {get_pmt_uid(mDOM_prod_id, int(channel), mdom_tt_dir)}", transform=ax.transAxes, ha='left', fontsize=10)
    ax.set_xlabel(f"{'tt [ns]'}", fontsize=22)
    ax.set_ylabel(f"{'count'}", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=12,ncols=1)
    ax.set_yscale('log')
    ax.set_xlim(0, 80)
    plot_name = f"{plotFolder}/transit_time_hist_after_refit"
    print(f"plot name: {plot_name}")
    # plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()

def get_transit_time_after_refit(MSU_mDOMs_list,mdom_tt_dir):
    refit_tt = []
    refit_tt2 = []
    original_tt = []
    for device in MSU_mDOMs_list:
        mDOM_prod_id, channel = device
        channel = channel.split("_")[-1]
        channel = int(channel)
        meas_files = glob.glob(f"{mdom_tt_dir}/{mDOM_prod_id}*/*.json")
        available_channels = [extract_channel(ifile) for ifile in meas_files]
        print(f"available channels for {mDOM_prod_id}: {list(set(available_channels))}")
        meas_files = [ifile for ifile in meas_files if extract_channel(ifile) == channel]
        print(f"meas files for {mDOM_prod_id} {meas_files}")
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
                initial_guess = [np.max(y_values), 55, 3, np.max(y_values)/1.5, 55, 3]
                # bounds = ([0, 30, 0.1,0, 30, 0.1], #lower limits
                #            [np.max(y_values)+100, 80, 6,np.max(y_values)+100, 80, 6])      # upper limits
                bounds = ([0, 45, 0.1,0, 45, 0.1], #lower limits
                        [np.max(y_values)+100, 65, 5,np.max(y_values)+100, 65, 5])      # upper limits
                
                poisson_errors = np.sqrt(np.asarray(y_values))
                poisson_errors[poisson_errors == 0] = 1
                try:
                    popt, pcov = curve_fit(double_gaussian,x_values,y_values,p0=initial_guess,bounds=bounds,sigma=poisson_errors,absolute_sigma=True)
                    perr = np.sqrt(np.diag(pcov))
                    perr_dof = perr / np.sqrt(len(x_values) - len(popt))
                    chi2 = np.sum(((y_values - double_gaussian(x_values, *popt))/poisson_errors)**2)
                    ndof = len(x_values) - len(popt)
                    reduced_chi2 = chi2 / ndof
                    A1, mu1, sigma1, A2, mu2, sigma2 = popt
                    refit_tt.append(mu1)
                    refit_tt2.append(mu2)

                    run = data["run_number"]
                    # if run in exclude_runs:
                    #     continue
                    # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f" tt {mu:.1f}\u00B1{sigma:.1f} ns PMT HV {hv:.1f} V {temperature:.1f} \u00b0C Run {run}",alpha=1)
                    tt, tt_spread, a, b, c, chi2, hv = extract_fit_params(ifile)
                    original_tt.append(b)
                except RuntimeError:
                    print(f"Failed to fit double gaussian for {mDOM_prod_id} channel {channel}")
                    continue
    return refit_tt,refit_tt2, original_tt     




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
    MSU_mDOMs_list = extract_json_mdom_pmt_list(transit_time_file,site="MSU")
    #msu mdoms
    print(f"msu mdoms: {len(MSU_mDOMs_list)}")
    # print(f"msu mdoms: {MSU_mDOMs_list}")
    print(f"msu mdoms unique: {len(set([pmt[0] for pmt in MSU_mDOMs_list]))}") #140 MSU mDOMs 


    # plot_single_transit_time_histogram_refit("mDOM_M134", 5, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[427,424])#Run 427 very off
    plot_single_transit_time_histogram_refit_double_gaussian("mDOM_M134", 5, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[427,424])#Run 427 very off
    plot_single_transit_time_histogram_refit_single_gaussian("mDOM_M134", 5, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[427,424])#Run 427 very off

    plot_single_transit_time_histogram_refit_double_gaussian("mDOM_M139", 9, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[])#Run 427 very off
    plot_single_transit_time_histogram_refit_single_gaussian("mDOM_M139", 9, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[])#Run 427 very off
    # for device in MSU_mDOMs_list:
    #     mdom, channel = device
    #     channel = channel.split("_")[-1]
    #     plot_single_transit_time_histogram_refit_double_gaussian(mdom, int(channel), mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[])#Run 427 very off
    refit_tt,refit_tt2,original_tt = get_transit_time_after_refit(MSU_mDOMs_list[:],mdom_tt_dir)
    plot_transit_time_hist_after_refit(refit_tt,refit_tt2, original_tt, plotFolder)

    


if __name__ == "__main__":
    main()
