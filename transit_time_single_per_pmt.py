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

            popt, pcov = curve_fit(gaussian,x_values,y_values,p0=initial_guess,bounds=bounds,sigma=poisson_errors,absolute_sigma=True)
            A1, mu1, sigma1 = popt

            chi2 = np.sum(((y_values - gaussian(x_values, *popt))/poisson_errors)**2)
            ndof = len(x_values) - len(popt)
            reduced_chi2 = chi2 / ndof

            run = data["run_number"]
            # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f" tt {mu:.1f}\u00B1{sigma:.1f} ns PMT HV {hv:.1f} V {temperature:.1f} \u00b0C Run {run}",alpha=1)
            tt, tt_spread, a, b, c, chi2, hv = extract_fit_params(ifile)
    return A1, mu1, sigma1, reduced_chi2



def extract_json_unique_tt(transit_time_file,mdom_tt_dir ,run_picks_file,need_refits_file,empties_file,site_str,filter_non_zero,check_outliers=None) -> list:
    '''
    Reads the json file with the list of DOMs and transit time measurements and returns a list of transit times
    and other data values
    for chi2, object key : "chi2"
    for transit time, object key:"mu"

    '''
    mu_list = []
    sigma_list = []
    chi2_list = []
    mdom_list_collection_checker = []
    mdom_collection_checker = []
    with open(run_picks_file, 'r') as f:
        run_picks_data = json.load(f)
    with open(need_refits_file, 'r') as f:
        need_refits_data = json.load(f)
    with open(empties_file, 'r') as f:
        empties_data = json.load(f)
    run_picks_data_mdoms = [ielt["mDOM"] for ielt in run_picks_data]
    need_refits_data_mdoms = [ielt["mDOM"] for ielt in need_refits_data]
    need_refits_data_mdoms_pmt = [ielt["mDOM"]+"_"+ielt["channel"] for ielt in need_refits_data]
    empties_data_mdoms = [ielt["mDOM"] for ielt in empties_data]
    empties_data_mdoms_pmt = [ielt["mDOM"]+"_"+ielt["channel"] for ielt in empties_data]
    # print(f"run picks data: {run_picks_data}")
    # print(f"need refits data: {need_refits_data}")
    # print(f"empties data: {empties_data}")
    print(f"run picks data mdoms: {need_refits_data_mdoms_pmt}")
    data_values = []    
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            if site_str in mdom:
                if mdom in run_picks_data_mdoms:
                    select_mdom_run_data = [ielt for ielt in run_picks_data if ielt["mDOM"] == mdom][0]
                    select_runs = select_mdom_run_data["select_run"]
                    for ichannel in tt_data["transit_times"]:
                        if len(tt_data["transit_times"][ichannel]) == 0:
                            print(f"{mdom} channel {ichannel} has no transit time data, skipping")
                        else:
                            # print(f"channel {ichannel.split('_')[-1]} transit time")
                            select_channel_run_data = [ielt for ielt in select_mdom_run_data["select_run"] if int(ielt["channel"]) == int(ichannel.split('_')[-1])][0]
                            select_run = select_channel_run_data["run"]
                            available_runs = [itt["run_number"] for itt in tt_data["transit_times"][ichannel]]

                            # print(f"select channel run data for {mdom} channel {ichannel}: {select_run} {available_runs}")
                            mu = [itt["mu"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                            sigma = [itt["sigma"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                            a = [itt["a"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                            b = [itt["b"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                            c = [itt["c"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                            chi2 = [itt["chi2"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                            applied_hv = [itt["applied HV"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                            temperature = [itt["temperature"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                            tt_info = [itt for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]

                            # print(f"tt data for {mdom} channel {ichannel} select runs: {tt}")
                            if len(b) > 0:
                                if abs(b[0])>150:
                                    print(f"transit time data for {mdom} channel {ichannel} {select_run} select runs: {b[0]}")
                                # print(f"transit time data for {mdom} channel {ichannel} select runs: {tt}")
                                mu_list.append(b[0])
                                sigma_list.append(c[0])
                                chi2_list.append(chi2[0])
                                if check_outliers:
                                    if b[0] < check_outliers[0]:
                                        print(f"run picks data")
                                        print(f"transit time data for {mdom} channel {int(ichannel.split('_')[-1])} {select_run} select runs: {b[0]:.1f} ns")
                                    if b[0]>check_outliers[1]:
                                        print(f"run picks data")
                                        print(f"large transit time data for {mdom} channel {int(ichannel.split('_')[-1])} {select_run} select runs: {b[0]:.1f} ns")
                                if mdom+ f"_{int(ichannel.split('_')[-1])}" not in mdom_list_collection_checker:
                                    mdom_list_collection_checker.append(mdom+ f"_{int(ichannel.split('_')[-1])}")
                                    mdom_collection_checker.append(mdom)
                                else:
                                    print(f"DUPLICATED ENTRY IN COLLECTION CHECKER, CHECK CODE {mdom} channel {ichannel}")
                else:
                    for ichannel in np.linspace(0,23,24):
                        print(f"tt info for {mdom} channel {ichannel} {tt}")
                        if mdom+ f"_{int(ichannel)}" in need_refits_data_mdoms_pmt:
                            A1, mu1, sigma1, reduced_chi2 = get_single_gaussian_fit(mdom, int(ichannel), mdom_tt_dir)
                            mu_list.append(mu1)
                            sigma_list.append(sigma1)
                            chi2_list.append(reduced_chi2)
                            if check_outliers:
                                if mu1 < check_outliers[0]:
                                    print(f"refit data")
                                    print(f"transit time data for {mdom} channel {ichannel} {itt['run_number']} runs: {mu1:.1f} ns")
                                if mu1>check_outliers[1]:
                                    print(f"refit data")
                                    print(f"large transit time data for {mdom} channel {ichannel} {itt['run_number']} runs: {mu1:.1f} ns")
                            if mdom+ f"_{int(ichannel)}" not in mdom_list_collection_checker:
                                mdom_list_collection_checker.append(mdom+ f"_{int(ichannel)}")
                                mdom_collection_checker.append(mdom)
                            else:
                                print(f"DUPLICATED ENTRY IN COLLECTION CHECKER, CHECK CODE {mdom} channel {ichannel}")
                        elif mdom+ f"_{int(ichannel)}" in empties_data_mdoms_pmt:
                            print(f"{mdom} channel {ichannel} has empty measurement, skipping")
                            continue                    
                        else:
                            for itt in tt_data["transit_times"][f"channel_{int(ichannel)}"]:
                                    if filter_non_zero:
                                        if abs(itt["mu"]) > 0:
                                            mu_list.append(itt["b"])
                                            sigma_list.append(itt["c"])
                                            chi2_list.append(itt["chi2"])
                                            if check_outliers:
                                                if itt["b"]< check_outliers[0]:
                                                    print(f"others data")
                                                    print(f"transit time data for {mdom} channel {ichannel} {itt['run_number']} runs: {itt['b']:.1f} ns")
                                                if itt["b"]>check_outliers[1]:
                                                    print(f"others data")
                                                    print(f"large transit time data for {mdom} channel {ichannel} {itt['run_number']} runs: {itt['b']:.1f} ns")
                                            if mdom+ f"_{int(ichannel)}" not in mdom_list_collection_checker:
                                                mdom_list_collection_checker.append(mdom+ f"_{int(ichannel)}")
                                                mdom_collection_checker.append(mdom)
                                            else:
                                                print(f"DUPLICATED ENTRY IN COLLECTION CHECKER, CHECK CODE {mdom} channel {ichannel}")
                                    else:
                                        mu_list.append(itt["b"])
                                        sigma_list.append(itt["c"])
                                        chi2_list.append(itt["chi2"])
                                        if check_outliers:
                                            if itt["b"]< check_outliers[0]:
                                                print(f"others data")
                                                print(f"transit time data for {mdom} channel {ichannel} {itt['run_number']} runs: {itt['b']:.1f} ns")
                                            if itt["b"]>check_outliers[1]:
                                                print(f"others data")
                                                print(f"large transit time data for {mdom} channel {ichannel} {itt['run_number']} runs: {itt['b']:.1f} ns")
                                        if mdom+ f"_{int(ichannel)}" not in mdom_list_collection_checker:
                                            mdom_list_collection_checker.append(mdom+ f"_{int(ichannel)}")
                                            mdom_collection_checker.append(mdom)
                                        else:
                                            print(f"here")
                                            print(f"DUPLICATED ENTRY IN COLLECTION CHECKER, CHECK CODE {mdom} channel {ichannel}")
    print(f"collected {len(mdom_list_collection_checker)} transit time values for {site_str} DOMs")
    print(f"collected {len(list(set(mdom_collection_checker)))} transit time values for {site_str} DOMs")
    for ichannel in range(24):
        for mdom in list(set(mdom_collection_checker)):
            if mdom+ f"_{int(ichannel)}" in mdom_list_collection_checker:
                pass
            else:
                print(f"missing transit time value for {site_str} {mdom} channel {ichannel}")
        

    return mu_list, sigma_list, chi2_list

def measurements_per_channel(transit_time_file) -> dict:
    masurement_number = []
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            for ichannel in tt_data["transit_times"]:
                n_measurements = len(tt_data["transit_times"][ichannel])
                masurement_number.append(n_measurements)
    return masurement_number

def plot_histogram(data, plot_name,bins=200) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(data, bins=bins,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=r'$tt_{FAT}$', alpha=0.7)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"runs per channel", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    ax.set_yscale('log')
    ax.legend(fontsize=8,ncols=1)
    ax.grid(True,alpha=0.6)
    print(f"plot name: {plot_name}")
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_transit_time_proxy_histogram(transit_times_msu, transit_times_desy, plot_name,bins=200) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(transit_times_msu, bins=bins,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=r'MSU $tt_{FAT}$ proxy', alpha=1)
    ax.hist(transit_times_desy, bins=bins,histtype='step',color=colorsCustom[1],linewidth=1.5 ,label=r'DESY $tt_{FAT}$ proxy', alpha=1)
    ax.hist(transit_times_desy+transit_times_msu, bins=bins,histtype='step',color=colorsCustom[2],linewidth=1.5 ,label=r'Both $tt_{FAT}$ proxy', alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"time [ns]", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    ax.set_yscale('log')
    ax.text(0.95, 0.95, fr"MSU $tt_{{{'FAT'}}}$: {np.mean(transit_times_msu):.0f}"+r"$\pm$" + fr" {np.std(transit_times_msu):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    ax.text(0.95, 0.85, fr"DESY $tt_{{{'FAT'}}}$: {np.mean(transit_times_desy):.0f}"+r"$\pm$" + fr" {np.std(transit_times_desy):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    ax.text(0.95, 0.75, fr"Both $tt_{{{'FAT'}}}$: {np.mean(transit_times_msu+transit_times_desy):.0f}"+r"$\pm$" + fr" {np.std(transit_times_msu+transit_times_desy):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    ax.grid(True,alpha=0.6)
    ax.set_xlim(45,65)
    ax.legend(fontsize=12,ncols=1)
    plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()


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
    # transit_times_single_per_pmt = extract_json_unique_tt(transit_time_file,run_picks_json,refit_json,empty_meas_json,obj_key="mu",filter_non_zero=True)
    transit_times_single_per_pmt_msu, sigma_list, chi2_list = extract_json_unique_tt(transit_time_file,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_M",filter_non_zero=False,check_outliers=[46,60])
    transit_times_single_per_pmt_desy, sigma_list, chi2_list = extract_json_unique_tt(transit_time_file,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_D",filter_non_zero=False,check_outliers=[46,60.5])
    # print(transit_times_single_per_pmt)
    print(f"number of transit times: {len(transit_times)}")
    print(f"number of transit times with single entry per pmt:MSU {len(transit_times_single_per_pmt_msu)} DESY {len(transit_times_single_per_pmt_desy)} Both {len(transit_times_single_per_pmt_msu+transit_times_single_per_pmt_desy)}")
    plot_transit_time_proxy_histogram(transit_times_single_per_pmt_msu,transit_times_single_per_pmt_desy,plotFolder+"/mDOM_transit_times_test",bins=np.linspace(45,65,41))
    measuerment_numbers = measurements_per_channel(transit_time_file)
    # print(f"measurement numbers per channel: {measuerment_numbers}")
    # plot_histogram(measuerment_numbers,plotFolder+"/mDOM_transit_time_measurement_numbers",bins=np.linspace(-0.5,20.5,22))
    plot_histogram(measuerment_numbers,plotFolder+"/mDOM_transit_time_measurement_numbers",bins=np.linspace(0,20,21))

if __name__ == "__main__":
    main()