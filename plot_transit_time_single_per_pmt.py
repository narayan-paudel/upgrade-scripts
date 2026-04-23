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
from transit_time_fit_functions import gaussian

import json

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})


# colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#1f78b4','#33a02c','#e31a1c','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#cab2d6','#ffff99','#6a3d9a','#b15928','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#1f78b4']


def extract_json_site(transit_time_file,site_str,obj_key,filter_non_zero) -> list:
    '''
    Reads the json file with the list of DOMs and transit time measurements and returns a list of transit times
    and other data values
    for chi2, object key : "chi2"
    for transit time, object key:"mu"

    '''
    data_values = []    
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            if site_str in mdom:
                for ichannel in tt_data["transit_times"]:
                    for itt in tt_data["transit_times"][ichannel]:
                        if len(itt) > 0:
                            if filter_non_zero:
                                if abs(itt["mu"]) > 0:
                                    data_values.append(itt[obj_key])
                            else:
                                data_values.append(itt[obj_key])
    return data_values

def extract_rse_site(transit_time_file, mdom_tt_dir, site_str, obj_key, filter_non_zero) -> list:
    '''
    Reads the json file with the list of DOMs and transit time measurements and returns a list of transit times
    and other data values
    for chi2, object key : "chi2"
    for transit time, object key:"mu"

    '''
    data_values = []    
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():            
            if site_str in mdom:
                for ichannel in tt_data["transit_times"]:                    
                    for itt in tt_data["transit_times"][ichannel]:
                        if len(itt) > 0:
                            run_number = itt["run_number"]
                            mu = itt["mu"]
                            sigma = itt["sigma"]
                            a = itt["a"]
                            b = itt["b"]
                            c = itt["c"]
                            chi2 = itt["chi2"]

                            meas_files = glob.glob(f"{mdom_tt_dir}/{mdom}*/*_{run_number}.json")
                            meas_files_channel = [ifile for ifile in meas_files if int(extract_channel(ifile)) == int(ichannel.split("_")[-1])]
                            # print(f"mdom {mdom} {ichannel} run {run_number} meas files: {meas_files_channel}")
                            if len(meas_files_channel) > 1:
                                print(f"more than one measurement file found for {mdom} channel {ichannel} run {run_number} {meas_files_channel}")
                            meas_files_channel = meas_files_channel[0]
                            with open(meas_files_channel, 'r') as f:
                                data = json.load(f)
                                transit_times = []
                                y_values = data["meas_data"][0]['y_values']
                                x_min = data["meas_data"][0]["x_min"]
                                x_max = data["meas_data"][0]["x_max"]
                                n_bins = data["meas_data"][0]["n_bins"]
                                x_values = np.linspace(x_min,x_max,n_bins)
                                x_label = data["meas_data"][0]["x_label"]
                                y_label = data["meas_data"][0]["y_label"]
                                y_fit = gaussian(x_values, a,b,c)
                                residuals = y_values - y_fit
                                RSS = np.sum(residuals**2)
                                if RSS > 10**6:
                                    print(f"large residual sum of squares for {mdom} channel {ichannel} run {run_number}: {RSS}")
                                else:
                                    data_values.append(RSS)
    return data_values


def plot_transit_time_proxy_histogram(transit_times_msu, transit_times_desy, plot_name,bins=200,xlim=[45,65],no_component=None) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    if no_component is None:
        ax.hist(transit_times_msu, bins=bins,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=r'MSU tt', alpha=1)
        ax.hist(transit_times_desy, bins=bins,histtype='step',color=colorsCustom[1],linewidth=1.5 ,label=r'DESY tt', alpha=1)
    ax.hist(transit_times_desy+transit_times_msu, bins=bins,histtype='step',color=colorsCustom[2],linewidth=1.5 ,label=r'Combined tt', alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"time [ns]", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    ax.set_yscale('log')
    if no_component is None:
        ax.text(0.95, 0.95, fr"MSU tt: {np.mean(transit_times_msu):.0f}"+r"$\pm$" + fr" {np.std(transit_times_msu):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
        ax.text(0.95, 0.85, fr"DESY tt: {np.mean(transit_times_desy):.0f}"+r"$\pm$" + fr" {np.std(transit_times_desy):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    ax.text(0.95, 0.75, fr"Combined tt: {np.mean(transit_times_msu+transit_times_desy):.0f}"+r"$\pm$" + fr" {np.std(transit_times_msu+transit_times_desy):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    ax.grid(True,alpha=0.6)
    ax.set_xlim(xlim[0],xlim[1])
    ax.legend(fontsize=12,ncols=1)
    plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_rse_histogram(rse_msu, rse_desy, plot_name,bins=200,xlim=[45,65],no_component=None) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    if no_component is None:
        ax.hist(rse_msu, bins=bins,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=r'MSU RSS', alpha=1)
        ax.hist(rse_desy, bins=bins,histtype='step',color=colorsCustom[1],linewidth=1.5 ,label=r'DESY RSS', alpha=1)
    ax.hist(rse_desy+rse_msu, bins=bins,histtype='step',color=colorsCustom[2],linewidth=1.5 ,label=r'Combined RSS', alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"Residual sum of squares", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    ax.set_yscale('log')
    if no_component is None:
        ax.text(0.95, 0.95, fr"MSU RSS: {np.mean(rse_msu):.0f}"+r"$\pm$" + fr" {np.std(rse_msu):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
        ax.text(0.95, 0.85, fr"DESY RSS: {np.mean(rse_desy):.0f}"+r"$\pm$" + fr" {np.std(rse_desy):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    # ax.text(0.95, 0.75, fr"Combined RSS: {np.mean(rse_msu+rse_desy):.0f}"+r"$\pm$" + fr" {np.std(rse_msu+rse_desy):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    ax.grid(True,alpha=0.6)
    # ax.set_xlim(xlim[0],xlim[1])
    ax.legend(fontsize=12,ncols=1)
    plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()



def main() -> None:
    C_msu = 23.9 #ns
    C_desy = 25.5 #ns
    upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"
    geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))
    mdom_tt_dir = home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    run_picks_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_run_picks.json"
    refit_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_tt_needing_refit.json"
    empty_meas_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_empty_meas.json"

    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    transit_times_single_per_pmt_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time_single_pick.json"
    # transit_times = extract_json(transit_time_file,obj_key="mu",filter_non_zero=False)
    # transit_times_single_per_pmt = extract_json_unique_tt(transit_time_file,run_picks_json,refit_json,empty_meas_json,obj_key="mu",filter_non_zero=True)
    transit_times_single_per_pmt_msu = extract_json_site(transit_times_single_per_pmt_file,site_str="_M",obj_key="b",filter_non_zero=False)
    transit_times_single_per_pmt_desy = extract_json_site(transit_times_single_per_pmt_file,site_str="_D",obj_key="b",filter_non_zero=False)

    transit_times_single_per_pmt_msu_after_correction = extract_json_site(transit_times_single_per_pmt_file,site_str="_M",obj_key="mu",filter_non_zero=False)
    transit_times_single_per_pmt_desy_after_correction = extract_json_site(transit_times_single_per_pmt_file,site_str="_D",obj_key="mu",filter_non_zero=False)

    # print(transit_times_single_per_pmt)
    #####################################
    #############plot histogram##########
    plot_transit_time_proxy_histogram(transit_times_single_per_pmt_msu, transit_times_single_per_pmt_desy, plotFolder+"/transit_time_proxy_histogram_final",bins=np.linspace(45,65,41),xlim=[45,65],no_component=True)
    plot_transit_time_proxy_histogram(transit_times_single_per_pmt_msu, transit_times_single_per_pmt_desy, plotFolder+"/transit_time_proxy_histogram_final_components",bins=np.linspace(45,65,41),xlim=[45,65],no_component=None)
    transit_times_single_per_pmt_desy_corrected = [tt - C_desy for tt in transit_times_single_per_pmt_desy]
    transit_times_single_per_pmt_msu_corrected = [tt - C_msu for tt in transit_times_single_per_pmt_msu]
    plot_transit_time_proxy_histogram(transit_times_single_per_pmt_msu_corrected, transit_times_single_per_pmt_desy_corrected,plotFolder+"/transit_time_proxy_histogram_final_corrected",bins=np.linspace(20,40,41),xlim=[20,40],no_component=True)
    plot_transit_time_proxy_histogram(transit_times_single_per_pmt_msu_corrected, transit_times_single_per_pmt_desy_corrected,plotFolder+"/transit_time_proxy_histogram_final_corrected_components",bins=np.linspace(20,40,41),xlim=[20,40],no_component=None)
    plot_transit_time_proxy_histogram(transit_times_single_per_pmt_msu_after_correction, transit_times_single_per_pmt_desy_after_correction,plotFolder+"/transit_time_proxy_histogram_final_after_corrected",bins=np.linspace(20,40,41),xlim=[20,40],no_component=True)
    plot_transit_time_proxy_histogram(transit_times_single_per_pmt_msu_after_correction, transit_times_single_per_pmt_desy_after_correction,plotFolder+"/transit_time_proxy_histogram_final_after_corrected_components",bins=np.linspace(20,40,41),xlim=[20,40],no_component=None)

    #######################################
    ###########plot residual###############
    rse_msu = extract_rse_site(transit_times_single_per_pmt_file, mdom_tt_dir, site_str="_M", obj_key="b", filter_non_zero=False)
    rse_desy = extract_rse_site(transit_times_single_per_pmt_file, mdom_tt_dir, site_str="_D", obj_key="b", filter_non_zero=False)
    # plot_rse_histogram(rse_msu, rse_desy, plotFolder+"/rse_histogram_final",bins=np.linspace(0,20,41),xlim=[0,20])
    plot_rse_histogram(rse_msu, rse_desy, plotFolder+"/rse_histogram_final",bins=np.linspace(-5,0.4*10**6,1000),xlim=[-5,0.4*10**6],no_component=True)
    plot_rse_histogram(rse_msu, rse_desy, plotFolder+"/rse_histogram_final_components",bins=np.linspace(-5,0.4*10**6,1000),xlim=[-5,0.4*10**6],no_component=None)
if __name__ == "__main__":
    main()