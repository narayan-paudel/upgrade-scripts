# !/usr/bin/env python
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import json

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})


colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']






def extract_json(transit_time_file,obj_key) -> list:
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
            # print(tt_data.keys())
            for ichannel in tt_data["transit_times"]:
                for itt in tt_data["transit_times"][ichannel]:
                    data_values.append(itt[obj_key])
    return data_values




def plot_transit_time_histogram(transit_times, plot_name) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(transit_times, bins=200,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=r'$tt_{FAT}$', alpha=0.7)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time [ns]", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    ax.set_yscale('log')
    ax.text(0.95, 0.95, fr"mean $tt_{{{'FAT'}}}$: {np.mean(transit_times):.0f}"+r"$\pm$" + fr" {np.std(transit_times):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=12,ncols=1)
    plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()
    
# transit_time_hist(total_mu_list)

def plot_transit_time_chi2(transit_times, chi2_values, plot_name) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.scatter(transit_times, chi2_values, color=colorsCustom[0], alpha=0.7, s=20)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time [ns]", fontsize=22)
    ax.set_ylabel(r"chi2", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    # ax.set_xlim(-10,70)
    plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')


def transit_time_chi2_hist2d(transit_times, chi2_values, plot_name) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    xbins = np.linspace(-400,1100,751)
    ybins = np.linspace(-10,250,66)
    h = ax.hist2d(transit_times, chi2_values, bins=[xbins, ybins], alpha=1,norm=plt.matplotlib.colors.LogNorm(0.3), cmap='viridis')
    # ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time [ns]", fontsize=22)
    ax.set_ylabel(r"chi2", fontsize=22)
    fig.colorbar(h[3], ax=ax)
    # ax.grid(True,alpha=0.6)
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    ax.set_xlim(-10,70)
    # ax.set_ylim(-10,50)
    plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')

# def individual_transit_time_histogram(mdom_prod_id,channel) -> None:

def get_unique_mdoms(transit_time_file) -> list:
    mdom_list = []
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom in data.keys():
            mdom_list.append(mdom)
    return len(list(set(mdom_list)))

def get_unique_mdoms_pmt(transit_time_file) -> list:
    data_values = []
    data_values_diff = []    
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            # print(tt_data.keys())
            for ichannel in tt_data["transit_times"]:
                for itt in tt_data["transit_times"][ichannel]:
                    # if abs(itt["mu"]) > 0:
                    data_values.append((mdom, ichannel))
                    if abs(itt["mu"]) > 0:
                        data_values_diff.append((mdom, ichannel))
    print(f"diff difference {len(list(set(data_values)))-len(list(set(data_values_diff)))}")
    print(f"{[ielt for ielt in list(set(data_values)) if ielt not in list(set(data_values_diff))]}")
    return len(list(set(data_values))), len(list(set(data_values_diff)))






def main() -> None:
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    transit_times = extract_json(transit_time_file,obj_key="mu")
    chi2_values = extract_json(transit_time_file,obj_key="chi2")
    plot_transit_time_histogram(transit_times,plotFolder+"/mDOM_transit_time_all_runs")
    # plot_transit_time_chi2(transit_times, chi2_values, plotFolder+"/mDOM_transit_time_chi2")
    # transit_time_chi2_hist2d(transit_times, chi2_values, plotFolder+"/mDOM_transit_time_chi2_hist2d")
    print(f"Unique mDOMs: {get_unique_mdoms(transit_time_file)}")
    print(f"Unique mDOMs-PMT combinations and non zero tt pmts: {get_unique_mdoms_pmt(transit_time_file)}")

if __name__ == "__main__":
    main()
