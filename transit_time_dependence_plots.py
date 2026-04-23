# !/usr/bin/env python
import glob
from subprocess import run

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


# colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#1f78b4','#33a02c','#e31a1c','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#1f78b4']


def extract_json(transit_time_file,obj_key,filter_non_zero) -> list:
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
                    if filter_non_zero:
                        if abs(itt["mu"]) > 0:
                            data_values.append(itt[obj_key])
                    else:
                        data_values.append(itt[obj_key])
    return data_values


def extract_json_tts_filter(transit_time_file,obj_key,filter_tts,filter_non_zero) -> list:
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
                    if itt["sigma"] < filter_tts:
                        if filter_non_zero:
                            if abs(itt["mu"]) > 0:
                                data_values.append(itt[obj_key])
                        else:
                            data_values.append(itt[obj_key])
    return data_values


def extract_json_min_tts_filter(transit_time_file,obj_key,filter_tts,filter_non_zero) -> list:
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
                this_channel_tt = []
                this_channel_tts = []
                for itt in tt_data["transit_times"][ichannel]:
                    if itt["sigma"] < filter_tts:
                        if filter_non_zero:
                            if abs(itt["mu"]) > 0:
                                this_channel_tt.append(itt[obj_key])
                                this_channel_tts.append(itt["sigma"])
                        else:
                            this_channel_tt.append(itt[obj_key])
                            this_channel_tts.append(itt["sigma"])
                if len(this_channel_tt) > 0:
                    min_tts = min(this_channel_tts)
                    for tt,tts in zip(this_channel_tt,this_channel_tts):
                        print()
                        if abs(tts-min_tts)<1e-3:
                            data_values.append(tt)
    return data_values



def extract_json_tts_chi2_filter(transit_time_file,obj_key,filter_tts,filter_chi2,filter_non_zero) -> list:
    '''
    Reads the json file with the list of DOMs and transit time measurements and returns a list of transit times
    and other data values
    for chi2, object key : "chi2"
    for transit time, object key:"mu"

    '''
    pmt_all = []
    pmt_filter = []
    data_values = []    
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            # print(tt_data.keys())
            for ichannel in tt_data["transit_times"]:
                for itt in tt_data["transit_times"][ichannel]:
                    run_number = itt["run_number"]
                    if itt["sigma"] < filter_tts:
                        if itt["chi2"] < filter_chi2:
                            if filter_non_zero:
                                if abs(itt["mu"]) > 0:
                                    data_values.append(itt[obj_key])
                                    pmt_filter.append((mdom,ichannel))
                                data_values.append(itt[obj_key])
                                pmt_filter.append((mdom,ichannel))
                    pmt_all.append((mdom,ichannel))
    print(f"diff difference chi2 tt filter {len(list(set(pmt_all)))-len(list(set(pmt_filter)))}")
    selected_pmts = [ielt for ielt in list(set(pmt_filter))]
    rejected_pmts = [ielt for ielt in list(set(pmt_all)) if ielt not in list(set(pmt_filter))]
    # print(f"rejected pmts: {rejected_pmts}")
    return data_values, selected_pmts,rejected_pmts

def extract_json_mdom_pmt_list(transit_time_file) -> list:
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
    return list(set(pmt_all))
    


def extract_fit_params(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    y_values = data["meas_data"][0]['y_values']
    # print(data["meas_data"])
    # print(len(data["meas_data"]))
    param_dict = {"Transit time":np.nan,"Transit time spread":np.nan,"a":np.nan,"b":np.nan,"c":np.nan,"chi2":np.nan,"applied HV":np.nan}
    for j,ielt in enumerate(data["meas_data"][1:]):
        param_dict[ielt["label"]] = ielt["value"]

    # for j,ielt in enumerate(data["meas_data"][1:]):
    #     if ielt["label"] == "Transit time":
    #         mu = ielt["value"]
    #     if ielt["label"] == "Transit time spread":
    #         sigma = ielt["value"]
    #     if ielt["label"] == "a":
    #         a = ielt["value"]
    #     if ielt["label"] == "b":
    #         mu = ielt["value"]


    return param_dict["Transit time"], param_dict["Transit time spread"],param_dict["a"]\
        ,param_dict["b"],param_dict["c"],param_dict["chi2"],param_dict["applied HV"]



def plot_transit_time_histogram(transit_times, plot_name,bins=200) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(transit_times, bins=bins,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=r'$tt_{FAT}$', alpha=0.7)
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
    ax.set_xlim(-10,70)
    ax.set_ylim(-10,50)
    plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')

def plot_transit_time_spread(transit_times, tts_values, plot_name) -> None:
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.scatter(transit_times, tts_values, color=colorsCustom[0], alpha=0.7, s=20)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time [ns]", fontsize=22)
    ax.set_ylabel(r"transit time spread [ns]", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.set_yscale('log')
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
    ax.set_ylim(-10,50)
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

def extract_channel(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)

    channel_dict = {"mb channel":np.nan}
    device_uid = data["device_uid"]
    pmt = data['subdevice_uid'].split('_')[-1]
    channel = data["meas_data"][2]["value"]
    channel_dict["mb channel"] = channel
    for j,ielt in enumerate(data["meas_data"][1:]):
        channel_dict[ielt["label"]] = ielt["value"]
    return channel_dict["mb channel"]

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

def get_unique_mdoms_pmt_tts_zero_filter(transit_time_file,filter_tts) -> list:
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
                    if itt["sigma"] < filter_tts and abs(itt["mu"]) > 0:
                        data_values_diff.append((mdom, ichannel))
    print(f"diff difference {len(list(set(data_values)))-len(list(set(data_values_diff)))}")
    print(f"{[ielt for ielt in list(set(data_values)) if ielt not in list(set(data_values_diff))]}")
    return len(list(set(data_values))), len(list(set(data_values_diff)))


def get_unique_mdoms_pmt_tts_zero_filter_28_68(transit_time_file,filter_tts) -> list:
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
                    if itt["sigma"] < filter_tts and abs(itt["mu"]) > 0 and 28 < itt["mu"] < 68:
                        data_values_diff.append((mdom, ichannel))
    print(f"diff difference {len(list(set(data_values)))-len(list(set(data_values_diff)))}")
    print(f"{[ielt for ielt in list(set(data_values)) if ielt not in list(set(data_values_diff))]}")
    return len(list(set(data_values))), len(list(set(data_values_diff)))

def extract_pmt_uid(json_file) -> str:
    with open(json_file, 'r') as f:
        data = json.load(f)
    device_uid = data["device_uid"]
    # print(f"pmt uid {data['pmt_uid']}")
    pmt = data['subdevice_uid'].split('_')[-1]
    return pmt


def get_pmt_uid(mdom_prod_id, channel,mdom_tt_dir) -> str:
    '''Returns the PMT UID for a given mDOM production ID and channel'''
    meas_files = glob.glob(f"{mdom_tt_dir}/{mdom_prod_id}*/*.json")
    # print([extract_channel(ifile) for ifile in meas_files])
    pmt_uid = list(set([extract_pmt_uid(ifile) for ifile in meas_files if extract_channel(ifile) == channel]))
    # print(f"PMT UID for {mdom_prod_id} channel {channel}: {pmt_uid}")
    if len(pmt_uid) > 1:
        print(f"multiple PMT UIDs found for {mdom_prod_id} channel {channel}: {pmt_uid}")
    elif len(pmt_uid) == 0:
        print(f"no PMT UID found for {mdom_prod_id} channel {channel}")
        return None
    return pmt_uid[-1]

def extract_run_number(json_file) -> str:
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["run_number"]

def get_run_number(mdom_prod_id, channel,mdom_tt_dir) -> str:
    '''Returns the PMT UID for a given mDOM production ID and channel'''
    meas_files = glob.glob(f"{mdom_tt_dir}/{mdom_prod_id}*/*.json")
    # print([extract_channel(ifile) for ifile in meas_files])
    run_number = list(set([extract_run_number(ifile) for ifile in meas_files if extract_channel(ifile) == channel]))
    print(run_number)
    if len(run_number) > 1:
        print(f"multiple run numbers found for {mdom_prod_id} channel {channel}: {run_number}")
    return run_number[-1]

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def plot_single_transit_time_histogram(mDOM_prod_id, channel, mdom_tt_dir, plotFolder,fit_line=False,fit_xlim=[-30,100],exclude_runs=[]) -> None:
    '''Plots the transit time histogram for a single mDOM and all its channels'''
    meas_files = glob.glob(f"{mdom_tt_dir}/{mDOM_prod_id}*/*.json")

    available_channels = [extract_channel(ifile) for ifile in meas_files]
    print(f"available channels for {mDOM_prod_id}: {list(set(available_channels))}")
    meas_files = [ifile for ifile in meas_files if extract_channel(ifile) == channel]
    if len(meas_files) == 0:
        print(f"meas files for {mDOM_prod_id} channel {channel}: {meas_files}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    x_label = ""
    y_label = ""
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
            run = data["run_number"]
            if run in exclude_runs:
                continue
            # ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f" tt {mu:.1f}\u00B1{sigma:.1f} ns PMT HV {hv:.1f} V {temperature:.1f} \u00b0C Run {run}",alpha=1)
            tt, tt_spread, a, b, c, chi2, hv = extract_fit_params(ifile)
            print(f"mdom {mDOM_prod_id} channel {channel} Run {run}: tt {tt:.1f} ns, tts {tt_spread:.1f} ns, a {a:.1f}, b {b:.1f}, c {c:.1f}, chi2 {chi2:.1f}, hv {hv:.1f} V")
            # if c < 20:
            ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom[i],label=f"Run {run}: tt {b:.0f} ns ({tt:.0f} ns) tts {tt_spread:.0f} ns chi2 {chi2:.1f}",alpha=1)
            if fit_line:
                x_values = np.linspace(fit_xlim[0],fit_xlim[1],1000)                
                ax.plot(x_values,gauss(x_values,a,b,c),ls='--',lw = 2.5,c=colorsCustom[i],alpha=1)

    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel}", transform=ax.transAxes, ha='left', fontsize=10)
    # ax.text(0.55, 0.95, f"{mDOM_prod_id} channel {channel} {get_pmt_uid(mDOM_prod_id, int(channel), mdom_tt_dir)}", transform=ax.transAxes, ha='left', fontsize=10)
    ax.set_xlabel(f"{x_label}", fontsize=22)
    ax.set_ylabel(f"{y_label}", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=8,ncols=1)
    plot_name = f"{plotFolder}/{mDOM_prod_id}_channel_{channel}"
    # plt.savefig(plot_name+".png",transparent=False,bbox_inches='tight')
    plt.savefig(plot_name+".pdf",transparent=False,bbox_inches='tight')
    plt.close()


def device_with_negative_tt(transit_time_file, mdom_tt_dir) -> list:
    '''
    Reads the json file with the list of DOMs and transit time measurements and returns a list of devices with negative transit time'''
    negative_tt_devices = []
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            for ichannel in tt_data["transit_times"]:
                for itt in tt_data["transit_times"][ichannel]:
                    if itt["mu"] < 0:
                        run_number = itt["run_number"]
                        print(f"{mdom} {ichannel} PMT {get_pmt_uid(mdom, int(ichannel.split('_')[-1]), mdom_tt_dir)} RUN {run_number}: {itt['mu']:.0f} ns")
                        negative_tt_devices.append(mdom+"_"+str(ichannel)+"_"+str(run_number))
                        break
    return list(set(negative_tt_devices))

def prod_id_to_icm_id(prod_id,geometry_files) -> str:
    '''
    production id
    '''
    icm_id_list = []
    for ifile in geometry_files:
        with open(ifile, 'r') as f:
            data = json.load(f)
        for idev in data[0]["devices"]:
            if idev["production_id"] == prod_id:
                icm_id_list.append(idev["icm_id"])
    if len(icm_id_list)>1:
        print(f"more than one icm id found for prod id {prod_id} {icm_id_list}")
    return icm_id_list[0]


def main() -> None:
    upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"
    geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))
    mdom_tt_dir = home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    transit_times = extract_json(transit_time_file,obj_key="mu",filter_non_zero=False)
    transit_times_nonzero = extract_json(transit_time_file,obj_key="mu",filter_non_zero=True)
    transit_times_with_tts_less_6 = extract_json_tts_filter(transit_time_file,obj_key="mu",filter_tts=6,filter_non_zero=True)
    chi2_with_tts_less_6 = extract_json_tts_filter(transit_time_file,obj_key="chi2",filter_tts=6,filter_non_zero=True)
    transit_times_with_minimum_nonzero_tts = extract_json_min_tts_filter(transit_time_file,obj_key="mu",filter_tts=6,filter_non_zero=True)
    chi2_values = extract_json(transit_time_file,obj_key="chi2",filter_non_zero=False)
    transit_time_spread = extract_json(transit_time_file,obj_key="sigma",filter_non_zero=False)
    #chi2 based filtering
    transit_time_chi2_tt_filter,chi2_selected_pmts,chi2_rejected_pmts = extract_json_tts_chi2_filter(transit_time_file,obj_key="mu",filter_tts=6,filter_chi2=20,filter_non_zero=True)
    plot_transit_time_histogram(transit_time_chi2_tt_filter,plotFolder+"/mDOM_transit_time_chi2_tt_filter_runs")
    plot_transit_time_histogram(transit_times,plotFolder+"/mDOM_transit_time_all_runs")
    plot_transit_time_histogram(transit_times_nonzero,plotFolder+"/mDOM_transit_time_nonzero_runs")
    plot_transit_time_histogram(transit_times_with_tts_less_6,plotFolder+"/mDOM_transit_time_tts_less_6_runs")
    plot_transit_time_histogram(transit_times_with_minimum_nonzero_tts,plotFolder+"/mDOM_transit_time_minimum_nonzero_tts_runs")
    plot_transit_time_spread(transit_times, transit_time_spread, plotFolder+"/mDOM_transit_time_tts")
    plot_transit_time_chi2(transit_times, chi2_values, plotFolder+"/mDOM_transit_time_chi2")
    transit_time_chi2_hist2d(transit_times, chi2_values, plotFolder+"/mDOM_transit_time_chi2_hist2d")
    transit_time_chi2_hist2d(transit_times_with_tts_less_6, chi2_with_tts_less_6, plotFolder+"/mDOM_transit_time_chi2_hist2d_tts_less_6")
    print(f"Unique mDOMs: {get_unique_mdoms(transit_time_file)}")
    print(f"Unique mDOMs-PMT combinations and non zero tt pmts: {get_unique_mdoms_pmt(transit_time_file)}")
    print(f"Unique mDOMs-PMT combinations and non zero tt pmts and tts < 6: {get_unique_mdoms_pmt_tts_zero_filter(transit_time_file,filter_tts=6)}")
    print(f"Unique mDOMs-PMT combinations and non zero tt pmts and tts < 6, 28 < mu < 68: {get_unique_mdoms_pmt_tts_zero_filter_28_68(transit_time_file,filter_tts=6)}")
    plot_single_transit_time_histogram("mDOM_D114",18,mdom_tt_dir, plotFolder)
    plot_single_transit_time_histogram("mDOM_D016",13,mdom_tt_dir, plotFolder)
    plot_single_transit_time_histogram("mDOM_D219",13,mdom_tt_dir, plotFolder)
    #transit time outside 28-68 ns range
    # plot_single_transit_time_histogram("mDOM_D054",17,mdom_tt_dir, plotFolder,fit_line=True)
    # plot_single_transit_time_histogram("mDOM_D090",14,mdom_tt_dir, plotFolder,fit_line=True)
    # plot_single_transit_time_histogram("mDOM_D092",1,mdom_tt_dir, plotFolder,fit_line=True)
    # plot_single_transit_time_histogram("mDOM_D090",6,mdom_tt_dir, plotFolder,fit_line=True )
    #

    # negative_tt_devices = device_with_negative_tt(transit_time_file, mdom_tt_dir)
    # print(f"Devices with negative transit time: {negative_tt_devices}")
    # #printing negative tt individual histograms
    # for device in negative_tt_devices:
    #     print(device.split("_channel_"))
    #     mdom, channel_run = device.split("_channel_")
    #     channel, run = channel_run.split("_")
    #     pmt_uid = get_pmt_uid(mdom, int(channel), mdom_tt_dir)
    #     print(f"Plotting transit time histogram for {mdom} channel {channel} PMT UID {pmt_uid} RUN {run} with negative transit time")
    #     plot_single_transit_time_histogram(mdom, int(channel), mdom_tt_dir, plotFolder)
    # plot_single_transit_time_histogram("mDOM_M162",13,mdom_tt_dir, plotFolder)
    #chi2 based filtering
    # for device in chi2_rejected_pmts:
    # #     print(device.split("_channel_"))
    #     mdom, channel = device
    #     channel = channel.split("_")[-1]
    #     # print(f"Plotting transit time histogram for {mdom} channel {channel} PMT UID {pmt_uid} RUN {run} with negative transit time")
    #     plot_single_transit_time_histogram(mdom, int(channel), mdom_tt_dir, plotFolder,fit_line=True)
    #getting ICM ID
    print(prod_id_to_icm_id("mDOM_D114",geometry_files))
    print(prod_id_to_icm_id("mDOM_D016",geometry_files))
    print(prod_id_to_icm_id("mDOM_D219",geometry_files))
    mdom_pmt_list = extract_json_mdom_pmt_list(transit_time_file)
    # recovered_mdom_list = ["mDOM_D074", "mDOM_D075", "mDOM_D032",
    #                         "mDOM_D076", "mDOM_D079", "mDOM_D041", "mDOM_D071",
    #                           "mDOM_D036", "mDOM_D047", "mDOM_D070", "mDOM_D035",
    #                             "mDOM_M168" ]
    recovered_mdom_list = ["mDOM_D074"]
    #recovered_pmt_mdom_list
    # plot_single_transit_time_histogram("mDOM_D074", 16, mdom_tt_dir, plotFolder,fit_line=True)
    # plot_single_transit_time_histogram("mDOM_D074", 18, mdom_tt_dir, plotFolder,fit_line=True)
    # plot_single_transit_time_histogram("mDOM_M030", 10, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    # plot_single_transit_time_histogram("mDOM_M030", 12, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    # plot_single_transit_time_histogram("mDOM_M049", 0, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[427])
    # plot_single_transit_time_histogram("mDOM_M049", 22, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[427])
    plot_single_transit_time_histogram("mDOM_M063", 7, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_M063", 5, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    #################################################################
    #fringe fit checking
    plot_single_transit_time_histogram("mDOM_D090",5,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D090",12,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D090",16,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D090",19,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D090",20,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D090",23,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D084",1,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D176",15,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D089",18,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D070",19,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_D211",10,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,70],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_M173",20,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,75],exclude_runs=[])
    plot_single_transit_time_histogram("mDOM_M093",5,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,75],exclude_runs=[658,629,649,562,128,647,135,657])
    plot_single_transit_time_histogram("mDOM_M033",17,mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,75],exclude_runs=[])



    
    # for device in mdom_pmt_list:
    #     mdom, channel = device
    #     if mdom in recovered_mdom_list:
    #         channel = channel.split("_")[-1]
    #         pmt_uid = get_pmt_uid(mdom, int(channel), mdom_tt_dir)
    #         print(f"mDOM {mdom} channel {channel} PMT UID {pmt_uid}")
    #         plot_single_transit_time_histogram(mdom, int(channel), mdom_tt_dir, plotFolder,fit_line=True)



if __name__ == "__main__":
    main()
