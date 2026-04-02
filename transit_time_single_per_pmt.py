# !/usr/bin/env python
import glob
from subprocess import run

import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from transit_time_dependence_plots import extract_json, extract_json_tts_filter, extract_json_min_tts_filter, extract_json_tts_chi2_filter,plot_single_transit_time_histogram
from transit_time_dependence_plots import prod_id_to_icm_id,plot_transit_time_histogram,extract_channel

import json

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})


# colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
colorsCustom = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#1f78b4','#33a02c','#e31a1c','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#cab2d6','#ffff99','#6a3d9a','#b15928','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#1f78b4']

def extract_json_unique_tt(transit_time_file,mdom_tt_dir ,run_picks_file,need_refits_file,empties_file,obj_key,filter_non_zero) -> list:
    '''
    Reads the json file with the list of DOMs and transit time measurements and returns a list of transit times
    and other data values
    for chi2, object key : "chi2"
    for transit time, object key:"mu"

    '''
    with open(run_picks_file, 'r') as f:
        run_picks_data = json.load(f)
    with open(need_refits_file, 'r') as f:
        need_refits_data = json.load(f)
    with open(empties_file, 'r') as f:
        empties_data = json.load(f)
    run_picks_data_mdoms = [ielt["mDOM"] for ielt in run_picks_data]
    need_refits_data_mdoms = [ielt["mDOM"] for ielt in need_refits_data]
    empties_data_mdoms = [ielt["mDOM"] for ielt in empties_data]
    # print(f"run picks data: {run_picks_data}")
    # print(f"need refits data: {need_refits_data}")
    # print(f"empties data: {empties_data}")
    print(f"run picks data mdoms: {run_picks_data_mdoms}")
    data_values = []    
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, tt_data in data.items():
            if mdom in need_refits_data_mdoms:
                channel_in_refit_data = [ielt for ielt in need_refits_data if ielt["mDOM"] == mdom][0]["channel"]
                run_in_refit_data = [ielt for ielt in need_refits_data if ielt["mDOM"] == mdom][0]["run"]
                print(f"{mdom} is in need refits data, channel in refit data: {channel_in_refit_data}, run in refit data: {run_in_refit_data}")
                for ichannel in tt_data["transit_times"]:
                    if int(channel_in_refit_data) == int(ichannel.split('_')[-1]):
                        hist_data_files = [ifile for ifile in glob.glob(f"{mdom_tt_dir}/{mdom}*/*{run_in_refit_data}.json") if int(extract_channel(ifile)) == int(channel_in_refit_data)][0]
                        print(f"hist data file for {mdom} channel {ichannel} run {run_in_refit_data}: {hist_data_files}")
                        with open(hist_data_files, 'r') as f:
                            data = json.load(f)
                        y_values = data["meas_data"][0]['y_values']
                        x_min = data["meas_data"][0]["x_min"]
                        x_max = data["meas_data"][0]["x_max"]
                        n_bins = data["meas_data"][0]["n_bins"]
                        x_values = np.linspace(x_min,x_max,n_bins)
                        y_sum = np.sum(y_values)
                        xy_sum = np.sum(ix_values*iy_values for ix_values,iy_values in zip(x_values,y_values))
                        mean = xy_sum/y_sum
                        print(f"mean transit time for {mdom} channel {ichannel} run {run_in_refit_data}: {mean}")
                        # data_values.append(mean)

                    # else:
                    #     for itt in tt_data["transit_times"][ichannel]:
                    #         if filter_non_zero:
                    #             if abs(itt["mu"]) > 0:
                    #                 data_values.append(itt[obj_key])
                    #         else:
                    #             data_values.append(itt[obj_key])
                # print(f"{mdom} is in need refits data, skipping")
                # continue
            elif mdom in empties_data_mdoms:
                print(f"{mdom} is in empties data, skipping")
                # continue
            elif mdom in run_picks_data_mdoms:
                select_mdom_run_data = [ielt for ielt in run_picks_data if ielt["mDOM"] == mdom][0]
                select_runs = select_mdom_run_data["select_run"]
                # print(f"{mdom} is in run picks data, select runs: {select_runs}")
                for ichannel in tt_data["transit_times"]:
                    if len(tt_data["transit_times"][ichannel]) == 0:
                        print(f"{mdom} channel {ichannel} has no transit time data, skipping")
                    else:
                        # print(f"channel {ichannel.split('_')[-1]} transit time")
                        select_channel_run_data = [ielt for ielt in select_mdom_run_data["select_run"] if int(ielt["channel"]) == int(ichannel.split('_')[-1])][0]
                        select_run = select_channel_run_data["run"]
                        available_runs = [itt["run_number"] for itt in tt_data["transit_times"][ichannel]]

                        # print(f"select channel run data for {mdom} channel {ichannel}: {select_run} {available_runs}")
                        tt = [itt["mu"] for itt in tt_data["transit_times"][ichannel] if int(itt["run_number"]) == int(select_run)]
                        # print(f"tt data for {mdom} channel {ichannel} select runs: {tt}")
                        if len(tt) > 0:
                            if abs(tt[0])>150:
                                print(f"transit time data for {mdom} channel {ichannel} {select_run} select runs: {tt}")
                            # print(f"transit time data for {mdom} channel {ichannel} select runs: {tt}")
                            data_values.append(tt[0])
            else:
                # print(tt_data.keys())
                for ichannel in tt_data["transit_times"]:
                    for itt in tt_data["transit_times"][ichannel]:
                        if filter_non_zero:
                            if abs(itt["mu"]) > 0:
                                data_values.append(itt[obj_key])
                        else:
                            data_values.append(itt[obj_key])



    return data_values

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
    transit_times_single_per_pmt = extract_json_unique_tt(transit_time_file,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,obj_key="mu",filter_non_zero=False)
    print(f"number of transit times: {len(transit_times)}")
    print(f"number of transit times with single entry per pmt: {len(transit_times_single_per_pmt)}")
    plot_transit_time_histogram(transit_times_single_per_pmt,plotFolder+"/mDOM_transit_times_test",bins=np.linspace(0,60,121))
    measuerment_numbers = measurements_per_channel(transit_time_file)
    print(f"measurement numbers per channel: {measuerment_numbers}")
    # plot_histogram(measuerment_numbers,plotFolder+"/mDOM_transit_time_measurement_numbers",bins=np.linspace(-0.5,20.5,22))
    plot_histogram(measuerment_numbers,plotFolder+"/mDOM_transit_time_measurement_numbers",bins=np.linspace(0,20,21))

if __name__ == "__main__":
    main()