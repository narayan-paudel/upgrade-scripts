# !/usr/bin/env python
import glob
from subprocess import run

import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from transit_time_dependence_plots import extract_json, extract_json_tts_filter, extract_json_min_tts_filter, extract_json_tts_chi2_filter,plot_single_transit_time_histogram

import json

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})


# colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
colorsCustom = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#1f78b4','#33a02c','#e31a1c','#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#ff7f00','#cab2d6','#ffff99','#6a3d9a','#b15928','#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666','#1f78b4']



def main() -> None:
    mdom_tt_dir = home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    ################################################################################
    #main file
    # transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    #selected files
    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time_select.json"
    ################################################################################
    # transit_times = extract_json(transit_time_file,obj_key="mu",filter_non_zero=False)
    # transit_times_nonzero = extract_json(transit_time_file,obj_key="mu",filter_non_zero=True)
    # transit_times_with_tts_less_6 = extract_json_tts_filter(transit_time_file,obj_key="mu",filter_tts=6,filter_non_zero=True)
    # chi2_with_tts_less_6 = extract_json_tts_filter(transit_time_file,obj_key="chi2",filter_tts=6,filter_non_zero=True)
    # transit_times_with_minimum_nonzero_tts = extract_json_min_tts_filter(transit_time_file,obj_key="mu",filter_tts=6,filter_non_zero=True)
    # chi2_values = extract_json(transit_time_file,obj_key="chi2",filter_non_zero=False)
    # transit_time_spread = extract_json(transit_time_file,obj_key="sigma",filter_non_zero=False)
    #chi2 based filtering
    transit_time_chi2_tt_filter,chi2_selected_pmts,chi2_rejected_pmts = extract_json_tts_chi2_filter(transit_time_file,obj_key="mu",filter_tts=6,filter_chi2=2,filter_non_zero=True)
    #transit time outside 28-68 ns range
    # for device in chi2_rejected_pmts:
    # #     print(device.split("_channel_"))
    #     mdom, channel = device
    #     channel = channel.split("_")[-1]
    #     # print(f"Plotting transit time histogram for {mdom} channel {channel} PMT UID {pmt_uid} RUN {run} with negative transit time")
    #     plot_single_transit_time_histogram(mdom, int(channel), mdom_tt_dir, plotFolder,fit_line=True)

    print(f"number of selected pmts with chi2 filter: {len(chi2_selected_pmts)}")
    ###########filtered pmts#########################
    # for device in chi2_selected_pmts:
    # #     print(device.split("_channel_"))
    #     mdom, channel = device
    #     channel = channel.split("_")[-1]
    #     # print(f"Plotting transit time histogram for {mdom} channel {channel} PMT UID {pmt_uid} RUN {run} with negative transit time")
    #     plot_single_transit_time_histogram(mdom, int(channel), mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80])
    ######################remaining pmts###################
    # for device in chi2_rejected_pmts:
    # #     print(device.split("_channel_"))
    #     mdom, channel = device
    #     channel = channel.split("_")[-1]
    #     # print(f"Plotting transit time histogram for {mdom} channel {channel} PMT UID {pmt_uid} RUN {run} with negative transit time")
    #     plot_single_transit_time_histogram(mdom, int(channel), mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80])


    inspected_mdoms = ["mDOM_D001","mDOM_D002","mDOM_D003","mDOM_D004","mDOM_D005","mDOM_D006","mDOM_D007","mDOM_D008",
                    "mDOM_D010","mDOM_D011","mDOM_D012","mDOM_D013","mDOM_D015","mDOM_D016"]
    mdoms_needing_refit = [{"mDOM": "mDOM_D007", "channel": "18", "run": "22"}, {"mDOM": "mDOM_D015", "channel": "18", "run": "2"}, {"mDOM": "mDOM_D015", "channel": "22", "run": "2"}]
    mdoms_with_runs_pick = [{"mDOM": "mDOM_D017", "channel": "0", "run": "152"}]
    empty_measurement = [{"mDOM": "mDOM_D016", "channel": "13", "run": "138"}]
    plot_single_transit_time_histogram("mDOM_D017", 0, mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[40,80],exclude_runs=[151,153,154])#Run 151 very off





if __name__ == "__main__":
    main()
