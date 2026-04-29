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

from transit_time_dependence_plots import (extract_json,extract_json_mdom_pmt_list,get_pmt_uid, extract_json_tts_filter, extract_json_min_tts_filter, extract_json_tts_chi2_filter,plot_single_transit_time_histogram)

def main() -> None:
    upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"
    geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))
    mdom_tt_dir = home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    mdom_pmt_list = extract_json_mdom_pmt_list(transit_time_file)
    which_site = "_M"    
    for device in mdom_pmt_list:
        mdom, channel = device
        if which_site in mdom:
            print(f"mDOM {mdom} channel {channel}")
            channel = channel.split("_")[-1]
            pmt_uid = get_pmt_uid(mdom, int(channel), mdom_tt_dir)
            print(f"mDOM {mdom} channel {channel} PMT UID {pmt_uid}")
            plot_single_transit_time_histogram(mdom, int(channel), mdom_tt_dir, plotFolder,fit_line=True,fit_xlim=[45,75])   

if __name__ == "__main__":
    main()
