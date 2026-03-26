#!/usr/bin/env python

import os
import glob
import re

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

plotFolder = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"


def inspect_transit_time(transit_time_file):
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        # print("data", data)
        flat_mDOM = []
        for mdom, channels in data.items():
            for channel, value in channels.items():
                if type(value) == list and len(value) == 1:
                    print(f"{mdom} {channel} {value}")

        flat_mDOM = list(set(flat_mDOM))
        print(f"mDOMs with zero transit time {len(flat_mDOM)} {flat_mDOM}")

inspect_transit_time('/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json')

def plot_transit_time_single_measurement(transit_time_file):
    single_measurements = []    
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, channels in data.items():
            for channel, value in channels.items():
                if type(value) == list and len(value) == 1:
                    if 100 > value[0] > 25:
                        single_measurements.append(value[0])
                    if value[0]> 100:
                        print(f"too large tt {mdom} {channel} {value}")
                    if value[0] < 25:
                        print(f"too small tt {mdom} {channel} {value}")
    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(single_measurements, bins=200,histtype='step', alpha=0.7, color='blue')
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time [ns]", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"mDOM_transit_time_single_measuremnet.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"mDOM_transit_time_single_measuremnet.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_transit_time_single_measurement('/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json')


def plot_transit_time_multiple_measurement(transit_time_file):
    single_measurements = []    
    with open(transit_time_file, 'r') as f:
        data = json.load(f)
        for mdom, channels in data.items():
            for channel, value in channels.items():
                if type(value) != list and np.isnan(value):
                    print(f"mDOM {mdom} channel  {channel} value {value}")
                    # if 100 > value[0] > 25:
                    #     single_measurements.append(value[0])
                    # if value[0]> 100:
                    #     print(f"too large tt {mdom} {channel} {value}")
                    # if value[0] < 25:
                    #     print(f"too small tt {mdom} {channel} {value}")
    
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(single_measurements, bins=200,histtype='step', alpha=0.7, color='blue')
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time [ns]", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"mDOM_transit_time_multiple_measuremnet.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"mDOM_transit_time_multiple_measuremnet.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_transit_time_multiple_measurement('/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json')


