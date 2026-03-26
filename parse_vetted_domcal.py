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

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']

plotFolder = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
domcal_data = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/research_notes/notes/vetted_domcal/vetted domcal JSON/"
# domcal_data = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/research_notes/notes/vetted_domcal/vetted domcal JSON/0a00000f17ae622d_degg.json"

mdom_domcal_files = sorted(glob.glob(domcal_data+"/*_mdom.json"))


mdom_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/") if f.is_dir()][:]
mdom_list_desy = []
mdom_list_desy_dvt = []
mdom_list_msu = []
for ifile in mdom_list:
    if ifile.split("/")[-1].split("_")[1][0]=="M":
        mdom_list_msu.append(ifile)
    elif "DVT" in ifile.split("/")[-1].split("_")[1]:
        mdom_list_desy_dvt.append(ifile)
    else:
        mdom_list_desy.append(ifile)


mdom_names = [istr.split("/")[-1] for istr in mdom_list]


upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"

geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))


string_list = ["88","89","90","91","92"]

def get_device_list(string,device):
    device_list = []
    for ifile in geometry_files:
        if string in ifile:
            # print(f"ifile {ifile}")
            with open(ifile, 'r') as f:
                data = json.load(f)
                # print(f"data {data[0]}")
            for idev in data[0]["devices"]:
                # print(f"this device #################################")
                # print(f"idev {idev}")                
                if device in idev["production_id"]:
                    device_list.append(idev["production_id"])
            # device_list.append(ifile.split("/")[-1].split("_")[0]+"_"+ifile.split("/")[-1].split("_")[1])
    # print(f"device list for {device} on string {string} {len(device_list)} {device_list}")
    return device_list

# deployed_device_list = get_device_list("87","mDOM") + get_device_list("88","mDOM") + get_device_list("89","mDOM") + get_device_list("90","mDOM") + get_device_list("91","mDOM") + get_device_list("92","mDOM")
deployed_device_list = get_device_list("88","mDOM") + get_device_list("89","mDOM") + get_device_list("90","mDOM") + get_device_list("91","mDOM") + get_device_list("92","mDOM")

combined_mdom_list = mdom_list_msu + mdom_list_desy
# combined_mean_plot(mdom_list_desy,mdom_list_msu)
def get_mdom_production_id(mdom_path_name):
    '''
    extract production id from mdom path name
    '''
    return mdom_path_name.split("/")[-1].split("_")[0]+"_"+mdom_path_name.split("/")[-1].split("_")[1]



combined_mdom_list_deployed = [imdom for imdom in combined_mdom_list if get_mdom_production_id(imdom) in deployed_device_list]


print(f"combined mdom list {len(combined_mdom_list_deployed)}")


def prod_id_to_icm_id(prod_id):
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

print(f"combined mdom list file {combined_mdom_list[0]} {get_mdom_production_id(combined_mdom_list[0])} {prod_id_to_icm_id(get_mdom_production_id(combined_mdom_list[0]))}")


def get_vetted_domcal_files(icm_id):
    '''
    get vetted domcal files for a given icm id
    '''
    return domcal_data+f"{icm_id}_mdom.json"

print(f"vetted domcal file {get_vetted_domcal_files(prod_id_to_icm_id(get_mdom_production_id(combined_mdom_list[0])))}")

nominal_gain = 5e6

def get_nominal_hv_icm_id(icm_id,channel):
    '''
    get nominal hv for a given icm id from geometry file
    '''
    domcal_file = get_vetted_domcal_files(icm_id)
    if not os.path.exists(domcal_file):
        print(f"vetted domcal file {domcal_file} does not exist")
        return np.nan
    else:
        with open(domcal_file, 'r') as f:
                data = json.load(f)
                slope = data['vetting'][str(channel)]['gain_slope']
                intercept = data['vetting'][str(channel)]['gain_intercept']
                nominal_hv = 10**((np.log10(nominal_gain) - intercept)/slope)
        return nominal_hv

def get_nominal_hv_prod_id(prod_id,channel):
    '''
    get nominal hv for a given production id and channel
    '''
    icm_id = prod_id_to_icm_id(prod_id)
    return get_nominal_hv_icm_id(icm_id,channel)

print(f"nominal hv for {prod_id_to_icm_id(get_mdom_production_id(combined_mdom_list[0]))} is {get_nominal_hv_icm_id(prod_id_to_icm_id(get_mdom_production_id(combined_mdom_list[0])),channel=0)}")

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


def get_fat_transit_time(prod_id,channel):
    '''
    get transit time for a given production id and channel
    '''
    transit_time_dir = [ifile for ifile in combined_mdom_list_deployed if prod_id in ifile][0]
    ifiles = sorted(glob.glob(transit_time_dir+"/m*_pmt_*.json"))
    ifiles = [ifile for ifile in ifiles if extract_channel(ifile) == channel]
    mu_list = []
    hv_list = []
    for ifile in ifiles:
        mu,sigma,a,b,c,chi2,ihv = extract_fit_params(ifile)
        if not np.isnan(ihv) and not np.isnan(mu) and 0<mu<100:
            # print(f"transit time for {prod_id} channel {channel} is {mu} ns at hv {ihv} V")
            mu_list.append(mu)
            hv_list.append(ihv)
    if len(mu_list)>0:
        return mu_list[0],hv_list[0]
    else:
        return np.nan, np.nan




def get_corrected_transit_time(prod_id,channel):
    '''
    get corrected transit time for a given production id and channel
    '''
    tt_fat, hv_fat = get_fat_transit_time(prod_id,channel)
    nominal_hv = get_nominal_hv_prod_id(prod_id,channel)
    tt_ice = tt_fat * (hv_fat/nominal_hv)**(0.5)
    return tt_ice

def get_both_transit_time(prod_id,channel):
    '''
    get corrected transit time for a given production id and channel
    '''
    tt_fat, hv_fat = get_fat_transit_time(prod_id,channel)
    nominal_hv = get_nominal_hv_prod_id(prod_id,channel)
    tt_ice = tt_fat * (hv_fat/nominal_hv)**(0.5)
    return tt_fat,tt_ice

tt_fat, hv_fat = get_both_transit_time(get_mdom_production_id(combined_mdom_list_deployed[0]),channel=0)
print(f"fat transit time {tt_fat} ns at hv {hv_fat}")


def transit_time_hist(mdom_list):
    transit_times_fat = []
    transit_times_ice = []

    for imdom in mdom_list:
        prod_id = get_mdom_production_id(imdom)
        for ichannel in range(24):
            tt_fat, hv_fat = get_fat_transit_time(prod_id,ichannel)
            nominal_hv = get_nominal_hv_prod_id(prod_id,ichannel)
            tt_ice = tt_fat * (hv_fat/nominal_hv)**(0.5)
            if not np.isnan(hv_fat) and not np.isnan(tt_fat) and not np.isnan(tt_ice) and 0<tt_fat<100:
                transit_times_fat.append(tt_fat)
                transit_times_ice.append(tt_ice)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(transit_times_fat, bins=200,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=r'$tt_{FAT}$', alpha=0.7)
    ax.hist(transit_times_ice, bins=200,histtype='step', color=colorsCustom[1],linewidth=1.5, label=r'$tt_{Ice}$ = $tt_{FAT}$ x $\sqrt{\frac{V_{FAT}}{V_{Ice}}}$', alpha=0.7)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time [ns]", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    ax.set_yscale('log')
    ax.text(0.95, 0.95, fr"mean $tt_{{{'FAT'}}}$: {np.mean(transit_times_fat):.0f}"+r"$\pm$" + fr" {np.std(transit_times_fat):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
    ax.text(0.95, 0.85, fr"mean $tt_{{{'Ice'}}}$: {np.mean(transit_times_ice):.0f}"+r"$\pm$" + fr" {np.std(transit_times_ice):.0f} ns", transform=ax.transAxes, ha='right', va='top')
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=12,ncols=1)
    plt.savefig(plotFolder+f"/fat_Ice_transit_time.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/fat_Ice_transit_time.pdf",transparent=False,bbox_inches='tight')
    plt.close()
    
# transit_time_hist(combined_mdom_list_deployed)


def transit_time_diff_hist(mdom_list):
    transit_times_diff = []

    for imdom in mdom_list:
        prod_id = get_mdom_production_id(imdom)
        for ichannel in range(24):
            tt_fat, hv_fat = get_fat_transit_time(prod_id,ichannel)
            nominal_hv = get_nominal_hv_prod_id(prod_id,ichannel)
            tt_ice = tt_fat * (hv_fat/nominal_hv)**(0.5)
            if not np.isnan(hv_fat) and not np.isnan(tt_fat) and not np.isnan(tt_ice) and 0<tt_fat<100:
                tt_diff = tt_ice - tt_fat
                transit_times_diff.append(tt_diff)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.hist(transit_times_diff, bins=200,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=f'diff', alpha=0.7)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(fr"$tt_{{{'Ice'}}}$ - $tt_{{{'FAT'}}}$", fontsize=22)
    ax.set_ylabel(r"count", fontsize=22)
    ax.set_yscale('log')
    ax.text(0.95, 0.95, fr"mean $tt_{{{'Diff'}}}$: {np.mean(transit_times_diff):.2f}"+f" \u00B1"+fr" {np.std(transit_times_diff):.2f} ns", transform=ax.transAxes, ha='right', va='top')
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=12,ncols=1)
    plt.savefig(plotFolder+f"/diff_transit_time.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/diff_transit_time.pdf",transparent=False,bbox_inches='tight')
    plt.close()
    
transit_time_diff_hist(combined_mdom_list_deployed)


