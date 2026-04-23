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





#############################################
from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsCustom4 = colorsCustom2 + colorsCustom2
colorsCustom16 = colorsCustom4 + colorsCustom4 + colorsCustom4 + colorsCustom4
colorsCustom64 = colorsCustom16+colorsCustom16+colorsCustom16+colorsCustom16
colorsIter = iter(colorsCustom)
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
colors = ['#1e4c7b','#498c9c','#89b3a2','#f1d27e','#d56b48','#1e1e3e',"#5a5a8c","#9c7bbc","#d45e7d",'#f7b1a1']
colornames = ["Blue, Bright, Coastal    Green    Modern    Orange    Spring    Teal    Vibrant    Yellow,    70s",
              "Aesthetic    Bright    Deep    Indigo    Purple    Red    Soft    Summer    Vintage Starman"]

R12199_tt = 43 #ns
R12199_tts = 1.92 #ns
R15458_02_tts = 1.4 #ns



def extract_histogram(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    meas_time = data["meas_time"]
    y_values = data["meas_data"][0]['y_values']
    x_min = data["meas_data"][0]["x_min"]
    x_max = data["meas_data"][0]["x_max"]
    n_bins = data["meas_data"][0]["n_bins"]
    x_values = np.linspace(x_min,x_max,n_bins)
    x_label = data["meas_data"][0]["x_label"]
    device_uid = data["device_uid"]
    pmt = data['subdevice_uid'].split('_')[-1]
    channel = data["meas_data"][2]["value"]
    run = data['run_number']
    temperature = data["meas_data"][0]["temperature"]
    # print(data["meas_data"][])

    # fit_x_values = np.linspace(x_min,x_max,n_bins)
    # fit_y_values = data["meas_data"][0]["fit_y_values"]
    return x_values,y_values,x_label,device_uid,pmt,channel,run,temperature,meas_time,


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
 

def extract_temperature(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["meas_data"][0]["temperature"]

def extract_HV(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    hv = -100
    if data["meas_data"][1]["label"] == "applied HV":
        hv = data["meas_data"][1]["value"]
    else:
        print(f"{data["meas_data"][1]["label"]} is not applied HV")
    return hv

def get_chi2(json_file):
    '''
    extract chi2 value for goodness of fit
    '''
    with open(json_file, 'r') as f:
        data = json.load(f)
    y_values = data["meas_data"][0]['y_values']
    fit_y_values = data["meas_data"][0]["fit_y_values"]
    chi2,pvalue = stats.chisquare(y_values,fit_y_values,ddof=3,sum_check=True)
    return chi2,pvalue

def gauss(x, A, mu, sigma):
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def gaussian_fit(x,y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[ max(y), mean, sigma])
    return popt


def mdom_transit_plot(file_list):
    n_files = len(file_list)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for i,ifile in enumerate(file_list):
        # print(f"no of file {i}")
        x_values,y_values,x_label,device_uid,pmt,channel,run,temperature,meas_time = extract_histogram(ifile)
        popt = gaussian_fit(x_values,y_values)
        # print(f"popt {popt}")
        ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom64[i],label=f"{str(device_uid)} {str(pmt)} R{run} T{temperature:.2f}",alpha=1)    
        # ax.plot(fit_x_values,fit_y_values,ls='--',lw = 2.5,c=colorsCustom2[i],alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"{}".format(x_label), fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.00),loc="lower center")
    plt.savefig(plotFolder+f"/mDOM_transit_time{device_uid}mDOM.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/mDOM_transit_time{device_uid}mDOM.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def make_transit_plots(mdom_list):
    for imdom in mdom_list[:]:
        ifiles = sorted(glob.glob(imdom+"/m*"))
        if len(ifiles)>0:
            mdom_transit_plot(ifiles)
        else:
            print(f"the files list is empty for {imdom}")

def check_meas_site(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data['meas_site']

def extract_run_number(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data['run_number']

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

def plot_mu(mu_list,b_list,meas_site):
    print(f"mu min {min(mu_list)}max{max(mu_list)}")
    print(f"b min {min(b_list)}max{max(b_list)}")
    print(np.asarray(mu_list)-np.asarray(b_list))
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(int(min(mu_list)-1),int(max(mu_list)+1),int((max(mu_list)-min(mu_list))/6)+1)
    # bins = np.linspace(0,100,401)
    bins = np.linspace(40,60,201)
    ax.hist(mu_list,bins=bins,histtype="step",color="orange",label=r"Transit time",lw=2.5,alpha=1)
    ax.hist(b_list,bins=bins,histtype="step",color="blue",label="b",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time peaks $\mu$ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=14)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_meansmDOM_{meas_site}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_meansmDOM_{meas_site}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_sigma(sigma_list,c_list,meas_site):
    print(f"sigma min {min(sigma_list)}max{max(sigma_list)}")
    print(f"c min {min(c_list)}max{max(c_list)}")
    print("diff", np.asarray(c_list)-np.asarray(sigma_list))
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(int(min(sigma_list)-1),int(max(sigma_list)+1),int((max(sigma_list)-min(sigma_list)))+1)
    # bins = np.linspace(-60,150,211)
    # bins = np.linspace(0,20,81)
    bins = np.linspace(0,5,51)
    print(bins)
    ax.hist(sigma_list,bins=bins,histtype="step",color="orange",lw=2.5,label="transit time spread",alpha=1)
    ax.hist(c_list,bins=bins,histtype="step",color="blue",lw=2.5,label="c",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time width $\sigma $ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=14)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_sigmamDOM_{meas_site}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_sigmamDOM_{meas_site}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_temp(temp_list,meas_site):
    print(f"temp min {min(temp_list)}max{max(temp_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(int(min(temp_list)-1),int(max(temp_list)+1),int((max(temp_list)-min(temp_list))/0.5)+1)
    bins = np.linspace(-30,40,141)
    # print(temp_list)
    print("bins")
    ax.hist(temp_list,bins=bins,histtype="step",color="orange",linewidth=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"temperature [$^{\circ}$C]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_temp_mDOM_{meas_site}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_temp_mDOM_{meas_site}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_hv(hv_list,meas_site):
    print(f"temp min {min(hv_list)}max{max(hv_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(int(min(hv_list)-1),int(max(hv_list)+1),int((max(hv_list)-min(hv_list))/0.5)+1)
    bins = np.linspace(75,105,31)
    # print(hv_list)
    print("bins")
    ax.hist(hv_list,bins=bins,histtype="step",color="orange",linewidth=2.5,alpha=1)
    # ax.hist(hv_list,histtype="step",color="orange",lw=4.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"applied HV [V]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_HV_mDOM_{meas_site}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_HV_mDOM_{meas_site}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_chi2(chi2_list,meas_site):
    print(f"min {min(chi2_list)}max{max(chi2_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(int(min(chi2_list)-1),int(max(chi2_list)+1),int((max(chi2_list)-min(chi2_list)))+1)
    ax.hist(chi2_list,bins=bins,histtype="step",color="orange",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" $\chi^{2} $ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../chi2_distmDOM_{meas_site}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_distmDOM_{meas_site}.pdf",transparent=False,bbox_inches='tight')
    plt.close()



def plot_chi2Sigma(chi2_list,sigma_list,meas_site):
    print(f"min {min(chi2_list)}max{max(chi2_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.plot(chi2_list,sigma_list,"o")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" $\chi^{2} $ [ns]", fontsize=22)
    ax.set_ylabel(r"$\sigma$", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../chi2_sigmamDOM_{meas_site}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_sigmamDOM_{meas_site}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_chi2mu(chi2_list,mu_list,meas_site):
    print(f"min {min(chi2_list)}max{max(chi2_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.plot(chi2_list,mu_list,"o")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" $\chi^{2} $ [ns]", fontsize=22)
    ax.set_ylabel(r"$\mu$", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../chi2_mumDOM_{meas_site}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_mumDOM_{meas_site}.pdf",transparent=False,bbox_inches='tight')
    plt.close()



def transit_params(mdom_list,meas_site):
    file_list = []
    for imdom in mdom_list[:]:
        ifiles = sorted(glob.glob(imdom+"/m*"))
        file_list += ifiles
    a_list = []
    b_list = []
    c_list = []
    mu_list = []
    sigma_list = []
    chi2_list = []
    pvalue_list = []
    zero_mu_files = []
    files_at_site = []
    temp_list = []
    hv_list = []
    for ifile in file_list:
        # print(f"running {ifile}")
        if check_meas_site(ifile) == meas_site:
            files_at_site.append(ifile)

            mu,sigma,a,b,c,chi2,ihv = extract_fit_params(ifile)
            itemp = extract_temperature(ifile)
            # if 88<ihv<91:
            #     print(f"HV close to 89 {ifile}") 

            # extract_fit_params(ifile)
            # print(chi2,pvalue)
            # chi2,pvalue  = get_chi2(ifile)
            if mu > -500 and sigma < 500 and not np.isnan(mu) and not np.isnan(sigma)\
                and not np.isnan(a) and not np.isnan(b) and not np.isnan(c) and not np.isnan(chi2):
                mu_list.append(mu)
                sigma_list.append(sigma)
                chi2_list.append(chi2)
                a_list.append(a)
                b_list.append(b)
                c_list.append(c)
                temp_list.append(itemp)
                hv_list.append(ihv)
                if c!=sigma:
                    print(f"{c}c and {sigma}")
            elif abs(mu-0.0)< 0.001 or abs(sigma-0.0)<0.001:
                zero_mu_files.append(ifile)
            else:
                print(f"please check {ifile}. It has {mu} transit time and {sigma} spread")
    print(f"There were {len(zero_mu_files)} empty TT meas out of {len(files_at_site)} files at {meas_site}")
    # print(mu_list)
    plot_mu(mu_list,b_list,meas_site)
    plot_sigma(sigma_list,c_list,meas_site)
    plot_chi2(chi2_list,meas_site)
    plot_chi2mu(chi2_list,mu_list,meas_site)
    plot_chi2Sigma(chi2_list,sigma_list,meas_site)
    plot_temp(temp_list,meas_site)
    plot_hv(hv_list,meas_site)
    return a_list,mu_list,sigma_list

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

def get_string_number(prod_id):
    string_list = []
    for ifile in geometry_files:
        with open(ifile, 'r') as f:
            data = json.load(f)
        string_number = data[0]["id"]
        for idev in data[0]["devices"]:
            if idev["production_id"] == prod_id:
                string_list.append(string_number)
    if len(string_list)>1:
        print(f"more than one string number found for prod id {prod_id} {string_list}")
    return string_list[0]

def extract_transit_time(mdom_list):
    '''
    gets mean and variance from json file in the given list of files.
    '''
    mdom_dict = {}
    for imdom in mdom_list[:]:
        mdom_name = imdom.split("/")[-1]
        mdom_name = mdom_name.split("_")[0]+"_"+mdom_name.split("_")[1]
        icm_id = prod_id_to_icm_id(mdom_name)
        # print(f"mDOM name {mdom_name} {icm_id}")
        ifiles = sorted(glob.glob(imdom+"/m*_pmt_*.json"))
        PMT_list = [re.search(r"(DM\d+)", ifilename).group(1) for ifilename in ifiles]
        PMT_list = list(set(PMT_list))
        # print(f"pmt list   {len(PMT_list)} {PMT_list}")
        pmt_dict = {}
        pmt_dict_new = {}
        for ipmt in PMT_list:
            # pmt_dict = {}
            ifiles_this_pmt = sorted(glob.glob(imdom+f"/m*_pmt_{ipmt}*.json"))
            # print(f"files for {ipmt} {len(ifiles_this_pmt)}")
            mu_list = []
            sigma_list = []
            hv_list = []
            tt_info_map = []
            for ifile in ifiles_this_pmt:
                mu,sigma,a,b,c,chi2,ihv = extract_fit_params(ifile)
                itemp = extract_temperature(ifile)
                channel = extract_channel(ifile)
                hv_list.append(ihv)
                run_number = extract_run_number(ifile)
                tt_info = {"run_number": run_number, "mu": mu, "sigma": sigma,"a": a, "b": b, "c": c, "chi2": chi2,"applied HV": ihv, "temperature": itemp}
                tt_info_map.append(tt_info)
                # if mu > -500 and sigma < 500 and not np.isnan(mu) and not np.isnan(sigma)\
                #     and not np.isnan(a) and not np.isnan(b) and not np.isnan(c) and not np.isnan(chi2):
                # if mu > -500 and sigma < 500 and not np.isnan(mu) and not np.isnan(sigma):
                mu_list.append(mu)
                sigma_list.append(sigma)
            # if len(mu_list)>0:
            #     if len(mu_list)==1:
            #         tt = mu_list[0]
            #         best_mu = tt
            #     else:
            #         tt = min(mu_list, key=lambda x: abs(x - R12199_tt))
            #         # print(f"mDOM {mdom_name} pmt {ipmt} has {len(mu_list)} measurements with mean {tt:.1f} and std {np.std(mu_list):.1f}")
            # best_mu = np.median(mu_list)
            # print(f"mDOM {mdom_name} pmt {ipmt} has {len(mu_list)} measurements")
            # print(f"mDOM {mdom_name} pmt {ipmt} channel {channel} has mean transit time {tt:.1f} and {mu_list}")
            pmt_dict[f"channel_{channel}"] = [mu_list,hv_list]
            pmt_dict_new[f"channel_{channel}"] = tt_info_map
        # mdom_dict[mdom_name] = pmt_dict
        # print(f"PMT dict new for mDOM {mdom_name} {pmt_dict_new}")
        channels = [f"channel_{i}" for i in range(0,24)]
        pmt_dict_ordered = {channel: pmt_dict_new.get(channel, []) for channel in channels}
        # print(f"pmt_dict_ordered for mDOM {mdom_name}: {pmt_dict_ordered}")
        mdom_dict[mdom_name] = {"icm_id": icm_id, "transit_times": pmt_dict_ordered}
    return mdom_dict


            


def get_transit_time_dict(mdom_list):
    mdom_dict = extract_transit_time(mdom_list)
    # print(mdom_dict)
    mDOMs = sorted(mdom_dict.keys())
    # print(mDOMs)
    channels = [f"channel_{i}" for i in range(0,24)]
    # print(channels)
    # ordered_mdom_dict = {mdom: {channel: mdom_dict[mdom].get(channel, np.nan) for channel in channels} for mdom in mDOMs}
    # print(ordered_mdom_dict)
    ordered_mdom_dict = mdom_dict
    with open('/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json', 'w') as f:
        json.dump(ordered_mdom_dict, f, indent=4)


def get_mdom_production_id(mdom_path_name):
    '''
    extract production id from mdom path name
    '''
    return mdom_path_name.split("/")[-1].split("_")[0]+"_"+mdom_path_name.split("/")[-1].split("_")[1]




def get_mDOM_missing_channels(mdom_list):
    '''
    gets mean and variance from json file in the given list of files.
    '''
    mdom_dict = {}
    for imdom in mdom_list[:]:
        mdom_name = imdom.split("/")[-1]
        mdom_name = mdom_name.split("_")[0]+"_"+mdom_name.split("_")[1]
        icm_id = prod_id_to_icm_id(mdom_name)
        # print(f"mDOM name {mdom_name} {icm_id}")
        ifiles = sorted(glob.glob(imdom+"/m*_pmt_*.json"))
        PMT_list = [re.search(r"(DM\d+)", ifilename).group(1) for ifilename in ifiles]
        PMT_list = list(set(PMT_list))
        # print(f"pmt list   {len(PMT_list)} {PMT_list}")
        pmt_dict = {}
        pmt_dict_new = {}
        channel_list = []
        for ipmt in PMT_list:
            # pmt_dict = {}
            ifiles_this_pmt = sorted(glob.glob(imdom+f"/m*_pmt_{ipmt}*.json"))
            # print(f"files for {ipmt} {len(ifiles_this_pmt)}")
            mu_list = []
            sigma_list = []
            hv_list = []
            tt_info_map = []
            for ifile in ifiles_this_pmt:
                mu,sigma,a,b,c,chi2,ihv = extract_fit_params(ifile)
                itemp = extract_temperature(ifile)
                channel = extract_channel(ifile)
                channel_list.append(channel)
                hv_list.append(ihv)
                run_number = extract_run_number(ifile)
                tt_info = {"run_number": run_number, "mu": mu, "sigma": sigma, "chi2": chi2,"applied HV": ihv, "temperature": itemp}
                tt_info_map.append(tt_info)
                # if mu > -500 and sigma < 500 and not np.isnan(mu) and not np.isnan(sigma)\
                #     and not np.isnan(a) and not np.isnan(b) and not np.isnan(c) and not np.isnan(chi2):
                # if mu > -500 and sigma < 500 and not np.isnan(mu) and not np.isnan(sigma):
                mu_list.append(mu)
                sigma_list.append(sigma)
            # if len(mu_list)>0:
            #     if len(mu_list)==1:
            #         tt = mu_list[0]
            #         best_mu = tt
            #     else:
            #         tt = min(mu_list, key=lambda x: abs(x - R12199_tt))
            #         # print(f"mDOM {mdom_name} pmt {ipmt} has {len(mu_list)} measurements with mean {tt:.1f} and std {np.std(mu_list):.1f}")
            # best_mu = np.median(mu_list)
            # print(f"mDOM {mdom_name} pmt {ipmt} has {len(mu_list)} measurements")
            # print(f"mDOM {mdom_name} pmt {ipmt} channel {channel} has mean transit time {tt:.1f} and {mu_list}")
            pmt_dict[f"channel_{channel}"] = [mu_list,hv_list]
            pmt_dict_new[f"channel_{channel}"] = tt_info_map
        channel_list = list(set(channel_list))
        if len(channel_list)<24:
            missing_channels = [f"channel_{i}" for i in range(0,24) if i not in channel_list]
            # print(f"mDOM {mdom_name} is missing transit time measurement for channels {len(missing_channels)} {missing_channels}")
            print(f"S{get_string_number(mdom_name)} {mdom_name} ({prod_id_to_icm_id(mdom_name)}) {missing_channels[0]}")
        # mdom_dict[mdom_name] = pmt_dict
        # print(f"PMT dict new for mDOM {mdom_name} {pmt_dict_new}")
        channels = [f"channel_{i}" for i in range(0,24)]
        pmt_dict_ordered = {channel: pmt_dict_new.get(channel, np.nan) for channel in channels}
        # print(f"pmt_dict_ordered for mDOM {mdom_name}: {pmt_dict_ordered}")
        mdom_dict[mdom_name] = {"icm_id": icm_id, "transit_times": pmt_dict_ordered}
    return mdom_dict


get_mDOM_missing_channels(combined_mdom_list_deployed)

def get_device_list(string,device,geometry_files):
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

def main():
    plotFolder = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
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

    # print(geometry_files)

    # string_list = ["87","88","89","90","91","92"]
    string_list = ["88","89","90","91","92"]


    # deployed_device_list = get_device_list("87","mDOM") + get_device_list("88","mDOM") + get_device_list("89","mDOM") + get_device_list("90","mDOM") + get_device_list("91","mDOM") + get_device_list("92","mDOM")
    deployed_device_list = get_device_list("88","mDOM",geometry_files) + get_device_list("89","mDOM",geometry_files) + get_device_list("90","mDOM",geometry_files) + get_device_list("91","mDOM",geometry_files) + get_device_list("92","mDOM",geometry_files)
    # print(f"device list {len(device_list)} {device_list}")

    print(f"deployed device in string 88 {len(get_device_list('88','mDOM',geometry_files))}")
    print(f"deployed device in string 89 {len(get_device_list('89','mDOM',geometry_files))}")
    print(f"deployed device in string 90 {len(get_device_list('90','mDOM',geometry_files))}")
    print(f"deployed device in string 91 {len(get_device_list('91','mDOM',geometry_files))}")
    print(f"deployed device in string 92 {len(get_device_list('92','mDOM',geometry_files))}")

    print(f"deployed mDOM list {len(deployed_device_list)}")
    combined_mdom_list = mdom_list_msu + mdom_list_desy
    # combined_mean_plot(mdom_list_desy,mdom_list_msu)
    combined_mdom_list_deployed = [imdom for imdom in combined_mdom_list if get_mdom_production_id(imdom) in deployed_device_list]

    deployed_mdom_missing_tt = [imdom for imdom in deployed_device_list if imdom not in [get_mdom_production_id(idom) for idom in combined_mdom_list]]
    deployed_mdom_missing_tt_icm_id = [prod_id_to_icm_id(imdom) for imdom in deployed_mdom_missing_tt]

    print(f"deployed mDOMs missing transit time measurement {len(deployed_mdom_missing_tt)} {deployed_mdom_missing_tt} {deployed_mdom_missing_tt_icm_id}")

    for prod_id, icm_id in zip(deployed_mdom_missing_tt, deployed_mdom_missing_tt_icm_id):
        print(f"{prod_id} ({icm_id})")

    # print(f"total number of deployed mDOMs {len(combined_mdom_list_deployed)}")

    print(f"FAT mDOM list of deployed mDOMs has {len(combined_mdom_list_deployed)} entries")
    get_transit_time_dict(combined_mdom_list_deployed) #gives all runs of transit time measurement





    print(f"Combined mDOM list has {len(combined_mdom_list)} entries")



    # for imdom in combined_mdom_list:
    #     print(f"{imdom} has production id {get_mdom_production_id(imdom)}")

    print(f"deployed mdoms {len(deployed_device_list)} {len(combined_mdom_list)} {len([get_mdom_production_id(imdom) for imdom in combined_mdom_list if get_mdom_production_id(imdom) in deployed_device_list])}")

    missing_FAT_entries = [imdom for imdom in deployed_device_list if get_mdom_production_id(imdom) not in [get_mdom_production_id(idom) for idom in combined_mdom_list]]
    print(f"missing FAT entries {len(missing_FAT_entries)} {missing_FAT_entries}")

if __name__ == "__main__":
    main()


# print(f" 87 {get_device_list("87","mDOM")}")
# print(f" 88 {get_device_list("88","mDOM")}")
# print(f" 89 {get_device_list("89","mDOM")}")
# print(f" 90 {get_device_list("90","mDOM")}")
# print(f" 91 {get_device_list("91","mDOM")}")
# print(f" 92 {get_device_list("92","mDOM")}")


# def list_mDOMs_missing_transit_time(string,mdom_list):
#     print(f"working on string {string}")
#     missing_mdoms = []
#     select_mdoms = []
#     deployed_device_list = get_device_list(f"{string}","mDOM")
#     for imdom in mdom_list:
#         production_id = get_mdom_production_id(imdom)
#         if production_id in deployed_device_list:
#             select_mdoms.append(imdom)
#     single_hv_missing = []

#     for imdom in select_mdoms[:]:
#         mdom_name = imdom.split("/")[-1]
#         mdom_name = mdom_name.split("_")[0]+"_"+mdom_name.split("_")[1]
#         ifiles = sorted(glob.glob(imdom+"/m*_pmt_*.json"))
#         PMT_list = [re.search(r"(DM\d+)", ifilename).group(1) for ifilename in ifiles]
#         PMT_list = list(set(PMT_list))

#         # print(f"pmt list   {len(PMT_list)} {PMT_list}")
#         pmt_dict = {}
#         channel_list = []
#         for ipmt in PMT_list:
#             # pmt_dict = {}
#             ifiles_this_pmt = sorted(glob.glob(imdom+f"/m*_pmt_{ipmt}*.json"))
#             # print(f"files for {ipmt} {len(ifiles_this_pmt)}")
#             mu_list = []
#             sigma_list = []
#             hv_list = []
#             missing_transit_mdoms = []
#             for ifile in ifiles_this_pmt:
#                 mu,sigma,a,b,c,chi2,ihv = extract_fit_params(ifile)
#                 itemp = extract_temperature(ifile)
#                 channel = extract_channel(ifile)
#                 channel_list.append(channel)
#                 hv_list.append(ihv)

#                 # if np.isnan(ihv):
#                     # single_hv_missing.append(mdom_name)


#                 # if mu > -500 and sigma < 500 and not np.isnan(mu) and not np.isnan(sigma)\
#                 #     and not np.isnan(a) and not np.isnan(b) and not np.isnan(c) and not np.isnan(chi2):
#                 # if mu > -500 and sigma < 500 and not np.isnan(mu) and not np.isnan(sigma):
#                 mu_list.append(mu)
#                 sigma_list.append(sigma)
#             # print(len(mu_list),mu_list[0])
#             # if mdom_name == "mDOM_M137":
#             #     print(f"mDOM {mdom_name} pmt {channel} check in nan {mu_list[0]}")
#             #     print(mu_list)
#             if len(mu_list)==1 and np.isnan(hv_list[0]):
#                 single_hv_missing.append(mdom_name)
#         # print(f"channel list {channel_list}")
#         for i in range(24):
#             if i not in channel_list:
#                 missing_transit_mdoms.append(mdom_name+"ch"+str(i))
#             #     print(f"mDOM {mdom_name} pmt {ipmt} has {len(mu_list)} measurements with mean {mu_list} and hv {hv_list}")
#         if len(missing_transit_mdoms)>0:
#             print(f"missing transit time for {mdom_name} {missing_transit_mdoms}")
#     single_hv_missing = list(set(single_hv_missing))
#     print(f"single hv missing {single_hv_missing}")
#     # if len(single_hv_missing)>0:
#     #     print(f'mDOMs with single measurement and missing HV {list(set(single_hv_missing))}')
#     missing_mdoms = list(set(missing_mdoms))
#     if len(missing_mdoms)>0:
#         print(f'mDOMs channel missing transit time {missing_mdoms}')

#     return missing_mdoms

# # list_mDOMs_missing_transit_time("87",combined_mdom_list)
# # list_mDOMs_missing_transit_time("88",combined_mdom_list)
# # list_mDOMs_missing_transit_time("89",combined_mdom_list)
# # list_mDOMs_missing_transit_time("90",combined_mdom_list)
# # list_mDOMs_missing_transit_time("91",combined_mdom_list)
# # list_mDOMs_missing_transit_time("92",combined_mdom_list)

# def mDOMs_with_all_channel_tt_hv(string,mdom_list):
#     print(f"working on string {string}")
#     missing_mdoms = []
#     select_mdoms = []
#     deployed_device_list = get_device_list(f"{string}","mDOM")
#     for imdom in mdom_list:
#         production_id = get_mdom_production_id(imdom)
#         if production_id in deployed_device_list:
#             select_mdoms.append(imdom)
#     single_hv_missing = []
#     all_channel_mdoms = []
#     mu_list_string = []

#     for imdom in select_mdoms[:]:
#         mdom_name = imdom.split("/")[-1]
#         mdom_name = mdom_name.split("_")[0]+"_"+mdom_name.split("_")[1]
#         ifiles = sorted(glob.glob(imdom+"/m*_pmt_*.json"))
#         PMT_list = [re.search(r"(DM\d+)", ifilename).group(1) for ifilename in ifiles]
#         PMT_list = list(set(PMT_list))

#         # print(f"pmt list   {len(PMT_list)} {PMT_list}")
#         pmt_dict = {}
#         channel_list = []

#         for ipmt in PMT_list:
#             # pmt_dict = {}
#             ifiles_this_pmt = sorted(glob.glob(imdom+f"/m*_pmt_{ipmt}*.json"))
#             # print(f"files for {ipmt} {len(ifiles_this_pmt)}")
#             mu_list = []
#             sigma_list = []
#             hv_list = []
#             missing_transit_mdoms_channel = []
#             for ifile in ifiles_this_pmt:
#                 mu,sigma,a,b,c,chi2,ihv = extract_fit_params(ifile)
#                 if not np.isnan(ihv) and not np.isnan(mu) and 0<mu<100:
#                     hv_list.append(ihv)
#                     mu_list.append(mu)
#                     mu_list_string.append(mu)
#                     itemp = extract_temperature(ifile)
#                     channel = extract_channel(ifile)
#                     channel_list.append(channel)
            
#         for i in range(24):
#             if i not in channel_list:
#                 missing_transit_mdoms_channel.append(mdom_name+"ch"+str(i))
#             #     print(f"mDOM {mdom_name} pmt {ipmt} has {len(mu_list)} measurements with mean {mu_list} and hv {hv_list}")
#         if len(missing_transit_mdoms_channel)>0:
#             # print(f"missing transit time for {mdom_name} channels {len(missing_transit_mdoms_channel)} {missing_transit_mdoms_channel}")
#             print(f"missing transit time for {mdom_name} channels {len(missing_transit_mdoms_channel)}")
#         if len(set(channel_list))==24:
#             all_channel_mdoms.append(mdom_name)
#     print(f"string {string} mDOMs with all channel transit time and hv {len(all_channel_mdoms)}")
#     # single_hv_missing = list(set(single_hv_missing))
#     # print(f"single hv missing {single_hv_missing}")
#     # # if len(single_hv_missing)>0:
#     # #     print(f'mDOMs with single measurement and missing HV {list(set(single_hv_missing))}')
#     # missing_mdoms = list(set(missing_mdoms))
#     # if len(missing_mdoms)>0:
#     #     print(f'mDOMs channel missing transit time {missing_mdoms}')

#     # return missing_mdoms
#     return mu_list_string


# # mu_list_88 = mDOMs_with_all_channel_tt_hv("88",combined_mdom_list)
# # mu_list_89 = mDOMs_with_all_channel_tt_hv("89",combined_mdom_list)
# # mu_list_90 = mDOMs_with_all_channel_tt_hv("90",combined_mdom_list)
# # mu_list_91 = mDOMs_with_all_channel_tt_hv("91",combined_mdom_list)
# # mu_list_92 = mDOMs_with_all_channel_tt_hv("92",combined_mdom_list)

# # total_mu_list = mu_list_88 + mu_list_89 + mu_list_90 + mu_list_91 + mu_list_92
# # print(f"total mu list {len(total_mu_list)} {np.min(total_mu_list)} {np.max(total_mu_list)}")


# def mDOMs_transit_time(string,mdom_list):
#     print(f"working on string {string}")
#     missing_mdoms = []
#     select_mdoms = []
#     deployed_device_list = get_device_list(f"{string}","mDOM")
#     for imdom in mdom_list:
#         production_id = get_mdom_production_id(imdom)
#         if production_id in deployed_device_list:
#             select_mdoms.append(imdom)
#     single_hv_missing = []
#     all_channel_mdoms = []
#     mu_list_string = []

#     for imdom in select_mdoms[:]:
#         mdom_name = imdom.split("/")[-1]
#         mdom_name = mdom_name.split("_")[0]+"_"+mdom_name.split("_")[1]
#         ifiles = sorted(glob.glob(imdom+"/m*_pmt_*.json"))
#         PMT_list = [re.search(r"(DM\d+)", ifilename).group(1) for ifilename in ifiles]
#         PMT_list = list(set(PMT_list))

#         # print(f"pmt list   {len(PMT_list)} {PMT_list}")
#         pmt_dict = {}
#         channel_list = []

#         for ipmt in PMT_list:
#             # pmt_dict = {}
#             ifiles_this_pmt = sorted(glob.glob(imdom+f"/m*_pmt_{ipmt}*.json"))
#             # print(f"files for {ipmt} {len(ifiles_this_pmt)}")
#             mu_list = []
#             sigma_list = []
#             hv_list = []
#             missing_transit_mdoms_channel = []
#             for ifile in ifiles_this_pmt:
#                 mu,sigma,a,b,c,chi2,ihv = extract_fit_params(ifile)
#                 if not np.isnan(mu) and 0<mu<100:
#                     hv_list.append(ihv)
#                     mu_list.append(mu)
#                     mu_list_string.append(mu)
#                     itemp = extract_temperature(ifile)
#                     channel = extract_channel(ifile)
#                     channel_list.append(channel)
            
#         for i in range(24):
#             if i not in channel_list:
#                 missing_transit_mdoms_channel.append(mdom_name+"ch"+str(i))
#             #     print(f"mDOM {mdom_name} pmt {ipmt} has {len(mu_list)} measurements with mean {mu_list} and hv {hv_list}")
#         if len(missing_transit_mdoms_channel)>0:
#             # print(f"missing transit time for {mdom_name} channels {len(missing_transit_mdoms_channel)} {missing_transit_mdoms_channel}")
#             print(f"missing transit time for {mdom_name} channels {len(missing_transit_mdoms_channel)}")
#         if len(set(channel_list))==24:
#             all_channel_mdoms.append(mdom_name)
#     print(f"string {string} mDOMs with all channel transit time and hv {len(all_channel_mdoms)}")
#     # single_hv_missing = list(set(single_hv_missing))
#     # print(f"single hv missing {single_hv_missing}")
#     # # if len(single_hv_missing)>0:
#     # #     print(f'mDOMs with single measurement and missing HV {list(set(single_hv_missing))}')
#     # missing_mdoms = list(set(missing_mdoms))
#     # if len(missing_mdoms)>0:
#     #     print(f'mDOMs channel missing transit time {missing_mdoms}')

#     # return missing_mdoms
#     return mu_list_string


# mu_list_88 = mDOMs_transit_time("88",combined_mdom_list)
# mu_list_89 = mDOMs_transit_time("89",combined_mdom_list)
# mu_list_90 = mDOMs_transit_time("90",combined_mdom_list)
# mu_list_91 = mDOMs_transit_time("91",combined_mdom_list)
# mu_list_92 = mDOMs_transit_time("92",combined_mdom_list)

# total_mu_list = mu_list_88 + mu_list_89 + mu_list_90 + mu_list_91 + mu_list_92
# print(f"total mu list {len(total_mu_list)} {np.min(total_mu_list)} {np.max(total_mu_list)}")



# def transit_time_hist(tt_list):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     ax.hist(tt_list, bins=200,histtype='step',color=colorsCustom[0],linewidth=1.5 ,label=r'$tt_{FAT}$', alpha=0.7)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#     ax.set_xlabel(r"transit time [ns]", fontsize=22)
#     ax.set_ylabel(r"count", fontsize=22)
#     ax.set_yscale('log')
#     ax.text(0.95, 0.95, fr"mean $tt_{{{'FAT'}}}$: {np.mean(tt_list):.0f}"+r"$\pm$" + fr" {np.std(tt_list):.0f}  ns", transform=ax.transAxes, ha='right', va='top')
#     ax.grid(True,alpha=0.6)
#     ax.legend(fontsize=12,ncols=1)
#     plt.savefig(plotFolder+f"/transit_time_no_hv_selection.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/transit_time_no_hv_selection.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
    
# # transit_time_hist(total_mu_list)

