#!/usr/bin/env python

import os
import glob

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import json

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/degg_transit"
degg_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/degg_transit/") if f.is_dir()]
# print(degg_list)


from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']

R5912_100_tt = 54#ns
R5912_100_tts = 1.01918616#ns

def extract_histogram(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    meas_time = data["meas_time"]
    y_values = data["meas_data"][0]['y_values']
    x_min = data["meas_data"][0]["x_min"]
    x_max = data["meas_data"][0]["x_max"]
    x_values = np.linspace(x_min,x_max,99)
    x_label = data["meas_data"][0]["x_label"]
    device_uid = data["device_uid"]
    pmt = data['subdevice_uid'].split('_')[-1]
    run = data['run_number']
    temperature = data["meas_data"][0]["temperature"]
    # print(data["meas_data"][])

    fit_x_values = np.linspace(x_min,x_max,99)
    fit_y_values = data["meas_data"][0]["fit_y_values"]
    return x_values,y_values,fit_x_values,fit_y_values,x_label,device_uid,pmt,run,temperature,meas_time,

def extract_temperature(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["meas_data"][0]["temperature"]

def extract_device(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["device_uid"],data["subdevice_uid"]


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




def degg_transit_plot(file_list):
    n_files = len(file_list)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for i,ifile in enumerate(file_list):
        x_values,y_values,fit_x_values,fit_y_values,x_label,device_uid,pmt,run,temperature,meas_time = extract_histogram(ifile)
        ax.step(x_values,y_values,ls='-',lw = 2.5,c=colorsCustom2[i],label=f"{str(device_uid)} {str(pmt)} R{run} T{temperature:.2f}",alpha=1)    
        ax.plot(fit_x_values,fit_y_values,ls='--',lw = 2.5,c=colorsCustom2[i],alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"{}".format(x_label), fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.00),loc="lower center")
    plt.savefig(plotFolder+f"/transit_time{device_uid}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/transit_time{device_uid}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def make_transit_plots(degg_list):
    for idegg in degg_list[:]:
        ifiles = sorted(glob.glob(idegg+"/DEgg*"))
        degg_transit_plot(ifiles)

def extract_fit_params(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    y_values = data["meas_data"][0]['y_values']
    fit_y_values = data["meas_data"][0]["fit_y_values"]
    # chi2,pvalue = stats.chisquare(y_values,fit_y_values,ddof=3,sum_check=True)
    chi2,pvalue = stats.chisquare(y_values,fit_y_values,ddof=3,sum_check=False)
    fit_params = data["meas_data"][0]["fit_params"]
    return *fit_params,chi2,pvalue

def plot_temp(temp_list,temp_limits):
    print(f"min {min(temp_list)}max{max(temp_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(int(min(temp_list)-1),int(max(temp_list)+1),int((max(temp_list)-min(temp_list))/0.5)+1)
    bins = np.linspace(-30,40,141)
    ax.hist(temp_list,bins=bins,histtype="step",color="orange",linewidth=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"temperature [$^{\circ}$C]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_temp{temp_limits[0]}_{temp_limits[1]}C.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_temp{temp_limits[0]}_{temp_limits[1]}C.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_mu(mu_list,temp_limits):
    shift = False
    if shift == True:
        mu_list = [imu-12 for imu in mu_list]
    print(f"min {min(mu_list)}max{max(mu_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(int(min(mu_list)-1),int(max(mu_list)+1),int((max(mu_list)-min(mu_list))/6)+1)
    # bins = np.linspace(50, 70,41)
    ax.hist(mu_list,bins=bins,histtype="step",color="orange",label=f"Chiba",linewidth=2.5,alpha=1)
    # ax.axvline(x=R5912_100_tt, ymin=0, ymax=1,ls="--",color="gray",label=f"{R5912_100_tt:.0f} ns",linewidth=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"PMT Time - Tabletop Time $\mu$ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=14)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_means{temp_limits[0]}_{temp_limits[1]}C.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_means{temp_limits[0]}_{temp_limits[1]}C.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_sigma(sigma_list,temp_limits):
    print(f"min {min(sigma_list)}max{max(sigma_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(int(min(sigma_list)-1),int(max(sigma_list)+1),int((max(sigma_list)-min(sigma_list)))+1)
    bins = np.linspace(0,5,51)
    ax.hist(sigma_list,bins=bins,histtype="step",color="orange",label=f"Chiba",linewidth=2.5,alpha=1)
    # ax.axvline(x=R5912_100_tts, ymin=0, ymax=1,ls="--",color="gray",label=f"{R5912_100_tts:.1f} ns",linewidth=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"PMT Time - Tabletop Time $\sigma $ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=14)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_sigma{temp_limits[0]}_{temp_limits[1]}C.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_sigma{temp_limits[0]}_{temp_limits[1]}C.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_chi2(chi2_list,temp_limits):
    print(f"min {min(chi2_list)}max{max(chi2_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(int(min(chi2_list)-1),int(max(chi2_list)+1),int((max(chi2_list)-min(chi2_list))/12)+1)
    ax.hist(chi2_list,bins=bins,histtype="step",color="orange",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" $\chi^{2} $", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../chi2_dist{temp_limits[0]}_{temp_limits[1]}C.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_dist{temp_limits[0]}_{temp_limits[1]}C.pdf",transparent=False,bbox_inches='tight')
    plt.close()



def plot_chi2Sigma(chi2_list,sigma_list,temp_limits):
    print(f"min {min(chi2_list)}max{max(chi2_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.plot(chi2_list,sigma_list,"o")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" $\chi^{2} $", fontsize=22)
    ax.set_ylabel(r"$\sigma$ [ns]", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../chi2_sigma{temp_limits[0]}_{temp_limits[1]}C.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_sigma{temp_limits[0]}_{temp_limits[1]}C.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_chi2mu(chi2_list,mu_list,temp_limits):
    print(f"min {min(chi2_list)}max{max(chi2_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.plot(chi2_list,mu_list,"o")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" $\chi^{2} $", fontsize=22)
    ax.set_ylabel(r"$\mu$ [ns]", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../chi2_mu{temp_limits[0]}_{temp_limits[1]}C.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_mu{temp_limits[0]}_{temp_limits[1]}C.pdf",transparent=False,bbox_inches='tight')
    plt.close()



def transit_params(degg_list,temp_limits):
    file_list = []
    for idegg in degg_list:
        ifiles = sorted(glob.glob(idegg+"/DEgg*"))
        file_list += ifiles
    a_list = []
    mu_list = []
    sigma_list = []
    chi2_list = []
    pvalue_list = []
    temp_list = []
    doms_with_temp = []
    for ifile in file_list:
        a,mu,sigma,chi2,pvalue = extract_fit_params(ifile)
        if sigma > 40:
            print(f"sigma {sigma} in {ifile}")


        itemp = extract_temperature(ifile)
        idevice,isubdevice = extract_device(ifile)
        if itemp < -10:
            doms_with_temp.append(isubdevice)
        

        # if -10 <itemp < 10:
        #     print(ifile,itemp)
        # print(chi2,pvalue)
        # chi2,pvalue  = get_chi2(ifile)
        if mu > -500 and sigma < 45 and temp_limits[0]<=itemp<temp_limits[1]:
            mu_list.append(mu)
            sigma_list.append(sigma)
            chi2_list.append(chi2)
            a_list.append(a)
            temp_list.append(itemp)
            # if chi2 < 200:
            #     print(f"chi2 {chi2} in {ifile}")
        # else:
        #     print(f"please check {ifile}")
    # print(mu_list)
    # doms_with_temp = [idom.split('_')[-1] for idom in doms_with_temp]
    # print(f"dom number {len(list(set(doms_with_temp)))}")
    # print(list(set(doms_with_temp)))
    plot_mu(mu_list,temp_limits)
    plot_sigma(sigma_list,temp_limits)
    plot_chi2(chi2_list,temp_limits)
    plot_chi2mu(chi2_list,mu_list,temp_limits)
    plot_chi2Sigma(chi2_list,sigma_list,temp_limits)
    plot_temp(temp_list,temp_limits)
    return a_list,mu_list,sigma_list


def count_DOM_PMTs(degg_list):
    file_list = []
    for idegg in degg_list[:]:
        ifiles = sorted(glob.glob(idegg+"/DEgg*"))
        ifiles = [jfile for jfile in ifiles if extract_temperature(jfile)<-10]
        ifiles = list(set([jfile.split("_")[-3] for jfile in ifiles]))
        if len(ifiles)<2:
            print(f"{idegg} has {len(ifiles)} DOM measurements")
            print(ifiles)

# count_DOM_PMTs(degg_list)

transit_params(degg_list,[-30,-10])
transit_params(degg_list,[-30,50])
    
    
# make_transit_plots(degg_list)