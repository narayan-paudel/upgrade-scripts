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

plotFolder = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
mdom_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/") if f.is_dir()][:]
print(mdom_list)


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
    run = data['run_number']
    temperature = data["meas_data"][0]["temperature"]
    # print(data["meas_data"][])

    # fit_x_values = np.linspace(x_min,x_max,n_bins)
    # fit_y_values = data["meas_data"][0]["fit_y_values"]
    return x_values,y_values,x_label,device_uid,pmt,run,temperature,meas_time,

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




def mdom_transit_plot(file_list):
    n_files = len(file_list)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for i,ifile in enumerate(file_list):
        # print(f"no of file {i}")
        x_values,y_values,x_label,device_uid,pmt,run,temperature,meas_time = extract_histogram(ifile)
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
    plt.savefig(plotFolder+f"/mDOM_transit_time{device_uid}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/mDOM_transit_time{device_uid}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def make_transit_plots(mdom_list):
    for imdom in mdom_list[:]:
        ifiles = sorted(glob.glob(imdom+"/m*"))
        if len(ifiles)>0:
            mdom_transit_plot(ifiles)
        else:
            print(f"the files list is empty for {imdom}")

def extract_fit_params(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    y_values = data["meas_data"][0]['y_values']
    fit_y_values = data["meas_data"][0]["fit_y_values"]
    # chi2,pvalue = stats.chisquare(y_values,fit_y_values,ddof=3,sum_check=True)
    chi2,pvalue = stats.chisquare(y_values,fit_y_values,ddof=3,sum_check=False)
    fit_params = data["meas_data"][0]["fit_params"]
    return *fit_params,chi2,pvalue

def plot_mu(mu_list):
    print(f"min {min(mu_list)}max{max(mu_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(int(min(mu_list)-1),int(max(mu_list)+1),int((max(mu_list)-min(mu_list))/6)+1)
    ax.hist(mu_list,bins=bins,histtype="step",color="orange",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time peaks $\mu$ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_means.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_means.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_sigma(sigma_list):
    print(f"min {min(sigma_list)}max{max(sigma_list)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(int(min(sigma_list)-1),int(max(sigma_list)+1),int((max(sigma_list)-min(sigma_list)))+1)
    ax.hist(sigma_list,bins=bins,histtype="step",color="orange",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time width $\sigma $ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_sigma.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_sigma.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_chi2(chi2_list):
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
    plt.savefig(plotFolder+f"/../chi2_dist.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_dist.pdf",transparent=False,bbox_inches='tight')
    plt.close()



def plot_chi2Sigma(chi2_list,sigma_list):
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
    plt.savefig(plotFolder+f"/../chi2_sigma.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_sigma.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_chi2mu(chi2_list,mu_list):
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
    plt.savefig(plotFolder+f"/../chi2_mu.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../chi2_mu.pdf",transparent=False,bbox_inches='tight')
    plt.close()



def transit_params(mdom_list):
    file_list = []
    for imdom in mdom_list:
        ifiles = sorted(glob.glob(imdom+"/*"))
        file_list += ifiles
    a_list = []
    mu_list = []
    sigma_list = []
    chi2_list = []
    pvalue_list = []
    for ifile in file_list:
        a,mu,sigma,chi2,pvalue = extract_fit_params(ifile)
        # print(chi2,pvalue)
        # chi2,pvalue  = get_chi2(ifile)
        if mu > -500 and sigma < 500:
            mu_list.append(mu)
            sigma_list.append(sigma)
            chi2_list.append(chi2)
            a_list.append(a)
        else:
            print(f"please check {ifile}")
    # print(mu_list)
    plot_mu(mu_list)
    plot_sigma(sigma_list)
    plot_chi2(chi2_list)
    plot_chi2mu(chi2_list,mu_list)
    plot_chi2Sigma(chi2_list,sigma_list)
    return a_list,mu_list,sigma_list





# transit_params(mdom_list)
    
    
make_transit_plots(mdom_list)