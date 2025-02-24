#!/usr/bin/env python

import os
import glob

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

# mdom_list_desy = [ifile for ifile in mdom_list if ifile.split("/")[-1].split("_")[1][0]=="D" and not "DVT" in ifile.split("/")[-1].split("_")[1]]
# mdom_list_desy_dvt = [ifile for ifile in mdom_list if "DVT" in ifile.split("/")[-1].split("_")[1]]
# # mdom_list_desy_dvt = [ifile for ifile in mdom_list if ifile.split("/")[-1].split("_")[1][0]=="D" and ifile.split("/")[-1].split("_")[1][1]=="V"]
# # test = [ifile.split("/")[-1].split("_")[1] for ifile in mdom_list_desy_dvt]
# # print(test)
# print(mdom_list_desy[:5])
# print(mdom_list_desy_dvt[0] in mdom_list_desy)
# print(mdom_list_msu)
# print(mdom_list_desy_dvt)
mdom_names = [istr.split("/")[-1] for istr in mdom_list]
# print(f"found {len(mdom_names)} mdoms in fatcat")
# print(mdom_names)


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
    run = data['run_number']
    temperature = data["meas_data"][0]["temperature"]
    # print(data["meas_data"][])

    # fit_x_values = np.linspace(x_min,x_max,n_bins)
    # fit_y_values = data["meas_data"][0]["fit_y_values"]
    return x_values,y_values,x_label,device_uid,pmt,run,temperature,meas_time,

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
        x_values,y_values,x_label,device_uid,pmt,run,temperature,meas_time = extract_histogram(ifile)
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

def extract_fit_params(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    y_values = data["meas_data"][0]['y_values']
    # print(data["meas_data"])
    # print(len(data["meas_data"]))
    param_dict = {"Transit time":np.nan,"Transit time spread":np.nan,"a":np.nan,"b":np.nan,"c":np.nan,"chi2":np.nan}
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
        ,param_dict["b"],param_dict["c"],param_dict["chi2"]

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
    ax.hist(temp_list,bins=bins,histtype="step",color="orange",lw=4.5,alpha=1)
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
    ax.hist(hv_list,bins=bins,histtype="step",color="orange",lw=4.5,alpha=1)
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

            mu,sigma,a,b,c,chi2 = extract_fit_params(ifile)
            itemp = extract_temperature(ifile)
            ihv = extract_HV(ifile)
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

def extract_params(mdom_list,meas_site):
    '''
    gets mean and variance from json file in the given list of files.
    '''
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
    for ifile in file_list:
        # print(f"running {ifile}")
        if check_meas_site(ifile) == meas_site:
            files_at_site.append(ifile)
        else:
            print(f"This file has mesurement from {check_meas_site(ifile)} instead of {meas_site}")
    for ifile in files_at_site:
            mu,sigma,a,b,c,chi2 = extract_fit_params(ifile)
            itemp = extract_temperature(ifile)
            if mu > -500 and sigma < 500 and not np.isnan(mu) and not np.isnan(sigma)\
                and not np.isnan(a) and not np.isnan(b) and not np.isnan(c) and not np.isnan(chi2):
                mu_list.append(mu)
                sigma_list.append(sigma)
                chi2_list.append(chi2)
                a_list.append(a)
                b_list.append(b)
                c_list.append(c)
                temp_list.append(itemp)
                # if c!=sigma:
                #     print(f"{c}c and {sigma}")
            elif abs(mu-0.0)< 0.001 or abs(sigma-0.0)<0.001:
                zero_mu_files.append(ifile)
            else:
                print(f"please check {ifile}. It has {mu} transit time and {sigma} spread")
    return mu_list,sigma_list




def combined_mean_plot(mdom_list_desy,mdom_list_msu):
    tt_desy,tts_desy = extract_params(mdom_list_desy,"desy")
    tt_msu,tts_msu = extract_params(mdom_list_msu,"msu")
    print(f"desy tt min {min(tt_desy)} max{max(tt_desy)}")
    print(f"desy tts min {min(tts_desy)} max{max(tts_desy)}")

    print(f"msu tt min {min(tt_msu)} max{max(tt_msu)}")
    print(f"msu tts min {min(tts_msu)} max{max(tts_msu)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # bins = np.linspace(int(min(mu_list)-1),int(max(mu_list)+1),int((max(mu_list)-min(mu_list))/6)+1)
    # bins = np.linspace(0,100,401)
    c1="#1C5E7D"#"#2b8cbe"
    c2="#4DC9CB"#"#c34a36"
    bins = np.linspace(-400,1100,3001)#full
    bins = np.linspace(35,55,41)#zoom

    ax.hist(tt_desy,bins=bins,histtype="step",color=c1,label=r"DESY",linewidth=2.5,alpha=1)
    ax.hist(tt_msu,bins=bins,histtype="step",color=c2,label="MSU",linewidth=2.5,alpha=1)
    ax.axvline(x=R12199_tt, ymin=0, ymax=1,ls="--",color="gray",label=f"{R12199_tt:.0f} ns",linewidth=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time peaks $\mu$ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=14)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_meansmDOM_combined.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_meansmDOM_combined.pdf",transparent=False,bbox_inches='tight')
    plt.close()
    
def combined_spread_plot(mdom_list_desy,mdom_list_msu):
    tt_desy,tts_desy = extract_params(mdom_list_desy,"desy")
    tt_msu,tts_msu = extract_params(mdom_list_msu,"msu")
    print(f"desy tt min {min(tt_desy)} max{max(tt_desy)}")
    print(f"desy tts min {min(tts_desy)} max{max(tts_desy)}")

    print(f"msu tt min {min(tt_msu)} max{max(tt_msu)}")
    print(f"msu tts min {min(tts_msu)} max{max(tts_msu)}")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    c1="#1C5E7D"#"#2b8cbe"
    c2="#4DC9CB"#"#c34a36"
    bins = np.linspace(0,400,801)
    bins = np.linspace(0,10,101)#zoom
    ax.hist(tts_desy,bins=bins,histtype="step",color=c1,linewidth=2.5,label="DESY",alpha=1)
    ax.hist(tts_msu,bins=bins,histtype="step",color=c2,linewidth=2.5,label="MSU",alpha=1)
    ax.axvline(x=R15458_02_tts, ymin=0, ymax=1,ls="--",color="gray",label=f"{R15458_02_tts:.1f} ns",linewidth=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"transit time width $\sigma $ [ns]", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=14)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/../transit_time_sigmamDOM_combined.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_sigmamDOM_combined.pdf",transparent=False,bbox_inches='tight')
    plt.close()
    


combined_mdom_list = mdom_list_msu + mdom_list_desy
# combined_mean_plot(mdom_list_desy,mdom_list_msu)
# combined_spread_plot(mdom_list_desy,mdom_list_msu)





transit_params(mdom_list,meas_site="desy")
transit_params(mdom_list,meas_site="msu")
    
    
# make_transit_plots(mdom_list[:1])