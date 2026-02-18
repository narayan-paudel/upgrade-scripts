#!/usr/bin/env python

import subprocess
import os
import glob

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from iminuit import Minuit
from iminuit.cost import LeastSquares

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

import json
import argparse

from json_tools import (extract_HV, extract_temperature, extract_device,
                         extract_tt, extract_meas_values,extract_run)

from pathlib import Path
home = str(Path.home())

'''
Reads the json file with the list of DOMs and downloads
transit measurements from fatcat database

usage: python mdom_transit_times.py -i
 /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_mdom.json

'''

print(home)
from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsCustom4 = colorsCustom2 + colorsCustom2
colorsCustom16 = colorsCustom4 + colorsCustom4 + colorsCustom4 + colorsCustom4
colorsCustom64 = colorsCustom16+colorsCustom16+colorsCustom16+colorsCustom16
colorsIter = iter(colorsCustom)
# colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
# colors = ['#1e4c7b','#498c9c','#89b3a2','#f1d27e','#d56b48','#1e1e3e',"#5a5a8c","#9c7bbc","#d45e7d",'#f7b1a1']
# colornames = ["Blue, Bright, Coastal    Green    Modern    Orange    Spring    Teal    Vibrant    Yellow,    70s",
#               "Aesthetic    Bright    Deep    Indigo    Purple    Red    Soft    Summer    Vintage Starman"]







cmdParser = argparse.ArgumentParser()
cmdParser.add_argument('-i', '--input',type=str, dest='input', 
                       default="/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_mdom.json"
                       ,help='input json file to read')
args = cmdParser.parse_args()

print("input file",args.input)
plotFolder = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"

with open(args.input,"r") as fh:
    mdom_list = json.load(fh)

data_path = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
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

print(f"mdom msu {len(mdom_list_msu)}")
print(f"mdom desy {len(mdom_list_desy)}")
print(f"mdom dvt {len(mdom_list_desy_dvt)}")


# folder_list = sorted(glob.glob(home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/*"))
folder_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/") if f.is_dir()]
folder_list = [ifolder.split("/")[-1] for ifolder in folder_list]
folder_list_msu = [ifolder for ifolder in folder_list if "_M" in ifolder]
folder_list_desy_dvt = [ifolder for ifolder in folder_list if "_D" in ifolder and "_DVT" in ifolder]
folder_list_desy = [ifolder for ifolder in folder_list if "_D" in ifolder and "_DVT" not in ifolder]

print(f"MSU folders {len(folder_list_msu)}")
print(f"desy folders {len(folder_list_desy)}")
# print(folder_list_desy)
missing_mdom = [imdom for imdom in mdom_list if imdom not in folder_list]
def check_tt(folder_list):
    mDOM_tt = 0
    mdom_tt_nonzero = 0
    for imDOM in folder_list_desy[:]:
        subdevice_list = []
        subdevice_dict = {}
        meas_list = glob.glob(data_path+f"{imDOM}/mDOM*")
        for ifile in meas_list:
            device, subdevice = extract_device(ifile)
            subdevice_list.append(subdevice)
        subdevice_list = list(set(subdevice_list))
        if len(subdevice_list) > 0:
            mDOM_tt += 1
        # print(f"subdevices {subdevice_list}")
        # tt_files = [extract_device(ifile)[1] for ifile in meas_list if extract_meas_values(ifile)[2] > 0
        #             and extract_meas_values(ifile)[0] is not np.nan]
        tt_files = [extract_device(ifile)[1] for ifile in meas_list if extract_meas_values(ifile)[2] > 0]
        tt_files = list(set(tt_files))
        if len(tt_files) < 24:
            print(f"has missing PMT measurement {len(tt_files)}")
            print(imDOM)
        else:
            mdom_tt_nonzero += 1
        print(len(tt_files)) 
        pmt_tt = []
        for isub in subdevice_list:
            isub_meas = [ifile for ifile in meas_list if extract_device(ifile)[1] == isub]
            # isub_meas = [ifile for ifile in isub_meas if extract_meas_values(ifile)[2] > 0 
            #              and extract_meas_values(ifile)[0] is not np.nan]  #select non zero tt and HV not nan 
            isub_meas = [ifile for ifile in isub_meas if extract_meas_values(ifile)[2] > 0 
                        ]  #select non zero tt and HV not nan      
            run_list = [extract_run(ifile) for ifile in isub_meas]
            # if len(run_list) == 0:
            #     print(f"{len(run_list)} runs with non zero tt")
            #     print(f"This is {imDOM} {isub}")
            # else:
            #     print(f"run list {run_list} max {max(run_list)}")
            if len(isub_meas) > 1:
                run_max = max(run_list)
                select_meas = [ifile for ifile in isub_meas if extract_run(ifile) == run_max][0]
                hv,mb_ch,mu,sigma,a,b,c,chi2 = extract_meas_values(select_meas)
                print(f"{isub} tt {mu} ns {extract_temperature(select_meas)}C {hv}V {extract_run(select_meas)}")
            elif len(isub_meas) == 1:
                hv,mb_ch,mu,sigma,a,b,c,chi2 = extract_meas_values(isub_meas[0])
                print(f"{isub} tt {mu} ns {extract_temperature(isub_meas[0])}C {hv}V {extract_run(isub_meas[0])}")

                # print(f"This is {imDOM}")
                # print(f"subdevice {isub} has {len(isub_meas)} meas\n")
                # for ifile in isub_meas:
                #     hv,mb_ch,mu,sigma,a,b,c,chi2 = extract_meas_values(ifile)
                    # print(f"{isub} tt {mu} ns {extract_temperature(ifile)}C {hv}V {extract_run(ifile)}")

        # for ifile in meas_list:
        #     temp = extract_temperature(ifile)
        #     hv = extract_HV(ifile)
        #     device, subdevice = extract_device(ifile)
        #     subdevice_dict[subdevice][]
        #     print(f"{subdevice} {temp} C {hv} V")
        #     subdevice_list.append(subdevice)
        # subdevice_list = list(set(subdevice_list))
        # print(len(subdevice_list))

        # print(meas_list)
        # print(f"mDOM {imDOM} has {len(meas_list)} meas")
    print(f"mDOM with at least one meas {mDOM_tt}")
    print(f"mDOM with at least one meas {mdom_tt_nonzero}")
    # print(f"There are {len(missing_mdom)} missing dom measurements")
    # print(missing_mdom)


def check_HV(folder_list):
    '''
    check to see if there are multiple measurement of tt at different HV
    of mDOMs in folder_list
    '''
    print("checking if multiple HV measurements per DOM")
    HV_diff_max = []
    # fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=1,ncols=1)
    # ax = fig.add_subplot(gs[0])
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    color_count = 0
    for imDOM in folder_list[:]:
        subdevice_list = []
        subdevice_dict = {}
        meas_list = glob.glob(data_path+f"{imDOM}/mDOM*")
        for ifile in meas_list:
            device, subdevice = extract_device(ifile)
            subdevice_list.append(subdevice)
        subdevice_list = list(set(subdevice_list))
        for isub in subdevice_list:
            tt_files = [ifile for ifile in meas_list if extract_device(ifile)[1] == isub and extract_meas_values(ifile)[2] > 0]
            HV_meas = [extract_meas_values(ifile)[0] for ifile in tt_files]
            tt_meas = [extract_meas_values(ifile)[2] for ifile in tt_files]
            tt_meas = [itt for itt,ihv in zip(tt_meas,HV_meas) if not np.isnan(ihv)]
            HV_meas = [ihv for ihv in HV_meas if not np.isnan(ihv)]

            if len(HV_meas)>1:
                HV_meas_diff = np.diff(HV_meas)
                HV_diff_max.append(max(abs(HV_meas_diff)))
                # print(HV_meas_diff)
                if max(abs(HV_meas_diff)) > 5 and len(HV_meas_diff)>3:
                    HV_meas_inv = [1.0/np.sqrt(iV) for iV in HV_meas]
                    ax.plot(HV_meas_inv,tt_meas,"o",color=colorsCustom[color_count],label=f"{imDOM}_{isub}")
                    z = np.polyfit(HV_meas_inv,tt_meas, 1)
                    p = np.poly1d(z)
                    x_fit = np.linspace(min(HV_meas_inv),max(HV_meas_inv),20)
                    y_fit = p(x_fit)
                    ax.plot(x_fit,y_fit,"-",color=colorsCustom[color_count])
                    print(HV_meas_diff)
                    print(f"{isub} in {imDOM} potentially has multiple HV meas, {max(abs(HV_meas_diff))}")
                    color_count += 1
                    print(f"color count {color_count}")
    print(f"maximum HV diff {max(HV_diff_max)}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"high voltage [1/sqrt(V)]", fontsize=22)
    ax.set_ylabel("transit time [ns]", fontsize=22)
    ax.set_ylim(0,50)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=12,ncols=2)
    plt.savefig(plotFolder+f"/../transit_time_vs_HV.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/../transit_time_vs_HV.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def line(x, α, β):
    return α + x * β

def check_HV_separate(folder_list):
    '''
    check to see if there are multiple measurement of tt at different HV
    of mDOMs in folder_list
    '''
    print("checking if multiple HV measurements per DOM")
    HV_diff_max = []

    color_count = 0
    for imDOM in folder_list[:]:
        subdevice_list = []
        subdevice_dict = {}
        meas_list = glob.glob(data_path+f"{imDOM}/mDOM*")
        for ifile in meas_list:
            device, subdevice = extract_device(ifile)
            subdevice_list.append(subdevice)
        subdevice_list = list(set(subdevice_list))
        for isub in subdevice_list:
            tt_files = [ifile for ifile in meas_list if extract_device(ifile)[1] == isub and extract_meas_values(ifile)[2] > 0]
            HV_meas = [extract_meas_values(ifile)[0] for ifile in tt_files]
            tt_meas = [extract_meas_values(ifile)[2] for ifile in tt_files]
            tt_sigma_meas = [extract_meas_values(ifile)[3] for ifile in tt_files]
            tt_meas = [itt for itt,ihv in zip(tt_meas,HV_meas) if not np.isnan(ihv)]
            tt_sigma_meas = [itt for itt,ihv in zip(tt_sigma_meas,HV_meas) if not np.isnan(ihv)]
            HV_meas = [ihv for ihv in HV_meas if not np.isnan(ihv)]

            if len(HV_meas)>1:
                HV_meas_diff = np.diff(HV_meas)
                HV_diff_max.append(max(abs(HV_meas_diff)))
                # print(HV_meas_diff)
                if max(abs(HV_meas_diff)) > 5 and len(HV_meas_diff)>3:
                    HV_meas_inv = [1.0/np.sqrt(iV) for iV in HV_meas]
                    fig = plt.figure(figsize=(8,5))
                    gs = gridspec.GridSpec(nrows=1,ncols=1)
                    ax = fig.add_subplot(gs[0])
                    ax.errorbar(HV_meas_inv,tt_meas,yerr=tt_sigma_meas,fmt="o",color=colorsCustom[color_count],label=f"transit time",alpha=1)
                    least_squares = LeastSquares(HV_meas_inv, tt_meas, tt_sigma_meas, line)
                    m = Minuit(least_squares, α=-40, β=1000)  # starting values for α and β
                    m.migrad()  # finds minimum of least_squares function
                    m.hesse()  # accurately computes uncertainties
                    ax.plot(HV_meas_inv, line(np.asarray(HV_meas_inv), *m.values), label="linear fit")
                    ax.text(0.78,0.1,s=f"{imDOM} {isub}",size=13,horizontalalignment='center',verticalalignment='center', transform=ax.transAxes)
 
                    # display legend with some fit info
                    fit_info = [
                        f"$\\chi^2$/$n_\\mathrm{{dof}}$ = {m.fval:.1f} / {m.ndof:.0f} = {m.fmin.reduced_chi2:.1f}",
                    ]
                    for p, v, e in zip(m.parameters, m.values, m.errors):
                        fit_info.append(f"{p} = ${v:.0f} \\pm {e:.0f}$")

                    # ax.plot(x_fit,y_fit,"-",color=colorsCustom[color_count])
                    print(HV_meas_diff)
                    print(f"{isub} in {imDOM} potentially has multiple HV meas, {max(abs(HV_meas_diff))}")
                    color_count += 1
                    print(f"color count {color_count}")                    
                    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
                    ax.set_xlabel(r"1/$\sqrt{\mathrm{high\, voltage\, [V]}}$", fontsize=22)
                    ax.set_ylabel("transit time [ns]", fontsize=22)
                    ax.set_ylim(0,50)
                    ax.grid(True,alpha=0.6)
                    ax.legend(title="\n".join(fit_info),fontsize=12)
                    plt.savefig(plotFolder+f"/../transit_time_vs_HV{imDOM}_{isub}.png",transparent=False,bbox_inches='tight')
                    plt.savefig(plotFolder+f"/../transit_time_vs_HV{imDOM}_{isub}.pdf",transparent=False,bbox_inches='tight')
                    plt.close()
    print(f"maximum HV diff {max(HV_diff_max)}")

check_HV(folder_list_desy)
check_HV(folder_list_msu)
# check_HV_separate(folder_list_msu)


# python mdom_meas_detail.py -i /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_mdom.json