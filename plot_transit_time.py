#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import json

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = "/home/narayan/research_ua/icecube/Upgrade/timing_calibration/plots/"

from customColors import qualitative_colors

colorsCustom = qualitative_colors(12)
colorsCustom2 = colorsCustom + colorsCustom
colorsIter = iter(colorsCustom)
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']

input_file = "/home/narayan/research_ua/icecube/Upgrade/timing_calibration/data/degg_transit/DEgg2020-1-001_v2/DEgg2020-1-001_v2_SQ0411_run138_n1.json"

with open(input_file, 'r') as f:
    data = json.load(f)

print(data["meas_data"][0]['y_values'])

y_values = data["meas_data"][0]['y_values']
x_min = data["meas_data"][0]["x_min"]
x_max = data["meas_data"][0]["x_max"]
x_values = np.linspace(x_min,x_max,99)
x_label = data["meas_data"][0]["x_label"]
device_uid = data["device_uid"]
# print(data["meas_data"][])

fit_x_values = np.linspace(x_min,x_max,99)
fit_y_values = data["meas_data"][0]["fit_y_values"]


def degg_transit_plot(x_values,y_values,fit_x_values,fit_y_values,x_label,device_uid):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.step(x_values,y_values,ls='-',lw = 2.5,c="blue",label=str(device_uid),alpha=1)
    ax.plot(fit_x_values,fit_y_values,ls='-',lw = 2.5,c="orange",label="fit",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"{}".format(x_label), fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    # ax.set_yscale("log")
    ax.legend(fontsize=16)
    plt.savefig(plotFolder+"/transit_time.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+"/transit_time.pdf",transparent=False,bbox_inches='tight')
    plt.close()
degg_transit_plot(x_values,y_values,fit_x_values,fit_y_values,x_label,device_uid)


