#!/usr/bin/env python

import json

import glob

from matplotlib import ticker,text

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


plotFolder = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/plots/"



moni_data_folder = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/moni-data/"

moni_files = sorted(glob.glob(moni_data_folder+"/*.json"))


moni_file = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/moni-data/fieldhub87_mDOM_12000030c2c6e62d.json"

def tilt_angle(nx,ny,nz):
    '''
    calculate tilt angle
    '''
    return np.arctan2(np.sqrt(nx*nx + ny*ny),nz)

def to_spherical(x,y,z):
    '''
    converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
    '''
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arctan2(np.sqrt(x*x + y*y),z)
    phi = np.arctan2(y,x)
    #####to ensure phi is positive
    if phi < 0.0:
        phi_pos = phi + 2*np.pi
    else:
        phi_pos = phi
    return r,theta,phi_pos

def plot_tilt_moni(tilt_list,time_list,board_name,hostname,port,icm_id):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    bins = np.linspace(8.8,10.1,40)
    time_list = pd.to_datetime(time_list, format='ISO8601')
    ax.plot(time_list,tilt_list,'o',markersize=2,alpha=0.5)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(0.5,0.7)
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    # print(ax.get_xticklabels(minor=False))
    # labels = [label.get_text()[0:10] for label in ax.get_xticklabels(minor=False)]
    # location = [label.get_position() for label in ax.get_xticklabels(minor=False)]
    # label_loc = [text.Text(iloc,0,ilabel) for iloc, ilabel in enumerate(zip(location,labels))]
    # ax.set_xticklabels(label_loc,rotation=45, ha="right")
    ###################
    # ax.set_xticklabels(labels,rotation=45, ha="right")
    ################
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    # print(f"xtick labels {labels}, {location}")
    ax.set_ylabel(f" tilt angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f"time [UTC]", fontsize=20)
    ax.tick_params("x", rotation=45) 
    plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    bins = np.linspace(8.8,10.1,40)
    time_list = pd.to_datetime(time_list, format='ISO8601')
    ax.plot(time_list,orientation_list,'o',markersize=2,alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(-180,180)
    # ax.set_ylim(0,360)
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.set_ylabel(f" orientation angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f"time [UTC]", fontsize=20)
    ax.tick_params("x", rotation=45) 
    plt.savefig(plotFolder+f"/orientation_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_vs_time_{board_name}_{hostname}_{port}_{icm_id}.pdf",transparent=False,bbox_inches='tight')
    plt.close()




for ifile in moni_files[:]:
    with open(ifile, "r") as fh:
        moni_data = json.load(fh)
        # print(f"moni data {moni_data}")
        count = 0
        tilt_list = []
        orientation_list = []
        time_list = []
        for ievent in moni_data[:]:
            # print(f"{ievent['moni_common']['accelerometer']}")
            gx,gy,gz = ievent['moni_common']['accelerometer']
            bx,by,bz = ievent['moni_common']['magnetometer']
            g,g_theta,g_phi = to_spherical(gx,gy,gz)
            B,B_theta,B_phi = to_spherical(bx,by,bz)
            g_tilt = tilt_angle(gx,gy,gz)
            # print(f"tilt {np.rad2deg(g_tilt)} degrees")
            # print(f"azimuth {np.rad2deg(g_phi)} degrees")
            time = ievent["header"]['utc_time']
            board_name = ievent['header']['board_name']
            hostname = ievent['header']['hostname']
            port = ievent['header']['port']
            icm_id = ievent['header']['icm_id']
            # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
            tilt_list.append(np.rad2deg(g_tilt))
            orientation_list.append(np.rad2deg(B_phi))
            time_list.append(time)
    print(f"plotting file {ifile}")
    # plot_tilt_moni(tilt_list,time_list,board_name,hostname,port,icm_id)
    # plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id)



def tilt_plot(string,device):
    tilt_list = []
    time_list = []
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id = ievent['header']['icm_id']
                    # print(f"{ievent['moni_common']['accelerometer']}")
                    gx,gy,gz = ievent['moni_common']['accelerometer']
                    bx,by,bz = ievent['moni_common']['magnetometer']
                    g,g_theta,g_phi = to_spherical(gx,gy,gz)
                    B,B_theta,B_phi = to_spherical(bx,by,bz)
                    g_tilt = tilt_angle(gx,gy,gz)
                    # print(f"tilt {np.rad2deg(g_tilt)} degrees")
                    # print(f"azimuth {np.rad2deg(g_phi)} degrees")
                    # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
                    tilt_list.append(np.rad2deg(g_tilt))
                    # orientation_list.append(np.rad2deg(B_phi))
                    time_list.append(time)
                print(f"plotting file {string} {device}")
    plot_tilt_moni(tilt_list,time_list,device,f"fieldhub{string}",f"all",f"all")
    # plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id)

# tilt_plot("87","mDOM")
# tilt_plot("88","mDOM")
# tilt_plot("89","mDOM")
# tilt_plot("90","mDOM")
# tilt_plot("91","mDOM")
# tilt_plot("92","mDOM")




def orientation_plot(string,device):
    orientation_list = []
    time_list = []
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id = ievent['header']['icm_id']
                    # print(f"{ievent['moni_common']['accelerometer']}")
                    gx,gy,gz = ievent['moni_common']['accelerometer']
                    bx,by,bz = ievent['moni_common']['magnetometer']
                    g,g_theta,g_phi = to_spherical(gx,gy,gz)
                    B,B_theta,B_phi = to_spherical(bx,by,bz)
                    g_tilt = tilt_angle(gx,gy,gz)
                    # print(f"tilt {np.rad2deg(g_tilt)} degrees")
                    # print(f"azimuth {np.rad2deg(g_phi)} degrees")
                    # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
                    orientation_list.append(np.rad2deg(B_phi))
                    time_list.append(time)
                print(f"plotting file {string} {device}")
    plot_orientation_moni(orientation_list,time_list,device,f"fieldhub{string}",f"all",f"all")
    # plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id)

# orientation_plot("87","mDOM")
# orientation_plot("88","mDOM")
# orientation_plot("89","mDOM")
# orientation_plot("90","mDOM")
# orientation_plot("91","mDOM")
# orientation_plot("92","mDOM")



def get_tilt(string,port,device,measurement_time):
    g_tilt_list = []
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id = ievent['header']['icm_id']
                    if time > measurement_time and port == port:
                        gx,gy,gz = ievent['moni_common']['accelerometer']
                        g_tilt = tilt_angle(gx,gy,gz)
                        g_tilt_list.append(np.rad2deg(g_tilt))
    return np.mean(g_tilt_list), np.std(g_tilt_list)

def get_calibrated_bxbybz(bx,by,bz,parameter):










def get_orientation(string,port,device,measurement_time):
    orientation_list = []
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id = ievent['header']['icm_id']
                    if time > measurement_time and port == port:
                        print(f"time {time} {measurement_time} {time > measurement_time}")
                        bx,by,bz = ievent['moni_common']['magnetometer']
                        B,B_theta,B_phi = to_spherical(bx,by,bz)
                        bx_cal, by_cal, bz_cal = get_calibrated_bxbybz(bx,by,bz,parameter=[-0.00273318, -0.00030861, -0.00270213, -0.04387554, 0.1290214])
                        orientation_list.append(np.rad2deg(B_phi))
    return np.mean(orientation_list), np.std(orientation_list)


def get_angle_between_devices(string1,string2,string3,port,device,measurement_time):
    tilt1, tilt1_std = get_tilt(string1,port,device,measurement_time)
    tilt2, tilt2_std = get_tilt(string2,port,device,measurement_time)
    tilt3, tilt3_std = get_tilt(string3,port,device,measurement_time)
    orientation1, orientation1_std = get_orientation(string1,port,device,measurement_time)
    orientation2, orientation2_std = get_orientation(string2,port,device,measurement_time)
    orientation3, orientation3_std = get_orientation(string3,port,device,measurement_time)
    print(f"tilt 1 {string1} {tilt1} +/- {tilt1_std} degrees")
    print(f"tilt 2 {string2} {tilt2} +/- {tilt2_std} degrees")
    print(f"tilt 3 {string3} {tilt3} +/- {tilt3_std} degrees")
    print(f"orientation 1 {string1} {orientation1} +/- {orientation1_std} degrees")
    print(f"orientation 2 {string2} {orientation2} +/- {orientation2_std} degrees")
    print(f"orientation 3 {string3} {orientation3} +/- {orientation3_std} degrees")

get_angle_between_devices("87","88","89","5170","mDOM","2026-02-21T00:00:00Z")


def get_calibrated_angles_between_devices(string1,string2,string3,port,device,measurement_time):
    tilt1, tilt1_std = get_tilt(string1,port,device,measurement_time)
    tilt2, tilt2_std = get_tilt(string2,port,device,measurement_time)
    tilt3, tilt3_std = get_tilt(string3,port,device,measurement_time)
    orientation1, orientation1_std = get_orientation(string1,port,device,measurement_time)
    orientation2, orientation2_std = get_orientation(string2,port,device,measurement_time)
    orientation3, orientation3_std = get_orientation(string3,port,device,measurement_time)
    print(f"calibrated tilt 1 {string1} {tilt1} +/- {tilt1_std} degrees")
    print(f"calibrated tilt 2 {string2} {tilt2} +/- {tilt2_std} degrees")
    print(f"calibrated tilt 3 {string3} {tilt3} +/- {tilt3_std} degrees")
    print(f"calibrated orientation 1 {string1} {orientation1} +/- {orientation1_std} degrees")
    print(f"calibrated orientation 2 {string2} {orientation2} +/- {orientation2_std} degrees")
    print(f"calibrated orientation 3 {string3} {orientation3} +/- {orientation3_std} degrees")
        
