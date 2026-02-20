#!/usr/bin/env python

import json

import glob

from matplotlib import ticker

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
    fig = plt.figure(figsize=(16,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    bins = np.linspace(8.8,10.1,40)
    ax.plot(time_list,tilt_list,'o',markersize=2,alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(0.5,0.7)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.set_ylabel(f" tilt angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f"time [UTC]", fontsize=20)
    plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id):
    fig = plt.figure(figsize=(16,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    bins = np.linspace(8.8,10.1,40)
    ax.plot(time_list,orientation_list,'o',markersize=2,alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(-180,180)
    # ax.set_ylim(0,360)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.set_ylabel(f" orientation angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f"time [UTC]", fontsize=20)
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
    plot_tilt_moni(tilt_list,time_list,board_name,hostname,port,icm_id)
    plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id)

    
            # if "2026-02-06" in time:
            #     print(f"event {count} time {time} tilt {np.rad2deg(g_tilt)} degrees rotation {np.rad2deg(g_phi)} degrees")
            #     count += 1
                # print(f" time {time} tilt {np.rad2deg(g_tilt)} degrees rotation {np.rad2deg(g_phi)} degrees")


        
