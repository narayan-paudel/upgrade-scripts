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

plotFolder = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_rotation/"

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']

mdom = "mDOM_M016_v1"
runs = {"PMT 8 North":776,"PMT 7 North": 777,"PMT 6 North": 778, "PMT 5 North": 779, "PMT 4 North": 780, "PMT 3 North": 781, "PMT 2 North": 782, "PMT 1 North": 783}
filemap = {}
for pmt, run in runs.items():
    filelist = glob.glob(home+f"/research_ua/icecube/Upgrade/timing_calibration/data/mdom_transit/{mdom}/{mdom}_mdom-mainboard_7f0000006732f842_{run}.json")
    if len(filelist)==0:
        print(f"missing {pmt} run {run} data")
    elif len(filelist)>1:
        print(f"multiple files found for {pmt} run {run} data")
    else:
        filemap[run] = filelist[0]

# print(filemap)
def extract_magnetometer_data(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    meas_data = data['meas_data']
    for idata in meas_data:
        times = idata['moni_times']
        for imoni in idata['monitoring']:
            if imoni["moni_name"] == 'Magnetometer Z':
                mag_z = imoni["moni_data"]
            elif imoni["moni_name"] == 'Magnetometer Y':
                mag_y = imoni["moni_data"]
            elif imoni["moni_name"] == 'Magnetometer X':
                mag_x = imoni["moni_data"]
    return mag_x, mag_y, mag_z,times

x_mean_list = []
x_allrun_list = []
y_mean_list = []
y_allrun_list = []
z_mean_list = []
z_allrun_list = []

r_mean_list = []
r_allrun_list = []
theta_mean_list = []
theta_allrun_list = []
phi_mean_list = []
phi_allrun_list = []



times_allrun_list = []

def plot_magnetometer_xyz(mag_x, mag_y, mag_z, times, run):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    times_sec = [(itime - times[0]) for itime in times]
    ax.plot(times_sec, np.asarray(mag_x)*10**6, "-o", c=colorsCustom[0], label=f"${'B_x'}$", alpha=1)
    ax.plot(times_sec, np.asarray(mag_y)*10**6, "-o", c=colorsCustom[1], label=f"${'B_y'}$", alpha=1)
    ax.plot(times_sec, np.asarray(mag_z)*10**6, "-o", c=colorsCustom[2], label=f"${'B_z'}$", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=16)
    ax.set_xlabel(r" $time$ [s]", fontsize=20)
    ax.set_ylabel(r" $B_i$ [$\mu$T]", fontsize=20)
    ax.set_ylim(-80, 50)
    plt.savefig(plotFolder+f"/magnetometer_xyz_{run}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometer_xyz_{run}.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_magnetometer_B(mag_x, mag_y, mag_z, times, run):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    times_sec = [(itime - times[0]) for itime in times]
    ax.plot(times_sec, np.sqrt(np.asarray(mag_x)**2+np.asarray(mag_y)**2+np.asarray(mag_z)**2)*10**6, "-o", c=colorsCustom[0], label=f"${'B_x'}$", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=16)
    ax.set_xlabel(r" $time$ [s]", fontsize=20)
    ax.set_ylabel(r" $B_i$ [$\mu$T]", fontsize=20)
    # ax.set_ylim(-35, 30)
    plt.savefig(plotFolder+f"/magnetometer_B_{run}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometer_B_{run}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_magnetometer_r(r, times, run):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    times_sec = [(itime - times[0]) for itime in times]
    ax.plot(times_sec, np.asarray(r)*10**6, "-o", c=colorsCustom[0], label=f"${'B_r'}$", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=16)
    ax.set_xlabel(r" $time$ [s]", fontsize=20)
    ax.set_ylabel(r" $B_i$ [$\mu$T]", fontsize=20)
    ax.set_ylim(0, 90)
    plt.savefig(plotFolder+f"/magnetometer_r_{run}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometer_r_{run}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

B_NTS = 53.06 #muT
inclination_NTS = 90 + 68.87 #degrees Down
declination_NTS = 90 - (-6.58) #degrees -West (+ve would be east)


def plot_magnetometer_r(r, times, run):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    times_sec = [(itime - times[0]) for itime in times]
    ax.plot(times_sec, np.asarray(r)*10**6, "-o", c=colorsCustom[0], label=f"${'B'}$", alpha=1)
    ax.axhline(B_NTS,0,1,ls="--",lw=2.5,label=f"B$_{{geo}}$ ({B_NTS:.1f} ${{\mu}}$T)",alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=16)
    ax.set_xlabel(r" $time$ [s]", fontsize=20)
    ax.set_ylabel(r" $B$ [$\mu$T]", fontsize=20)
    ax.set_ylim(0, 90)
    plt.savefig(plotFolder+f"/magnetometer_r_{run}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometer_r_{run}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_magnetometer_theta(theta, times, run):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    times_sec = [(itime - times[0]) for itime in times]
    ax.plot(times_sec, np.rad2deg(np.asarray(theta)), "-o", c=colorsCustom[0], label=f"${'Theta'}$", alpha=1)
    ax.axhline(inclination_NTS,0,1,ls="--",lw=2.5,label=f"$\delta$ ({inclination_NTS:.1f}\u00b0)",alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=16)
    ax.set_xlabel(r" $time$ [s]", fontsize=20)
    ax.set_ylabel(r" $\theta$ [deg]", fontsize=20)
    ax.set_ylim(0, 180)
    plt.savefig(plotFolder+f"/magnetometer_theta_{run}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometer_theta_{run}.pdf",transparent=False,bbox_inches='tight')
    plt.close()
def plot_magnetometer_phi(phi, times, run):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    times_sec = [(itime - times[0]) for itime in times]
    ax.plot(times_sec, np.rad2deg(np.asarray(phi)), "-o", c=colorsCustom[0], label=f"${'Phi'}$", alpha=1)
    # ax.axhline(declination_NTS,0,1,ls="--",lw=2.5,label=f"$I$ ({declination_NTS:.1f}\u00b0)",alpha=1.0)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(fontsize=16)
    ax.set_xlabel(r" $time$ [s]", fontsize=20)
    ax.set_ylabel(r" $\phi$ [deg]", fontsize=20)
    ax.set_ylim(-180, 180)
    plt.savefig(plotFolder+f"/magnetometer_phi_{run}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometer_phi_{run}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def to_spherical(x,y,z):
    '''
    converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
    '''
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arctan2(np.sqrt(x*x + y*y),z)
    phi = np.arctan2(y,x)
    return r,theta,phi

def spherical_lists(x_list,y_list,z_list):
    r_list = []
    θ_list = []
    φ_list = []
    for x,y,z in zip(x_list,y_list,z_list):
        r,theta,phi = to_spherical(x,y,z)
        r_list.append(r)
        θ_list.append(theta)
        φ_list.append(phi)
    return r_list,θ_list,φ_list

for run, json_file in filemap.items():
    mag_x, mag_y, mag_z,times = extract_magnetometer_data(json_file)
    r,theta,phi = spherical_lists(mag_x, mag_y, mag_z)
    plot_magnetometer_xyz(mag_x, mag_y, mag_z,times, run)
    # plot_magnetometer_r(r, times, run)
    # plot_magnetometer_theta(theta, times, run)
    # plot_magnetometer_phi(phi, times, run)
    x_allrun_list.extend(mag_x)
    y_allrun_list.extend(mag_y)
    z_allrun_list.extend(mag_z)
    r_allrun_list.extend(r)
    theta_allrun_list.extend(theta)
    phi_allrun_list.extend(phi)
    times_allrun_list.extend(times)
    print(f"Run {run}:")
    print(f"  mag x: {np.mean(mag_x)}")
    print(f"  mag y: {np.mean(mag_y)}")
    print(f"  mag z: {np.mean(mag_z)}")
    x_mean_list.append(np.mean(mag_x))
    y_mean_list.append(np.mean(mag_y))
    z_mean_list.append(np.mean(mag_z))
    r_mean_list.append(np.mean(r))
    theta_mean_list.append(np.mean(theta))
    phi_mean_list.append(np.mean(phi))

plot_magnetometer_xyz(x_allrun_list, y_allrun_list, z_allrun_list, times_allrun_list, run="allruns")
plot_magnetometer_B(x_allrun_list, y_allrun_list, z_allrun_list, times_allrun_list, run="allruns")
plot_magnetometer_r(r_allrun_list, times_allrun_list, run="allruns")
plot_magnetometer_theta(theta_allrun_list, times_allrun_list, run="allruns")
plot_magnetometer_phi(phi_allrun_list, times_allrun_list, run="allruns")


print(f"length of x_mean_list: {len(x_mean_list)}")
print(f"length of y_mean_list: {len(y_mean_list)}")
print(f"length of z_mean_list: {len(z_mean_list)}")


#ellipse: AX2+BXY+CY2+DX+EY+F=0

def get_ellipse(h,k,a,b):
  x = []
  y = []
  # for iangle in np.linspace(0,2*np.pi,181):
  for iangle in np.linspace(0,2*np.pi,37):
    x.append(h + a * np.cos(iangle))
    y.append(k + b*np.sin(iangle))
  return np.asarray(x),np.asarray(y)

def rotate(x,y,theta):
  return x*np.cos(theta)-y*np.sin(theta),x*np.sin(theta)+y*np.cos(theta)


def get_rotated_ellipse(h,k,a,b,theta):
  x,y = get_ellipse(0,0,a,b)
  x,y = rotate(x,y,theta)
  return x+h,y+k

def get_noisy_ellipse(h,k,a,b,theta,noise_scale):
  x,y = get_rotated_ellipse(h,k,a,b,theta)
  x_noisy = [i+np.random.normal(0,noise_scale*abs(np.mean(x))) for i in x]
  y_noisy = [i+np.random.normal(0,noise_scale*abs(np.mean(y))) for i in y]
  return np.asarray(x_noisy), np.asarray(y_noisy)

x,y = get_noisy_ellipse(5,5,8,12,2*np.pi/7,0.03)


def fit_ellipse(x,y):
  x = np.asarray(x)
  y = np.asarray(y)
  A = np.stack([x**2,x*y,y**2,x,y]).T
  b = np.ones_like(x)
  w = np.linalg.lstsq(A, b, rcond=None)[0].squeeze()
  return w

def get_origin(A,B,C,D,E):
  return (2*C*D - B*E)/(B**2-4*A*C), (2*A*E-B*D)/(B**2-4*A*C)


def get_rotation(A,B,C,D,E):
  return 0.5*np.arctan2(-B,C-A)

def get_axes(A,B,C,D,E,F):
  t1 = B**2-4*A*C
  t2 = 2*(A*E**2+C*D**2-B*D*E+(B**2-4*A*C)*F)
  t3 = ((A+C)+np.sqrt((A-C)**2+B**2))
  t4 = ((A+C)-np.sqrt((A-C)**2+B**2))
  print(t1,t2,t3,t4)
  print(t2*t3/t1)
  return -np.sqrt(t2*t3)/t1, -np.sqrt(t2*t4)/t1

def get_eclipse_parameters(w):
  a,b = get_axes(*w,-1)
  theta = get_rotation(*w)
  x0,y0 = get_origin(*w)
  return np.asarray([x0,y0]),a,b,theta

def circle_from_ellipse(x,y,w):
  c0,a,b,theta = get_eclipse_parameters(w)
  r = np.sqrt(a*b)
#   ct, st = np.cos(theta), np.sin(theta)
  ct, st = np.cos(np.pi/2-theta), np.sin(np.pi/2-theta) #this works for multiple rotation
  R = np.array([[ct, -st],
                  [st,  ct]]) #rotation matrix for ellipse
  S = np.diag([a/r, b/r]) #scale for semi-major and semiminor axes
  A = (R.T)@S@R #to rotate the ellipse to align axes,scale to circle and rotate back
  x_circ = []
  y_circ = []
  for ix,iy in zip(x,y):
    ix_circ,iy_circ = (A@ (np.asarray([ix,iy])-c0).T).T
    x_circ.append(ix_circ)
    y_circ.append(iy_circ)
  return x_circ,y_circ


# utah_ellipse_params = [-0.00119216, 0.00016461, -0.00135584, 0.07516923, -0.03476016] #if using scaled values 10**6
utah_ellipse_params = [-1.19216116e+09,  1.64605814e+08, -1.35583796e+09,  7.51692271e+04, -3.47601605e+04]

def corrected_ellipse(x,y):
  w = fit_ellipse(x,y)
#   w = utah_ellipse_params
  x_circ,y_circ = circle_from_ellipse(x,y,w)
  return x_circ,y_circ



colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']

def plot_corrected_heading_360(bx_mean_list,by_mean_list,label=""):
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    bx_calibrated, by_calibrated = corrected_ellipse(bx_mean_list, by_mean_list)
    headings_original = [np.rad2deg(np.arctan2(by, bx)) for bx, by in zip(bx_mean_list, by_mean_list)]
    print(f"headings original: {headings_original}")
    headings_original = [i+360 if i<-1.0 else i for i in headings_original]
    print(f"headings original: {headings_original}")
    headings_corrected = [np.rad2deg(np.arctan2(byc, bxc)) for bxc,byc in zip(bx_calibrated, by_calibrated)]
    print(f"headings corrected: {headings_corrected}")
    headings_corrected = [i+360 if i<-1.0 else i for i in headings_corrected]
    # headings_corrected = [i+360 if 0<i<42 else i for i in headings_corrected]
    # headings_corrected = [i+360 if 0<i<34.0 else i for i in headings_corrected]
    # print(f"headings original: {headings_original}")
    print(f"headings corrected: {headings_corrected}")
    if label == "partial":
        ax.plot(np.asarray([0,1,2,3,4])*45, headings_original, "-o", c=colorsCustom[ncolor], label=f"{'raw'}", alpha=1)
        ax.plot(np.asarray([0,1,2,3,4])*45, headings_corrected, "--o", c=colorsCustom[ncolor+2], label=f"{'calibrated'}", alpha=1)
    else:
        ax.plot(np.asarray([0,1,2,3,4,5,6,7])*45, headings_original, "-o", c=colorsCustom[ncolor], label=f"{'raw'}", alpha=1)
        ax.plot(np.asarray([0,1,2,3,4,5,6,7])*45, headings_corrected, "--o", c=colorsCustom[ncolor+2], label=f"{'calibrated'}", alpha=1)
    ncolor += 1
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xlim(0,360)
    ###############################################
    ##############################################
    # ax.set_xticks(np.linspace(0,360,9)-roll)
    ax.set_xticks(np.linspace(0,360,9))
    # ax.set_xticks(np.linspace(0,360,9))
    ax.set_yticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(-180,180,9))
    #############################################
    #############################################
    # ax.set_aspect('equal')
    ax.legend(ncols=2,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" mDOM rotation [$^{\circ}$]", fontsize=20)
    ax.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_B_NTS{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_NTS{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()
 

print(f"length of x_mean_list: {len(x_mean_list)}")
print(f"length of y_mean_list: {len(y_mean_list)}")
print(f"length of z_mean_list: {len(z_mean_list)}")

plot_corrected_heading_360(x_mean_list[:5], y_mean_list[:5],label="partial")

PMT_labels = [8,7,6,5,4,11,10,9]
def plot_xy_360(mag_x, mag_y):
    loop_mag_x = mag_x
    loop_mag_y = mag_y
    # loop_mag_x.append(mag_x[0])
    # loop_mag_y.append(mag_y[0])
    loop_mag_x = np.array(loop_mag_x)*10**6 #convert to microTesla
    loop_mag_y = np.array(loop_mag_y)*10**6 #convert to microTesla
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.plot(loop_mag_x[:], loop_mag_y[:], "-o", c="b", label=f"{""}", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    for pmt,bx,by in zip(PMT_labels,loop_mag_x,loop_mag_y):
        ax.text(bx, by, f"{pmt}", fontsize=10, ha='right', va='bottom')
    ax.grid(True,alpha=0.6)
    ax.set_aspect('equal')
    # ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=20)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xy.pdf",transparent=False,bbox_inches='tight')
    plt.close()


# plot_xy_360(x_mean_list[:5], y_mean_list[:5]) #removing weird ponts for partial fits
plot_xy_360(x_mean_list[:], y_mean_list[:])

def plot_xyz_360(mag_x, mag_y, mag_z):
    loop_mag_x = mag_x
    loop_mag_y = mag_y
    loop_mag_z = mag_z
    # loop_mag_x.append(mag_x[0])
    # loop_mag_y.append(mag_y[0])
    loop_mag_x = np.array(loop_mag_x)*10**6 #convert to microTesla
    loop_mag_y = np.array(loop_mag_y)*10**6 #convert to microTesla
    loop_mag_z = np.array(loop_mag_z)*10**6 #convert to microTesla
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0],projection='3d')
    ax.plot(loop_mag_x[:], loop_mag_y[:], loop_mag_z[:], "-o", c="b", label=f"{""}", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=12)
    for pmt,bx,by,bz in zip(PMT_labels,loop_mag_x,loop_mag_y,loop_mag_z):
        ax.text(bx, by, bz, f"{pmt}", fontsize=10, ha='right', va='bottom')
    ax.grid(True,alpha=0.6)
    ax.set_aspect('equal')
    # ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=12)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=12)
    ax.set_zlabel(r" $B_z$ [$\mu$T]", fontsize=12)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xyz.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xyz.pdf",transparent=False,bbox_inches='tight')
    plt.show()
    # plt.close()


plot_xyz_360(x_mean_list[:], y_mean_list[:], z_mean_list[:])





def plot_xy_calibrated(mag_x, mag_y):
    
    # mag_x.append(mag_x[0])
    # mag_y.append(mag_y[0])
    mag_x = np.array(mag_x)*10**6 #convert to microTesla
    mag_y = np.array(mag_y)*10**6 #convert to microTesla
    mag_x_cal, mag_y_cal = corrected_ellipse(mag_x, mag_y)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.plot(mag_x, mag_y, "-o", c="b", label=f"{"raw"}", alpha=1)
    ax.plot(mag_x_cal, mag_y_cal, "-o", c="r", label=f"{"cal"}", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.set_aspect('equal')
    # ax.legend(loc="upper right",ncols=2,fontsize=16)
    ax.legend(fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=20)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xy_calibrated.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xy_calibrated.pdf",transparent=False,bbox_inches='tight')
    plt.close()


plot_xy_calibrated(x_mean_list[:5], y_mean_list[:5])