#!/usr/bin/env python

import subprocess
import os
import glob

import json
import argparse

from json_tools import extract_HV, extract_temperature, extract_device

from pathlib import Path
home = str(Path.home())

'''
Reads the json file with the list of DOMs and downloads
transit measurements from fatcat database

usage: python mdom_transit_times.py -i
 /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_mdom.json

'''

print(home)

cmdParser = argparse.ArgumentParser()
cmdParser.add_argument('-i', '--input',type=str, dest='input', 
                       help='input json file to read')
args = cmdParser.parse_args()

print("input file",args.input)

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
for imDOM in folder_list_desy[:2]:
    subdevice_list = []
    subdevice_dict = {}
    meas_list = glob.glob(data_path+f"{imDOM}/mDOM*")
    for ifile in meas_list:
        device, subdevice = extract_device(ifile)
        subdevice_list.append(subdevice)
    subdevice_list = list(set(subdevice_list))
    print(f"subdevices {subdevice_list}")
    for isub in subdevice_list:
        isub_meas = [ifile for ifile in meas_list if extract_device(ifile)[1] == isub]
        print(f"subdevice {isub} has {len(isub_meas)} meas\n")
    # for ifile in meas_list:
    #     temp = extract_temperature(ifile)
    #     hv = extract_HV(ifile)
    #     device, subdevice = extract_device(ifile)
    #     subdevice_dict[subdevice][]
    #     print(f"{subdevice} {temp} C {hv} V")
    #     subdevice_list.append(subdevice)
    # subdevice_list = list(set(subdevice_list))
    print(len(subdevice_list))

    # print(meas_list)
    print(f"mDOM {imDOM} has {len(meas_list)} meas")
# print(f"There are {len(missing_mdom)} missing dom measurements")
# print(missing_mdom)

# python mdom_meas_detail.py -i /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_mdom.json