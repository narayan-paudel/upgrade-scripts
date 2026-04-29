#!/usr/bin/env python

import subprocess
import os
import glob

import json
import argparse

from pathlib import Path
home = str(Path.home())

print(home)

cmdParser = argparse.ArgumentParser()
cmdParser.add_argument('-i', '--input',type=str, dest='input', 
                       help='input json file to read')
args = cmdParser.parse_args()

print("input file",args.input)

with open(args.input,"r") as fh:
    degg_list = json.load(fh)

def get_device_list(string,device,geometry_files):
    device_list = []
    for ifile in geometry_files:
        if string in ifile:
            with open(ifile, 'r') as f:
                data = json.load(f)
            for idev in data[0]["devices"]:
                if device in idev["production_id"]:
                    device_list.append(idev["production_id"])
    return device_list


upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"

geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))
deployed_device_list = (get_device_list("88","DEgg",geometry_files) + 
get_device_list("89","DEgg",geometry_files) + 
get_device_list("90","DEgg",geometry_files) + 
get_device_list("91","DEgg",geometry_files) + 
get_device_list("92","DEgg",geometry_files))

print(f"{len(get_device_list('88','DEgg',geometry_files))} DEggs on string 88")
print(f"{len(get_device_list('89','DEgg',geometry_files))} DEggs on string 89")
print(f"{len(get_device_list('90','DEgg',geometry_files))} DEggs on string 90")
print(f"{len(get_device_list('91','DEgg',geometry_files))} DEggs on string 91")
print(f"{len(get_device_list('92','DEgg',geometry_files))} DEggs on string 92")

# print(degg_list)
# degg_exclude_list = ["DEgg2020-1-015_v2","DEgg2020-2-017_v1","DEgg2020-2-021_v1",
#                      "DEgg2020-2-062_v1","DEgg2020-2-064_v1","DEgg2020-2-066_v1",
#                      "DEgg2020-2-067_v1","DEgg2020-2-068_v1","DEgg2020-2-069_v1",
#                      "DEgg2020-2-072_v1","DEgg2020-2-074_v1","DEgg2020-2-075_v1"]
degg_exclude_list = []

deployed_degg_list = [idegg for idegg in degg_list if idegg.split('_v')[0] in deployed_device_list]
missing_deggs = [idegg for idegg in deployed_device_list if idegg not in [jdegg.split('_v')[0] for jdegg in degg_list]]
#missing_deggs = ['DEgg2021-0-004', 'DEgg2021-0-001']



print(f"{len(deployed_degg_list)} out of {len(degg_list)} are deployed missing {len(missing_deggs)} deployed DEggs do not have measurements")
print("missing deployed DEggs", missing_deggs)

for idom in deployed_degg_list[:]:
    if idom in degg_exclude_list:
        continue
    print(f"reading {idom} dom")
    subprocess.call(f"python get_degg_linearity_meas.py {idom}",shell=True)

DEggs_pass_list = []

# folder_list = sorted(glob.glob(home+"/research_ua/icecube/upgrade/timing_calibration/data/degg_transit/*"))
folder_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/degg_linearity/") if f.is_dir()]
folder_list = [ifolder.split("/")[-1] for ifolder in folder_list]
missing_degg = [idegg for idegg in deployed_degg_list if idegg not in folder_list]
print(f"There are additional {len(missing_degg)} missing dom linearity measurements")
print(missing_degg)









# python degg_transit_times.py -i /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_degg.json