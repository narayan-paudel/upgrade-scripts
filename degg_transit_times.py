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

# print(degg_list)
degg_exclude_list = ["DEgg2020-1-015_v2","DEgg2020-2-017_v1","DEgg2020-2-021_v1",
                     "DEgg2020-2-062_v1","DEgg2020-2-064_v1","DEgg2020-2-066_v1",
                     "DEgg2020-2-067_v1","DEgg2020-2-068_v1","DEgg2020-2-069_v1",
                     "DEgg2020-2-072_v1","DEgg2020-2-074_v1","DEgg2020-2-075_v1"]

for idom in degg_list[:]:
    if idom in degg_exclude_list:
        continue
    print(f"reading {idom} dom")
    subprocess.call(f"python get_degg_transit_times.py {idom}",shell=True)

DEggs_pass_list = []

# folder_list = sorted(glob.glob(home+"/research_ua/icecube/upgrade/timing_calibration/data/degg_transit/*"))
folder_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/degg_transit/") if f.is_dir()]
folder_list = [ifolder.split("/")[-1] for ifolder in folder_list]
missing_degg = [idegg for idegg in degg_list if idegg not in folder_list]
print(f"There are {len(missing_degg)} missing dom measurements")
print(missing_degg)









# python degg_transit_times.py -i /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_degg.json