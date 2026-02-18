#!/usr/bin/env python

import subprocess
import os
import glob

import json
import argparse

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

# print(mdom_list)

mdom_exclude_list = ["mDOM_D032_v1","mDOM_D035_v1","mDOM_D036_v1",
                     "mDOM_D041_v1","mDOM_D047_v1","mDOM_D070_v1","mDOM_D071_v1",
                     "mDOM_D072_v1","mDOM_D074_v1","mDOM_D075_v1","mDOM_D076_v1",
                     "mDOM_D079_v1"," mDOM_D124_v2","mDOM_D147_v1","mDOM_D148_v1",
                     "mDOM_D178_v1","mDOM_M015_v1","mDOM_M081_v1","mDOM_M168_v1",
                     "mdom_DVT_01_v2","mdom_DVT_02_v2","mdom_DVT_03_v2","mdom_DVT_04_v3",
                     "mdom_DVT_05_v2","mdom_DVT_07_v1","mdom_DVT_08_v1","mdom_DVT_09_v1",
                     "mdom_DVT_10_v1","mdom_DVT_11_v1"]

for idom in mdom_list[:]:
    if idom in mdom_exclude_list:
        continue
    print(f"reading {idom} dom")
    subprocess.call(f"python get_mdom_transit_times.py {idom}",shell=True)

# folder_list = sorted(glob.glob(home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/*"))
folder_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/") if f.is_dir()]
folder_list = [ifolder.split("/")[-1] for ifolder in folder_list]
missing_mdom = [imdom for imdom in mdom_list if imdom not in folder_list]
print(f"There are {len(missing_mdom)} missing dom measurements")
print(missing_mdom)









# python degg_transit_times.py -i /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_degg.json