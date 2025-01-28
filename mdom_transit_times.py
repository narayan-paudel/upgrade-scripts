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
    mdom_list = json.load(fh)

# print(mdom_list)
for idom in mdom_list[:]:
    print(f"reading {idom} dom")
    subprocess.call(f"python get_mdom_transit_times.py {idom}",shell=True)

# folder_list = sorted(glob.glob(home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/*"))
folder_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/") if f.is_dir()]
folder_list = [ifolder.split("/")[-1] for ifolder in folder_list]
missing_mdom = [imdom for imdom in mdom_list if imdom not in folder_list]
print(f"There are {len(missing_mdom)} missing dom measurements")
print(missing_mdom)









# python degg_transit_times.py -i /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_degg.json