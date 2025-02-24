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

def extract_temperature(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["meas_data"][0]["temperature"]

def extract_device(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["device_uid"]
def extract_subdevice(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["subdevice_uid"]
def extract_yvalues(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["meas_data"][0]['y_values']

degg_list_pass = ["DEgg2020-1-001","DEgg2020-1-002","DEgg2020-1-003","DEgg2020-1-004",
                  "DEgg2020-1-005","DEgg2020-1-006","DEgg2020-1-007","DEgg2020-1-008",
                  "DEgg2020-1-009","DEgg2020-1-010","DEgg2020-1-011","DEgg2020-1-012",
                  "DEgg2020-1-013","DEgg2020-1-014","DEgg2020-1-016","DEgg2020-1-017",
                  "DEgg2020-1-018","DEgg2020-1-019","DEgg2020-1-020","DEgg2020-1-021",
                  "DEgg2020-1-022","DEgg2020-1-023","DEgg2020-1-024","DEgg2020-1-025",
                  "DEgg2020-1-026","DEgg2020-1-027","DEgg2020-1-028","DEgg2020-1-030",
                  "DEgg2020-1-031","DEgg2020-1-032","DEgg2020-1-033","DEgg2020-1-034",
                  "DEgg2020-1-035","DEgg2020-1-036","DEgg2020-1-037","DEgg2020-1-038",
                  "DEgg2020-1-039","DEgg2020-1-040","DEgg2020-1-041","DEgg2020-1-042",
                  "DEgg2020-1-043","DEgg2020-1-044","Degg2020-1-045","DEgg2020-1-046",
                  "Degg2020-1-047","DEgg2020-1-048","Degg2020-1-049","Degg2020-1-050",
                  "DEgg2020-2-001","DEgg2020-2-002","DEgg2020-2-003","DEgg2020-2-004",
                  "DEgg2020-2-005","DEgg2020-2-006","DEgg2020-2-007","DEgg2020-2-008",
                  "DEgg2020-2-010","DEgg2020-2-012","DEgg2020-2-013","DEgg2020-2-014",
                  "DEgg2020-2-015","DEgg2020-2-016","DEgg2020-2-017","DEgg2020-2-018",
                  "DEgg2020-2-019","DEgg2020-2-022","DEgg2020-2-023","DEgg2020-2-024",
                  "Degg2020-2-025","DEgg2020-2-026","DEgg2020-2-027","DEgg2020-2-028",
                  "DEgg2020-2-029","Degg2020-2-030","DEgg2020-2-031","DEgg2020-2-033",
                  "DEgg2020-2-034","DEgg2020-2-035","DEgg2020-2-036","DEgg2020-2-037",
                  "DEgg2020-2-038","DEgg2020-2-039","DEgg2020-2-040","DEgg2020-2-041",
                  "DEgg2020-2-042","DEgg2020-2-043","DEgg2020-2-044","DEgg2020-2-045",
                  "DEgg2020-2-047","DEgg2020-2-048","DEgg2020-2-049","DEgg2020-2-050",
                  "DEgg2020-2-051","DEgg2020-2-052","DEgg2020-2-053","DEgg2020-2-055",
                  "DEgg2020-2-056","DEgg2020-2-057","DEgg2020-2-059","DEgg2020-2-060",
                  "DEgg2020-2-061","DEgg2020-2-062","DEgg2020-2-063","DEgg2020-2-064",
                  "DEgg2020-2-065","DEgg2020-2-066","DEgg2020-2-067","DEgg2020-2-068",
                  "DEgg2020-2-070","DEgg2020-2-071","DEgg2020-2-072","DEgg2020-2-073",
                  "DEgg2020-2-074","DEgg2020-2-075","DEgg2020-2-076","DEgg2020-2-077",
                  "DEgg2020-2-078","DEgg2020-2-079","DEgg2020-2-080","DEgg2020-2-081",
                  "DEgg2020-2-084","DEgg2020-2-085","DEgg2020-2-086","DEgg2020-2-087",
                  "DEgg2020-2-088","DEgg2020-2-090","DEgg2020-2-091","DEgg2020-2-092",
                  "DEgg2020-2-093","DEgg2020-2-094","DEgg2020-2-095","DEgg2020-2-096",
                  "DEgg2020-2-097","DEgg2020-2-098","DEgg2020-2-099","DEgg2020-2-100",
                  "DEgg2021-0-001","DEgg2021-0-003","DEgg2021-0-004","DEgg2021-0-005",
                  "DEgg2021-3-001","DEgg2021-3-002","DEgg2021-3-004","DEgg2021-3-005",
                  "DEgg2021-3-006","DEgg2021-3-007","DEgg2021-3-008","DEgg2021-3-009",
                  "DEgg2021-3-010","DEgg2021-3-011","DEgg2021-3-012","DEgg2021-3-013",
                  "DEgg2021-3-014","DEgg2021-3-015","DEgg2021-3-016","DEgg2021-3-018",
                  "DEgg2021-3-019","DEgg2021-3-020","DEgg2021-3-021","DEgg2021-3-022",
                  "DEgg2021-3-023","Degg2021-3-024","DEgg2021-3-025","DEgg2021-3-027",
                  "DEgg2021-3-028","DEgg2021-3-029","DEgg2021-3-030","DEgg2021-3-031",
                  "DEgg2021-3-032","DEgg2021-3-033","DEgg2021-3-034","Degg2021-3-035",
                  "Degg2021-3-036","DEgg2021-3-037","DEgg2021-3-038","DEgg2021-3-039",
                  "DEgg2021-3-040","DEgg2021-3-041","DEgg2021-3-042","DEgg2021-3-043",
                  "DEgg2021-3-044","DEgg2021-3-045","DEgg2021-3-046","DEgg2021-3-047",
                  "DEgg2021-3-048","DEgg2021-3-049","DEgg2021-3-050","DEgg2021-3-051",
                  "DEgg2021-3-052","DEgg2021-3-053","DEgg2021-3-054","DEgg2021-3-055",
                  "DEgg2021-3-056","DEgg2021-3-057","DEgg2021-3-058","DEgg2021-3-059",
                  "Degg2021-3-060","DEgg2021-3-061","DEgg2021-3-062","DEgg2021-3-063",
                  "Degg2021-3-064","DEgg2021-3-065","DEgg2021-3-066","DEgg2021-3-068",
                  "DEgg2021-3-069","DEgg2021-3-070","Degg2021-3-071","DEgg2021-3-072",
                  "DEgg2021-3-073","Degg2021-3-074","DEgg2021-3-075","DEgg2021-3-076",
                  "DEgg2021-3-077","DEgg2021-3-078","DEgg2021-3-079","DEgg2021-3-080",
                  "DEgg2021-3-081","DEgg2021-3-084",
                  "DEgg2021-3-085","DEgg2021-3-087","DEgg2021-3-088","DEgg2021-3-089",
                  "DEgg2021-3-090","DEgg2021-3-091","DEgg2021-3-092","DEgg2021-3-093",
                  "DEgg2021-3-095","DEgg2021-3-096","DEgg2021-3-097","DEgg2021-3-098",
                  "DEgg2021-3-099","DEgg2021-3-101","DEgg2021-3-102","DEgg2021-3-103",
                  "DEgg2021-3-104","DEgg2021-3-105","DEgg2021-3-106","DEgg2021-3-108",
                  "DEgg2021-3-109","Degg2021-3-110","DEgg2021-3-111","DEgg2021-3-112",
                  "DEgg2021-3-113","DEgg2021-3-114","DEgg2021-3-115","DEgg2021-3-116",
                  "Degg2021-3-117","DEgg2021-3-118","DEgg2021-3-119","Degg2021-3-120",
                  "DEgg2021-3-121","DEgg2021-3-122","DEgg2021-3-123","DEgg2021-3-124",
                  "DEgg2021-3-125","DEgg2021-3-126","DEgg2021-3-127","DEgg2021-3-128",
                  "DEgg2021-3-129","DEgg2021-3-130","DEgg2021-3-131","DEgg2021-3-132",
                  "DEgg2021-3-133","DEgg2021-3-134","DEgg2021-3-136","DEgg2021-3-137",
                  "DEgg2021-3-138","DEgg2021-3-140",
                  "DEgg2021-3-141","DEgg2021-3-142","DEgg2021-3-143","DEgg2021-3-144",
                  "DEgg2021-3-145","DEgg2021-3-146","DEgg2021-3-147","DEgg2021-3-148",
                  "DEgg2021-3-149","DEgg2021-3-151","DEgg2021-3-152","DEgg2021-3-153",
                  "DEgg2021-3-154","DEgg2021-3-155","DEgg2021-3-156","DEgg2021-3-157",
                  "DEgg2021-3-158","DEgg2021-3-159","DEgg2021-3-160","DEgg2024-4-001",
                  "DEgg2024-4-002","DEgg2024-4-003","DEgg2024-4-004","DEgg2024-4-006"]

print(f"DEGG pass list {len(degg_list_pass)}")

# folder_list = sorted(glob.glob(home+"/research_ua/icecube/upgrade/timing_calibration/data/degg_transit/*"))
folder_list = [f.path for f in os.scandir(home+"/research_ua/icecube/upgrade/timing_calibration/data/degg_transit/") if f.is_dir()]
folder_list = [ifolder.split("/")[-1] for ifolder in folder_list]
missing_degg = [idegg for idegg in degg_list if idegg not in folder_list]
folder_list_name = [idegg.split("_v")[0] for idegg in folder_list]
folder_list_name = [idegg.lower() for idegg in folder_list_name]
# print(folder_list_name)
missing_degg_pass = [idegg for idegg in degg_list_pass if idegg.lower() not in folder_list_name]
found_degg_pass = [idegg for idegg in degg_list_pass if idegg.lower() in folder_list_name]
print(f"found pass degg {len(found_degg_pass)} missing {missing_degg_pass}")
print(f"There are {len(missing_degg)} missing dom measurements")
print(missing_degg)
n = 0
for idegg in folder_list[:2]:
    file_list = sorted(glob.glob(home+f"/research_ua/icecube/upgrade/timing_calibration/data/degg_transit/{idegg}/*"))
    file_list_20 = [ifile for ifile in file_list if extract_temperature(ifile) < -10] #seperates meas taken at ~ -20 C
    device_list = [extract_subdevice(ifile) for ifile in file_list_20]
    device_sets = list(set(device_list))
    for idevice in device_sets:
        meas_file = [ifile for ifile in file_list_20 if extract_subdevice(ifile) == idevice]
        latest_meas_list = []
        latest_meas = meas_file[0]
        ymeas_latest = extract_yvalues(latest_meas)
        print(f"first file {meas_file[1]}")
        for ifile in meas_file[1:]:
            iymeas = extract_yvalues(ifile)
            if iymeas == ymeas_latest:
                print(f"duplicate measurements in {idegg}")
                print(f"{iymeas} and {ymeas_latest}")
            else:
                print(f"unique measurements in {idegg} at similar temp")
                print(f"{extract_temperature(ifile)} and {extract_temperature(meas_file[0])}")

        if len(meas_file) > 1:
            print(f"{idevice}")
    if len(device_sets)>2:
        print(f"DOM {idegg} has {len(device_sets)} devices and {len(device_list)} meas")
        if len(file_list_20) == 4:
            n+=1
            print(f"DEgg {idegg} has {len(file_list_20)} meas files")
            for ifile in file_list_20:
                print(f"Meas file {ifile.split("/")[-1]} has temp {extract_temperature(ifile)}")
print(n)




# python degg_meas_detail.py -i /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_degg.json