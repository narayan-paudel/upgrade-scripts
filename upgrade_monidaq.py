#!/usr/bin/env python

import json

import subprocess
import os

monidaq_devices_path = "./monidaq_device_list.json"

BASE_URL="https://verical.icecube.wisc.edu"

data_dir = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/moni-data/"

VERICAL_USER="icecube"
VERICAL_PASS="skua"

with open(monidaq_devices_path, "r") as fh:
    device_list = json.load(fh)
    for fh in device_list[1:2]:
        print(fh["fieldhub_id"])
        print(len(fh["devices"]))
        for idevice in fh["devices"]:
            print(idevice["icm_id"], idevice["board_name"])
            command = f"curl -s -u {VERICAL_USER}:{VERICAL_PASS} -o {data_dir}{fh['fieldhub_id']}_{idevice['board_name']}_{idevice['icm_id']}.json {BASE_URL}/monidaq/in-ice/{fh['fieldhub_id']}/{idevice['icm_id']}/raw"
            subprocess.run(command, shell=True, capture_output=True, text=True)
        # curl -s -u $VERICAL_USER:$VERICAL_PASS -o monidaq_device_list.json "$BASE_URL/monidaq/device-list/raw?download=1"


#fieldhub92 {'icm_id': 'df0000290cb0ae2d', 'board_name': 'DEgg'}, {'icm_id': 'e2000030c1ebf72d', 'board_name': 'mDOM'}

for idevice in device_list[0]["devices"]:
    print(idevice["icm_id"], idevice["board_name"])


# https://verical.icecube.wisc.edu/monidaq/in-ice/fieldhub87/0300000f170f952d/raw
