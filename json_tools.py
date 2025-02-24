#!/usr/bin/env python

import json



def extract_temperature(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["meas_data"][0]["temperature"]

def extract_HV(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    hv = -100
    if data["meas_data"][1]["label"] == "applied HV":
        hv = data["meas_data"][1]["value"]
    else:
        print(f"{data["meas_data"][1]["label"]} is not applied HV")
    return hv

def extract_device(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    device_uid = data["device_uid"]
    pmt = data['subdevice_uid'].split('_')[-1]
    return device_uid, pmt