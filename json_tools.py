#!/usr/bin/env python

import json

import numpy as np



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

def extract_run(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data["run_number"]


def extract_tt(json_file):
    with open(json_file,"r") as f:
        data = json.load(f)
    tt = -100
    if data["meas_data"][3]["label"] == "Transit time":
        tt = data["meas_data"][3]["value"]
    else:
        print(f"{data["meas_data"][3]["label"]} is not transit time for {json_file}")
    return tt

def extract_meas_values(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    param_dict = {"applied HV":np.nan,"mb channel":np.nan,
                  "Transit time":np.nan,"Transit time spread":np.nan,"a":np.nan,"b":np.nan,"c":np.nan,"chi2":np.nan}
    for j,ielt in enumerate(data["meas_data"][1:]):
        param_dict[ielt["label"]] = ielt["value"]
    return param_dict["applied HV"],param_dict["mb channel"],param_dict["Transit time"], param_dict["Transit time spread"],param_dict["a"]\
        ,param_dict["b"],param_dict["c"],param_dict["chi2"]