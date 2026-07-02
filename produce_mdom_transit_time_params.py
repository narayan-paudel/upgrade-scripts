
import pandas as pd
import sys
import os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import json
from IPython.display import display, JSON
import glob
from matplotlib.ticker import MaxNLocator





string_maps = glob.glob('/Users/epaudel/research_ua/icecube/pole_calibration/string_map/*json')
from transit_time_file_paths import geometry_files
from mdom_transit_time_constants import mean_tt

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("-i",'--input', type=str, default='/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_times_v1.json', help='Input file path for the transit time parameters (default: /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_times_v1.json)')
parser.add_argument("-v",'--version', type=int, default=1, help='Version of the transit time parameters to produce (default: 1)')
parser.add_argument("-o",'--output', type=str, default='/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/transit_time_params.json', help='Output file path for the transit time parameters (default: /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/transit_time_params.json)')
args = parser.parse_args()

json_data = {}
# json_data['version'] = 0
json_data['version'] = args.version

defaultval = False

with open(args.input, "r") as f:
    tt_data = json.load(f)

ctr = 0 
for item in geometry_files:
    with open(item, "r") as f:
        data = json.load(f)
        for i in data[0]['devices']:
            if i['type'] == 'mdom':
                found = False
                nchan = 24
                #print(i['icm_id'])
                for tt_item in tt_data:
                    if tt_data[tt_item]['icm_id'] == i['icm_id']:
                        found = True
                        ctr += 1
                        #print(tt_data[tt_item]['transit_times']['channel_0'][0]['sigma'])
                        #print(len(tt_data[tt_item]['transit_times']['channel_0']))
                        json_data[i['icm_id']] = {}
                        json_data[i['icm_id']]['string'] = data[0]['id']
                        json_data[i['icm_id']]['position'] = i['position']
                        json_data[i['icm_id']]['transit_time_data'] = []
                        for c in range(0,nchan,1):
                            chanfound = False
                            chan_id = 'channel_'+str(c)
                            if len(tt_data[tt_item]['transit_times'][chan_id]) < 1:
                                print(i['icm_id'],c)
                                json_data[i['icm_id']]['transit_time_data'].append({
                                    'channel': c,
                                    'transit_time': mean_tt['mu'],
                                    'transit_time_spread': mean_tt['sigma'],
                                    'default' : True,
                                    }) 
                                continue
                            else:
                                if not tt_data[tt_item]['transit_times'][chan_id][0]:
                                    continue
                            #print(c)
                            if len(tt_data[tt_item]['transit_times'][chan_id]) > 0:
                                if tt_data[tt_item]['transit_times'][chan_id][0]:
                                    defaultval = tt_data[tt_item]['transit_times'][chan_id][0]['default']
                                    json_data[i['icm_id']]['transit_time_data'].append({
                                        'channel': c,
                                        'transit_time': tt_data[tt_item]['transit_times'][chan_id][0]['mu'],
                                        'transit_time_spread': tt_data[tt_item]['transit_times'][chan_id][0]['sigma'],
                                        'default' : defaultval,
                                        })
                              
                if not found:
                    #print(i['icm_id'],posfinder(i['icm_id']),'not found')
                    nchan = 24
                    json_data[i['icm_id']] = {}
                    json_data[i['icm_id']]['string'] = data[0]['id']
                    json_data[i['icm_id']]['position'] = i['position']
                    json_data[i['icm_id']]['transit_time_data'] = []
                    for c in range(0,nchan,1):
                        json_data[i['icm_id']]['transit_time_data'].append({
                                    'channel': c,
                                    'transit_time': mean_tt['mu'],
                                    'transit_time_spread': mean_tt['sigma'],
                                    'default' : True,
                                    }) 

#print(ctr)

with open(args.output, 'w') as file:
            json.dump(json_data, file, indent=4)
print(f'Transit time parameters saved to {args.output}')
