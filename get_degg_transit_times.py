#!/usr/bin/env python

import os
import argparse
import json
from bson.objectid import ObjectId

from fatcat_db.forwarder import Tunnel
from fatcat_db.mongoreader import MongoReader

from pathlib import Path
home = str(Path.home())


# 2025-01-15 - example script to grab degg pmt transit time measurements
#Ubuntu or mac
# output_dir = "/home/narayan/research_ua/icecube/Upgrade/timing_calibration/data/degg_transit/"
output_dir = home+"/research_ua/icecube/Upgrade/timing_calibration/data/degg_transit/"


def main():

    cmdParser = argparse.ArgumentParser()
    cmdParser.add_argument('-nt', '--no-tunnel', dest='tunnel', action='store_false',
                           help='Do not port forward mongodb server')
    cmdParser.add_argument(dest='uid',
                           help='The DEgg device UID')
    args = cmdParser.parse_args()
    
    # open ssh tunnel to mongo port
    if args.tunnel:
        tunnel = Tunnel()
    
    # connect to mongo
    mongo = MongoReader(database='production_calibration')
    if not mongo.isConnected:
        return
    
    docs = list(mongo.db.measurements.find({'device_uid': args.uid,
                                            'meas_stage': 'fat',
                                            'meas_name': 'pmt-timing-resolution-info'}))
    
    print('{0} measurements found'.format(len(docs)))
    no_measurements = []
    
    if len(docs)==0:
        print(f"no measurement found for {args.uid}")

    pmts = []
    for n, doc in enumerate(docs):
        del doc['_id']
        del doc['insert_meta']
        
        degg = doc['device_uid']
        pmt = doc['subdevice_uid'].split('_')[-1]
        run = doc['run_number']

        deggdir = output_dir + degg
        if not os.path.exists(deggdir):
            os.makedirs(deggdir)
        filename = (deggdir+'/{0}_{1}_run{2}_n{3}.json'.format(degg, pmt, run, n+1))
        print('writing', filename)
        with open(filename, 'w') as jfile:
            json.dump(doc, jfile, separators=(', ', ': '), indent=4)

    
            
    del mongo

    
if __name__ == "__main__":
    main()

