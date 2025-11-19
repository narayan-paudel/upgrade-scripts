#!/usr/bin/env python

import os
import argparse
import json
from bson.objectid import ObjectId

from fatcat_db.forwarder import Tunnel
from fatcat_db.mongoreader import MongoReader

from pathlib import Path
home = str(Path.home())


# 2025-01-15 - example script to grab mdom pmt transit time measurements
output_dir = home+"/research_ua/icecube/Upgrade/timing_calibration/data/mdom_transit/"


def main():

    cmdParser = argparse.ArgumentParser()
    cmdParser.add_argument('-nt', '--no-tunnel', dest='tunnel', action='store_false',
                           help='Do not port forward mongodb server')
    cmdParser.add_argument(dest='uid',
                           help='The mDOM device UID')
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
                                            'meas_name': "mb_time_series"}))
    
    print('{0} measurements found'.format(len(docs)))

    pmts = []
    for doc in docs:
        del doc['_id']
        del doc['insert_meta']
        
        mdom = doc['device_uid']
        pmt = doc['subdevice_uid']
        pmts.append(pmt)
        run = doc['run_number']

        mdomdir = output_dir+mdom
        if not os.path.exists(mdomdir):
            os.makedirs(mdomdir)
        filename = (mdomdir+'/{0}_{1}_{2}.json'.format(mdom, pmt, run))
        print('writing', filename)
        with open(filename, 'w') as jfile:
            json.dump(doc, jfile, separators=(', ', ': '), indent=4)

    
    for pmt in pmts:
        dm = pmt.split('_')[-1]
        docs = list(mongo.db.measurements.find({'subdevice_uid': pmt,
                                                'meas_site': 'aachen',
                                                'meas_stage': 'calibration',
                                                'meas_name': 'arrivaltime-spectrum'}))
        for n, doc in enumerate(docs):
            del doc['_id']
            del doc['insert_meta']
            del doc['daq_config']

            pmtdir = mdomdir+'/'+pmt
            if not os.path.exists(pmtdir):
                os.makedirs(pmtdir)
            filename = (pmtdir+'/{0}_{1}.json'.format(pmt, n+1))
            print('writing', filename)
            with open(filename, 'w') as jfile:
                json.dump(doc, jfile, separators=(', ', ': '), indent=4)
            
            
    del mongo

    
if __name__ == "__main__":
    main()

