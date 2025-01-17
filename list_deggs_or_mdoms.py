#!/usr/bin/env python

import os
import argparse
import json

from fatcat_db.forwarder import Tunnel
from fatcat_db.mongoreader import MongoReader


def main():

    cmdParser = argparse.ArgumentParser()
    cmdParser.add_argument('-nt', '--no-tunnel', dest='tunnel', action='store_false',
                           help='Do not port forward mongodb server')
    cmdParser.add_argument(dest='device',
                           help='Specify either degg or mdom')
    args = cmdParser.parse_args()
    
    # open ssh tunnel to mongo port
    if args.tunnel:
        tunnel = Tunnel()
    
    # connect to mongo
    mongo = MongoReader(database='production_calibration')
    if not mongo.isConnected:
        return

    xdevice = (args.device).lower()
    if 'degg' in xdevice:
        xdevice = 'degg'
    if 'mdom' in xdevice:
        xdevice = 'mdom'
    
    docs = list(mongo.db.devices.find({'device_type': xdevice}, {'uid': 1}))
    uids = sorted([doc['uid'] for doc in docs])

    # remove old v1 versions
    tmp_uids = uids
    for uid in tmp_uids:
        if uid[-2:] == 'v2':
            oldid = uid[:-2]+'v1'
            print('removing', oldid)
            uids.remove(oldid)
    # remove old v2 versions
    tmp_uids = uids
    for uid in tmp_uids:
        if uid[-2:] == 'v3':
            oldid = uid[:-2]+'v2'
            print('removing', oldid)
            uids.remove(oldid)
    
    filename = 'uids_'+xdevice+'.json'
    print('writing', filename)
    with open(filename, 'w') as jfile:
        json.dump(uids, jfile, separators=(', ', ': '), indent=4)
    
    
    
    del mongo

    
if __name__ == "__main__":
    main()

