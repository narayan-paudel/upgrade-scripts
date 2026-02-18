# upgrade-scripts

scripts for transit time calibration for IceCube Upgrade PMTs


To get list of mDOMs and DEggs from FAT database:


* Start the python environment (eg. venv_upgrade)
* cd /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts

To get updated list of FAT DEggs and mDOMs

* python list_deggs_or_mdoms.py degg
* python list_deggs_or_mdoms.py mdom

List saved as 
* /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/domlist/uids_mdom/degg.json

Getting transit time data for devices in the list
* python degg_transit_times.py -i ../data/domlist/uids_degg.json
* python degg_transit_times.py -i ../data/domlist/uids_mdom.json

Plot transit time
*  python plot_transit_time_mDOMs.py





