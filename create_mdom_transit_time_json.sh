
echo "mdom transit time table of all runs"
python mdom_transit_time_table.py
echo "single transit time per pmt" 
python write_transit_time_single_per_pmt.py
echo "single transit time per pmt with gaussian refit"
python write_transit_time_single_gaussian_refit.py
echo "Copying the new json file to mdom_transit_times_v1.json"
VERSIONED_TT_FILE=/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_times_v1.json
cp /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time_gaussian_refit.json $VERSIONED_TT_FILE
python produce_mdom_transit_time_params.py -i $VERSIONED_TT_FILE -v 1 -o /Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/transit_time_params.json