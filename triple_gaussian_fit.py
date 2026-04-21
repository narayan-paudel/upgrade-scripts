import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erfc


from pathlib import Path
home = str(Path.home())


import glob
import json

from transit_time_dependence_plots import (extract_json, extract_json_tts_filter, extract_json_min_tts_filter, extract_json_tts_chi2_filter,plot_single_transit_time_histogram,
                                            prod_id_to_icm_id,extract_channel,extract_fit_params,gauss,get_pmt_uid,extract_run_number)


from refit_bad_tt_fits import merge_bins_reduceat
from transit_time_fit_functions import (plot_model_fit,plot_triple_gaussian_fit_residual,
                                        get_transit_time_data,triple_gaussian,gaussian,
                                        gaussian_rise_exp_decay,
                                        plot_double_gaussian_fit_residual,plot_single_gaussian_fit_residual,
                                       plot_gaussian_rise_exponential_decay_fit)



def main():
    upgrade_commissioning_scripts = home+"/research_ua/icecube/software/upgrade_commissioning_scripts/"
    geometry_files = sorted(glob.glob(upgrade_commissioning_scripts+"/geometry/string_*geometry*.json"))
    mdom_tt_dir = home+"/research_ua/icecube/upgrade/timing_calibration/data/mdom_transit/"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    ################################################################################
    #main file
    # transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time.json"
    #selected files
    transit_time_file = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_transit_time_select.json"
    ################################################################################
    ##########run picks file###########
    run_picks_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_run_picks.json"
    refit_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mdom_tt_needing_refit.json"
    empty_meas_json = home+"/research_ua/icecube/upgrade/timing_calibration/scripts/mDOM_tt_empty_meas.json"
    plotFolder: str = home+"/research_ua/icecube/Upgrade/timing_calibration/plots/mdom_transit"
    #####################################
    # x_values, y_values = get_transit_time_data("mDOM_M055", 13,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_M",filter_non_zero=False,check_outliers=None)
    # plot_model_fit(x_values, y_values,plotFolder)
    # plot_triple_gaussian_fit_residual(x_values, y_values,plotFolder)
    # plot_double_gaussian_fit_residual(x_values, y_values,plotFolder)
    # plot_single_gaussian_fit_residual(x_values, y_values,plotFolder)
    # plot_gaussian_rise_exponential_decay_fit(x_values, y_values,plotFolder)
    ############another example###########
    x_values, y_values = get_transit_time_data("mDOM_M139", 9,mdom_tt_dir,run_picks_json,refit_json,empty_meas_json,site_str="_M",filter_non_zero=False,check_outliers=None)
    # plot_model_fit(x_values, y_values,plotFolder)
    plot_gaussian_rise_exponential_decay_fit(x_values, y_values,plotFolder)
    plot_triple_gaussian_fit_residual(x_values, y_values,plotFolder)
    plot_double_gaussian_fit_residual(x_values, y_values,plotFolder)
    plot_single_gaussian_fit_residual(x_values, y_values,plotFolder)


if __name__ == "__main__":
    main()