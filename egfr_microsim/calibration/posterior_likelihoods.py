from argparse import ArgumentParser

import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.calibration.posterior as posterior

"""
Calculate log likelihoods for a specified experiment, save
a file for each calibration target separately.
"""    
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest="experiment_string",
        help="unique experiment name",
    )
    parser.add_argument(
        "-R", 
        "--parameter_set_number", 
        dest = "parameter_set_number", 
        help = "number of parameter sets in the experiment",
        type = int
    )    
    parser.add_argument(
        "--params_base_path",
        dest = "params_base_path",
        help = "location of parameter dictionary",
        default = "parameter_files/params_base.json",
    )
    parser.add_argument(
        "--calib_targets",
        dest="calib_targets",
        default=[
           "calib_targ_DM.csv",
           "calib_targ_HT.csv",
           "calib_targ_sex.csv",
           "calib_targ_ICE.csv",
           "calib_targ_SDI.csv"],
        nargs="+",
    )

    args = parser.parse_args()

    ### Setup, reading in parameters
    exp_dir, params_dir, sum_stats_dir, loglik_dir, coverage_dir = path_helpers.get_dirs(
        args, ["param_samples", "sum_stats", "loglik", "coverage_analysis"]
    )
    params_base = path_helpers.get_params_base(args, allow_except=True)
    all_target_loss = posterior.calculate_likelihoods(
        params_base, 
        args, 
        sum_stats_dir, 
        loglik_dir
    )

