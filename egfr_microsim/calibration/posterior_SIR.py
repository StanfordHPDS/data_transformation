from argparse import ArgumentParser
import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.helpers.param_helpers as param_helpers
import egfr_microsim.calibration.posterior as posterior
import egfr_microsim.calibration.aggregate as aggregate

"""
Calculate posterior with sample importance resampling for a specified experiment, 
based on specified calibration targets. Assumed that posterior_likelihoods completed
successfully.
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
            "calib_targ_DM_combined.csv",
            "calib_targ_HT_combined.csv",
            "calib_targ_ICE.csv",
            "calib_targ_SDI.csv"
        ],
        nargs="+",
        help = "a list of calibration target files"
    )
    parser.add_argument(
        "--resample_size",
        dest = "resample_size", 
        default = 1e4, 
        type = int,
        help = "size of the posterior sample"
    )

    args = parser.parse_args()

    ### Setup, reading in parameters
    exp_dir, params_dir, sum_stats_dir, loglik_dir, coverage_dir = path_helpers.get_dirs(
        args, ["param_samples", "sum_stats", "loglik", "coverage_analysis"]
    )
    params_base = path_helpers.get_params_base(args, allow_except=True)
    min_age = params_base['cohort']['init_age']

    ### Reading in calibration-specific parameters
    calib_targets = path_helpers.calib_targets(args.calib_targets, params_base)
    calib_targets = [el.query("age>@min_age") for el in calib_targets]
    
    # get 95%CI from prior distribution
    prior = posterior.get_prior(params_base)

    # q (to calculate p)
    q_counts_cats_all = aggregate.combine(sum_stats_dir, int(args.parameter_set_number/1000)).query("age>@min_age")
    all_target_loss = posterior.read_likelihoods(params_base, args, loglik_dir)

    # calculate a posterior over a simple mean of likelihoods across the 3 targets
    target_name = "sum_lik_param"
    all_target_loss_mean = all_target_loss.sum(axis=1)

    ### Calculating the posterior using sample importance resampling (SIR)
    try:
        param_set_df = param_helpers.get_param_df_from_1k(params_dir, args.parameter_set_number)
    except:
        param_set_df = param_helpers.get_param_df(params_dir)

    posterior_distr, sampled_ids = posterior.calculate_SIR(all_target_loss_mean, param_set_df, params_base, args)

    joint_frame = posterior.get_joint_frame(
        prior, 
        posterior_distr, 
        save = True, 
        target_name = "sum_lik_param", 
        loglik_dir = loglik_dir
        )

    posterior_plots = posterior.posterior_coverage_plots(
        q_counts_cats_all, 
        sampled_ids, 
        calib_targets, 
        args, 
        coverage_dir
    )
