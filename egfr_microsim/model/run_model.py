import os
from argparse import ArgumentParser, BooleanOptionalAction

import egfr_microsim.model.egfr_model as egfr_model
import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.helpers.param_helpers as param_helpers
import egfr_microsim.model.helpers.sampling as sampling
import egfr_microsim.model.sample_param_sets as sample_param_sets

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest = "experiment_string",
        help = "unique experiment name",
    )
    parser.add_argument(
        "--sample", 
        dest = "sample", 
        help = "should cohorts and param_sets be sampled?", 
        action = BooleanOptionalAction
    )
    parser.add_argument(
        "-N", 
        "--cohort_size", 
        dest = "cohort_size", 
        help = "number of people in a single cohort", 
        type = int, 
        default = 1e4
    )
    parser.add_argument(
        "-eq",
        "--equation", 
        dest = "interv_under_eq", 
        help = "which equation (09 or 21) should be used to assign interventions?",
        type = str,
        default = "09"
        )
    parser.add_argument(
        "-i", 
        "--cohort_id", 
        dest = "cohort_id", 
        help = "if providing pre-sampled cohorts - unique id of cohort to run", 
        type = int, 
        default = 2137
        )
    parser.add_argument(
        "--param_set_id", 
        dest = "param_set_id", 
        help = "if providing pre-sampled parameter sets - unique id of parameter set to run", 
        type = int, 
        default = 2137
        )
    parser.add_argument(
        "--run_two_scenarios", 
        dest="run_two_scenarios", 
        help = "run 2 scenarios with common random numbers", 
        action=BooleanOptionalAction
        )
    parser.add_argument(
        "--print_timing", 
        dest="print_timing", 
        help = "print how long trajectory generation took", 
        action=BooleanOptionalAction
    )
    parser.add_argument(
        "--params_base_path",
        dest = "params_base_path",
        help = "location of parameter dictionary",
        default = "parameter_files/params_base.json",
    )
    parser.add_argument(
        "-r",
        "--rejection_sampling",
        dest = "rejection_sampling",
        help = "apply rejection sampling using provided parameter value ordering",
        default = True,
    )
    parser.add_argument(
        "-R", 
        "--parameter_set_number", 
        dest = "parameter_set_number", 
        help = "number of parameter sets to run",
        type=int
    )
    parser.add_argument(
        "--multiplier",
        dest="multiplier",
        help = "upsampling multiplier applied to parameter_set_number \
            - 16 needed for rejection sampling",
        type=int,
        default=16,
    )
    args = parser.parse_args()
    params_base = path_helpers.get_params_base(args)

    """
    we use helper functions to get standard path names where files will be
    read from and written to.
    
    exp_dir is the experiment directory, named the same as args.experiment_string
    cohorts_dir is a sub-directory of the exp_dir where sampled cohorts are stored
    traj_dir is a sub-directory of the exp_dir where sampled trajectories are stored
    """

    exp_dir, cohorts_dir, traj_dir, param_dir = path_helpers.get_dirs(
        args, ["cohort_samples", "trajectories", "param_samples"]
    )

    traj_string = "_".join(("traj", str(args.cohort_id), str(args.param_set_id)))
    traj_filename = os.path.join(traj_dir, "".join((traj_string, ".csv")))

    df_slopes = param_helpers.get_slopes(params_base)

    if args.sample:
        df_slopes = param_helpers.get_slopes(params_base)
        df_cohort = sampling.sample_cohort(params_base)
        #slopes_sample = sample_param_sets.sample_simple_slopes(df_slopes).slope.values
        slopes_sample = sample_param_sets.sample_all_parameter_sets(
            df_slopes, 
            args, 
            params_base, 
            param_dir = param_dir
        ).values[0] 

    else:
        df_cohort = param_helpers.get_cohort(cohorts_dir, args.cohort_id)
        df_param = param_helpers.get_param_df(param_dir)
        slopes_sample = df_param.loc[args.param_set_id, :]

    cohort_size = df_cohort.reset_index().pid.unique().shape[0]

    params = params_base.copy()
    params = param_helpers.update_slope_params(params, slopes_sample)

    if args.run_two_scenarios:
        print("in two scenarios")
        df_bl, df_cf = egfr_model.get_trajectories_two_scenarios(
            df_cohort, params
        )
        traj_filename_bl = os.path.join(traj_dir, "".join((traj_string, "_bl", ".csv")))
        traj_filename_cf = os.path.join(traj_dir, "".join((traj_string, "_cf", ".csv")))
        df_bl.to_csv(traj_filename_bl)
        df_cf.to_csv(traj_filename_cf)
        print(
            "generated a new trajectory file under", traj_filename_bl, traj_filename_cf
        )
    else:
        df_new = egfr_model.get_trajectories(
            df_cohort, params, print_timing=args.print_timing
        )
        df_new.to_csv(traj_filename)
        print("generated a new trajectory file under", traj_filename)