import pandas as pd
import functools
from copy import deepcopy
import time
from tqdm import tqdm
from argparse import ArgumentParser, BooleanOptionalAction
from concurrent.futures import ProcessPoolExecutor

import egfr_microsim.analysis.helpers.aggregate_trajectory_stats as aggregate_trajectory_stats
import egfr_microsim.model.egfr_model as egfr_model
import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.helpers.param_helpers as param_helpers

def run_single_param(
    df_init,
    params_base,
    df_param,
    agg_index,
    max_attempts,
    param_set_id
):
    """
    Run simulation experiments for a single cohort and single parameter set.
    Repeat for up to max_attempts if it fails.

    Args:
        df_init (pandas.DataFrame): initial sampled cohort frame
        params_base (dict): parameters dictionary, before updating slope parameters
        df_param (pandas.DataFrame): a data frame with multiple parameter sets
        agg_index (pandas.DataFrame.index): index to be used for aggregated result frames
        max_attempts (int): Maximum number of attempts re-running a single 
            experiment if it fails
        param_set_id (int): the index by which df_param will be queried to identify
            the desired parameter set
    Returns:
        An aggregated summary of the simulation
    """
    params = deepcopy(params_base)
    params = param_helpers.update_slope_params(
        params, df_param.loc[param_set_id, :]
    )

    traj_complete = False
    i = 0
    while not traj_complete:
        try:
            df_new = egfr_model.get_trajectories(df_init, params)
            traj_complete = True
        except:
            if i >= (max_attempts - 1):
                raise RuntimeError("Tried running trajectory %i times, failed, param id %i" % (max_attempts,param_set_id))
            else:
                i += 1
                print("Tried running trajectory and failed for param id %i, attempt #%i" % (param_set_id, i))

    df_new = df_new.assign(male = lambda x: x.sex.map({"F": 0, "M": 1}))

    df_agg_all = aggregate_trajectory_stats.get_counts_strata(
        df_new.reset_index()
    ).reindex(agg_index)

    return df_agg_all


def run_parallel(func, n_cores, arguments, show_progress = False):
    """
    Run the provided function for each parameter with parallelization

    Args:
        func (functools.partial): a function to run
        n_cores (int): number of cores for parallelization
        arguments (iterable): an interable of arguments to pass to func
        show_progress (optional, bool): if True, show a progress bar
    Returns:
        A list of outputs (one per element in arguments)
    """
    with ProcessPoolExecutor(n_cores) as exe:
        if show_progress:
            df_aggs = list(tqdm(exe.map(func, arguments), total=len(arguments)))
        else:
            df_aggs = exe.map(func, arguments)
    return df_aggs


def run_non_parallel(func, arguments, show_progress = False):
    """
    Run the provided function for each parameter, without
    parallelization

    Args:
        func (functools.partial): a function to run
        arguments (iterable): an interable of arguments to pass to func
        show_progress (optional, bool): if True, show a progress bar
    Returns:
        A list of outputs (one per element in arguments)
    """
    df_aggs = []
    iterable = tqdm(arguments) if show_progress else arguments
    for param_set_id in iterable:
        df_aggs.append(func(param_set_id))

    return df_aggs


def one_cohort_many_params(
    cohort_id, params, args, cohorts_dir, df_param
):
    """
    Run simulation experiments for a single cohort for 
    each provided parameter set

    Args:
        cohort_id (int): unique cohort id
        params (dict): parameters dictionary
        args: command line arguments
        cohorts_dir (str): cohort samples directory
        df_param (pandas.DataFrame): parameters frame (one row per parameter set)
    Returns:
        a pandas.DataFrame with columns containing summaries of
            simulation experiments (one per parameter set)
    """
    df_init = param_helpers.get_cohort(cohorts_dir, cohort_id)
    cohort_size = df_init.reset_index().pid.unique().shape[0]
    params["cohort"]["n"] = cohort_size

    # prepare frames for saving aggregate values
    counts_agg_all = get_agg_frames(args.parameter_set_number)

    # transformation of the run_single_param function, with all but the last function
    # argument filled in - will be used below for parallelizing function calls
    partial_run_single_param = functools.partial(
        run_single_param,
        df_init,
        params,
        df_param,
        counts_agg_all.index,
        args.max_attempts
    )
    
    if int(args.n_cores) > 1:
        df_aggs = run_parallel(
            partial_run_single_param, 
            args.n_cores, 
            range(args.parameter_set_number),
            args.show_progress
            )

    else:
        df_aggs = run_non_parallel(
            partial_run_single_param, 
            range(args.parameter_set_number),
            args.show_progress
            )

    # move elements from list to a single frame
    for i, df_agg_i in enumerate(df_aggs):
        counts_agg_all[i] = df_agg_i

    return counts_agg_all


def get_agg_frames(parameter_set_number):
    """
    Generate an empty frame for storing experiment summaries,
    with index corresponding to aggregation strata and one column
    per parameter set
    
    Args:
        parameter_set_number: number of parameter sets
    Returns:
        an empty pandas.DataFrame
    """

    counts_agg_all = pd.DataFrame(
        index = aggregate_trajectory_stats.get_agg_all_index(),
        columns = list(range(parameter_set_number)),
    )

    return counts_agg_all


def save_agg_frame(frame, sum_stats_dir, cohort_id, filename):
    """
    Save a summary aggregate file, with NaN rows removed, 
    and other NaNs filled with 0 (for storage efficiency)
    
    Args:
        frame: frame to save
        sum_stats_dir (str): summary statistics directory
        cohort_id: unique cohort id
        filename: name of file to save
    Returns:
        None
    """

    counts_agg_filename = path_helpers.agg_counts_path(
        sum_stats_dir, cohort_id, filename = filename
    )
    frame_updated = frame.dropna(how="all").fillna(0).astype("int")
    frame_updated.to_feather(counts_agg_filename)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest="experiment_string",
        help="unique experiment name",
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
        "-i", 
        "--cohort_id", 
        dest = "cohort_id", 
        help = "if providing pre-sampled cohorts - unique id of cohort to run", 
        type = int
        )
    parser.add_argument(
        "--parameter_set_id", 
        dest = "parameter_set_id", 
        help = "id of parameter set file", 
        type = int
        )
    parser.add_argument(
        "--parameter_set_number", 
        dest = "parameter_set_number", 
        help = "Number of parameter sets in a file", 
        type = int 
        )
    parser.add_argument(
        "-eq",
        "--interv_under_eq", 
        dest = "interv_under_eq", 
        help = "which equation (09 or 21) should be used to assign interventions?",
        type = str,
        default = "09"
        )
    parser.add_argument(
        "--params_base_path",
        dest = "params_base_path",
        help = "location of parameter dictionary",
        default = "parameter_files/params_base.json",
    )
    parser.add_argument(
        "--n_cores",
        dest = "n_cores",
        help = "Number of cores for parallelization",
        default = 1,
        type = int
    )
    parser.add_argument(
        "--show_progress",
        dest = "show_progress",
        help = "Display progress bar",
        action = BooleanOptionalAction
    )
    parser.add_argument(
        "--max_attempts",
        dest = "max_attempts",
        help = "Maximum number of attempts re-running a single experiment if it fails",
        default = 5,
        type = int
    )
    args = parser.parse_args()

    t_beg = time.time()
    print("Started cohort", args.cohort_id)

    exp_dir, cohorts_dir, params_dir, sum_stats_dir = path_helpers.get_dirs(
        args, ["cohort_samples", "param_samples", "sum_stats"]
    )
    params_base = path_helpers.get_params_base(args)
    df_param = param_helpers.get_param_df(params_dir, i=args.parameter_set_id)

    counts_agg_all = one_cohort_many_params(
        args.cohort_id, 
        params_base, 
        args, 
        cohorts_dir, 
        df_param
    )

    if args.parameter_set_id is not None:
        filename =  "df_agg_all_%i_" % args.parameter_set_id
    else:
        filename = "df_agg_all_"

    save_agg_frame(
        counts_agg_all, 
        sum_stats_dir, 
        args.cohort_id, 
        filename = filename
    )

    t_end = time.time()

    print("completed a calibration job for cohort", args.cohort_id)

    print("calibration took", round((t_end - t_beg) / 60, 2), "minutes")
    print(
        "time per parameter set:",
        round((t_end - t_beg) / (args.parameter_set_number), 2),
        "seconds",
    )

