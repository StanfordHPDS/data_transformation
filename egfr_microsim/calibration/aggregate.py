import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser

import egfr_microsim.analysis.helpers.aggregate_trajectory_stats as aggregate_trajectory_stats
import egfr_microsim.model.helpers.path_helpers as path_helpers
# generate NaNs when dividing by zero
np.seterr(divide="ignore")

def generate_agg(args, sum_stats_dir):
    """
    Calculates a single aggregate count frame per parameter set  
    by averaging counts across all cohorts (rounding averages
    to integers)
    
    Args:
        args: command line arguments
        sum_stats_dir (str): summary statistics directory
    Returns:
        a pandas.DataFrame with one column
    """
    # empty frame for aggregation across files
    counts_agg_all = pd.DataFrame(
        index = aggregate_trajectory_stats.get_agg_all_index(),
        columns = list(range(args.cohort_number)),
    )
    # empty final frame
    counts_cats = pd.DataFrame(
        index = aggregate_trajectory_stats.get_agg_all_index(),
        columns = list(range(args.parameter_set_number)),
    )
    try:
        if args.parameter_set_id == None:
            agg_type = "df_agg_all_"
        else:
            agg_type = "df_agg_all_%i_" % args.parameter_set_id
    except: 
        agg_type = "df_agg_all_"

    # read in all aggregate count frames for a given experiment
    all_count_files = path_helpers.read_all_agg_counts(sum_stats_dir, args.cohort_number, agg_type)

    for param_id in range(args.parameter_set_number):
        counts_agg = counts_agg_all.copy()
        for cohort_id in range(args.cohort_number):
            counts_agg[cohort_id] = (
                all_count_files[cohort_id]
                .iloc[:, param_id]
                .reindex(counts_agg_all.index)
            )
        counts_cat = (
            counts_agg.sum(axis=1)
            .div(args.cohort_number)
            .round(0)
            .astype(int)
        )

        counts_cats.loc[:, param_id] = counts_cat

    return counts_cats

def combine(sum_stats_dir, parameter_set_group_number, overwrite = False):
    """
    Generates a single file combining all parameter set-level frames (single
    column) into a single frame, with one column per parameter set
    
    Args:
        sum_stats_dir (str): summary statistics directory
        parameter_set_group_number (int): parameter set file index
        overwrite (bool): if True, generate aggregated frame again, even if
            it already exists
    Returns:
        a pandas.DataFrame with one column per parameter set
    """
    final_path = os.path.join(sum_stats_dir, "counts_cat_all.feather")
    if os.path.exists(final_path) and not overwrite:
        df_full = pd.read_feather(final_path)
    else:
        agg_frames = []
        for i in range(parameter_set_group_number):
            df = pd.read_feather(
                os.path.join(sum_stats_dir, "counts_cat_all_%i.feather" % i)
            )
            df.columns = [el + i * 1000 for el in df.columns.astype(int)]
            agg_frames.append(df)
        df_full = pd.concat(agg_frames, axis=1)
        df_full.to_feather(final_path)
    return df_full

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
        default=1000,
        help = "number of parameter sets to run",
        type=int
    )
    parser.add_argument(
        "--parameter_set_id", 
        dest = "parameter_set_id", 
        help = "unique id of parameter set to run", 
        type = int
        )
    parser.add_argument(
        "-M",
        "--cohort_number", 
        dest = "cohort_number", 
        help = "number of cohorts to run", 
        type = int
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
        default = 1, 
        type = int,
        help = "number of cores to use (if parallelizing)"
    )

    args = parser.parse_args()

    exp_dir, sum_stats_dir = path_helpers.get_dirs(args, ["sum_stats"])

    counts_cats = generate_agg(args, sum_stats_dir)

    if args.parameter_set_id is not None:
        file_name = "counts_cat_all_%i" % args.parameter_set_id
    else:
        file_name = "counts_cat_all"

    counts_cats.to_feather(
        os.path.join(sum_stats_dir, "%s.feather" % file_name)
    )