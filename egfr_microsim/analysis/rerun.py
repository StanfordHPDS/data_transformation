import os
import pandas as pd
from argparse import ArgumentParser
import egfr_microsim.analysis.helpers.aggregate_trajectory_stats as aggregate_trajectory_stats
import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.egfr_model as egfr_model
import egfr_microsim.model.helpers.param_helpers as param_helpers
import warnings
from pathlib import Path

# to suppress dtype warning triggered by
# egfr_microsim/results/aggregate_trajectory_stats.py:83
warnings.simplefilter(action="ignore", category=FutureWarning)

def get_posterior_params(params_base, posterior_param_sample):
    """
    Read in posterior parameters from the posterior_param_sample
    file
    """
    params = params_base.copy()
    params = param_helpers.update_slope_params(
        params, slopes_sample=posterior_param_sample
    )
    return params


def process_output(df_output, cohort_id):
    """
    Generate summary statistics for a trajectory data frame
    # TODO in the future: reuse code from calibration_batch
    """

    df_output = (
        df_output.assign(
            male = lambda x: x.sex.map({"F": 0, "M": 1}),
            cohort = cohort_id
        )
        .set_index("cohort", append = True)
    )

    df_agg_all = aggregate_trajectory_stats.get_counts_strata(
        df_output.reset_index(), strata=["age", "dm", "ht", "male", "race"]
    )

    return df_output, df_agg_all

def rerun(args, output_dir = None):
    """
    This script simulate trajectories for a pre-defined number of cohorts 
    based on the posterior calculated for the specified experiment.
    """
    exp_dir, loglik_dir, cohorts_dir, results_dir = path_helpers.get_dirs(
        args, ["loglik", "cohort_samples", "results"]
    )
    params_base = path_helpers.get_params_base(args)

    if output_dir is None:
        output_dir = results_dir
    else:
        output_dir = os.path.join(path_helpers.repo_path, output_dir)
        Path(output_dir).mkdir(parents=False, exist_ok=True)
        #os.mkdir(output_dir, exist_ok=True)
        print(output_dir)

    # TODO: document where this came from
    if args.posterior_file is None:
        try:
            df = pd.read_feather(
                os.path.join(loglik_dir, "sum_lik_param_posterior.feather")
            )
        except FileNotFoundError as e:
            e.strerror = "Posterior file does not exist; make sure to run posterior calculation or provide custom file path through command line arguments"
            raise e
        posterior = df.dropna(axis=1).query("param=='posterior'")
    else:
        posterior = pd.read_csv(args.posterior_file)

    # TODO: adapt this to sample from posterior rather than taking mean values
    params = get_posterior_params(
        params_base, posterior.mu.values
    )
    df_agg_all_bl_list = []
    df_agg_all_cf_list = []

    trajectories_bl = []
    trajectories_cf = []

    
    for cohort_id in range(args.cohort_number):
        df_cohort = param_helpers.get_cohort(cohorts_dir, cohort_id)

        df_bl, df_cf = egfr_model.get_trajectories_two_scenarios(
            df_cohort, params, crn_death = True
        )

        df_bl, df_agg_all_bl = process_output(df_bl, cohort_id)
        df_cf, df_agg_all_cf = process_output(df_cf, cohort_id)

        df_agg_all_bl_list.append(df_agg_all_bl)
        df_agg_all_cf_list.append(df_agg_all_cf)

        trajectories_bl.append(df_bl)
        trajectories_cf.append(df_cf)

    df_agg_all_bl = pd.concat(df_agg_all_bl_list, axis=1)
    df_agg_all_bl.columns = list(range(args.cohort_number))

    df_agg_all_cf = pd.concat(df_agg_all_cf_list, axis=1)
    df_agg_all_cf.columns = list(range(args.cohort_number))

    save_dir = os.path.join(output_dir, "data")
    os.makedirs(save_dir, exist_ok=True)
    df_agg_all_bl.to_feather(os.path.join(save_dir, "counts_cat_all_bl.feather"))
    df_agg_all_cf.to_feather(os.path.join(save_dir, "counts_cat_all_cf.feather"))

    pd.concat(trajectories_bl).to_feather(
        os.path.join(output_dir, "data", "trajectories_bl.feather")
    )
    pd.concat(trajectories_cf).to_feather(
        os.path.join(output_dir, "data", "trajectories_cf.feather")
    )

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest="experiment_string",
        help="unique experiment name",
    )
    parser.add_argument("--cohort_number", dest="cohort_number", type=int, default=100)
    parser.add_argument("--cohort_size", dest="cohort_size", type=int, default=10000)
    parser.add_argument("--interv_under_eq", dest="interv_under_eq", default="09", type=str)

    parser.add_argument(
        "--params_base_path",
        dest="params_base_path",
        default="parameter_files/params_base.json",
    )
    # this file must be a csv file with a column named mu and 8 rows corresponding to covariate combinations
    parser.add_argument(
        "--posterior_file", dest="posterior_file", type=str, default=None
    )
    parser.add_argument(
        "--output_dir", dest="output_dir", type=str, default=None
    )
    args = parser.parse_args()
    rerun(args, args.output_dir)
