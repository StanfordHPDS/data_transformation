import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser, BooleanOptionalAction

import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.analysis.plot_cover as plot_cover
import egfr_microsim.calibration.aggregate as aggregate


def get_calibration_strata(target):
    """
    Get a list of variables defining the calibration stratum
    Args:
        target (pandas.DataFrame): a target frame

    Returns:
        a list<str> of stratum variables
    """
    calib_strata = target.drop(
            columns=["x", "N", "frac"]
        ).columns.to_list()
    return calib_strata
    
def get_stderr(df):
    """
    Calculate binomial standard error for each row of the target table,
    as well as min and max values corresponding to the 95% CI

    Args:
        df (pandas.DataFrame): a frame containing the following columns
            frac: a fraction of counts in a category
            N: a total number of trials
    Returns:
        the same frame with three additional columns 
        (stderr, min and max corresponding to 95% CI)
    """

    return df.assign(
        stderr = lambda x: np.sqrt((x.frac * (1 - x.frac)) / x.N),
        min = lambda x: x.frac - 1.96 * x.stderr,
        max = lambda x: x.frac + 1.96 * x.stderr,
    )

def get_target_distributions(target, calib_strata):
    """
    Calculate multinomial standard deviation for each stratum
    to obtain the distribution of the calibration targets
    
    Args:
        target: x, N and frac (x/N) coming from target data
        calib_strata: variables defining a calibration stratum (column names)
    Returns:
        Distribution summaries (mean and 95% confidence interval)
    """    
    target = get_stderr(target)

    target_distr = (
        target.set_index(calib_strata)
        .rename(columns = {"frac": "mean"})
        .filter(["mean", "min", "max"])
        .assign(param = "target")
    )
    return target_distr

def get_coverage_distribution(df):
    """
    Calculate distribution of probabilities within each stratum,
    across all parameter sets (mean and 95% CI)

    Args:
        df: A frame with fractions corresponding to the distributions of
            experiment counts within a particular calibration target
    Returns:
        Distribution summaries (mean and 95% confidence interval) per row
    """
    distribution = (
        df.assign(
            mean = df.mean(axis = 1),
            min = df.quantile(q = 0.025, axis = 1),
            max = df.quantile(q = 0.975, axis = 1),
        )
        .filter(["mean", "min", "max"])
        .assign(param="sampled")
    )
    return distribution


def combine_prior_post(prior, posterior, calib_strata_stage):
    """
    Make a summary table with prior and posterior parameter distributions
    (mean and 95% CI)
    """
    df_coverage = (
        pd.concat((prior, posterior))
        .reset_index()
        .set_index(calib_strata_stage + ["param"])
    )
    return df_coverage


def get_target_fractions(calib_target, counts_cats, calib_strata, args):
    """
    Takes a frame of experiment counts (counts_cats) and 
    calculate corresponding fractions

    Args:
        calib_target: a table containing the calibration target 
        counts_cats: a table containing experiment summary counts
        calib_strata: calibration strata
        args: command line arguments

    Returns:
        a table derived from counts_cats that contains fractions 
        corresponding to the distribution of simulated individuals
        across calibration strata 
    """
    # generate a frame of probabilities (calibration targets p)
    calib_target = calib_target.set_index(calib_strata)
    
    # for coverage analysis
    # just to compare against calibration targets
    calib_frac_frame = pd.DataFrame(
        index=calib_target.index, columns=list(range(args.parameter_set_number))
    )

    # counts_cats come from calibration (X)
    # iterate across all parameter sets
    # pick the counts for that set
    calib_frac_frames = []
    for j in range(counts_cats.shape[1]):
        calib_target_count = (
            counts_cats[j]
            .reset_index()
            .rename(columns={j: "counts"})
            .groupby(calib_strata)
            .sum()
            .filter(["counts"])
            .rename(columns={"counts": "x_counts"})
        )
        # coverage analysis
        calib_target_frac = (
            calib_target_count / calib_target_count.groupby(calib_strata[:-1]).sum()
        )
        calib_frac_frames.append(calib_target_frac)
    calib_frac_frame = pd.concat(calib_frac_frames, axis=1)
    calib_frac_frame.columns = list(range(counts_cats.shape[1]))
    return calib_frac_frame



def get_coverage(target, calib_frac_frame, calib_strata):
    """
    Calculates coverage of a specified calibration targets by simulation results.

    Args:
        calib_strata: variables defining a calibration stratum (column names)
        target: x, N and frac (x/N) coming from target data
        calib_frac_frame: p coming from experiments

    Returns:
        a summary table with target and simulated coverage distributions
        mean and 95% CI)
    """
    # prior distribution of prevalence across stages
    target_distr = get_target_distributions(target, calib_strata)
 
    # coverage distribution of prevalence across stages
    experiment_distr = get_coverage_distribution(
        calib_frac_frame
    )
    
    df_coverage = combine_prior_post(
        target_distr, 
        experiment_distr, 
        calib_strata
    ).astype(float)

    replace_dict = {"dm": "diabetes", "ht": "hypertension"}
    index_names = df_coverage.index.names
    index_names_updated = [
        replace_dict.get(a) if replace_dict.get(a) else a for a in index_names
    ]

    df_coverage.index.names = index_names_updated

    return df_coverage

def update_target_combined_g3(target, calib_strata):
    """
    Aggregate stages G3a and G3b into G3 in the target table

    Args:
        target: target table
        calib_strata: a list of variables defining the target strata
    Returns:
        An updated target frame
    """
    target =  target.assign(
        stage = lambda x: x.stage.map(
            {"G1G2": "G1G2",
             "G1": "G1",
             "G2": "G2",
             "G3a": "G3", 
             "G3b": "G3",
             "G4": "G4",
             "G5": "G5"
             }
             ))
    target = target.groupby(calib_strata).agg({"x": sum, "N": max, "frac": sum})
    return target.reset_index()

def update_exp_agg_combined_g1g2(calib_frac_frame, calib_strata):
    """
    Aggregate stages G1 and G2 into a combined stage G1G2 in the provided frame

    Args:
        calib_frac_frame: a table derived from counts_cats that contains
            fractions corresponding to the distribution of simulated 
            individuals across calibration strata 
        calib_strata: a list of variables defining the target strata
    Returns:
        An updated target frame
    """
    orig_columns = calib_frac_frame.columns
    calib_frac_frame =  (calib_frac_frame
                         .reset_index()
                         .assign(
                             stage = lambda x: x.stage.map(
                                 {"G1G2": "G1G2",
                                  "G1": "G1G2",
                                  "G2": "G1G2",
                                  "G3": "G3", 
                                  "G3a": "G3a", 
                                  "G3b": "G3b",
                                  "G4": "G4",
                                  "G5": "G5"}
                             )))
                         
    calib_frac_frame = calib_frac_frame.groupby(calib_strata)[orig_columns].sum()
    return calib_frac_frame

def update_strata_stages(q_counts_cats_all, target, calib_strata):
    """
    Update stages (G3a, G3b --> G3 or G1, G2 --> G1G2) to make them consistent
    between calibration target and the aggregate counts frame

    Args:
        q_counts_cats_all: aggregate counts frame
        target: target frame
        calib_strata: a list of variables defining the target strata
    Returns:
        An updated target frame and aggregate counts frame
    """

    if "G3" in q_counts_cats_all.reset_index().stage.unique():
        target = update_target_combined_g3(target, calib_strata)
    
    if "G1G2" in target.reset_index().stage.unique():
        q_counts_cats_all = update_exp_agg_combined_g1g2(q_counts_cats_all, calib_strata)
    return target, q_counts_cats_all

def update_strata(q_counts_cats_all, target, calib_strata):
    """
    Update strata to make them consistent
    between calibration target and the aggregate counts frame

    Args:
        q_counts_cats_all: ggregate counts frame
        target: target frame
        calib_strata: a list of variables defining the target strata
    Returns:
        An updated target frame and aggregate counts frame
    """
    target, q_counts_cats_all = update_strata_stages(q_counts_cats_all, target, calib_strata)
    try:
        calib_strata_addl = [el for el in calib_strata if el not in ["age", "stage"]][0]
    except:
        calib_strata_addl = []
    if calib_strata_addl in ["ICE_q3", "SDI_q3"]:    
        sdvi_strata = []
        for i in range(1,4):
            sdvi_strata.append(
                q_counts_cats_all
                .assign(**{calib_strata_addl: i})
                .set_index(calib_strata_addl, append=True)
            )
        q_counts_cats_all = pd.concat(sdvi_strata)
        target = target.query("%s != 0" % calib_strata_addl)
    
    return target, q_counts_cats_all
    
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest = "experiment_string",
        help = "unique experiment name",
    )
    parser.add_argument(
        "-R", 
        "--parameter_set_number", 
        dest = "parameter_set_number", 
        help = "number of parameter sets in the experiment",
        type = int
    )
    parser.add_argument(
        "-M",
        "--cohort_number", 
        dest = "cohort_number", 
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
        "--recalculate", 
        dest = "recalculate", 
        action = BooleanOptionalAction
    )

    args = parser.parse_args()

    exp_dir, sum_stats_dir, coverage_dir = path_helpers.get_dirs(
        args, ["sum_stats", "coverage_analysis"]
    )
    params_base = path_helpers.get_params_base(args, allow_except=True)

    min_age = params_base['cohort']['init_age']    
    q_counts_cats_all = aggregate.combine(sum_stats_dir, int(args.parameter_set_number/1000)).query("age>@min_age")

    calib_targets = path_helpers.calib_targets(args.calib_targets, params_base)
    calib_targets = [el.query("age>@min_age") for el in calib_targets]

    # for each parameter set, calculate a frame of probabilities corresponding to
    # the calibration target's categories
    for k in range(len(args.calib_targets)):
        target = calib_targets[k]
        target_name = args.calib_targets[k].split(".")[0]
        calib_strata = get_calibration_strata(target)

        target, q_counts_cats_all1 = update_strata(q_counts_cats_all, target, calib_strata)
        calib_strata = get_calibration_strata(target)
        # if previously saved, don't re-calculate
        calib_frac_frame_path = os.path.join(
                    coverage_dir, "tables", "_".join((target_name, "frac.feather"))
                )    
        try:
            assert not args.recalculate
            calib_frac_frame = pd.read_feather(calib_frac_frame_path)
            calib_strata = get_calibration_strata(target)
        except:
            calib_frac_frame = get_target_fractions(target, q_counts_cats_all1, calib_strata, args)
            calib_frac_frame.to_feather(calib_frac_frame_path)
        
        df_coverage = get_coverage(target, calib_frac_frame, calib_strata)
        col_var = [el for el in df_coverage.index.names if el not in ["age", "stage"]][0]
        p = plot_cover.plot_coverage(
            df_coverage, x_var="age", col_var=col_var
        )
        p.save(
            os.path.join(
                coverage_dir, "plots", "_".join((target_name, "coverage.png"))
            ),
            bbox_inches="tight",
        )
        df_coverage.to_feather(
            os.path.join(
                coverage_dir, "tables", "_".join((target_name, "coverage.feather"))
            )
        )
