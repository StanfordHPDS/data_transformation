import pandas as pd
import os
from argparse import ArgumentParser, BooleanOptionalAction
import seaborn.objects as so
import matplotlib.pyplot as plt

import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.analysis.helpers.plots_read_agg as plots_read_agg
import egfr_microsim.calibration.coverage_analysis as coverage_analysis

def get_counts_cat_path(results_dir, traj_type):
    """
    Generate the path to experiment aggregate counts for a specific trajectory
    type (bl for reference and cf for counterfactual)
    """
    return os.path.join(results_dir, "data", "counts_cat_all_"+traj_type+".feather")

def update_count_frame(counts_cats, ages):
    """
    Format a provided aggregate counts frame for plotting
    """
    counts_cats_updated = (
        plots_read_agg.update_var(
            counts_cats.reset_index(),
            var_name = "male", 
            update_dict = {"Female": 0, "Male": 1}
        )
        .query("age in @ages")
        .assign(age = lambda x: x.age.astype(int))
        .set_index(["age", "dm", "ht", "stage", "male", "race"])
    )
    return counts_cats_updated
    
def read_agg(result_dir, ages):
    """
    Read all aggregate experiment results and include them in a list
    """

    counts = [] 
    
    for i,traj_type in enumerate(["bl", "cf"]):

        counts_cats = pd.read_feather(get_counts_cat_path(result_dir, traj_type))

        counts_cats = plots_read_agg.update_var(counts_cats, var_name='male')
        counts_cats = plots_read_agg.update_var(counts_cats, var_name='race')
        counts_cats = update_count_frame(counts_cats, ages)
        counts.append(counts_cats)
    return counts

def get_cov_frame(counts_cats, target, calib_strata, label, args):
    """
    Calculate calibration target coverage for a given table of 
    experiment counts
    """

    target, counts_cats_updated = coverage_analysis.update_strata(
        counts_cats,
        target, 
        calib_strata
    )
    
    calib_frac_frame = coverage_analysis.get_target_fractions(
        target, 
        counts_cats_updated, 
        calib_strata, 
        args
    )
    
    df_coverage = coverage_analysis.get_coverage(
        target, 
        calib_frac_frame, 
        calib_strata
    
    )
    df_coverage = plots_read_agg.update_var(
        df_coverage, 
        var_name = 'param',
        update_dict = {'target': 'data', 'sampled': label}
    )

    return df_coverage


def plot_coverage_paper(df, col_var, x_var="age"):
    """
    Generate manuscript-ready coverage plots with 
    updated labels for a given calibration target
    """

    strip_g = {
        "G1": "1",
        "G2": "2",
        "G3a": "3a",
        "G3b": "3b",
        "G4": "4",
        "G5": "5"
    }
    col_var_map = {
        "diabetes": {0: "Diabetes", 1: "No diabetes"},
        "male": {0:"Female", 1:"Male"},
        "hypertension": {0: "HTN", 1: "No HTN"},
        "ICE_q3": {1: "ICE Q1", 2: "ICE Q2", 3: "ICE Q3"},
        "SDI_q3": {1: "SDI Q1", 2: "SDI Q2", 3: "SDI Q3"},
        "param": {
            "data": "Cohort", 
            "counterfactual": "Counterfactual", 
            "reference": "Reference"}
    }
    updated_df = (
        df
        .reset_index()
        .replace({
            col_var:col_var_map[col_var], 
            "param": col_var_map["param"],
            "stage": strip_g})
    )
    if col_var in ["ICE_q3", "SDI_q3"]:
        size = (10,6)
    else:
        size = (10,4)
    p = (
        so.Plot(
            data=updated_df,#df.reset_index(),
            x=x_var,
            y="mean",
            color="param",
            ymin="min",
            ymax="max",
        )
        .facet(row=col_var, col="stage")
        .share(y="row")
        .add(so.Line())
        .add(so.Band())
        .add(so.Range())
        .layout(size=size)
        .label(
            x=x_var.capitalize(),
            y="Prevalence",
            color=""
        )
        .theme(plt.style.library["seaborn-v0_8-whitegrid"])
    )
    return p

def plot_coverage(df, col_var, x_var="age"):
    """
    Generate coverage plots for a given calibration target
    (for internal analyses)
    """
    p = (
        so.Plot(
            data = df.reset_index(),
            x = x_var,
            y = "mean",
            color = "param",
            ymin = "min",
            ymax = "max",
        )
        .facet(col = col_var, row="stage")
        .share(y = "row")
        .add(so.Line())
        .add(so.Band())
        .add(so.Range())
        .layout(size=(10, 10))
        .label(
            x = x_var.capitalize(),
            y = "Prevalence",
            color = "Parameter",
            col = "".join((col_var.capitalize(), ":")),
            row = "Stage:",
        )
    )
    return p
    
def generate(params, args, result_dir, exp_result_dir):
    """
    Generate coverage plots for the reference and counterfactual experiments
    compared to calibration targets. Save resulting plots.
    """

    min_age = params['cohort']['init_age'] 
    
    calib_targets = path_helpers.calib_targets(args.calib_targets, params)
    calib_targets = [el.query("age>@min_age") for el in calib_targets]
    ages = calib_targets[0].age.unique()
    
    counts_cats_bl, counts_cats_cf = read_agg(result_dir, ages)

    for k in range(len(args.calib_targets)):

        target = calib_targets[k]
        target_name = args.calib_targets[k].split(".")[0]
        calib_strata = coverage_analysis.get_calibration_strata(target)
        df_coverage_cf = get_cov_frame(counts_cats_cf, target, calib_strata, "counterfactual", args)
        df_coverage_bl = get_cov_frame(counts_cats_bl, target, calib_strata, "reference", args)
        df_coverage = pd.concat([df_coverage_bl, df_coverage_cf])

        if target_name in ["calib_targ_ICE", "calib_targ_SDI"]:
            # remove no location elements
            target = target.query("%s != 0" % calib_strata[1])

        p = plot_coverage_paper(df_coverage, col_var = df_coverage.index.names[-3])
    
        p.share(y = "col").save(
            os.path.join(
                exp_result_dir, 
                "plots", 
                "".join(("stage_distr_", target_name, "_share.png"))
            ), 
            bbox_inches = "tight",
            dpi = 800
        )

        p.save(
            os.path.join(
                exp_result_dir, 
                "plots", 
                "".join(("stage_distr_", target_name, "_noshare.png"))
            ),
            bbox_inches="tight",
            dpi = 800
        )


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest = "experiment_string",
        help = "unique experiment name",
    )
    parser.add_argument(
        "--params_base_path",
        dest = "params_base_path",
        help = "location of parameter dictionary",
        default = "parameter_files/params_base.json",
    )
    parser.add_argument(
        "-o",
        "--output_filename",
        dest = "output_filename",
        help = "name of output png file containing plot",
        default = "stage_diff_scenarios_time_egfr",
    )
    parser.add_argument(
        "--recalculate",
        dest = "recalculate",
        help = "recalculate data needed for plot directly from trajectories file",
        action = BooleanOptionalAction
    )
    parser.add_argument(
        "-t",
        "--trajectories_path",
        dest = "trajectories_path",
        help = "path to trajectories file (relative to repository top), in feather format",
        default = None
    )
    parser.add_argument(
        "--parameter_set_number",
        dest = "parameter_set_number",
        help = "Number of parameter sets in an experiment",
        default = 1
    )
    parser.add_argument(
        "--calib_targets",
        dest="calib_targets",
        default=[
           "calib_targ_DM.csv",
           "calib_targ_HT.csv",
           "calib_targ_sex.csv",
           "calib_targ_ICE.csv",
           "calib_targ_SDI.csv"
        ],
        nargs="+",
        help = "a list of calibration target files"
    )
    
    args = parser.parse_args()

    exp_dir, result_dir = path_helpers.get_dirs(
        args, ["results"]
    )
    exp_result_dir = path_helpers.get_exp_result_dir(args)
    params = path_helpers.get_params_base(args, allow_except=True)
    generate(params, args, result_dir, exp_result_dir)
