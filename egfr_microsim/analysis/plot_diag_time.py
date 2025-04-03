import numpy as np
import pandas as pd
import os
from argparse import ArgumentParser, BooleanOptionalAction

import egfr_microsim.model.egfr_model as egfr_model
import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.analysis.helpers.plots_read_agg as plots_read_agg
from tqdm import tqdm

import matplotlib.pyplot as plt
import seaborn as sns

def calculate_earliest_diag(result_dir, trajectories, cutoffs = None):
    """
    Calculate the earliest time of diagnosis 
    and eGFR value at the earliest time of diagnosis under eGFR21-eGFR09
    """
    if cutoffs == None:
        cutoffs = [90, 60, 45, 30, 15]
        
    print("Calculating age and eGFR values at earliest \
            possible diagnoses into stages G2-G5")
    
    all_diff = []
    cohorts = trajectories.reset_index().cohort.unique()
    for i in tqdm(cohorts):
        bl_at_cutoff_09 = egfr_model.add_age_at_cutoff(
            trajectories.query("cohort==@i"), cutoffs, interv_under_eq="09"
        )
        bl_at_cutoff_21 = egfr_model.add_age_at_cutoff(
            trajectories.query("cohort==@i"), cutoffs, interv_under_eq="21"
        )
        diff = (
            pd.merge(
                bl_at_cutoff_09.filter(
                    ["egfr", "sex", "race", "stage_diag"]
                ).reset_index(),
                bl_at_cutoff_21.filter(
                    ["egfr",  "sex", "race", "stage_diag"]
                ).reset_index(),
                how="inner", 
                on = ["pid",  "sex", "race", "stage_diag"], 
                suffixes = ("_09", "_21")
            )
            .assign(cohort = i)
        )
        all_diff.append(diff)
    
    stage_diff = (pd.concat(all_diff)
                  .assign(
                      age_diff = lambda x: x.age_21-x.age_09, 
                      egfr_diff = lambda x: x.egfr_21-x.egfr_09
                  ) 
                 ).set_index([ "stage_diag", "sex", "race"])
    
    
    stage_diff.to_feather(
        os.path.join(result_dir, "data", "stage_diff_scenarios.feather")
    )
    return stage_diff

def save_diff_agg_summary(stage_diff, exp_result_dir):
    """
    Save summary statistics corresponding to the difference in 
    earliest time of diagnosis  and eGFR value at the earliest time 
    of diagnosis under eGFR21-eGFR09
    """
    diff_agg = (stage_diff
     .groupby([ "stage_diag",  "sex", "race"])
     .agg({
         "age_diff": ["mean", np.std, "min", "max"], 
         "egfr_diff": ["mean", np.std,  "min", "max"]
     })
    )

    diff_agg.to_csv(
        os.path.join(
            exp_result_dir, 
            "tables", 
            "stage_diff_scenarios_agg.csv"
        )
    )
    
    diff_agg_age = diff_agg.age_diff.reset_index()
    # needed for manuscript
    diff_agg_age.to_csv(
        os.path.join(
            exp_result_dir, 
            "tables", 
            "stage_diff_scenarios_agg_age.csv"
        )
    )

    diff_agg_egfr = diff_agg.egfr_diff.reset_index()
    # needed for manuscript
    diff_agg_egfr.to_csv(
        os.path.join(
            exp_result_dir, 
            "tables", 
            "stage_diff_scenarios_agg_egfr.csv"
        )
    )

def get_stage_diff(args, result_dir, exp_result_dir):
    """
    Calculate the difference in earliest time of diagnosis 
    and eGFR value at the earliest time of diagnosis under eGFR21-eGFR09
    """
    if args.recalculate:   
        trajectories_bl = plots_read_agg.read_trajectory(
            result_dir, 
            traj_path = args.trajectories_path
        )
        stage_diff = calculate_earliest_diag(result_dir, trajectories_bl)
        stage_diff.to_feather(
            os.path.join(
                result_dir, 
                "data", 
                "stage_diff.feather"
            )
        )
        save_diff_agg_summary(stage_diff, exp_result_dir)
    else:
        stage_diff = pd.read_feather(
            os.path.join(
                result_dir, 
                "data", 
                "stage_diff.feather"
            )
        )
    return stage_diff
        
def get_plot_df(stage_diff):
    """
    Format the table to plot
    """
    stage_diff_toplot = (stage_diff
                         .filter(["age_diff", "egfr_diff"])
                         .reset_index()
                         .assign(
                             stage_diag = lambda x: x.stage_diag.astype("category")
                         )
                        )
    stage_diff_toplot = plots_read_agg.update_var(
        stage_diff_toplot, 
        var_name = 'race'
    )
    stage_diff_toplot = plots_read_agg.update_var(
        stage_diff_toplot, 
        var_name = 'sex', 
        update_dict = {"M": "Male", "F": "Female"}
    )    
    stage_diff_toplot = plots_read_agg.update_var(
        stage_diff_toplot, 
        var_name = "stage_diag"
    )
    
    stage_diff_toplot_stacked = (
        pd.DataFrame(
            stage_diff_toplot
            .set_index(["stage_diag", "sex", "race"])
            .stack()
        )
        .reset_index()
        .rename(columns = {"level_3": "diff_type", 0: "diff_value"})
    )
    return stage_diff_toplot_stacked

def gen_save_plot(to_plot, exp_result_dir, filename):
    """
    Save bar plot of the difference in earliest time of diagnosis 
    and eGFR value at the earliest time of diagnosis under eGFR21-eGFR09
    """
    g = sns.catplot(
        data = to_plot, 
        y = "diff_value", 
        x = "stage_diag", 
        hue = "race",
        kind = "box", 
        col = "diff_type", 
        sharey = False, 
        showfliers = False
    )
    
    ax1, ax2 = g.axes[0]
    
    ax1.axhline(0, ls = '--', color = "grey")
    ax1.set_ylim(-17, 17)
    ax1.set_title("")
    ax1.set_ylabel('Time difference (years)', fontsize=12, position=(5, 0.5))
    ax1.grid(axis = "y")
    
    ax2.axhline(0, ls='--',  color="grey")
    ax2.set_ylim(-10, 10)
    ax2.set_title("")
    ax2.set_ylabel("\neGFR difference", fontsize=12, )
    ax2.grid(axis = "y")
    
    g.set_xlabels("CKD stage", fontsize=12)
    plt.tight_layout()
    plt.legend(title = "Race",facecolor='white', frameon=True)
    g.legend.remove()
    
    g.savefig(
        os.path.join(exp_result_dir, "plots", filename + ".png"),
        dpi = 800
    )

def generate(result_dir, exp_result_dir, args):
    """
    Generate and save bar plot of the difference in earliest time of diagnosis 
    and eGFR value at the earliest time of diagnosis under eGFR21-eGFR09
    """

    stage_diff = get_stage_diff(args, result_dir, exp_result_dir)

    stage_diff_toplot_stacked = get_plot_df(stage_diff)

    gen_save_plot(
        stage_diff_toplot_stacked, 
        exp_result_dir, 
        args.output_filename
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
    args = parser.parse_args()

    exp_dir, result_dir = path_helpers.get_dirs(
        args, ["results"]
    )
    exp_result_dir = path_helpers.get_exp_result_dir(args)
    generate(result_dir, exp_result_dir, args)



