from argparse import ArgumentParser
import numpy as np
import pandas as pd
import os
import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.analysis.helpers.plots_read_agg as plots_read_agg

import matplotlib.pyplot as plt
import seaborn as sns
import seaborn.objects as so

def get_both_death(trajectories_bl, trajectories_cf):
    """
    For two trajectory frames, identify death times and merge them
    so that each person has a single row and both death times and eGFR values
    """

    death_times_bl = pd.DataFrame(
        trajectories_bl
        .query("death==1")
        .reset_index(level="age")
        .set_index(["male", "race"], append=True)
    )
    death_times_cf = pd.DataFrame(
        trajectories_cf
        .query("death==1")
        .reset_index(level="age")
        .set_index(["male", "race"], append=True)
    )
    
    both_death = pd.merge(
        death_times_bl.filter(["age"]).reset_index(),
        death_times_cf.filter(["age"]).reset_index(),
        how = "inner",
        on = ["pid", "cohort", "male", "race"],
        suffixes = ("_09", "_21")
    )
    
    return both_death

def get_death_frames(result_dir):
    """
    Identify death times for individuals under both scenarios,
    separately including all individuals and those who developed CKD stage 3a+
    """

    trajectories_bl = plots_read_agg.read_trajectory(
        result_dir, mode = "reference"
    )
    trajectories_cf = plots_read_agg.read_trajectory(
        result_dir, mode = "counterfactual"
    )
    both_death = get_both_death(trajectories_bl, trajectories_cf)
    
    # people who reached G3a before death
    reached_g3_bl = trajectories_bl.query(
        "death==1 & g3==1"
    ).reset_index("age").index
    
    both_death = get_both_death(trajectories_bl, trajectories_cf)
    both_death_CKD = both_death.set_index(["pid", "cohort"]).loc[reached_g3_bl]
    return both_death, both_death_CKD

def percentile(n):
    """
    From Stack Overflow:
    https://stackoverflow.com/questions/19894939/calculate-arbitrary-percentile-on-pandas-groupby
    
    Applies the np.percentile function to pandas.groupby().agg()
    """
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'perc_%s' % n
    return percentile_

def plot_diff_le(LEs_df_diff, row_var, color_var=None, 
                 ranges=False, top_text_loc = 0.6,
                 bottom_text_loc = 0.5
                ):
    """
    Generate the final difference in life expectancy plot
    for the manuscript
    """
    if ranges:
        y, ymin, ymax = "mean", "perc_5", "perc_95"
    else:
        y, ymin, ymax = "diff", None, None
    fig = plt.Figure(figsize=(8,5))
    p = (
        so.Plot(
            data=LEs_df_diff,
            x="age",
            y=y,
            ymin=ymin,
            ymax=ymax,
            color=color_var,
        )
        .add(so.Line())
        .add(so.Band())
        .label(y="Additional life expectancy (years)", x="Age")
        .theme(plt.style.library["seaborn-v0_8-whitegrid"])
        .facet(col="sex", row = row_var)
        .scale(color={"Black": sns.color_palette("deep")[0], "Non-Black":  sns.color_palette("deep")[1]})
        .scale(
            x=so.Continuous().tick(every=10),
        )
        .on(fig)
        .plot()
    )
    legend = fig.legends.pop(0)
    fig.text(0.15, top_text_loc, "Black", fontsize=12, color=sns.color_palette("deep")[0],  fontweight="bold")
    fig.text(0.15, bottom_text_loc, "Non-Black", fontsize=12, color=sns.color_palette("deep")[1], fontweight="bold")
    fig.text(0.58, top_text_loc, "Black", fontsize=12, color=sns.color_palette("deep")[0], fontweight="bold")
    fig.text(0.58, bottom_text_loc, "Non-Black", fontsize=12, color=sns.color_palette("deep")[1], fontweight="bold")
    
    return p

def read_agg(name, results_dir):
    """
    Read experiment summary (aggregate) file
    """
    df = pd.read_feather(os.path.join(results_dir, "data", name+".feather"))
    return df

def summarize_le_diff(death_df, CI_type = "normal_CI"):
    """
    Generate and save summary tables describing differences in life expectancy
    and eGFR values based on the death data frame
    """
    
    if CI_type == "percentile":
        agg_dict = [np.mean, 
                  "min", 
                  "max", 
                  percentile(5), 
                  percentile(95)]
        cols = ["sex", "race", "mean", "min", "max", "perc_5", "perc_95", "age"]
    
    if CI_type == "normal_CI":
        agg_dict = [
             np.mean, 
              "min", 
              "max", 
              np.std
         ]
        cols = ["sex", "race", "mean", "min", "max", "std", "age"]

    all_age_summs = []

    for age in range(30,101,5):
        age_summ = (death_df
         .query("age_09 >= @age")
         .assign(surv_diff = lambda x: x.age_21-x.age_09)
         .groupby(["male", "race"])
         .agg(
             {"surv_diff" : agg_dict}
         ).reset_index()
         .assign(age = age)
        )
        age_summ.columns = cols
        all_age_summs.append(age_summ)
    all_age_summs1 = pd.concat(all_age_summs).reset_index(drop=True)


    if CI_type =="normal_CI":
        all_age_summs1 = all_age_summs1.assign(
            perc_5 = lambda x: x["mean"] - 1.96* x["std"],
            perc_95 = lambda x: x["mean"] + 1.96* x["std"],
            male = lambda x: x.sex ).drop(columns="sex") 
    all_age_summs1 = plots_read_agg.update_var(all_age_summs1, "race")
    all_age_summs1 = plots_read_agg.update_var(all_age_summs1, var_name='male').rename(columns={"male": "sex"})
    return all_age_summs1

def plot_and_save(df, filename, exp_result_dir,
                  top_text_loc = 0.55,
                  bottom_text_loc = 0.45
                 ):
    """
    Generate and save plots
    """

    p_diff = plot_diff_le(
        df.reset_index(drop=True),
        row_var = None, 
        color_var = "race", 
        ranges = True,
        top_text_loc = top_text_loc,
        bottom_text_loc = bottom_text_loc
    )
    p_diff.save(
        os.path.join(exp_result_dir, "plots", filename + ".png"),
        dpi = 800
    )
    
    df.to_csv(
        os.path.join(exp_result_dir, "tables", filename + ".csv")
    )
    return p_diff

def generate(result_dir, exp_result_dir):
    """
    Generate plots of the difference in life expectancy across age comparing 
    the two scenarios (reference and counterfactual), stratified by race, sex
    and CKD status
    """

    both_death, both_death_CKD = get_death_frames(result_dir)

    all_age_summs_CKD = summarize_le_diff(both_death_CKD, CI_type = "normal_CI")
    all_age_summs_all = summarize_le_diff(both_death, CI_type = "normal_CI")

    p_diff_all = plot_and_save(
        all_age_summs_all, 
        "life_exp_all", 
        exp_result_dir
    )
    p_diff_CKD = plot_and_save(
        all_age_summs_CKD, 
        "life_exp_CKD", 
        exp_result_dir
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
        "-t",
        "--trajectories_path",
        dest = "trajectories_path",
        help = "path to trajectories file (relative to repository top), in feather format",
        default = None
    )
    args = parser.parse_args()

    exp_dir, sum_stats_dir, result_dir, cohorts_dir = path_helpers.get_dirs(
        args, ["sum_stats", 'results', "cohort_samples"]
    )
    
    exp_result_dir = path_helpers.get_exp_result_dir(args)
    
    generate(result_dir, exp_result_dir)
