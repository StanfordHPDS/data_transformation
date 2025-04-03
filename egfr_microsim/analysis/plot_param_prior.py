import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from argparse import ArgumentParser
import os

import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.helpers.param_helpers as param_helpers

"""
A dictionary mapping unique parameter ids to their descriptions
for plotting
"""
name_dict = {
    0: "Healthy",
    1: "CKD stage 3a+",
    2: "Hypertension",
    3: "Hypertension,\n CKD stage 3a+",
    4: "Diabetes",
    5: "Diabetes,\n CKD stage 3a+",
    6: "Diabetes,\n Hypertension",
    7: "Diabetes,\n Hypertension,\n CKD stage 3a+"
}

def get_prior(params, original = True):
    """
    Read the prior (mean) file, format
    """
    prior = (
        param_helpers.get_slopes(params, original = original)
        .rename(columns = {"mu": "mean_prior"})
        .assign(param_id = lambda x: x.idx.map(name_dict))
        .filter(["param_id", "mean_prior"])
    )
    return prior

def get_posterior(exp_result_dir, loglik_dir = None):
    """
    Read the posterior file, format
    """

    try:
        post_file = pd.read_csv(
            os.path.join(
                exp_result_dir, 
                "tables", 
                "sum_lik_param_posterior.csv"
            )
        )
    except:
        post_file = pd.read_feather(
            os.path.join(
                loglik_dir, 
                "sum_lik_param_posterior.feather"
            )
        )
        post_file.to_csv(
            os.path.join(
                exp_result_dir, 
                "tables", 
                "sum_lik_param_posterior.csv"
            )
        )
    posterior = (
        post_file
        .query("param=='posterior'")
        .rename(columns = {"mu": "mean_posterior", "min": "min_posterior", 
                          "max": "max_posterior"})
        .reset_index()
        .reset_index()
        .assign(param_id = lambda x: x.index.map(name_dict))
        .filter(["mean_posterior", "min_posterior", "max_posterior", "param_id"])
    )
    return posterior

def get_prior_sample_df(param_set_df):
    """
    Format the frame containing sampled parameters,
    update names corresponding to parameters
    """

    to_plot = (
        param_set_df
        .stack()
        .reset_index()
    
    )
    to_plot.columns = ["param_set_id", "param_id", "value"]
    to_plot = (
        to_plot
        .assign(
            param_id = lambda x:  x.param_id.astype(int).map(name_dict)
        )
    )
    return to_plot

def get_handles():
    """
    Formatting handles for plot legends
    """
    black_diam = mlines.Line2D(
        [], [], color='black', marker='d',
        linestyle='None', markersize=10, label='Prior mean'
    )
    red_diam = mlines.Line2D(
        [], [], color='red', marker='d', 
        linestyle='None', markersize=10, label='Posterior mean'
    )
    
    return [black_diam, red_diam]

def gen_save_plot(to_plot, exp_result_dir, args):
    """
    Generate and save the parameter distribution plot
    """

    max_y = round(to_plot.value.max(),0)
    
    g = sns.catplot(
        data=to_plot, y="value", x="param_id",
        kind="box", showfliers = False,
        height=5, aspect=1.7,
        width=0.6
    )
    g.map(
        sns.scatterplot,
        'param_id',
        'mean_posterior', 
        color='red', 
        marker = "d", 
        s = 100)
    g.map(
        sns.scatterplot, 
        'param_id', 
        'mean_prior', 
        color='black',
        marker = "d", 
        s = 100)
    
    for patch in g.axes.flat[0].patches:
        fc = patch.get_facecolor()
        patch.set_facecolor(plt.matplotlib.colors.to_rgba(fc, alpha=0.3))
    
    g.set_xticklabels(rotation = 45, size = 12)
    g.set_yticklabels(size = 11)
    g.set_axis_labels("Parameter", "Value", fontsize=12)
    plt.grid(axis = 'y')
    plt.ylim(0, max_y)
    plt.legend(
        handles = get_handles(), 
        loc = "upper left", 
        labelspacing = 1,
        fontsize = 12
    )

    g.savefig(
        os.path.join(exp_result_dir, "plots", args.output_filename + ".png"),
        dpi = 800
    )

def generate(params_dir, exp_result_dir, args, params, loglik_dir):
    """
    Generate a plot of the sample from prior, Hoerger et al prior mean
    and calibrated posterior and save as a table and a figure
    """
    try:
        param_set_df = param_helpers.get_param_df_from_1k(
            params_dir, 
            args.parameter_set_number
        )
    except:
        param_set_df = param_helpers.get_param_df(params_dir)

    prior_sample = get_prior_sample_df(param_set_df)
    posterior = get_posterior(exp_result_dir, loglik_dir)
    prior_mean = get_prior(params)

    plot_df = (prior_sample
                .merge(posterior, how="left", on="param_id")
                .merge(prior_mean, how="left", on="param_id")
               )
    
    gen_save_plot(plot_df, exp_result_dir, args)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest = "experiment_string",
        help = "unique experiment name",
    )
    parser.add_argument(
        "-r", 
        "--parameter_set_number", 
        dest="parameter_set_number", 
        help = "number of parameter sets in experiment",
        type=int
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
        default = "param_distrib",
    )
    args = parser.parse_args()

    exp_dir, params_dir, loglik_dir = path_helpers.get_dirs(
        args, ["param_samples", "loglik"]
    )
    exp_result_dir = path_helpers.get_exp_result_dir(args)

    params_base = path_helpers.get_params_base(args, allow_except=True)
    generate(params_dir, exp_result_dir, args, params_base, loglik_dir)

