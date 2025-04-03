import os
import pandas as pd
from sklearn.linear_model import LinearRegression
from argparse import ArgumentParser
import seaborn.objects as so

import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.helpers.param_helpers as param_helpers
import egfr_microsim.calibration.aggregate as aggregate
import egfr_microsim.calibration.coverage_analysis as coverage_analysis

def get_coef_target(reg, y):
    """
    Format coefficient frame for plotting
    """
    coef_target = pd.DataFrame(data = reg.coef_, index = y.columns)
    coef_target_vertical = pd.DataFrame(coef_target.stack(future_stack=True))
    coef_target_vertical.index.names = coef_target_vertical.index.names[:-1]+["param_id"]
    coef_target_vertical = (
        coef_target_vertical
        .rename(columns={0:"coeff"})
        .reset_index()
         .assign(
             param_id = lambda x: x.param_id.astype(dtype="category")
        )
    )
    return coef_target_vertical

def save_coefficients_plot(p, target_name, coverage_dir):
    """
    Save generated plot
    """
    p.save(
        os.path.join(
            coverage_dir, "plots", "_".join((target_name, "coeffs.png"))
        ),
        bbox_inches="tight",
    )

def save_coefficients_table(df, target_name, coverage_dir):
    """
    Save generated table
    """
    path = os.path.join(coverage_dir, "tables", "_".join((target_name, "coeffs.csv")))
    df.to_feather(path)


def plot_coefficients(coef_target_vertical, col_var):
    """
    Generate coefficients plots
    """
    p = (
        so.Plot(
            data = coef_target_vertical.reset_index(),
            x = "age",
            y = "coeff",
            color = "param_id"
        )
        .facet(col = col_var, row = "stage")
        .share(y = "row")
        .layout(size = (7, 20))
        .add(so.Dot())
        .label(
            x = "Age",
            y = "Coefficient",
            color = "Parameter id",
            col = "".join((col_var.capitalize(), ":")),
            row = "Stage:",
        )
    )
    return p

def generate_plots(args, calib_targets):
    """
    Run parameter regressions and generate plots
    """
    plots = []
    for k in range(len(args.calib_targets)):
        target = calib_targets[k]
        
        target_name = args.calib_targets[k].split(".")[0]
        calib_strata = coverage_analysis.get_calibration_strata(target)
        col_var = [el for el in calib_strata if el not in ["age", "stage"]][0]
        
        calib_frac_frame_path = os.path.join(
                    coverage_dir, "tables", "_".join((target_name, "frac.feather"))
                )
        calib_frac_frame = pd.read_feather(calib_frac_frame_path)
        target, calib_frac_frame = coverage_analysis.update_strata_stages(calib_frac_frame, target, calib_strata)

        X = param_set_df_norm
        y = calib_frac_frame.T
        
        reg = LinearRegression().fit(X, y)

        coef_target_vertical = get_coef_target(reg, y)
        p = plot_coefficients(coef_target_vertical, col_var)
        save_coefficients_plot(p, target_name, coverage_dir)
        save_coefficients_table(coef_target_vertical, target_name, coverage_dir)
        plots.append(p)
    return plots
    

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest = "experiment_string",
        help = "unique experiment name"
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
        "--cohort_number", 
        dest = "M", 
        type = int, 
        default = 100,
        help = "Number of cohorts in the experiment"
    )
    parser.add_argument(
        "--parameter_set_number", 
        dest = "parameter_set_number", 
        help = "Number of parameter sets in a file", 
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

    args = parser.parse_args()

    # 1. get parameters
    repo_path = path_helpers.repo_path
    exp_dir, params_dir, sum_stats_dir, loglik_dir, coverage_dir = path_helpers.get_dirs(
        args, ["param_samples", "sum_stats", "loglik", "coverage_analysis"]
    )
    params_base = path_helpers.get_params_base(args)
    params_base = param_helpers.update_params(params_base)

    df_slopes = param_helpers.get_slopes(params_base)

    min_age = params_base['cohort']['init_age']  

    q_counts_cats_all = aggregate.combine(sum_stats_dir, int(args.parameter_set_number/1000)).query("age>@min_age")
    calib_targets = path_helpers.calib_targets(args.calib_targets, params_base)
    calib_targets = [el.query("age>@min_age") for el in calib_targets]
    try:
        param_set_df = param_helpers.get_param_df_from_1k(params_dir, args.parameter_set_number)
    except:
        param_set_df = param_helpers.get_param_df(params_dir)

    param_set_df_norm = (param_set_df - param_set_df.mean()) / param_set_df.std()

    plots = generate_plots(args, calib_targets)