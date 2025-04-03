import os
import pandas as pd
import numpy as np
import scipy as sp

import egfr_microsim.calibration.coverage_analysis as coverage_analysis
import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.helpers.param_helpers as param_helpers
import egfr_microsim.analysis.plot_cover as plot_cover
import egfr_microsim.calibration.aggregate as aggregate

# allow division by 0, return NaN
np.seterr(divide='ignore')

def calc_multinomial_loss(x):
    """
    Calculate multinomial log loss for a frame x
    Args:
        x: a pandas.DataFrame, containing two columns
            x_counts: multinomial counts
            p_experiments: multinomial parameters
    Returns:
        Log loss per frame
    """
    return sp.stats.multinomial.logpmf(
        x = x.x.values,
        p = x.p.values,
        n = x.x.values.sum()
)


def get_ith_param_loglik(loglik_frame, i, calib_strata):
    """
    Calculate log likelihoods for all calibration strata for
    parameter set i (within a single calibration target defined by calib_strata).DS_Store
    The function reshapes the loglik_frame to correspond to the ith parameter set only,
    and applies calc_multinomial_loss to each stratum.

    Args:
        loglik_frame: a data frame containing x_counts and <parameter_set_number>
            number of columns corresponding to multinomial distribution frequencies
            associated with each parameter-specific experiment 
        i: indexing parameters from <parameter_set_number>
        calib_strata: list of column names defining the calibration stratum
    Returns:
        a pandas.DataFrame containing log likelihoods for all calibration strata for
        parameter set i 
    """
    p_col_name = "_".join(("p", str(i)))

    # select the p_experiments column corresponding to parameter i
    loglik_frame_i = (loglik_frame
                    .loc[:, [p_col_name, "x"]]
                    .rename(columns = {p_col_name: "p"})
                    # remove rows with p==0 to prevent inf likelihood
                    .query("p != 0")
    )

    loglik_i = (loglik_frame_i
              .groupby(calib_strata)
              .apply(calc_multinomial_loss)
    )
    return loglik_i

def get_target_strata(columns):
    """
    Given a target frame X_target, return the columns by which 
    it has been stratified

    Args:
        columns (list): list of column names
    Returns
        a list of column names
    """
    return [el for el in columns if el not in ["stage", "n_ppl"]]

def get_p_per_experiment(q_counts_cats_all, calib_strata):
    """
    Calculates distribution frequencies corresponding to counts
    across calibration strata for an experiment summary table q_counts_cats_all
    Args:
        q_counts_cats_all: counts across calibration strata for an experiment summary table
        calib_strata: a list of indices defining the calibration strata
    Returns:
        a table of distribution frequencies
    """    
    calib_strata_stage = calib_strata + ["stage"]

    # q_counts_cats_all has more granular categories (defined in index) than 
    # needed for a given calibration target - requires summation
    q_counts_cats = q_counts_cats_all.groupby(calib_strata_stage).sum()
    
    # calculate counts within each strata (across stages)
    M_counts_cats = (q_counts_cats
                     .groupby(calib_strata)
                     .sum()
                    )
    
    # calculate frequencies per calibration stratum
    p_counts_cats = q_counts_cats.div(M_counts_cats)

    # update column names to p_i for each i in <parameter_set_number>
    p_counts_cats.columns = ['_'.join(('p', str(i))) for i in p_counts_cats.columns]

    return p_counts_cats

def get_X_target(x_counts, calib_strata):
    """
    Format the x_counts table: set index to calib_strata_stage, 
    update column name

    Args:
        x_counts: table containing the calibration target
        calib_strata: a list of indices defining the calibration strata
    Returns:
        A formatted table
    """
    X_target = (x_counts
                .set_index(calib_strata)
                .rename(columns={"n_ppl": "x_counts"})
    )
    return X_target


def single_calibration_target_loss(single_target_counts, q_counts_cats_all, parameter_set_number):
    """
    Calculate loss for a single calibration target

    Args:
        single_target_counts: table containing the calibration target (x)
        q_counts_cats_all: table containing experiment counts corresponding
            to the calibration target (q)
        parameter_set_number: number of parameter sets in experiment
    """
    calib_strata = coverage_analysis.get_calibration_strata(single_target_counts)
    calib_strata_nostage = [el for el in calib_strata if el != "stage"]

    # update column names to be consistent between targets and experimental results
    single_target_counts, q_counts_cats_all = coverage_analysis.update_strata(
        q_counts_cats_all, single_target_counts, calib_strata
        )
    
    P_exper = get_p_per_experiment(q_counts_cats_all, calib_strata_nostage) 
    X_target = get_X_target(single_target_counts, calib_strata) 
    
    multinomial_frame = P_exper.join(X_target, on = calib_strata)

    log_liks_k = []

    # for each parameter set, calculate log likelihood
    for i in range(parameter_set_number):
        loglik_i = get_ith_param_loglik(multinomial_frame, i, calib_strata_nostage)
        log_liks_k.append(loglik_i.sum())
    return log_liks_k

def get_joint_frame(prior, posterior, save = True, 
                    target_name = None, loglik_dir = None):
    """
    Displays prior and posterior in a single table
    """
    joint_frame = pd.concat((prior, posterior))

    if save:
        assert (target_name is not None) and (loglik_dir is not None), \
             "target_name and loglik_dir required to save"
        filename = os.path.join(loglik_dir, "_".join((target_name, "posterior")))
        joint_frame.to_feather(".".join((filename, "feather")))
    return joint_frame


def calc_importance_sample_weights(param_set_loss):
    """
    Calculate importance weights per each parameter set based on the 
    log loss associated with a single parameter set

    Args:
        param_set_loss: a frame with a single log loss value per parameter set

    Returns:
        importance weights
    """  
    weights = sp.special.softmax(param_set_loss, axis = None)

    return pd.Series(weights, index=param_set_loss.index)

def get_weighted_sample(weights, resample_size = 1e4):
    """
    For Sample Importance Resampling: using an indexed array of weights,
    generate a weighted sample and return a list of indices indicating 
    selected samples. The size of returned list is resample_size.
    
    Args:
        weights: a vector of values and weights to use for weighted sampling
        resample_size: number of values to resample
    Returns:
        a weighted sample of size resample_size
    """
    param_ids = weights.index

    samples = np.random.choice(
        param_ids, 
        size = int(resample_size), 
        replace = True, 
        p = weights
    )
    return samples

def calc_SIR(samples, param_set_df, params):
    """
    Given a sample of parameters, calculate the mean and 95% CI 
    of the resulting distribution

    Args:
        samples: a list of indices pointing to param_set_df entries
        param_set_df: all parameter sets considered in the experiment 
            (sampled from prior)
        params: a dictionary of experiment parameters
    Returns:
        mean and 95% CI for each parameter
    """
    df_slopes_index = param_helpers.get_slopes_index(params)
    resampled_params = param_set_df.to_numpy()[samples]
    posteriors = pd.DataFrame(
        index = df_slopes_index,
        data = {
            "mu": resampled_params.mean(axis=0),
            "min": np.quantile(resampled_params, q=0.025, axis=0),
            "max": np.quantile(resampled_params, q=0.975, axis=0)
            }
        ).assign(param = "posterior")
    return posteriors

def get_prior(params):
    """
    Calculate the mean, standard deviation and 95% CI of the sampled
    prior

    Args:
        params: dictionary of experiment parameters

    """
    df = param_helpers.get_slopes(params)
    df_slopes_index = param_helpers.get_slopes_index(params)
    prior = (
        df
        .set_index(df_slopes_index.names)
        .filter(["mu", "stddev"])
        .assign(
            min=lambda x: x.mu - 1.96 * x.stddev, 
            max=lambda x: x.mu + 1.96 * x.stddev
        )
        .drop(columns="stddev")
        .assign(param="prior")
    )
    return prior

def save_log_liks(log_liks, loglik_dir, args, k):
    """
    Save the provided log likelihood table

    Args:
        log_liks: a frame of log likelihoods
        loglik_dir: a directory path to save log likelihood tables in
        args: command line arguments
        k: the index corresponding to the calibration target corresponding 
            to the log likelihood table

    Returns:
        None
    """
    (
        pd.DataFrame(log_liks)
        .to_feather(
            os.path.join(loglik_dir, args.calib_targets[k].split(".")[0] + ".feather")
        )
    )


def calculate_log_likelihoods_all(calib_targets, q_counts_cats_all, args, loglik_dir):
    """
    Calculates log likelihoods across calibration targets, saves in 
    a single table

    Args:
        calib_targets: a list of calibration target tables
        q_counts_cats_all: counts across calibration strata for an 
            experiment summary table
        args: command line arguments
        loglik_dir: a directory path to save log likelihood tables in

    Returns:
        A table containing all calcylated log likelihood loss values
    """
    log_liks = []
    for k in range(len(calib_targets)):
        target = calib_targets[k]
        calib_strata = coverage_analysis.get_calibration_strata(target)
        calib_strata_addl = [el for el in calib_strata if el not in ["age", "stage"]][0]

        # for SDVI-based targets, calculate loss separately for each stratum, 
        # considering the entire cohort to correspond to the stratum
        if calib_strata_addl in ["ICE_q3", "SDI_q3"]:
            # remove 0th strata - corresponds to null values
            target = target.query("%s != 0" % calib_strata_addl)
            log_liks_sdvi = []
            for i in range(1,4):
                target_i = target.query("%s == @i" % calib_strata_addl).drop(columns = [calib_strata_addl])
                log_liks_ki = single_calibration_target_loss(
                                target_i,
                                q_counts_cats_all,
                                args.parameter_set_number
                            )
                log_liks_sdvi.append(log_liks_ki)
            log_liks_k = np.array(log_liks_sdvi).sum(axis=0)
        else:
            log_liks_k = single_calibration_target_loss(
                                target,
                                q_counts_cats_all,
                                args.parameter_set_number
                            )
        save_log_liks(log_liks_k, loglik_dir, args, k)
        log_liks.append(log_liks_k)
        
    all_target_loss = pd.DataFrame(log_liks).T
    all_target_loss.index.name = "param_set"
    all_target_loss.columns.name = "calib_target"   
    return all_target_loss

def update_resampled_plot_strata(calib_strata, resampled_counts_cats, target):
    """
    Update names of SDVI strata before plotting

    Args:
        calib_strata: list of column names defining the calibration stratum
        resampled_counts_cats: counts across calibration strata for an 
            experiment summary table corresponding to the posterior (after resampling) 
        target: a calibration target table

    Returns:
        Counts and target tables with updated strata names
    """
    calib_strata_addl = [el for el in calib_strata if el not in ["age", "stage"]][0]

    # for SDVI-based targets, calculate loss separately for each stratum, 
    # considering the entire cohort to correspond to the stratum
    if calib_strata_addl in ["ICE_q3", "SDI_q3"]:    
        sdvi_strata = []
        for i in range(1,4):
            sdvi_strata.append(
                resampled_counts_cats
                .assign(**{calib_strata_addl: i})
                .set_index(calib_strata_addl, append=True)
            )
        resampled_counts_cats = pd.concat(sdvi_strata)
        target = target.query("%s != 0" % calib_strata_addl)
    return resampled_counts_cats, target

def save_posterior_coverage_summary(calib_strata, target, resampled_counts_cats,
                                    args, target_name, save_dir = None):
    """
    Save plots and tables corresponding to posterior coverage

    Args:
        calib_strata: list of column names defining the calibration stratum
        target: a calibration target table
        resampled_counts_cats: counts across calibration strata for an 
            experiment summary table corresponding to the posterior (after resampling) 
        args: command line arguments
        target_name (str): name of the calibration target 
        save_dir: directory in which to save summary

    Returns:
        a seaborn.objects plot
    """

    calib_frac_frame = coverage_analysis.get_target_fractions(
        target, resampled_counts_cats, calib_strata, args
    )

    target, calib_frac_frame = coverage_analysis.update_strata(calib_frac_frame, target, calib_strata)
    calib_strata = coverage_analysis.get_calibration_strata(target)
    
    df_coverage = coverage_analysis.get_coverage(target, calib_frac_frame, calib_strata)

    df_coverage = (
        df_coverage.reset_index(level=-1)
        .assign(
            param=lambda x: x.param.map(
                {"target": "target", "sampled": "posterior"}
            )
        )
        .set_index("param", append=True)
    )
    p = plot_cover.plot_coverage(
        df_coverage, x_var="age", col_var=df_coverage.index.names[-3]
    )

    if save_dir is not None:
        df_coverage.to_feather(
            os.path.join(
                save_dir, "tables", "_".join((target_name, "posterior.feather"))
            )
        )
        p.save(
            os.path.join(
                save_dir, "plots", "_".join((target_name, "posterior.png"))
            ),
            bbox_inches="tight",
        )

    return p

def calculate_likelihoods(params, args, sum_stats_dir, loglik_dir):
    """
    A wrapper function around calculate_log_likelihoods_all

    Args:
        params: a dictionary of experiment parameters
        args: command line arguments
        sum_stats_dir: path to the directory containing experiment summary statistics
        loglik_dir: a directory path to save log likelihood tables in

    Returns:
        Calculated loss for all targets
    """
    calib_targets = path_helpers.calib_targets(args.calib_targets, params)

    min_age = params['cohort']['init_age']
    q_counts_cats_all = aggregate.combine(
        sum_stats_dir, 
        int(args.parameter_set_number / 1000)
        ).query("age > @min_age")

    all_target_loss = calculate_log_likelihoods_all(
        calib_targets, 
        q_counts_cats_all, 
        args, 
        loglik_dir
        )
    
    return all_target_loss

def read_likelihoods(params, args, loglik_dir):
    """
    Read files with previously calculated log likelihoods, return
    as a single frame

    Args:
        params: dictionary of experiment parameters
        args: command line arguments
        loglik_dir: a directory path containing log likelihood tables

    Returns:
        A frame containing all log likelihoods for the experiment
    """
    calib_targets = path_helpers.calib_targets(args.calib_targets, params)
    log_liks = []
    for k in range(len(calib_targets)):
        log_liks_k = pd.read_feather(
            os.path.join(loglik_dir, args.calib_targets[k].split(".")[0] + ".feather")
            )
        log_liks.append(log_liks_k.values.reshape(-1))
    
    all_target_loss = pd.DataFrame(log_liks).T
    all_target_loss.index.name = "param_set"
    all_target_loss.columns.name = "calib_target"
    
    return all_target_loss

def calculate_SIR(all_target_loss_mean, param_set_df, params, args):
    """
    Sample importance resampling: calculate weights for an experiment
    based on log likelihoods, resample from the prior sample using weights
    to obtain the posterior and return a summary of the posterior 
    (mean and normal 95% CI)

    Args:
        all_target_loss_mean: 
        param_set_df: all parameter sets considered in the experiment 
            (sampled from prior)
        params: a dictionary of experiment parameters
        args: command line arguments

    Returns:
        posterior summary and ids of sampled parameters
    """
    SIR_weights = calc_importance_sample_weights(all_target_loss_mean)
    sampled_ids = get_weighted_sample(SIR_weights, resample_size = args.resample_size)
    posterior_distr = calc_SIR(
        sampled_ids, 
        param_set_df.dropna(axis=1), 
        params
        ).assign(param = "posterior")
    print((SIR_weights > 0).sum(), "number of non-zero weights")
    return posterior_distr, sampled_ids

def posterior_coverage_plots(q_counts_cats_all, sampled_ids, calib_targets, 
                             args, coverage_dir = None):
    """
    Generate summary tables and plots of posterior coverage,
    comparing the distribution across calibration targets
    based on the posterior with the calibration targets themselves

    Args:
        q_counts_cats_all: counts across calibration strata for an 
            experiment summary table
        sampled_ids: prior parameter sample ids selected during
            sample importance resampling
        calib_targets: a list of calibration target tables
        args: command line arguments
        coverage_dir: A directory to save coverage plots in 

    Returns:
        A list of plots
    """
    resampled_counts_cats = q_counts_cats_all.iloc[:, sampled_ids].copy()
    resampled_counts_cats.columns = list(range(resampled_counts_cats.columns.shape[0]))

    posterior_plots = []

    # calculate coverage of the posterior across calibration targets
    for k in range(len(calib_targets)):
        target = calib_targets[k]
        calib_strata = coverage_analysis.get_calibration_strata(target)
        resampled_counts_cats1, target = update_resampled_plot_strata(calib_strata, resampled_counts_cats, target)

        p = save_posterior_coverage_summary(
        
            calib_strata = calib_strata, 
            target = target,
            resampled_counts_cats = resampled_counts_cats1, 
            args=args,
            target_name = args.calib_targets[k].split(".")[0],
            save_dir=coverage_dir
            )
        posterior_plots.append(p)
    return posterior_plots