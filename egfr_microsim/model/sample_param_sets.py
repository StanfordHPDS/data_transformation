import os
import numpy as np
import pandas as pd
from functools import reduce
from argparse import ArgumentParser, BooleanOptionalAction
from scipy.stats import qmc, truncnorm

import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.helpers.param_helpers as param_helpers

def regular_lhs(n_variables, n_param_sets):
    """
    Function:
        Latin Hypercube Sampling (LHS), generating a random sample that covers
        the entire multidimensional parameter space
    Args:
        n_variables: dimension of the parameter space (L)
        n_param_sets: number of parameter sets to sample (R*multiplier)
    Returns:
        An L by (R*multiplier) array in [0,1).
    """
    sampler = qmc.LatinHypercube(
        d = n_variables, optimization="random-cd"
        )
    sample = sampler.random(n = n_param_sets)
    return sample

def truncated_normal(df_slopes, lhs_sample = None, parameter_set_number = None):
    """
    Function: Generating a truncated normal sample (using a pre-sampled LHS, if provided)
    Args:
        df_slopes: a frame with L rows containing the following columns
            mu: mean value of a parameter
            stddev: standard deviation of a parameter
            min: minimum value of a parameter
            max: maximum value of a parameter
        lhs_sample: (optional) an L by (R*multiplier) array containing an LHS sample that covers 
            the parameter space
        parameter_set_number: number of parameter sets
    Returns:
        an (R*multiplier) by L (if lhs_sample provided) OR R by L array (otherwise)
        containing sampled parameter sets 
    """
    n_params = df_slopes.shape[0]

    try:
        l_trunc = df_slopes["min"]
        r_trunc = df_slopes["max"].fillna(np.inf)
    except:    
        l_trunc = np.zeros(n_params) 
        r_trunc = np.ones(n_params) * (np.inf)

    a = (l_trunc - df_slopes.mu) / df_slopes.stddev
    b = (r_trunc - df_slopes.mu) / df_slopes.stddev

    if lhs_sample is not None:
        df_param = np.empty(lhs_sample.shape)
    else:
        df_param = np.empty((parameter_set_number, n_params))
        random_uniform_sample = np.random.uniform(0, 1, (parameter_set_number, n_params))
        
    for j in range(df_slopes.shape[0]):
        x = df_slopes.iloc[j]
        truncated_normal = truncnorm(a[j], b[j], loc=x.mu, scale=x.stddev)
        if lhs_sample is not None:
            df_param[:, j] = truncated_normal.ppf(lhs_sample.loc[:,j])
        else:
            df_param[:, j] = truncated_normal.ppf(random_uniform_sample[:,j])

    df_param = pd.DataFrame(df_param)
    return df_param

def apply_constraint(df_param, constraint):
    """
    Function:
        Rejects parameters in df_param which don't fulfill a constraint
    Args:
        df_param: a frame of parameter sets
        constraint: in the form (target_idx, others_idx), where:
            target_idx: int, parameter index value
            others_idx: list or int, indices of parameters that must be lower
                        than target_idx value
            in order to enforce the rule that more comorbidities -> higher slope
    Returns:
        a frame of parameter sets with rows violating the constraint removed
    """
    target_idx, others_idx = constraint
    target = df_param.iloc[:, target_idx]
    others = df_param.iloc[:, others_idx]
    filter = target.ge(others.max(axis=1))
    return df_param[filter]

def rejection_sampling(df_param, args, params, constraints=None):
    """
    Function:
        Filters out parameter sets which don't fulfill constraints and then
        sub-samples remaining parameter sets to get a total of args.parameter_set_number (R)
    Args:
        df_param: a frame of parameter sets
        constraints: a list in the form [(target_idx, others_idx), ...], where:
            target_idx: int, parameter index value
            others_idx: list or int, indices of parameters that must be lower
                        than target_idx value
            in order to enforce the rule that more comorbidities -> higher slope
    Returns:
        an L by R frame of parameter sets that satisfy constraints
    """
    if constraints is None:
        constraints = param_helpers.get_slopes_order_list(params)

    df_param_constrained = reduce(
        apply_constraint,
        constraints,
        df_param[df_param.ge(df_param.iloc[:, 0], axis=0)].dropna(axis=0),
    )

    try:
        subset_index = np.random.choice(
            df_param_constrained.index, args.parameter_set_number, replace=False
        )
    except:
        if args.allow_replacement:
            print("Not enough paramters left after applying constraints; parameters sampled with replacement in rejection_sampling")
            subset_index = np.random.choice(
                df_param_constrained.index, args.parameter_set_number, replace=True
            )
        else:
            raise ValueError("Not enough paramters left after applying constraints; \
                 pass --allow_replacement or calculate an LHS sample with a larger multiplier")
    
    df_param_constrained_subset = df_param_constrained.loc[
        subset_index
    ].reset_index(drop=True)

    return df_param_constrained_subset

def param_sample_exists(args, param_dir):
    """
    Function:
        Check if parameter sets have been sampled before
    Args:
        args: command line arguments
        param_dir: the location of parameter sets
    Returns:
        a boolean indicating whether param_dir containts parameter samples
    """
    if args.index is not None:
        return os.path.exists(os.path.join(param_dir, "params%i.feather" % args.index))
    else:
        return os.path.exists(os.path.join(param_dir, "params.feather"))

def sample_all_parameter_sets(
        df_slopes, args, params, presampled_lhs_file = None, param_dir = None
        ):
    """
    Function:
        Sample R parameter sets from a truncated normal distribution 
    Args:
        df_slopes: frame containing covariate combinations, and prior values
            defined by mean (mu), standard deviation, as well as 
            min and max values of a parameter (optional)
        args: command line arguments
        params: parameters dictionary, containing the relative order of slopes
        presampled_lhs_file: if specified, LHS will not be recalculated
        param_dir: the location of parameter sets
    Returns:
        an L by R frame of parameter sets
    """
    if args.rejection_sampling:
        sampled_params_number = args.parameter_set_number * args.multiplier
    else:
        sampled_params_number = args.parameter_set_number

    if presampled_lhs_file is not None:
        lhs_sample = pd.read_feather(os.path.join(param_dir, presampled_lhs_file))
        lhs_sample = lhs_sample.loc[(args.index*1000) : (args.index*1000 + sampled_params_number-1)]
    else:
        lhs_sample = regular_lhs(
            n_variables = df_slopes.shape[0],
            n_param_sets = sampled_params_number)

    df_param = truncated_normal(
        df_slopes,
        lhs_sample = pd.DataFrame(lhs_sample),
    )
    
    if args.rejection_sampling:
        df_param = rejection_sampling(df_param, args, params)

    return df_param

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
        type=int
    )
    parser.add_argument(
        "-L",
        "--parameter_set_size",
        dest = "parameter_set_size", 
        type = int, 
        default = 8
    )
    parser.add_argument(
        "--presampled_lhs_file", 
        dest="presampled_lhs_file", 
        type = str, 
        default = None
    )
    parser.add_argument(
        "-i",
        "--index", 
        dest = "index", 
        type = int, 
        default = None
    )
    parser.add_argument(
        "--params_base_path",
        dest = "params_base_path",
        help = "location of parameter dictionary",
        default = "parameter_files/params_base.json",
    )
    parser.add_argument(
        "-r",
        "--rejection_sampling",
        dest = "rejection_sampling",
        help = "apply rejection sampling using provided parameter value ordering",
        default = True,
    )
    parser.add_argument(
        "--overwrite_existing",
        dest="overwrite_existing",
        help = "overwrite existing sampled parameter sets",
        action=BooleanOptionalAction
    )
    parser.add_argument(
        "--allow_replacement",
        dest="allow_replacement",
        help = "allow sampling with replacement if the sample is too small \
            after applying rejection sampling",
        action=BooleanOptionalAction
    )
    parser.add_argument(
        "--multiplier",
        dest="multiplier",
        help = "upsampling multiplier applied to parameter_set_number \
            - 16 needed for rejection sampling",
        type=int,
        default=16,
    )
    args = parser.parse_args()
    repo_path = path_helpers.repo_path
    exp_dir, param_dir = path_helpers.get_dirs(args, ["param_samples"])
    params_base = path_helpers.get_params_base(args, allow_except=True)

    if not args.overwrite_existing:
        assert not param_sample_exists(args, param_dir), "Parameters already sampled - pass --overwrite_existing to overwrite"

    df_slopes = param_helpers.get_slopes(params_base)
    assert df_slopes.shape[0] == args.parameter_set_size

    df_param = sample_all_parameter_sets(
        df_slopes, 
        args, 
        params_base, 
        presampled_lhs_file = args.presampled_lhs_file,
        param_dir = param_dir
        )

    if args.index is not None:
        pd.DataFrame(df_param).to_feather(
            os.path.join(param_dir, "params%i.feather" % args.index)
        )
    else:
        pd.DataFrame(df_param).to_feather(os.path.join(param_dir, "params.feather"))
