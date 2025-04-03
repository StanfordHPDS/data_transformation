import os
import pandas as pd
from argparse import ArgumentParser

import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.sample_param_sets as sample_param_sets

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
        "--parameter_set_size", 
        dest="parameter_set_size", 
        help = "number of unique parameters in each parameter set",
        type = int, 
        default = 8
    )
    parser.add_argument(
        "--multiplier",
        dest="multiplier",
        help = "upsampling multiplier applied to parameter_set_number - 16 needed for rejection sampling",
        type=int,
        default=1,
    )
    
    args = parser.parse_args()
    repo_path = path_helpers.repo_path
    exp_dir, param_dir = path_helpers.get_dirs(args, ["param_samples"])
    sampled_params_number = args.parameter_set_number * args.multiplier

    lhs_sample = sample_param_sets.regular_lhs(
        n_variables = args.parameter_set_size,
        n_param_sets = sampled_params_number
    )
    if (sampled_params_number % 1000) == 0:
        filename = "lhs_sample_{0}k.feather".format(int(sampled_params_number / 1000))
    else:
        filename = "lhs_sample_{0}.feather".format(sampled_params_number)
    save_filename = os.path.join(param_dir, filename)
    print(save_filename)
    pd.DataFrame(lhs_sample).to_feather(save_filename)
