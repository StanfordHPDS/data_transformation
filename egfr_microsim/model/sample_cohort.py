import os
from argparse import ArgumentParser

import egfr_microsim.model.helpers.sampling as sampling
import egfr_microsim.model.helpers.path_helpers as path_helpers

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest = "experiment_string",
        help = "unique experiment name",
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
        "-M",
        "--cohort_number", 
        dest = "cohort_number", 
        help = "number of cohorts to sample", 
        type = int
    )
    parser.add_argument(
        "--params_base_path",
        dest = "params_base_path",
        help = "location of parameter dictionary",
        default = "parameter_files/params_base.json",
    )
    args = parser.parse_args()

    exp_dir, cohorts_dir = path_helpers.get_dirs(args, ["cohort_samples"])
    params_base = path_helpers.get_params_base(args, allow_except=True)

    for j in range(args.cohort_number):
        df_init = sampling.sample_cohort(params_base)
        df_init.to_csv(os.path.join(cohorts_dir, "".join(("cohort", str(j), ".csv"))))

