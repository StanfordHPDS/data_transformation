import os
from argparse import ArgumentParser
from egfr_microsim.model.helpers.path_helpers import repo_path

"""
For a specified experiment, identifies jobs not yet completed and generates
a file containing all the remaining commands necessary to run.
Assumes that calibration_batch is ran with --parameter_set_number = 1000.
"""
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-e",
        "--experiment_string",
        dest = "experiment_string",
        help = "unique experiment name"
    )
    parser.add_argument(
        "-o", 
        "--output_filename", 
        dest = "output_filename",
        help = "Name of text file to save"
    )
    parser.add_argument(
        "-r",
        "--param_number", 
        dest = "R", 
        type = int,
        help = "Number of parameter sets in the experiment"
    )
    parser.add_argument(
        "--cohort_number", 
        dest = "M", 
        type = int, 
        default = 100,
        help = "Number of cohorts in the experiment"
    )
    parser.add_argument(
        "-b",
        "--begin_from_param_n", 
        dest = "b", 
        type = int, 
        default = 0,
        help = "Only consider parameter files indexed above this value"
    )
    args = parser.parse_args()
    exp = args.experiment_string

    agg_path0 = "outputs/%s/sum_stats/" % exp
    agg_path = os.path.join(repo_path, agg_path0)
    filename = "df_agg_all_%i_%i.feather"
    command = "python egfr_microsim/calibration/calibration_batch.py --experiment_string %s \
        --cohort_id %i --parameter_set_id %i --parameter_set_number 1000 --interv_under_eq $equation \
        --cohort_size $N --params_base_path $param_file"
    
    files = []
    commands = []
    for cohort_id in range(args.M):
        for parameter_set_id in range(int(args.b/1000), int(args.R/1000)): 
            path0 = filename % ((parameter_set_id, cohort_id))
            if not os.path.exists(os.path.join(agg_path, path0)):
                files.append(os.path.join(agg_path, path0))
                commands.append(command % ((exp, cohort_id, parameter_set_id)))

    print(len(files))
    with open(os.path.join(agg_path, args.output_filename + ".txt"), "w+") as f:
        f.write("\n".join(commands))

