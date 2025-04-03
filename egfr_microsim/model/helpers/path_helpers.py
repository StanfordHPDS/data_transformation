from pathlib import Path
import os
import git
import json
import pandas as pd

repo_path = os.path.join(
    git.Repo(".", search_parent_directories=True).working_tree_dir, "egfr_microsim"
)

def get_dirs(args, dirs = []):
    """
    Generates requested standard paths for an experiment defined in args.
    Create directories if they have not been created before.
    Args:
        args: command line arguments which will be used to update the parameters dictionary
        dirs: a list of directory names to create 
    Returns:
        a list of directory paths
    """
    exp_dir = os.path.join(repo_path, "outputs", args.experiment_string)
    Path(exp_dir).mkdir(parents=True, exist_ok=True)

    allowed_dirs = ["cohort_samples", "param_samples", "sum_stats", 
                    "coverage_analysis", "loglik", "trajectories", "results"]
    
    assert set(dirs).issubset(allowed_dirs), "dirs must be one or more of %s, but is %s" % (str(allowed_dirs), str(dirs))

    dir_paths = []
    dir_paths.append(exp_dir)
    for dir_name in dirs:
        dir = os.path.join(exp_dir, dir_name)
        dir_paths.append(dir)
        Path(dir).mkdir(parents = True, exist_ok = True)
        if dir_name in ["coverage_analysis", "analysis", "results"]:
            Path(os.path.join(dir, "plots")).mkdir(parents=True, exist_ok=True)
            Path(os.path.join(dir, "tables")).mkdir(parents=True, exist_ok=True)
        if dir_name in ["results"]:
            Path(os.path.join(dir, "data")).mkdir(parents=True, exist_ok=True)
    if len(dir_paths) == 1:
        dir_paths = dir_paths[0]
    return dir_paths

def get_exp_result_dir(args):
    """
    Generates requested standard paths for experiment results defined in args.
    Create directories if they have not been created before.
    Args:
        args: command line arguments which will be used to update the parameters dictionary
    Returns:
        a directory path
    """
    results_path = os.path.join(repo_path, "exp_results", args.experiment_string)
    Path(results_path).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(results_path, "plots")).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(results_path, "tables")).mkdir(parents=True, exist_ok=True)
    return results_path
    
def from_repo(path, rel_path = None):
    """
    Generates a complete path using the path relative to the git repo and the path
    to the git repo

    Args:
        path (str): path relative to the git repo
        (optional) rel_path (str): if provided, path to the git repo
    Returns:
        a complete path
    """
    if rel_path:
        return os.path.join(rel_path, path)
    else:
        return os.path.join(repo_path, path)


def get_params_base(args, rel_path = None, allow_except = False):
    """
    Generates a baseline parameters dictionary based on the baseline parameters file,
    expands relative paths, updates dictionary with values provided in args.

    Args:
        args: command line arguments which will be used to update the parameters dictionary
        rel_path: if provided, use this instead of the path of the current git repository to
            expand file paths (useful when function is called from outside of the package)
        allow_except: allow missing n or interv_under_eq parameters (allowed only
            when sampling cohorts or parameter sets)
    Returns:
        an updated dictionary
    """
    params_base = json.load(open(from_repo(args.params_base_path), "r"))

    # expand all relative paths
    for key, value in params_base["model"]["file_paths"].items():
        params_base["model"]["file_paths"][key] = from_repo(value, rel_path)
        
    # not required when sampling parameters
    if not allow_except:
        params_base["cohort"]["n"] = int(args.cohort_size)
        params_base["model"]["interv_under_eq"] = args.interv_under_eq

    else:
        try:
            params_base["cohort"]["n"] = int(args.cohort_size)
        except:
            "cohort_size not provided"
        # not required when sampling cohorts
        try: 
            params_base["model"]["interv_under_eq"] = args.interv_under_eq
        except:
            "interv_under_eq not provided"
    return params_base

def rename_calib_targets(calib_targets):
    """
    Re-formats calibration targets (updates column names, age values, 
    sex/male column)
    Args:
        calib_targets: a list<pandas.DataFrame> of calibration target frames
    Returns:
        a list<pandas.DataFrame> of re-fromatted calibration target frames
    """
    map_cols = {
        "g_stage": "stage",
        "age_group": "age",
        "diabetes": "dm",
        "hypertension": "ht",
        "n": "x",
        "group_sum": "N"
    }

    for i in range(len(calib_targets)):
        calib_targets[i] = (
            calib_targets[i]
            .rename(columns=map_cols)
            .assign(age=lambda x: x.age + 5)
        )

        if "sex" in calib_targets[i].columns:
            calib_targets[i] = (
                calib_targets[i]
                .rename(columns={"sex": "male"})
                .assign(male=lambda x: x.male.map({"F": 0, "M": 1}))
            )
    return calib_targets

def calib_targets(filenames, params):
    """
    Reads in and re-formats calibration targets.
    Args:
        filenames (list<str>): a list of calibration target filenames (without
            relative paths)
        params (dict): parameters dictionary
    Returns:
        a list<pandas.DataFrame> of calibration target frames
    """
    targets = []
    for el in filenames:
        targets.append(
            pd.read_csv(
                os.path.join(
                    repo_path,
                    params["model"]["file_paths"]["calib_targets_path"],
                    el,
                ),
                index_col = 0
            )
        )
    targets_renamed = rename_calib_targets(targets)
    return targets_renamed

def agg_counts_path(sum_stats_dir, cohort_id, filename="df_agg_all_"):
    """
    Identifies the path to an aggregate counts file corresponding to
    cohort_id
    Args:
        sum_stats_dir (str): summary statistics path
        cohort_id (int): cohort index  
        filename (str): name of aggregate count files
    Returns:
        (str) path string
    """
    counts_agg_filename = os.path.join(
        sum_stats_dir, "".join((filename, str(cohort_id), ".feather"))
    )
    return counts_agg_filename

def read_all_agg_counts(sum_stats_dir, cohort_number, filename="df_agg_all_"):
    """
    Reads in all aggregate count frames for a given experiment
    (one per cohort)
    Args:
        sum_stats_dir (str): summary statistics path
        cohort_number (int): number of cohorts in the experiment 
        filename (str): name of aggregate count files
    Returns:
        list<pandas.DataFrame> of all aggregate count frames for a given experiment
    """
    all_count_files = [
        pd.read_feather(
            agg_counts_path(
                sum_stats_dir, cohort_id, filename=filename
            )
        )
        for cohort_id in range(cohort_number)
    ]
    return all_count_files