import pandas as pd
import os

import egfr_microsim.model.helpers.path_helpers as path_helpers
import egfr_microsim.model.helpers.param_helpers as param_helpers

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

update_dicts = {
    'sex': {'F': 'Female', 'M': 'Male'},
    'male': {0: 'Female', 1: 'Male'},
    'race': {'B': 'Black', 'NB': 'Non-Black'},
    'stage_diag': {'G1': '1', 'G2': '2',
                  'G3a': '3a', 'G3b': '3b',
                  'G4': '4', 'G5': '5'}
}

def update_var(df, col_name='race', update_dict=None):
    """
    Update values in a specified column using a dictionary
    """
    if update_dict == None:
        update_dict = update_dicts[col_name]
    try:
        return df.assign(
            **{col_name: lambda x: x[col_name].map(update_dict)}
        )
    except:
        new_var = df.index.get_level_values(col_name).map(update_dict)
        df = (
            df
            .droplevel(level=col_name)
            .assign(**{col_name: new_var})
            .set_index(col_name, append=True)
        )
        return df

def get_initial_counts(args, cohorts_dir):
    """
    Calculate the number of people in each sex and age group at the beginning of the simulation
    (to use as denominator in survival calculations)
    """
    group_sizes = []
    for cohort_id in range(args.cohort_number):
        df_cohort = param_helpers.get_cohort(cohorts_dir, cohort_id).reset_index()
        numbers = df_cohort.groupby(["sex", "race"]).nunique().assign(cohort = cohort_id).filter(["pid", "cohort"])
        group_sizes.append(numbers)

    group_sizes_df = (pd.concat(group_sizes)
                .set_index("cohort", append=True)
                .rename(columns={"pid": "group_size"})
    )
    group_sizes_df = update_var(group_sizes_df, var_name='sex')
    group_sizes_df = update_var(group_sizes_df, var_name='race')
    initial_counts = group_sizes_df.reset_index().pivot(index = ["sex", "race"], columns="cohort",values="group_size")
    return initial_counts

def read_trajectory(result_dir, traj_path = None, mode = "reference"):
    """
    Read reference/counterfactual trajectory files
    """
    filename = "trajectories_bl.feather" if mode == "reference" else "trajectories_cf.feather"
    if traj_path is None:
        traj_path = os.path.join(
            result_dir, 
            "data", 
            filename
        )
    else:
        traj_path = os.path.join(
            path_helpers.repo_path, 
            traj_path
        )

    trajectories_bl = pd.read_feather(traj_path)
    return trajectories_bl

