import pandas as pd
import itertools
import egfr_microsim.model.egfr_formulas as egfr_formulas

"""
suppressing the following FutureWarning:

FutureWarning: the convert_dtype parameter is deprecated and will 
be removed in a future version. Do ``ser.astype(object).apply()`` instead 
if you want ``convert_dtype=False``

(Warning appears despite following instructions in text above)
"""

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)

def get_agg_all_index():
    """
    Generate an index for the aggregate experiment summary
    """

    ages = list(range(25, 100, 10))
    stages = ["G1", "G2", "G3a", "G3b", "G4", "G5"]
    dm = [0, 1]
    ht = [0, 1]
    male = [0, 1]

    age_stage_index = pd.MultiIndex.from_tuples(
        itertools.product(ages, dm, ht, male, stages)
    ).set_names(["age", "dm", "ht", "male", "stage"])

    return age_stage_index

def egfr_at_age(trajectories_df, at_age, strata = None):
    """
    For each individual in the trajectory frame, identify an eGFR value
    at the specified age
    """

    if strata == None:
        strata=["dm", "ht", "male", "race"]
    
    # Identify the last entry in the trajectory table for an individual that precedes
    # an age of interest, and use slope at that entry to calculate eGFR at age of interest
    df_age = (
        trajectories_df.query("age < @at_age")
        .groupby("pid")
        .last()
        .assign(
            egfr=lambda x: (x.egfr - (at_age - x.age) * x.slope).round(2),
            stage = lambda x: (
                x.egfr
                .astype(float)
                .apply(egfr_formulas.egfr_to_stages)
            ),
            age = at_age,
        )
        .filter(strata + ["egfr", "age", "stage", "death"])
    )

    # for initial age entries
    if df_age.shape[0] == 0:
        df_age = (trajectories_df
                  .set_index("pid")
                  .query("age == @at_age")
                  .assign(
                      stage = lambda x: (
                          x.egfr
                          .astype(float)
                          .apply(
                              egfr_formulas.egfr_to_stages, 
                              convert_dtype = False
                              )
                          )
                      )
                  .filter(strata + ["egfr", "age", "stage", "death"])
                 )
    
    return df_age

def get_counts_strata(
    trajectory_df, ages=None, strata=["age", "dm", "ht", "male"]
):
    """
    Counting the number of people within each stratum and stage, by age
    """
    
    if ages == None:
        ages = range(30, 101, 5)
        
    age_dfs = []
    ## Get eGFR values at specific ages for each individual, match stages
    for age in ages:
        age_df = egfr_at_age(trajectory_df, age)
        age_dfs.append(age_df)

    stage_counts = (
        pd.concat(age_dfs)
        .query("stage != 'D'")
        .groupby(strata + ["stage"])
        .count()
        .rename(columns={"egfr": "counts"})
        .filter(["counts"])
    )
    return stage_counts