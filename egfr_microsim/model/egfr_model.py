import os
import time

import numpy as np
import pandas as pd
from copy import deepcopy

import egfr_microsim.model.helpers.survival_helpers as survival_helpers
import egfr_microsim.model.egfr_formulas as egfr_formulas
import egfr_microsim.model.helpers.param_helpers as param_helpers
import egfr_microsim.model.interventions as interventions

"""
Dictionary of types for each variable
"""
dtype_dict = {
    "ht": int,
    "dm": int,
    "death": int,
    "g3": int,
    "diag": int,
    "nephr": int,
    "sex": str,
    "race": str,
    "slope": float,
    "egfr": float,
}

def agg_max_stage(x):
    """ 
    Function:
        This custom aggregate function handles the mixed-type stage_diag column.
        Columns where a stage is reached are of type string, and the rest are 0s.
    Args:
        x: stage_diag column
    Returns:
        a single column value
    """
    if len(set(x.values)) == 1:
        ret = x.values[0]
    else:
        ret = max([el for el in x.fillna(0) if el != 0])
    return ret
    
"""
Used to combine duplicate entries - when two values exist for the same 
individual at a given age, defines which value should be kept
"""
agg_dictionary_one_age = {
    'egfr': 'mean',
    'ht': 'max', 
    'dm': 'max', 
    'death': 'max', 
    'g3': 'max', 
    'diag': 'max',
    'nephr': 'max', 
    # those two won't differ within a single person
    'sex': 'first', 
    'race': 'first',
    # take the last slope, since it will consider all updates
    'slope': 'last',
    # since we're combining at_cutoff (first) with frame (second), 
    # at_cutoff has the most up-to-date stage
    'stage_diag': 'first'}
    
def integrate_duplicate_rows(df, keep_var = "egfr", ignore_stages=False):
    """
    Integrate rows in df with duplicate pid and keep_var values
    """
    assert (keep_var in ["age", "egfr"]), "keep_var must be one of age, egfr"

    if keep_var == "age":
        agg_var = "egfr"
    elif keep_var == "egfr":
        agg_var = "age"

    agg_dictionary = {
        agg_var: 'mean',
        'ht': 'max', 
        'dm': 'max', 
        'death': 'max', 
        'g3': 'max', 
        'diag': 'max',
        'nephr': 'max', 
        'sex': 'first', 
        'race': 'first',
        'slope': 'last'
    }

    if df.reset_index().set_index(['pid','egfr']).index.duplicated().sum() == 0:
        return df

    if 'stage_diag' in df.columns and not ignore_stages:
        agg_dictionary['stage_diag'] = agg_max_stage
    
    integrated_df = (
        df
        .reset_index()
        .groupby(["pid", keep_var])
        .agg(agg_dictionary)
        .assign(**{agg_var: lambda x: x[agg_var].round(2)})
    )

    if keep_var == "egfr":
        integrated_df = (
            integrated_df
            .reset_index()
            .set_index(["pid", "age"])
            .sort_index()
        )

    return integrated_df

def if_print_timing(print_timing, prev_time, i):
    """
    Function:
        a helper function tor timing portions of the code
    Args:
        print_timing: boolean indicator for generating a print statement
        prev_time: float, indicating previous time to compare to current time
        i: int counter
    Returns:
        An increased counter and current time
    """
    this_time = time.time()
    if print_timing:
        try:
            print("part %i: %.2fs" % (i, this_time - prev_time))
        except:
            pass
    return i+1, this_time

def get_init_egfr(df):
    """
    Function:
        For each person (pid), identifies the initial eGFR value 
        (at initiation, which is the youngest age)
    Args:
        df: a trajectory pandas DataFrame indexed by pid and age
    Returns:
        DataFrame with updated eGFR values, indexed by pid
    """
    return (df
            .reset_index(level="age")
            .assign(init_egfr=lambda x: x.groupby(level="pid")["egfr"].first().ffill())
            )

def get_decline(df):
    """
    Function:
        Updates eGFR values based on the slope and intervention indicators.
    Args:
        df: a trajectory pandas DataFrame indexed by pid and age
    Returns:
        DataFrame with updated eGFR values, indexed by pid and age
    """
    assert df.index.names == ['pid', 'age']

    return (
        get_init_egfr(df)
        .assign(
            dt=lambda x: x.groupby(level="pid")["age"].diff().shift(-1),
            degfr=lambda x: (x.dt * x.slope).round(2),
            degfr_cum=lambda x: x.groupby(level="pid")["degfr"].cumsum().shift(1).fillna(0),
            egfr=lambda x: x.init_egfr - x.degfr_cum,
        )
        .drop(columns=["dt", "degfr", "degfr_cum", "init_egfr"])
        .set_index("age", append=True)
    )

def entries_at_cutoffs(df, cutoffs):
    """
    Function:
        Identifies elements of df including eGFR values corresponding 
        to cutoff values
    Args:
        df: a trajectory pandas DataFrame, indexed by pid and age
        cutoffs: an integer or list of integers corresponding to eGFR values
    Returns:
        DataFrame with selected eGFR values, indexed by pid and age    
    """
    if type(cutoffs) != type([]):
        cutoffs = [cutoffs]
    return df.query("egfr.round(2) in @cutoffs")

def identify_age_at_cutoff(df, cutoff):
    """
    Function:
        For a specified eGFR value, identifies individuals in the
        trajectory df who reach the value before death (or age 101),
        and calculates an age at which the value is reached.
        It inherits all indicators and slopes from the last recorded 
        observation preceding cutoff.
    Args:
        df: a trajectory pandas DataFrame, indexed by pid and age
        cutoff: an integer corresponding to an eGFR value
    Returns:
        a frame with one or zero entries per individual 
    """

    if 'pid' or 'age' in df.index.names:
        df = df.reset_index()

    # identify people who already have an entry exactly at cutoff 
    # and don't re-calculate their values
    ignore_indices = entries_at_cutoffs(df, cutoff).pid.values.tolist()

    ret_df = (
        df.query("egfr > @cutoff & pid not in @ignore_indices")
        .groupby("pid")
        .last()
        .query("death != 1")
        .assign(
            age=lambda x: (
                x.age + (x.egfr - cutoff) / (x.slope)
            ).round(2),
            egfr=cutoff,
        )
        .reset_index()
        .query("age <= 101")
    )

    # in case ret_df ends up finding duplicate entries by "pid", "age"
    # (possible under rounding errors), use corresponding entry from df
    # updating egfr to cutoff (old egfr would have been +- 0.02 off)
    df =  df.set_index(["pid", "age"])
    ret_df = ret_df.set_index(["pid", "age"])
    overlap = df.index.intersection(ret_df.index)
    if overlap.shape[0] > 0:
        ret_df.loc[overlap] = df.loc[overlap].assign(egfr=cutoff)

    return ret_df.reset_index()

def identify_age_at_cutoff_many(df, cutoff_list=None):
    """
    Function:
        Runs identify_age_at_cutoff for each cutoff in cutoff_list and returns 
        a resulting, concatenated frame
    Args:
        df: a trajectory pandas DataFrame, indexed by pid and age
        cutoff: a list of integers corresponding to eGFR values
    Returns:
        a frame with one or zero entries per individual per cutoff
    """
    if cutoff_list is None:
        cutoff_list = [60, 45, 30, 15]
    at_cutoffs = []
    for cutoff in cutoff_list:
        at_cutoffs.append(
            identify_age_at_cutoff(df, cutoff=cutoff)
            )

    calculated_at_cutoff = pd.concat(at_cutoffs)
    return calculated_at_cutoff

def update_slopes(df, params):
    """
    Function:
        Deterministically updates the slope variable in each row to the
        corresponding provided slope values conditional on the set of covariates
        (dm, ht, g3), and then adjusts them based on intervention indicators
        (reducing the speed of decline if an intervention is applied)
    Args:
        df: a trajectory pandas DataFrame, indexed by pid and age
        params: parameters dictionary, containing the slopes corresponding to
            binary indicators (params["model"]["slopes"]) and multiplicative
            slope reductions related to interventions 
            (params["model"]["interv_params"])
    Returns:
        a trajectory frame with updated slope values
    """

    slope_df = params["model"]["slopes"]

    slopes_from_table = (
        df
        .drop(columns=["slope"], errors="ignore")
        .merge(
            slope_df.filter(["dm", "ht", "g3", "slope"]),
            how="left",
            on=["dm", "ht", "g3"],
        )
        .set_index(df.index)
    )
    
    slopes_add_interv = (
            slopes_from_table.assign(
                slope = lambda y: y.apply(
                    lambda x: interventions.get_slope_reduction(
                        x.slope, x.nephr, x.diag, params["model"]["interv_params"]
                    ),
                    axis = 1
                )
            )
        )

    return slopes_add_interv

def update_frame_with_cutoff_event(
    df, at_cutoff, indicator_name, params, multiple_indic=False
):
    """
    Function:
        Updates a trajectory frame with events that occur at a cutoff, by merging
        them in and updating indicators and progression slopes following the events.
        Drops stage_diag.
    Args:
        df: a trajectory pandas DataFrame, indexed by pid and age
        at_cutoff: events identified at a cutoff value
        params: parameters dictionary, containing slope values conditional on indicators
            and reductions in progression speed related to interventions
        multiple_indic: a boolean indicating if at_cutoff containts
            updates multiple columns in df (e.g. two interventions: diag and nephr)
    Returns:
        an updated trajectory frame
    """
    full_at_cutoff = (
        pd.concat([at_cutoff, df])
        .sort_index()
        # if two entries at a single age, aggregate using the agg_dictionary_one_age dictionary
        .groupby(level = ["pid", "age"])
        .agg(agg_dictionary_one_age)
        .fillna(0)
    )

    if multiple_indic:
        assert(indicator_name == ["diag", "nephr"]), "only implemented for indicator_name == ['diag', 'nephr'])"
        to_update = []
        full_at_cutoff = full_at_cutoff.assign(
            # propagate the indicator for the rest of each person's trajectory
            diag = lambda x: x.groupby(level='pid')["diag"].cummax(),
            nephr = lambda x: x.groupby(level='pid')["nephr"].cummax()
        )
        to_update = full_at_cutoff.query("(diag == 1) or (nephr == 1)")

    else:
        full_at_cutoff = full_at_cutoff.assign(
            # propagate the indicator for the rest of each person's trajectory
            **{indicator_name: lambda x: x.groupby(level='pid')[indicator_name].cummax()}
        )
        to_update = full_at_cutoff.query("{0}==1".format(indicator_name))

    if to_update.shape[0] == 0:
        return full_at_cutoff 

    df_ev_subset = update_slopes(to_update, params)
    updated_subset = get_decline(df_ev_subset)
    full_at_cutoff.update(updated_subset)
    full_at_cutoff = full_at_cutoff.drop(columns="stage_diag")
    full_at_cutoff = integrate_duplicate_rows(full_at_cutoff, keep_var = "egfr")

    return full_at_cutoff

def add_age_at_cutoff(df, cutoff_list, interv_under_eq):
    """
    Function:
        Generates a row per person corresponding to the time they would be considered
        to have reached a specific eGFR cutoff, under the specified eGFR equation (09 or 21)
    Args:
        df: a trajectory dataframe
        cutoff_list: an integer or list of integers corresponding to eGFR values
        interv_under_eq: "09" or "21, determines which equation is used to establish whether a threshold
        has been reached
    Returns:
        a frame with one or zero entries per individual per cutoff 
        (zero if death would occur before reaching cutoff)
    """
    if cutoff_list == None:
        cutoff_list = [60, 30, 15]

    # if only a single cutoff is passed
    if type(cutoff_list) != type([]):
        cutoff_list = [cutoff_list]

    if interv_under_eq == "09":

        df_as_09 = egfr_formulas.translate_trajectory(
            df.reset_index(), change_type = "21to09"
        )

        # run this function again, on 2009 trajectory
        # interv_under_eq="21" is passed to indicate that no translation
        # of the trajectory is required
        at_cutoff_as_09 = add_age_at_cutoff(
            df_as_09, cutoff_list, interv_under_eq = "21"
        )

        # take the resulting trajectory and update it back to 2021
        if at_cutoff_as_09.shape[0] > 0:
            at_cutoff = egfr_formulas.update_egfr_equation(
                at_cutoff_as_09.reset_index(), change_type = "09to21"
            ).set_index(["pid", "age"])

        else:
            at_cutoff = at_cutoff_as_09

    # eGFR21 already used by model - no need to translate trajectory
    elif interv_under_eq == "21":
        
        calculated_at_cutoff = (
            identify_age_at_cutoff_many(df, cutoff_list)
            .assign(stage_diag = lambda x: x.egfr.map(egfr_formulas.egfr_to_g_stages))
        )

        already_at_cutoff = (
                entries_at_cutoffs(
                    df.reset_index(), cutoff_list
                )
                .assign(
                    stage_diag = lambda x: x.egfr.round(2).map(egfr_formulas.egfr_to_g_stages)
            )
        )
        
        at_cutoff = (
            pd.concat([calculated_at_cutoff, already_at_cutoff])
            .set_index(["pid", "age"])
        )

    return at_cutoff

def combine_with_death(df_events, df_death):
    """
    Function:
        Takes in the final events data frame before death is sampled, and merges with the
        frame containing each person's time of death, filtering out previously sampled 
        events proceeding death time.
    Args:
        df_events: a trajectory dataframe with all events (g3, interventions) already identified
        df_death: a frame with each person's sampled time of death
    Returns:
        A frame with final trajectories from intiation age to death
    """

    df_pre_death = (
        df_events
        .merge(df_death.filter(["pid", "death_age"]), on="pid")
        .query("age < death_age")
        .drop(columns=["death_age"])
        .set_index(["pid", "age"])
        .sort_index()
    )

    df_with_death = (
        pd.concat(
            [
                df_pre_death,
                (
                    df_death.rename(columns={"death_age": "age"}).set_index(
                        ["pid", "age"]
                    )
                ),
            ]
        )
        .sort_index()
        .fillna({"death": 0})
    )
    return df_with_death

def add_event_at_cutoff(df, cutoff, event_name, params):
    """
    Function:
        Identify individuals who will reach an event triggered by a
        specific cutoff (under eGFR21) before death, and add 
        corresponding values to the trajectory table, updating all consecutive
        entries accordingly.
    Args:
        df: a trajectory dataframe
        cutoff: an integer corresponding to an eGFR value
        event_name: string value corresponding to the event added at cutoff
        params: parameters dictionary, containing the slopes corresponding to
            binary indicators (params["model"]["slopes"])
    Returns:
        an updated trajectory frame
    """
    at_cutoff = add_age_at_cutoff(
        df, cutoff, interv_under_eq="21"
        ).assign(**{event_name: 1})

    df_ret = update_frame_with_cutoff_event(
        df, at_cutoff, event_name, params
    )

    return df_ret


def model_part1(df_cohort, params):
    """
    PART 1

    Function: 
        Generating initial trajectories, from initial age to max age,
        with slope values adjusted for diabetes and hypertension only
    Args:
        df_cohort: a cohort frame containing pre-sampled times of diabetes and
            hypertension initiation, indexed by pid and age
        params: parameters dictionary, containing the slopes corresponding to
            binary indicators (params["model"]["slopes"])
    Returns:
        DataFrame with updated eGFR and slope values, indexed by pid and age
    """
    df_cohort_slopes = update_slopes(df_cohort, params)
    df_init = get_decline(df_cohort_slopes)
    return df_init

def model_part2(df_init, params, common_rn_dict = None):
    """
    PART 2
     
    Function: 
        Adding events: reaching cutoff=60 under eGFR21 (g3),
        time of diagnosis and nephrology initiation
    Args:
        df_init: a trajectory frame including times of initiation, death (assumed at 101),
            diabetes and hypertension
        params: parameters dictionary, containing the slopes corresponding to
            binary indicators (params["model"]["slopes"])
        common_rn_dict: pre-sampled common random number values
    Returns:
        DataFrame with updated eGFR, slope, g3, diag and nephr values,
        indexed by pid and age
    """

    # first, add time of reaching G3a (based on the eGFR21 equation)
    df_g3 = add_event_at_cutoff(df_init, 60, "g3", params)

    df_temp = df_g3

    for cutoff in [60, 45, 30, 15]:

        # for each person, if they reach egfr value=cutoff before death, 
        # find the age of reaching cutoff
        at_cutoff = add_age_at_cutoff(df_temp, cutoff, params["model"]["interv_under_eq"])
        pid_at_cutoff = at_cutoff.reset_index().pid.values

        # for those at cutoff, sample diagnosis
        at_cutoff_subsamp_diag = interventions.sample_interventions(
            at_cutoff, 
            "diag", 
            params,
            common_rn = common_rn_dict["diag"][cutoff][pid_at_cutoff] if common_rn_dict is not None else None
        )
        at_cutoff_subsamp_nephr = interventions.sample_interventions(
            at_cutoff_subsamp_diag, 
            "nephr", 
            params,
            common_rn = common_rn_dict["nephr"][cutoff][pid_at_cutoff] if common_rn_dict is not None else None
        )        
        # update frame with the newly sampled values, recalculate the trajectory following cutoff
        df_temp = update_frame_with_cutoff_event(
            df_temp, at_cutoff_subsamp_nephr, ["diag", "nephr"], params, multiple_indic=True
        )

    df_events = df_temp
    
    return df_events

def model_part3(df_events, params):
    """
    PART 3
     
    Function: 
        Generates individual hazard functions corresponding to individual trajectories,
        based on hazard ratios specific to the time of diabetes incidence and the 
        time of reaching eGFR ranges
    Args:
        df_events: a trajectory dataframe with all events (g3, interventions)
            previously identified
        params: parameters dictionary, containing a hazards ratios table (params["model"]["death_HRs"])
    Returns:
        A frame containing individual hazard functions
    """
    assert df_events.query("age > 101").shape[0]==0, "Ages go above the maximum age of 101"

    calculated_at_cutoff = identify_age_at_cutoff_many(
        df_events, 
        cutoff_list = [60, 45, 30, 15]
        )
    
    df_events_with_cutoffs = pd.concat([
        df_events, 
        calculated_at_cutoff.set_index(["pid", "age"])
        ]).sort_index()

    df_events_with_ranges = (df_events_with_cutoffs.reset_index(level="age")
                                .assign(
                                    egfr_range_all = lambda x: egfr_formulas.egfr_to_ranges(x.egfr),
                                    egfr_range = lambda x: x.egfr_range_all.map(pd.Series({1:2, 3:2, 4:2,
                                                                                           2:2, 5:5, 6:6,
                                                                                           7:7, 8:8})),
                                    age = lambda x: np.minimum(x.age, 100.99).astype(int)
                                    )
                                )

    mortality_HRs = survival_helpers.get_conditional_mortality_HRs(
        df_events_with_ranges,
        HRs = param_helpers.get_HRs(params)
        )

    hazard_func_sex = param_helpers.get_hazard_functions_pivot(params)

    final_individual_hazards = survival_helpers.update_hazards(
        mortality_HRs,
        hazard_func_sex,
    )

    return final_individual_hazards

def model_part4(final_individual_hazards, df_events, params, sample_year=None,
                sample_cont=None):
    """
    PART 4
    
    Function:
        sampling death based on individual hazard functions 
        (final_individual_hazards), censoring all events following sampled death
        from the events frame (df_events)
    Args:
        final_individual_hazards: a frame containing individual hazard functions
        df_events: a frame containing complete trajectorite from initiation age
            until age 101
        params: parameters dictionary
        sample_year: (optional) a pre-sampled array of values [0,1] for
            common random numbers mortality sampling
        sample_cont: (optional) a pre-sampled array of values [0,1] for
            making discrete years continuous
    Returns:
        Final trajectories data frame from intiation age until death
    """
    
    if "pid" or "age" not in df_events.columns:
        df_events = df_events.reset_index()

    death_ages = survival_helpers.get_death_ages(
        final_individual_hazards, 
        params, 
        sample_year,
        sample_cont
        )
    pids = df_events.pid.drop_duplicates().values

    df_death = pd.DataFrame(
        {
            "pid": pids,
            "death_age": death_ages.round(2),
            "death": 1,
        }
    )
    df_with_death = combine_with_death(df_events, df_death)

    # update eGFR value at death
    df_final_to_update = (
        # grab last two entries per person (death and the event preceding death)
        df_with_death.groupby(level=["pid"])
        .nth([-1, -2])
    )
    updated_final_subset = get_decline(df_final_to_update)
    df_with_death.update(updated_final_subset)

    return df_with_death.ffill(axis=0).astype(dtype_dict)

def get_trajectories(
    df_cohort, params, save_path_id=None, print_timing=False
):    
    """
    Function:
        This function takes in an initial cohort and a set of parameters, and samples 
        eGFR trajectories from intiation until death
    Args:
        df_cohort: a cohort frame containing pre-sampled times of diabetes and
            hypertension initiation, indexed by pid and age
        params: parameters dictionary
        save_path_id: (optional) filename under which the trajectory will be saved
        print_timing: (optional) if True, print how much time the calculation took
    Returns:
        Final trajectories data frame from intiation age until death
    """
    i, end0 = if_print_timing(print_timing, None, 0)

    df_init = model_part1(df_cohort, params)

    df_events = model_part2(df_init, params)
    
    final_hazards = model_part3(df_events, params)
    
    df_with_death = model_part4(final_hazards, df_events, params)

    i, end1 = if_print_timing(print_timing, end0, i)

    if save_path_id:
        save_path = os.path.join("output", save_path_id)
        df_with_death.to_csv(save_path)

    return df_with_death.astype(dtype_dict)


def get_trajectories_two_scenarios(df_cohort, params, 
                                   crn_death = True, 
                                   crn_interv = False):
    """
    Function:
        This function takes in an initial cohort and a set of parameters, and samples 
        eGFR trajectories from intiation until death under both reference and
        counterfactual scenarios
    Args:
        df_cohort: a cohort frame containing pre-sampled times of diabetes and
            hypertension initiation, indexed by pid and age
        params: parameters dictionary
        crn_death: (optional) if True, apply common random sampling for death times
        crn_interv: (optional) if True, apply common random sampling for intervention assignment
    Returns:
        Final trajectories data frame from intiation age until death
    """
    df_init = model_part1(df_cohort, params)

    if crn_interv:
        sample_interv = {
            "diag": {
                60: np.random.uniform(low=0, high=1, size=params["cohort"]["n"]),
                45: np.random.uniform(low=0, high=1, size=params["cohort"]["n"]),
                30: np.random.uniform(low=0, high=1, size=params["cohort"]["n"]),
                15: np.random.uniform(low=0, high=1, size=params["cohort"]["n"]),
            },
            "nephr": {
                60: np.random.uniform(low=0, high=1, size=params["cohort"]["n"]),
                45: np.random.uniform(low=0, high=1, size=params["cohort"]["n"]),
                30: np.random.uniform(low=0, high=1, size=params["cohort"]["n"]),
                15: np.random.uniform(low=0, high=1, size=params["cohort"]["n"]),
            }
        }
    else:
        sample_interv = None

    # reference
    df_events_bl = model_part2(df_init, params, common_rn_dict = sample_interv)

    final_hazards_bl = model_part3(
        df_events_bl, params
    )

    # counterfactual
    params_cf = deepcopy(params)
    params_cf["model"]["interv_under_eq"] = "21"
    df_events_cf = model_part2(df_init, params_cf, common_rn_dict = sample_interv)

    final_hazards_cf = model_part3(
        df_events_cf, params_cf
    )

    # common random numbers
    if crn_death:
        sample_year = np.random.uniform(low=0, high=1, size=params["cohort"]["n"])
        sample_cont = np.random.uniform(low=0, high=1, size=params["cohort"]["n"])
    else:
        sample_year = None
        sample_cont = None
        
    df_with_death_bl = model_part4(
        final_hazards_bl, df_events_bl, params_cf, sample_year, sample_cont
    )
    df_with_death_cf = model_part4(
        final_hazards_cf, df_events_cf, params, sample_year, sample_cont
    )

    return df_with_death_bl.astype(dtype_dict), df_with_death_cf.astype(dtype_dict)