import numpy as np
import pandas as pd
import egfr_microsim.model.helpers.survival_helpers as survival_helpers
from itertools import chain

def get_comorb_hazard(comorb_incid, index_cols=[]):
    """
    Function:
        Calculates an age- dependent hazard function based on incidence statistics
        in age buckets
    Args: 
        comorb_incid: a frame with comorbidity incidence values
        index_cols: additional variables defining comorbidity incidence hazard functions
    Returns:
        A hazard function, indexed by age and index_cols
    """
    comorb_incid = comorb_incid.assign(
        hazard = lambda x: x.inc_per_1k/1000,
        age_diff = lambda x: (x.age_max - x.age_min + 1).astype(int),
    )

    ages = [
        range(el1, el2 + 1) for el1, el2 in zip(comorb_incid.age_min.values, comorb_incid.age_max.values)
    ]
    ages = list(chain.from_iterable(ages))
    index = comorb_incid.index.repeat(comorb_incid["age_diff"])

    hazard = (
        comorb_incid.reindex(index)
        .assign(age=ages)
        .set_index(["age"] + index_cols)
        .filter(["hazard"])
    )

    return hazard

def get_cdfs(hazard, index_cols=[]):
    """
    Function:
        Calculates a CDF of comorbidity incidence corresponding to the hazard
    Args: 
        comorb_incid: a frame with comorbidity incidence values
        index_cols: additional variables defining comorbidity incidence hazard functions
    Returns:
        cdf: A CDF function, indexed by age and index_cols, corresponding to the probability of
            developing a comorbidity among the whole population (does not sum to 1)
        lifetime_risk: of developing a comorbidity
        cdf_develops: a CDF function corresponding to the probability of developing a comorbidity
            among those who will
    """
    if len(index_cols) > 0:
        cumulative_hazard = hazard.groupby(index_cols).cumsum()
    else:
        cumulative_hazard = hazard.cumsum()

    cdf = cumulative_hazard.apply(lambda x: 1 - np.exp(-x)).rename(
        columns={"hazard": "prob"}
    )
    if len(index_cols) > 0:
        lifetime_risk = cdf.groupby(index_cols).max()
    else:
        lifetime_risk = cdf.iloc[-1]

    cdf_develops = cdf / lifetime_risk

    return cdf, lifetime_risk, cdf_develops

def get_cdfs_from_incid(comorb_incid, index_cols=[]):
    """
    Function:
        Calculates a CDF of comorbidity incidence based on incidence values
    Args: 
        comorb_incid: a frame with comorbidity incidence values
        index_cols: additional variables defining comorbidity incidence hazard functions
    Returns:
        cdf: A CDF function, indexed by age and index_cols, corresponding to the probability of
            developing a comorbidity among the whole population (does not sum to 1)
        lifetime_risk: of developing a comorbidity
        cdf_develops: a CDF function corresponding to the probability of developing a comorbidity
            among those who will
    """
    hazard = get_comorb_hazard(comorb_incid, index_cols)
    cdf, lifetime_risk, cdf_develops = get_cdfs(hazard, index_cols)
    return cdf, lifetime_risk, cdf_develops


def get_comorb_age(cohort_df, cdf, lifetime_risk, comorb):
    """
    Function:
        Sample ages of comorbidity incidence for individuals in cohort_df based on 
        the piecewise exponential frailty model defined by the CDF and a lifetime risk
    Args: 
        cohort_df: Initial frame with demographic values (sex, race) and eGFR
        cdf: CDF defining the frailty model
        lifetime_risk: proportion of individuals who ever develop comorbidity, defining
            the frailty model
        comorb: name of comorbidity
    Returns:
        a cohort frame with ages of 
    """

    coh_size = cohort_df.shape[0]

    # step 1: sampling indicator functions (gets comorbidity)
    lifetime_prob = dict(
        lifetime_risk.reset_index()
        .assign(sex = lambda x: x.male.map({0: "F", 1: "M"}))
        .filter(["sex", "prob"])
        .values
    )
    
    gets_comorb = np.random.binomial(
        size=coh_size, n=1, p = cohort_df.sex.map(lifetime_prob)
    )

    cohort_df = cohort_df.assign(gets_comorb=gets_comorb)

    # step 2: sampling time to event for those with comorbidity
    comorb_cdf = (
        cdf.reset_index()
        .assign(sex = lambda x: x.male.map({0: "F", 1: "M"}))
        .pivot(columns="age", index="sex", values="prob")
    )

    all_cdf = (
        cohort_df.query("gets_comorb==1")
        .reset_index()
        .merge(comorb_cdf, on="sex", how="left")
        .set_index("pid")
        .drop(columns=["egfr", "sex", "gets_comorb", "race"])
    )

    age_of_comorb = survival_helpers.sample_cat_vectorized(
        all_cdf, distr="cdf"
    ) + np.random.uniform(low=0.0, high=1.0, size=all_cdf.shape[0]).round(2)

    return pd.DataFrame(
        data={"pid": all_cdf.index, "age": age_of_comorb, comorb: int(1)}
    )
