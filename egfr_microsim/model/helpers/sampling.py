import numpy as np
import pandas as pd

import egfr_microsim.model.helpers.param_helpers as param_helpers
import egfr_microsim.model.helpers.comorb_helpers as comorb_helpers
from egfr_microsim.model.egfr_model import dtype_dict
dtype_dict.pop("slope")

def sample_dm_ht_times(df_demo_egfr, params):
    """
    Function:
        Sampling times of hypertension (ht) and diabetes (dm) incidence
    Args: 
        df_demo_egfr: Initial frame with demographic values (sex, race) and eGFR
        params: a parameters dictionary defining the cohort
    Returns:
        A cohort frame with diabetes and hypertension events
    """
    # read in parameters (annual hazards, grouped by age)
    df_dm = param_helpers.get_comorb_incidence(params, "diabetes")
    df_ht = param_helpers.get_comorb_incidence(params, "hypertension")

    # get continuous cumulative hazard functions and corresponding CDFs
    cdf_dm, lifetime_risk_dm, cdf_develops_dm = comorb_helpers.get_cdfs_from_incid(df_dm)
    cdf_ht, lifetime_risk_ht, cdf_develops_ht = comorb_helpers.get_cdfs_from_incid(df_ht, index_cols=["male"])

    # DM depends only on age, while HT depends on sex as well. We reshape data
    # corresponding to DM so it has the same format as HT (duplicating for male and female)
    lifetime_risk_dm_sex = pd.DataFrame(
        {"male": [0,1], "prob": [lifetime_risk_dm.values[0]]*2}
    ).set_index("male")

    cdf_develops_dm_sex = pd.concat(
        [cdf_develops_dm.assign(male=1), cdf_develops_dm.assign(male=0)]
    ).set_index(["male"], append=True)

    # for people in the df_demo_egfr frame, sample the age of incidence
    # from a piece exponential frailty model
    ht_age = comorb_helpers.get_comorb_age(
        df_demo_egfr, cdf_develops_ht, lifetime_risk_ht, comorb="ht"
    )
    diab_age = comorb_helpers.get_comorb_age(
        df_demo_egfr, cdf_develops_dm_sex, lifetime_risk_dm_sex, comorb="dm"
    )

    # generate final frame
    df_events = (
        pd.concat([ht_age, diab_age])
        .set_index(["pid", "age"])
        .sort_index()
        .fillna(0)
        .astype(int)
    )

    return df_events


def sample_cohort(params):
    """
    Function:
        Samples a cohort based on parameters, assigning initial ages,
        eGFR values sampled from a normal distribution, and incidence times 
        of diabetes and hypertension sampled from national incidence functions,
        and a death event at age 101
    Args: 
        params: a parameters dictionary defining the cohort
    Returns:
        A cohort frame
    """
    pids = range(params["cohort"]["n"])

    df_demo_egfr = pd.DataFrame(
        {
            "pid": pids,
            "egfr": np.random.normal(
                params["cohort"]["init_egfr"]["mu"],
                params["cohort"]["init_egfr"]["sigma"],
                params["cohort"]["n"],
            ),
            "sex": np.vectorize({0: "F", 1: "M"}.get)(
                np.random.binomial(1, params["cohort"]["male"], params["cohort"]["n"])
            ),
            "race": np.vectorize({0: "NB", 1: "B"}.get)(
                np.random.binomial(1, params["cohort"]["black"], params["cohort"]["n"])
            ),
        }
    ).set_index("pid")

    df_events = sample_dm_ht_times(df_demo_egfr, params)

    init_ev = pd.DataFrame({"age": params["cohort"]["init_age"], "pid": pids})
    death_ev = pd.DataFrame(
        {"age": params["cohort"]["max_age"], "death": 1, "pid": pids}
    )
    df_start_end = pd.concat([init_ev, death_ev]).set_index(["pid", "age"])
    min_age = params["cohort"]["init_age"]

    df_init = (
        pd.concat([df_events, df_start_end])
        .sort_index()
        .fillna(0)
        # propagate indicators to later ages
        .groupby(["pid"])
        .cummax()
        .query("age >= @min_age")
        .assign(g3=0, diag=0, nephr=0)
        .join(df_demo_egfr)
        # address duplicates
        .groupby(["pid", "age"])
        .last()
        .assign(egfr = lambda x: x.egfr.round(2))
        .astype(dtype_dict)
    )

    return df_init 