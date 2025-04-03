import numpy as np

"""
A dictionary mapping cutoff values to stages. Cutoff values
are the highest eGFR values that can be used to classify an individual
to a given stage.
"""
egfr_to_g_stages = {
    90: 'G2',
    60: 'G3a',
    45: 'G3b',
    30: 'G4',
    15: 'G5'
}

"""
Parameter sources:

2009:
Levey AS, Stevens LA, Schmid CH, et al. A new equation to estimate glomerular filtration rate. 
Annals of internal medicine 2009; 150: 604–612.

2021:
Inker LA, Eneanya ND, Coresh J, et al. New creatinine- and cystatin C–based equations to 
estimate GFR without race. N Engl J Med. 2021;385:1737-1749.
"""
CKD_EPI_params = {
    "21": {
        "kappa": {"F": 0.7, "M": 0.9},
        "alpha": {"F": -0.241, "M": -0.302},
        "sex_adj": {"F": 1.012, "M": 1},
        "race_adj": {"B": 1, "NB": 1},
        "max_exp": -1.2,
        "age_base": 0.9938,
        "beta": 142,
    },
    "09": {
        "kappa": {"F": 0.7, "M": 0.9},
        "alpha": {"F": -0.329, "M": -0.411},
        "sex_adj": {"F": 1.018, "M": 1},
        "race_adj": {"B": 1.159, "NB": 1},
        "max_exp": -1.209,
        "age_base": 0.993,
        "beta": 141,
    },
}

def CKD_EPI(scr, sex, age, race, formula):
    """
    Function:
        Calculates the eGFR using the CKD-EPI 2021 or 2009 Creatinine formula
    Args:
        scr: serum creatinine value (float)
        sex: in ["M", "F"]
        age: in years (float)
        race: in ["B", "NB"] (needed for 2009 formula)
        formula: in ["09", "21"], corresponding to CKD-EPI Creatinine-based 2009 or 
            CKD-EPI Creatinine-based 2021 formula
    Returns:
        eGFR value
    """
    params = CKD_EPI_params[formula]
    eGFR = (
        params["beta"]
        * pow(min(scr / params["kappa"][sex], 1), params["alpha"][sex])
        * pow(max(scr / params["kappa"][sex], 1), params["max_exp"])
        * pow(params["age_base"], age)
        * params["sex_adj"][sex]
        * params["race_adj"][race]
    )

    return eGFR

def creat_from_eGFR(egfr, sex, age, race, formula):
    """
    Function:
        Calculates serum creatinine value at a given level of eGFR,
        using the CKD-EPI 2021 or 2009 Creatinine formula
    Args:
        egfr: egfr value (float)
        sex: in ["M", "F"]
        age: in years (float)
        race: in ["B", "NB"] (needed for 2009 formula)
        formula: in ["09", "21"], corresponding to CKD-EPI Creatinine-based 2009 or 
            CKD-EPI Creatinine-based 2021 formula
    Returns:
        serum creatinine value
    """
    params = CKD_EPI_params[formula]

    a = egfr / (
        params["beta"]
        * pow(params["age_base"], age)
        * params["sex_adj"][sex]
        * params["race_adj"][race]
    )
    if a == 0:
        print(
            "a==0 in creat_from_eGFR, egfr", egfr, "sex", sex, "age", age, "race", race
        )
        return 0
    
    # case when Scr <= kappa, corresponding to min(Scr/kappa, 1)^alpha part of the eq
    if a >= 1:
        SeCr = pow(a, 1 / params["alpha"][sex]) * params["kappa"][sex]
    # case when Scr > kappa, corresponding to max(Scr/kappa, 1)^math_exp part of the eq
    else:
        SeCr = pow(a, 1 / params["max_exp"]) * params["kappa"][sex]
    return SeCr

def eGFR_2021_to_2009(egfr_2021, sex, age, race):
    """
    Function:
        Calculates CKD-EPI 2009 eGFR value corresponding to the CKD-EPI 2021 eGFR value
        at provided values of sex, age and race
    Args:
        egfr_2021: egfr value expressed in terms of CKD-EPI 2021 (float)
        sex: in ["M", "F"]
        age: in years (float)
        race: in ["B", "NB"] (needed for 2009 formula)
    Returns:
        egfr value expressed in terms of CKD-EPI 2009
    """
    if egfr_2021 == 0:
        return 0
    elif egfr_2021 < 0:
        return None
    SeCr = creat_from_eGFR(egfr_2021, sex, age, race, formula = "21")
    if SeCr == 0:
        return 0
    egfr_2009 = CKD_EPI(SeCr, sex, age, race, formula = "09")
    return egfr_2009

def eGFR_2009_to_2021(egfr_2009, sex, age, race):
    """
    Function:
        Calculates CKD-EPI 2021 eGFR value corresponding to the CKD-EPI 2009 eGFR value
        at provided values of sex, age and race
    Args:
        egfr_2009: egfr value expressed in terms of CKD-EPI 2009 (float)
        sex: in ["M", "F"]
        age: in years (float)
        race: in ["B", "NB"] (needed for 2009 formula)
    Returns:
        egfr value expressed in terms of CKD-EPI 2021
    """
    if egfr_2009 == 0:
        return 0
    elif egfr_2009 < 0:
        return None
    SeCr = creat_from_eGFR(egfr_2009, sex, age, race, formula = "09")
    if SeCr == 0:
        return 0
    egfr_2021 = CKD_EPI(SeCr, sex, age, race, formula = "21")
    return egfr_2021

def update_egfr_equation(df_traj, change_type="21to09"):
    """
    Function:
        Updates eGFR values in a trajectory frame from eGFR 2021 to 2009 or 2009 to 2021
    Args:
        df_traj: a trajectory dataframe
        change_type: in ["21to09", "09to21"] - direction of eGFR value updated
    Returns:
        trajectory frame with eGFR values expressed in terms of the selected eGFR equation
    """
    if change_type == "21to09":
        change_fun = eGFR_2021_to_2009
    elif change_type == "09to21":
        change_fun = eGFR_2009_to_2021

    df_traj_both = df_traj.assign(
        egfr=lambda y: y.apply(
            lambda x: change_fun(
                x.egfr,
                x.sex,
                x.age,
                x.race,
            ),
            axis=1,
        )
        .round(2)
        .drop(columns="slope")
    )

    return df_traj_both

def translate_trajectory(df_traj, change_type="21to09"):
    """
    Function:
        Translates a trajectory frame from eGFR 2021 to 2009 or 2009 to 2021
        and calculates updated slope values
    Args:
        df_traj: a trajectory dataframe
        change_type: in ["21to09", "09to21"] - direction of eGFR value updated
    Returns:
        trajectory frame with eGFR values expressed in terms of the selected eGFR equation
        and updated slopes
    """

    df_traj_translated = update_egfr_equation(df_traj, change_type).fillna(0)

    df_traj_updated = (
        df_traj_translated.set_index(["pid"])
        .assign(
            degfr=lambda x: x.groupby(level=0)["egfr"].diff().shift(-1),
            dt=lambda x: x.groupby(level=0)["age"].diff().shift(-1),
            slope=lambda x: (-x.degfr / x.dt).round(3),
        )
        .drop(columns=["dt", "degfr"])
        .groupby("pid")
        .ffill()
    )
    return df_traj_updated


def egfr_to_stages(egfr):   
    """
    Assigns a CKD stage corresponding to an eGFR value
    G1:  > 90
    G2:  60 < eGFR <= 90
    G3a: 45 < eGFR <= 60
    G3b: 30 < eGFR <= 45
    G4:  15 < eGFR <= 30
    G5:  0  < eGFR <= 15

    Args:
        egfr (pandas.Series, np.array or int): a column of eGFR values
    Returns:
        a series of CKD stages
    """
    conditions = [
        egfr > 90,
        (egfr > 60) & (egfr <= 90),
        (egfr > 45) & (egfr <= 60),
        (egfr > 30) & (egfr <= 45),
        (egfr > 15) & (egfr <= 30),
        (egfr > 0) & (egfr <= 15)
    ]
    
    choices = ["G1", "G2", "G3a", "G3b", "G4", "G5"]
    
    stages = np.select(conditions, choices, default=None)

    return stages

def egfr_to_ranges(egfr): 
    """
    Assigns a number corresponding to a range that each eGFR value 
    falls into:
    1: > 105
    2: 90 < eGFR <= 105
    3: 75 < eGFR <= 90
    4: 60 < eGFR <= 75
    5: 45 < eGFR <= 60
    6: 30 < eGFR <= 45
    7: 15 < eGFR <= 30
    8: eGFR <= 15
    These will used to identify appropriate hazard ratios from Fox et al

    Args:
        egfr (pandas.Series): a column of eGFR values
    Returns:
        a series of range ids corresponding to each eGFR value
    """
    conditions = [
        egfr > 105,
        (egfr > 90) & (egfr <= 105),
        (egfr > 75) & (egfr <= 90),
        (egfr > 60) & (egfr <= 75),
        (egfr > 45) & (egfr <= 60),
        (egfr > 30) & (egfr <= 45),
        (egfr > 15) & (egfr <= 30),
        (egfr > 0) & (egfr <= 15),
        egfr <= 0
    ]
    
    choices = [1, 2, 3, 4, 5, 6, 7, 8, 8]
    
    stages = np.select(conditions, choices, default=None)

    return stages