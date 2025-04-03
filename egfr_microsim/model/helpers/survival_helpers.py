import numpy as np
import pandas as pd

def get_conditional_mortality_HRs(df, HRs):
    """
    Assigns a mortality hazard ratio corresponding to each entry in the
    trajectory table df (depending on egfr_range and diabetes status)

    Args:
        df (pandas.DataFrame): a trajectory frame indexed by pid and age
        HRs (pandas.DataFrame): a frame containing hazards ratios corresponding
            to egfr ranges and diabetes status
    Returns:
        df with a column corresponding to mortality hazard ratios for each entry
    """
    mortality_HRs = (
        df.reset_index()
        .merge(
            HRs.filter(["egfr_range", "dm", "mort_HR"]),
            how="left",
            on=["egfr_range", "dm"],
        )
    )
    return mortality_HRs

def update_hazards(mortality_HRs, hazard_func_sex):
    """
    Generates an age-dependent hazard function for individual, based on their
    sex- and age-dependent baseline hazards and trajectory-dependent 
    mortality hazard ratios.

    Args:
        mortality_HRs: a table of year-by-year mortality_HRs for each person, 
            conditional on their covariates and trajectories at each age
        hazard_func_sex: an age- and sex-dependent baseline hazard from life tables 
    """
    min_age = mortality_HRs.age.min()
    
    hazard_before_init = pd.pivot_table(
        hazard_func_sex.query("age < @min_age"),
        values="mort_HR",
        index=["pid", "sex"],
        columns=["age"],
    )

    hazard_func_sex_post_init = hazard_func_sex.query("age >= @min_age")

    final_hazards = (
        pd.pivot_table(
            mortality_HRs,
            values="mort_HR",
            index=["pid", "sex"],
            columns=["age"],
        )
        .ffill(axis=1)
    )

    for sex in ["F", "M"]:
        final_hazards.loc[:, sex, :] = (
            final_hazards.loc[:, sex, :]
            .mul(hazard_func_sex_post_init.query("sex==@sex").mort_HR.values)
            .to_numpy()
        )

    final_hazards_from_birth = (
        final_hazards
        .reset_index()
        .merge(hazard_before_init, on="sex")
        .set_index(["pid", "sex"])
        .sort_index(axis=1)
    )

    return final_hazards_from_birth

def cdf_from_hazards(hazards, axis=None):
    """
    Calculates the CDF from hazards, using a piecewise exponential function
    discretized by age. Adjusts last entry in CDF to equal 1.
    
    Args:
        hazards (pandas.DataFrame or numpy.ndarray): a frame containing the hazard function,
            indexed by age
        axis (int): indicating orientation of the table (horizontal with axis=1
            or vertical with axis=0)
    Returns:
        cdf (pandas.DataFrame or numpy.ndarray): the CDF function
    """
    cum_hazard = np.cumsum(hazards, axis=axis)
    cdf = cum_hazard.apply(lambda x: 1 - np.exp(-x))

    # force setting last entry to 1,
    # since lim x-->inf (1-exp(-x)) = 1
    if axis==1:
        cdf.loc[:, 100] = 1
    else:
        cdf.loc[100] = 1
    return cdf

def pdf_from_hazards(hazards, axis=None):
    """
    Calculates the PDF from hazards, using a piecewise exponential function
    discretized by age.
    
    Args:
        hazards (pandas.DataFrame or numpy.ndarray): a frame containing the hazard function,
            indexed by age
        axis (int): indicating orientation of the table (horizontal with axis=1
            or vertical with axis=0)
    Returns:
        pdf (pandas.DataFrame or numpy.ndarray): the PDF function
    """
    cum_hazard = np.cumsum(hazards, axis=axis)
    exp_part = cum_hazard.apply(lambda x: np.exp(-x))
    pdf = exp_part*hazards

    return pdf

def survival_from_hazards(hazards, axis=None):
    """
    Calculates the survival function from hazards, using a piecewise exponential function
    discretized by age. Adjusts first entry in the survival function to equal 1.
    
    Args:
        hazards (pandas.DataFrame or numpy.ndarray): a frame containing the hazard function,
            indexed by age
        axis (int): indicating orientation of the table (horizontal with axis=1
            or vertical with axis=0)
    Returns:
        surv (pandas.DataFrame or numpy.ndarray): the survival function
    """
    cum_hazard = np.cumsum(hazards, axis=axis)
    surv = cum_hazard.apply(lambda x: np.exp(-x))

    if axis==1:
        surv.loc[:, 0] = 1
    else:
        surv.loc[0] = 1
    
    return surv


def pdf_from_cdf(cdf, axis=None):
    """
    Calculates PDF from CDF, discretized by age.
    
    Args:
        cdf (pandas.DataFrame or numpy.ndarray): a frame containing the CDF function,
            indexed by age
        axis (int): indicating orientation of the table (horizontal with axis=1
            or vertical with axis=0)
    Returns:
        pdf (pandas.DataFrame or numpy.ndarray): the PDF function
    """

    pdf = cdf - cdf.shift(1, axis=axis).fillna(0)
    assert np.all(pdf.sum(axis=axis).round(10) == 1), "PDF must sum to 1"
    return pdf

def cdf_from_pdf(pdf, axis=None):
    """
    Calculates CDF from PDF, discretized by age. Adjusts last entry in CDF to equal 1.
    
    Args:
        pdf (pandas.DataFrame or numpy.ndarray): a frame containing the PDF function,
            indexed by age
        axis (int): indicating orientation of the table (horizontal with axis=1
            or vertical with axis=0)
    Returns:
        cdf (pandas.DataFrame or numpy.ndarray): the CDF function
    """
    
    # assert that the pdf has the right form (sums to 1)
    # (with precision up to 10 floating points)
    assert np.all(pdf.sum(axis=axis).round(10) == 1)

    cdf = np.cumsum(pdf)
    if axis==1:
        cdf.loc[:, 100] = 1
    else:
        cdf.loc[100] = 1

    return cdf


def conditional_on_survival(age_surv, distr, distr_type, axis=None):
    """
    Adjusts a provided distribution `distr` to reflect conditioning on survival 
    until age `age_surv`

    For CDF: turns P(T <= t) --> P(T <= t | T > age_surv)
    For SURV: turns P(T > t) --> P(T > t | T > age_surv)
    
    Args:
        age_surv (int): the age until which one survived
        distr (numpy.ndarray): containing either CDF, survival or PDF functions,
            indexed by age
        distr_type (string): indicating the type of distribution provided (CDF, survival or PDF)
        axis (int): indicating orientation of the table (horizontal with axis=1
            or vertical with axis=0)
    Returns:
        distr_t (numpy.ndarray of floats): a truncated function of provided type
    """
    assert (int(age_surv) == age_surv), "age_surv must be an int"
    assert (age_surv > 0 & age_surv < 100), "age_surv must be in range [1,99]"
    assert (distr_type in ["cdf", "surv"]), "func_type must be one of ['cdf', 'surv']"

    if distr_type == "cdf":
        surv = 1-distr
    elif distr_type == "surv":
        surv = distr

    if axis == 0:
        surv_till_age = surv.loc[age_surv-1]
        surv_t = surv.loc[age_surv:]/surv_till_age
        surv_t.loc[100] = 0

    elif axis == 1:
        surv_till_age = surv.loc[:,age_surv-1]
        surv_t = surv.loc[:, age_surv:].div(surv_till_age, axis=0)
        surv_t.loc[:, 100] = 0

    if distr_type == "cdf":
        cdf_t = 1-surv_t
        return cdf_t
    elif distr_type == "surv":
        return surv_t

def sample_cat_vectorized(df, surv_sample=None, distr="pdf", axis=1):
    """
    vectorized multivariate sampling, following
    https://github.com/DARTH-git/darthtools/blob/main/R/samplev.R

    Optionally, accepts a pre-sampled array of values [0,1] for
    common random numbers mortality sampling
    
    Args:
        df (pandas DataFrame): containing either CDF or PDF functions,
            indexed by ['pid', 'sex'] with columns corresponding to ages
        surv_sample (numpy.ndarray): an optional pre-sampled vector of probabilities
            which will be used for sampling from the provided CDF or PDF
        distr (string): indicating the type of distribution provided (pdf or cdf)
        axis (int): indicating orientation of the table (horizontal with axis=1
            or vertical with axis=0)
    Returns:
        death_ages (numpy.ndarray of floats): an array of sampled death ages
    """
    assert (distr in ["cdf", "pdf"]), "func_type must be one of ['cdf', 'pdf']"

    cat_names = df.columns
    dim = df.shape[0]
    if distr == "pdf":
        cdf = df.cumsum(axis=axis).to_numpy()
    elif distr == "cdf":
        cdf = df.to_numpy()

    if surv_sample is None:
        surv_sample = np.random.uniform(low=0, high=1, size=dim)
    # returns the index of the first True value
    surv_age = (cdf >= surv_sample.reshape(-1, 1)).argmax(axis=axis)
    death_ages = cat_names[surv_age].to_list()

    return death_ages

def get_death_ages(final_hazards, params, sample_year=None, sample_cont=None):
    """
    Samples age of death based on individual hazards functions
    
    Args:
        final_hazards (pandas.DataFrame): an N (cohort size) by 100 (max age) frame 
            of hazard functions
        params (dict): parameters dictionary
        (optional) sample_year (np.array): a pre-sampled array of values [0,1] for
        common random numbers mortality sampling (N (cohort size) by 100 (max age))
        (optional) sample_cont (np.array): a pre-sampled array of values [0,1] for
            making discrete years continuous
    Returns:
        an N-length array with age of death per each individual
    """

    if sample_cont is None:
        sample_cont = np.random.uniform(low=0.0, high=1.0, size=final_hazards.shape[0])

    surv_final = survival_from_hazards(final_hazards, axis=1)

    surv_updated = conditional_on_survival(
        params["cohort"]["init_age"], surv_final, distr_type = "surv", axis = 1
        )    

    cdf_updated = 1-surv_updated

    death_ages = sample_cat_vectorized(
        cdf_updated, sample_year, distr="cdf"
    ) + sample_cont

    return death_ages