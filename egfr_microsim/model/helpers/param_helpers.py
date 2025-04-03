import os
import numpy as np
import pandas as pd
import types

def dict_to_sns(d):
    """
    Transforms a dictionary into a SimpleNamespace object, to allow references of type
    dictionary.key rather than dictionary[key]
    Args:
        d (dict)
    Returns:
        a SimpleNamespace object
    """
    return types.SimpleNamespace(**d)

def get_cohort(cohorts_dir, cohort_id):
    """
    Reads in a cohort file corresponding to cohort_id in cohorts_dir,
    using standard file names

    Args:
        cohorts_dir (str): directory path
        cohort_id (int): cohort index  
    Returns:
        a pandas.DataFrame defining relative ordering of slopes
    """
    return pd.read_csv(
        os.path.join(cohorts_dir, "".join(("cohort", str(cohort_id), ".csv")))
    ).set_index(["pid", "age"])

def get_slopes(params, original = False):
    """
    Reads in a file with slope parameters

    Args:
        params (dict): parameters dictionary
        original (bool): get the original Hoerger parameters
    Returns:
        a pandas.DataFrame object containing slope parameters
    """
    if original:
        return pd.read_csv(params["model"]["file_paths"]["slopes_original_path"])
    else:
        return pd.read_csv(params["model"]["file_paths"]["slopes_path"])

def get_comorb_incidence(params, comorb):
    """
    Reads in a file with comorbidity incidence values

    Args:
        params (dict): parameters dictionary
        comorb (str): comorbidity ("diabetes" or "hypertension")
    Returns:
        a pandas.DataFrame object containing comorbidity incidence values
    """
    if comorb == "diabetes":
        return pd.read_csv(params["model"]["file_paths"]["diabetes_incid_path"])
    elif comorb == "hypertension":
        return pd.read_csv(params["model"]["file_paths"]["hypertension_incid_path"])

def get_slopes_index(params):
    """
    Provides the index of the slopes frame (all variable combinations, in the same
    order as in the slopes input data frame)

    Args:
        params (dict): parameters dictionary
    Returns:
        a pandas.DataFrame.index object
    """
    df_slopes = get_slopes(params)
    df_slopes_index = df_slopes.set_index(["dm", "ht", "g3"]).index
    return df_slopes_index

def get_HRs(params):
    """
    Reads in a file with mortality hazards ratios based on eGFR range and 
    diabetes status. Ranges are described in 
    egfr_microsim.model.egfr_formulas.egfr_to_ranges

    Args:
        params (dict): parameters dictionary
    Returns:
        a pandas.DataFrame object containing mortality hazards ratios
    """
    return (
        pd.read_csv(params["model"]["file_paths"]["HRs_path"])
        .filter(["dm", "egfr_range", "mort_HR"])
    )

def get_slopes_order(params):
    """
    Reads in a file defining the desired relative order between sampled parameters.
    It takes form of an LxL lower triangular matrix (with L = number of parameters 
    in a single param set), indexed by parameter ids, where 
        array[i,j] = 0 enforces no ordering between the values i and j
        array[i,j] = 1 enforces that value of parameter i > value of parameter j
    in all sampled parameter sets

    Args:
        params (dict): parameters dictionary
    Returns:
        a pandas.DataFrame defining relative ordering of slopes
    """
    return (
        pd.read_csv(params["model"]["file_paths"]["slopes_order_path"],
                    index_col=0)
    )

def get_slopes_order_list(params):
    """
    Reads in a matrix defining relative ordering of parameter values and
    transforms it into a list.

    Args:
        params (dict): parameters dictionary
    Returns:
        a list in the form [(target_idx, others_idx), ...], where:
            target_idx: int, parameter index value
            others_idx: list or int, indices of parameters that must be lower
                        than target_idx value
            in order to enforce the rule that more comorbidities -> higher slope
    """
    slopes_relative_order = get_slopes_order(params)

    return_list = []
    for i in slopes_relative_order.index:
        slope_i = slopes_relative_order.iloc[i]
        larger_than = slope_i[slope_i==1].index.astype(int).to_list()
        if len(larger_than) > 0:
            return_list.append((i, larger_than))

    return return_list

def get_param_df(params_dir, i = None):
    """
    Reads in pre-sampled parameter sets from standard file locations.

    Args:
        params_dir (str): directory path
        (optional) i (int): parameter set id. If provided, selects
            a parameter file indexed by i.
    Returns:
        a pandas.DataFrame containing pre-sampled parameter sets
    """
    if i != None:
        df_param = pd.read_feather(os.path.join(params_dir, "params%i.feather" % i))
    else:
        df_param = pd.read_feather(os.path.join(params_dir, "params.feather"))
    return df_param

def get_param_df_from_1k(params_dir, parameter_set_number):
    """
    Reads in pre-sampled parameter sets distributed across 1k-parameter
    set files, and combines them into a single file, using standard 
    file locations.

    Args:
        params_dir (str): directory path
        parameter_set_number: number of parameter sets in the experiment
    Returns:
        a pandas.DataFrame containing a single file containing all
        pre-sampled parameter sets
    """
    all_params = []
    for i in range(int(parameter_set_number / 1000)):
        df_param = get_param_df(params_dir, i)
        df_param.index = range(i * 1000, (i + 1) * 1000)
        all_params.append(df_param)
    param_set_df = pd.concat(all_params)
    return param_set_df

def update_slope_params(params, slopes_sample):
    """
    Reads in the file with baseline slopes (for formatting only) 
    and updates their values with slopes_sample

    Args:
        params (dict): parameters dictionary
        slopes_sample (int): values of slopes (ordered according to 
        the index of params["model"]["slopes"])
    Returns:
        a pandas.DataFrame containing a single file containing all
        pre-sampled parameter sets
    """ 
    df_slopes = get_slopes(params)
    params["model"]["slopes"] = df_slopes.filter(["dm", "ht", "g3"]).assign(
        slope=slopes_sample
    )
    return params

def transform_lifetables(life_table):
    """
    Formats life table for later processing

    Args:
        life_table (pandas DataFrame)

    Returns:
        life_table (pandas DataFrame)
    """

    life_table.columns = ["Age"] + list(life_table[0:1].values[0][1:])
    life_table = life_table.iloc[1:103].reset_index(drop=True)
    life_table.index.name = "age"

    life_table = life_table.assign(
        qx=lambda x: x.qx.astype(float),
        lx=lambda x: x.lx.str.replace(",", "").astype(int),
        dx=lambda x: x.dx.str.replace(",", "").astype(int),
        Lx=lambda x: x.Lx.str.replace(",", "").astype(int),
        Tx=lambda x: x.Tx.str.replace(",", "").astype(int),
        ex=lambda x: x.ex.astype(float),
    )
    return life_table

def get_LT_functions(params, function_name, rel_path = None):
    """
    Reads in data tables from CSV files, using provided paths,
    and returns hazard functions (using the qx column)

    Args:
        params (dict): parameters dictionary
        function_name (string): defining the column in CDC life tables
        (optional) rel_path (str): if provided, used to identify
            path to the egfr_microsim repo in a local file system
    Returns:
        the selected life tables function (pandas DataFrame):
        indexed by ages 0-100, with 2 columns (1 per sex)
    """
    female_LT_path = params["model"]["file_paths"]["female_LT_path"]
    male_LT_path = params["model"]["file_paths"]["male_LT_path"]

    if rel_path:
        female_LT_path = os.path.join(rel_path, female_LT_path)
        male_LT_path = os.path.join(rel_path, male_LT_path)

    male_LT = pd.read_csv(male_LT_path)
    female_LT = pd.read_csv(female_LT_path)

    male_LTs = transform_lifetables(male_LT)
    female_LTs = transform_lifetables(female_LT)

    LT_func = pd.DataFrame()

    for sex, life_table in zip(("M", "F"), (male_LTs, female_LTs)):
        LT_func[sex] = life_table[function_name]
    return LT_func

def get_hazard_functions(params, rel_path=None):
    """
    A wrapper around get_LT_functions for hazard functions.

    Args:
        params (dict): parameters dictionary
        (optional) rel_path (str): if provided, used to identify
            path to the egfr_microsim repo in a local file system
    Returns:
        the hazards function (pandas DataFrame):
        indexed by ages 0-100, with 2 columns (1 per sex)
    """ 
    hazard_func = get_LT_functions(params, "qx", rel_path)

    return hazard_func

def get_hazard_functions_pivot(params, rel_path = None):
    """
    Reshapes hazard functions frame from 100x2 to 2x100

    Args:
        params (dict): parameters dictionary
        (optional) rel_path (str): if provided, used to identify
            path to the egfr_microsim repo in a local file system
    Returns:
        the hazards function (pandas DataFrame):
        indexed by sex, with 100 columns (1 per age)
    """ 
    hazard_func = get_hazard_functions(params, rel_path)
    hazard_func_sex = (
        hazard_func.reset_index()
        .melt(
            value_vars=["M", "F"], id_vars="age", var_name="sex", value_name="mort_HR"
        )
        .assign(pid="all")
    )
    return hazard_func_sex