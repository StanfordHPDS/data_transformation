import numpy as np

def frac_already_intervened(df, event_name):
    """
    Function:
        Calculates the fraction of individuals in df who already
        have been assigned a given intervention
    Args:
        df: a trajectory dataframe
        event_name: name of intervention (in ["nephr", "diag"])
    Returns:
        A float value corresponding to intervened fraction
    """
    return df.query("%s==1" % event_name).shape[0]/df.shape[0]

def intervention_prob(stage_diag, params, event_name):
    """
    Function:
        Identifies probability of an intervention in a given stage
    Args:
        stage_diag: CKD stage (in ["G3a", "G3b", "G4", "G5"])
        event_name: name of intervention (in ["nephr", "diag"])
        params: the parameters dictionary
    Returns:
        A float value corresponding to intervention probability
    """
    interv_params = params["model"]["interv_params"]
    return interv_params[event_name]["freq"][stage_diag]
        
def conditional_prob(interv_prob, frac_intervened):
    """
    Function:
        Calculates the probability of being assigned a new intervention
        in order to reach the overall desired intervention probability,
        given that a fraction of individuals already has an intervention
    Args:
        interv_prob: overall desired intervention probability
        frac_intervened: fraction of individuals that already has an intervention
    Returns:
        A float value corresponding to the conditional intervention probability
    """
    return (interv_prob-frac_intervened)/(1-frac_intervened)

def sample_interventions(df, event_name, params, common_rn = None):
    """
    Function:
        Sample and assign intervention indicators to individuals in
        a trajectory dataframe. The dataframe must already come with 
        entries corresponding to eGFR values and ages at which an 
        intervention could be applied
    Args:
        df: a trajectory dataframe
        event_name: name of intervention (in ["nephr", "diag"])
        params: the parameters dictionary
        common_rn: pre-sampled common random number value for random sampling
    Returns:
        trajectory dataframe with updated intervention indicators for the 
        selected intervention
    """
    stage_diags_in_frame = df.stage_diag.unique()

    assert (stage_diags_in_frame.shape[0] == 1), "Intervention sampling for \
        frames with multiple distinct values of stage_diag not implemented"

    frac_intervened = frac_already_intervened(df, event_name)
    interv_prob = intervention_prob(stage_diags_in_frame[0], params, event_name)
    
    prob = conditional_prob(
        interv_prob = interv_prob,
        frac_intervened = frac_intervened,
        )
    
    if event_name == "nephr":
        # adjust probability to reflect sampling only from among diagnosed
        prob = prob/frac_already_intervened(df, "diag")

    if common_rn is not None:
        samples = (common_rn <= prob)*1
    else:
        samples = [np.random.binomial(1, prob) for i in range(df.shape[0])]

    final_df = df.assign(
        **{event_name: samples}
        )
        
    if event_name == "nephr":
        # sampling only among those with diagnosis
        # zero out probability for those not diagnosed
        final_df = final_df.assign(nephr = lambda x: x.nephr * x.diag)

    return final_df

def get_slope_reduction(slope, interv_nephr, interv_diag, interv_params):
    """
    Function:
        Get speed of progression adjusted for the multiplicative reduction
        associated with applied interventions
    Args:
        slope: a float value of yearly eGFR reduction speed
        interv_nephr: a boolean indicator (nephrology intervention)
        interv_diag: a boolean indicator (diagnosis intervention)
        interv_params: parameters dictionary containting multiplicative reductions in 
            progression speed
    Returns:
        an adjusted slope
    """
    if interv_nephr == 1:
        reduction = max(interv_params["nephr"]["eff"], interv_params["diag"]["eff"])
    elif interv_diag == 1:
        reduction = interv_params["diag"]["eff"]
    else:
        return slope
    
    return slope * (1-reduction)
