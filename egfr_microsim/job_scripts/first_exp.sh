# cohort_number
M=11
# cohort_size
N=10000 
# parameter_set_number
R=100 
exp_name="calibration_experiment1"
## Which equation should be used to assign interventions
eq="09"
# number of available cores
C=10
multiplier=16
cov_targets=("calib_targ_DM.csv" "calib_targ_HT.csv" "calib_targ_sex.csv" "calib_targ_DM_combined.csv" "calib_targ_HT_combined.csv" "calib_targ_ICE.csv" "calib_targ_SDI.csv")

python egfr_microsim/model/presample_lhs.py -e $exp_name --parameter_set_number $R --multiplier $multiplier
python egfr_microsim/model/sample_param_sets.py -e $exp_name --parameter_set_number $R --overwrite_existing
python egfr_microsim/model/sample_cohort.py -e $exp_name --cohort_size $N --cohort_number $M

for i in $(seq 0 $((M - 1)));
do
    python egfr_microsim/calibration/calibration_batch.py --experiment_string $exp_name --cohort_id $i --parameter_set_number $R --cohort_size $N --interv_under_eq $eq --n_cores $C 
done

# aggregate results from calibration using 
python egfr_microsim/calibration/aggregate.py -e $exp_name --parameter_set_number $R --cohort_number $M --n_cores $C 

# Calculate coverage analysis

python egfr_microsim/calibration/coverage_analysis.py -e $exp_name --parameter_set_number $R --cohort_number $M --calib_targets $cov_targets

# Calculate a posterior using sample importance resampling
python egfr_microsim/calibration/posterior_likelihoods.py -e $exp_name --parameter_set_number $R --calib_targets $cov_targets
python egfr_microsim/calibration/posterior_SIR.py -e $exp_name --parameter_set_number $R --calib_targets calib_targ_DM.csv calib_targ_HT.csv calib_targ_sex.csv

