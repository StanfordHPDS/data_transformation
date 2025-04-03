# cohort_number
M=100
# cohort_size
N=10000 
# parameter_set_number
R=100000
exp_name="full_100k_local"
## Which equation should be used to assign interventions
eq="09"
multiplier=20
cov_targets=("calib_targ_DM.csv" "calib_targ_HT.csv" "calib_targ_sex.csv" "calib_targ_DM_combined.csv" "calib_targ_HT_combined.csv" "calib_targ_ICE.csv" "calib_targ_SDI.csv")
n_thousands=$(expr $R / 1000)
n_lhs_sampled=$(expr $n_thousands \* $multiplier)

# sample parameter sets
python egfr_microsim/model/presample_lhs.py -e $exp_name --parameter_set_number $R --multiplier $multiplier
for i in $(seq 0 $(expr $n_thousands - 1));
do
    python egfr_microsim/model/sample_param_sets.py -e $exp_name --parameter_set_number 1000 -i $i --presampled_lhs_file "lhs_sample_${n_lhs_sampled}k.feather" 
done

# sample cohorts
python egfr_microsim/model/sample_cohort.py -e $exp_name --cohort_size $N --cohort_number $M

# run experiment (batches of 1k parameter sets at a time)
for i in $(seq 0 $((M*R/1000 - 1)));
do
    python egfr_microsim/calibration/calibration_batch.py --experiment_string $exp_name --cohort_id $(expr $i % $M) --parameter_set_id $(expr $i / $M) --parameter_set_number 1000 --interv_under_eq $equation --cohort_size $N --params_base_path $param_file
done

# aggregate results from calibration using 
python egfr_microsim/calibration/aggregate.py -e $exp_name --parameter_set_number $R --cohort_number $M 

# Calculate coverage analysis
python egfr_microsim/calibration/coverage_analysis.py -e $exp_name --parameter_set_number $R --cohort_number $M --calib_targets $cov_targets

# Calculate a posterior using sample importance resampling
python egfr_microsim/calibration/posterior_likelihoods.py -e $exp_name --parameter_set_number $R --calib_targets $cov_targets
python egfr_microsim/calibration/posterior_SIR.py -e $exp_name --parameter_set_number $R --calib_targets calib_targ_DM.csv calib_targ_HT.csv calib_targ_sex.csv

