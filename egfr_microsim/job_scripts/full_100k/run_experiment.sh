#!/bin/bash
group_home="__________"
egfr_microsim_path=$group_home/egfr_microsim
cd $egfr_microsim_path
exp_name="full_100k"
exp_dir="egfr_microsim/job_scripts/"$exp_name 
output_dir=$group_home/slurm_outputs/$exp_name/

## PART 1 - sampling
# parameter sampling
jid_lhs=$(sbatch --job-name=lhs --output=$output_dir/lhs/%j.out $exp_dir/01_sampling/presample_lhs.sh "$exp_name" "$egfr_microsim_path")
jid_lhs=$(echo "$jid_lhs" | awk '{print $4}')

jid_par=$(sbatch --job-name=params --output $output_dir/params/%j.out --dependency=afterok:$jid_lhs $exp_dir/01_sampling/sample_params.sh "$exp_name" "$egfr_microsim_path")
jid_par=$(echo "$jid_par" | awk '{print $4}')

# cohort sampling
sbatch --job-name=cohorts --output=$output_dir/cohorts/%j.out $exp_dir/01_sampling/sample_cohorts.sh "$exp_name" "$egfr_microsim_path"

## PART 2 - running the experiment matrix
# total number of arrays should equal M*R/1000
# each job runs 10 batches
batch_per_job=10
jid_run1=$(sbatch --array 0-999 --dependency=afterok:$jid_par --job-name="run3" --time=40:00:00 --output=$output_dir/run3/%j.out --error=$output_dir/run3/%j.err $exp_dir/02_run/run_calib_mult.sh "$exp_name" "$egfr_microsim_path" $batch_per_job)
jid_run1=$(echo "$jid_run1" | awk '{print $4}')

# run failed jobs
missing_filename="missing_experiments"
n_missing=$(cat $egfr_microsim_path/egfr_microsim/outputs/$exp_name/sum_stats/$missing_filename.txt | wc -l)
jid_failed=$(sbatch --array 1-$((n_missing+1)) --job-name="failed" --output=$output_dir/failed/%j.out --error=$output_dir/failed/%j.err $exp_dir/02_run/run_failed.sh $exp_name $egfr_microsim_path $missing_filename)
echo $jid_failed