#!/bin/bash
group_home="__________"
egfr_microsim_path=$group_home/egfr_microsim
cd $egfr_microsim_path
exp_name="full_100k"
exp_dir="egfr_microsim/job_scripts/"$exp_name 
output_dir=$group_home/slurm_outputs/$exp_name/

## PART 3 - aggregation and analysis
# total array number should equal R/1000
jid_agg=$(sbatch --job-name=agg --array=0-99 --output=$output_dir/agg/%j.out $exp_dir/03_summarize/aggregate.sh "$exp_name" "$egfr_microsim_path")
jid_agg=$(echo "$jid_agg" | awk '{print $4}')

sbatch --job-name=cov --dependency=afterok:$jid_agg --output=$output_dir/cov/%j.out --error=$output_dir/cov/%j.err $exp_dir/03_summarize/coverage_separate.sh "$exp_name"

# posterior calculation
jid_post=$(sbatch --dependency=afterok:$jid_agg --job-name=post --output=$output_dir/post/%j.out --error=$output_dir/post/%j.err $exp_dir/04_posterior/posterior_separate.sh "$exp_name" "$egfr_microsim_path")
jid_post=$(echo "$jid_post" | awk '{print $4}')

sbatch --job-name=post --output=$output_dir/post1/%j.out --error=$output_dir/post1/%j.err $exp_dir/04_posterior/posterior_separate.sh "$exp_name" "$egfr_microsim_path"

# SIR calculation
jid_sir=$(sbatch --dependency=afterok:$jid_post --job-name=sir --output=$output_dir/sir/%j.out --error=$output_dir/sir/%j.err $exp_dir/04_posterior/posterior_SIR.sh "$exp_name" "$egfr_microsim_path")
jid_sir=$(echo "$jid_sir" | awk '{print $4}')

# rerun
sbatch --dependency=afterok:$jid_sir --job-name=rerun --output=$output_dir/rerun/%j.out --error=$output_dir/rerun/%j.err $exp_dir/04_posterior/rerun.sh $exp_name "$egfr_microsim_path"


sbatch --job-name=rerun4 --output=$output_dir/rerun/%j.out --error=$output_dir/rerun/%j.err $exp_dir/04_posterior/rerun.sh $exp_name "$egfr_microsim_path"

