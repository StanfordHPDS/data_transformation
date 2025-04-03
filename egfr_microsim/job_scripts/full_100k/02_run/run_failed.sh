#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=400M
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL

# array size: # of failed jobs (start at 1)

exp_name=$1
egfr_microsim_path=$2
missing_filename=$3
source $HOME/.bashrc
conda activate egfr_microsim
cd $egfr_microsim_path

#missing_filename="missing_experiments"
json_file='egfr_microsim/job_scripts/'$exp_name'/exp_params.json' 
param_file='job_scripts/'$exp_name'/params_base.json' 

equation=$(jq -r '.interv_under_eq' $json_file)
exp_name=$(jq -r '.exp_name' $json_file)
#N=$(jq -r '.N' $json_file)
N=$(jq -r '.cohort_size' $json_file) 

echo $(sed -n "$SLURM_ARRAY_TASK_ID p" egfr_microsim/outputs/$exp_name/sum_stats/${missing_filename}.txt)
eval $(sed -n "$SLURM_ARRAY_TASK_ID p" egfr_microsim/outputs/$exp_name/sum_stats/${missing_filename}.txt)
