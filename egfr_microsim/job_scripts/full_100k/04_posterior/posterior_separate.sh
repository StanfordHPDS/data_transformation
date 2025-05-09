#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --array=5-6 
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1500M
#SBATCH --time=00:40:00
#SBATCH --mail-type=ALL
#SBATCH -p sherrir

# array size: 1
source $HOME/.bashrc
conda activate egfr_microsim
cd /home/groups/sherrir/agataf/egfr_microsim

exp_name=$1
json_file='egfr_microsim/job_scripts/'$exp_name'/exp_params.json'  
param_file='job_scripts/'$exp_name'/params_base.json' 
cov_targets=("calib_targ_DM.csv" "calib_targ_HT.csv" "calib_targ_sex.csv" "calib_targ_DM_combined.csv" "calib_targ_HT_combined.csv" "calib_targ_ICE.csv" "calib_targ_SDI.csv")

M=$(jq -r '.cohort_number' $json_file) 
N=$(jq -r '.cohort_size' $json_file)  
R=$(jq -r '.param_sets' $json_file)
equation=$(jq -r '.interv_under_eq' $json_file)
exp_name=$(jq -r '.exp_name' $json_file)

echo "initiated"
echo "task id" $SLURM_ARRAY_TASK_ID
cov_target=${cov_targets[$SLURM_ARRAY_TASK_ID]}
echo "cov_target" $cov_target

python egfr_microsim/calibration/posterior_likelihoods.py -e $exp_name --parameter_set_number $R --calib_targets $cov_target --params_base_path $param_file

echo "completed likelihoods"


python egfr_microsim/calibration/posterior_likelihoods.py -e $exp_name --parameter_set_number 100000 --calib_targets calib_targ_ICE.csv --params_base_path job_scripts/'$exp_name'/params_base.json