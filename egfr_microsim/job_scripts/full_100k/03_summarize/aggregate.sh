#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH -p sherrir

# array size: R/1000
source $HOME/.bashrc
conda activate egfr_microsim
cd /home/groups/sherrir/agataf/egfr_microsim

exp_name=$1
json_file='egfr_microsim/job_scripts/'$exp_name'/exp_params.json'  
param_file='job_scripts/'$exp_name'/params_base.json' 

M=$(jq -r '.cohort_number' $json_file) 
N=$(jq -r '.cohort_size' $json_file)  
R=$(jq -r '.param_sets' $json_file)
equation=$(jq -r '.interv_under_eq' $json_file)
exp_name=$(jq -r '.exp_name' $json_file)

echo "initiated"
echo "experiment $exp_name"
echo "param_set_number $R"

python egfr_microsim/calibration/aggregate.py --experiment_string $exp_name --parameter_set_id $SLURM_ARRAY_TASK_ID --cohort_number $M --params_base_path $param_file

