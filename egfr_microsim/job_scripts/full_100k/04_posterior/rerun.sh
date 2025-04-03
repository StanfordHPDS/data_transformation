#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH -p sherrir
# array size: 1
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

python egfr_microsim/analysis/rerun.py -e $exp_name --cohort_number 100

echo "completed rerun"
