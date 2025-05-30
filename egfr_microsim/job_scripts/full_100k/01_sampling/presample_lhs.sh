#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0
#SBATCH --mem=1G
#SBATCH --time=5:00:00
#SBATCH --mail-type=ALL
#SBATCH -p sherrir

exp_name=$1
egfr_microsim_path=$2
source $HOME/.bashrc
conda activate egfr_microsim
cd $egfr_microsim_path

json_file='egfr_microsim/job_scripts/'$exp_name'/exp_params.json' 

# parameter_set_number
R=$(jq -r '.param_sets' $json_file) 
exp_name=$(jq -r '.exp_name' $json_file)
multiplier=$(jq -r '.multiplier' $json_file)

echo "initiated LHS calculation, multiplier $multiplier"

python egfr_microsim/model/presample_lhs.py -e $exp_name --parameter_set_number $R --multiplier $multiplier

echo "completed"
