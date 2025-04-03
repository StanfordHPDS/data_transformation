#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH -p sherrir

# array size: R/1000

exp_name=$1
egfr_microsim_path=$2
source $HOME/.bashrc
conda activate egfr_microsim
cd $egfr_microsim_path

json_file='egfr_microsim/job_scripts/'$exp_name'/exp_params.json' 
param_file='job_scripts/'$exp_name'/params_base.json' 

R=$(jq -r '.param_sets' $json_file) 
exp_name=$(jq -r '.exp_name' $json_file)
multiplier=$(jq -r '.multiplier' $json_file)

echo "initiated experiment $exp_name, task id $SLURM_ARRAY_TASK_ID"

n_thousands=$(expr $R / 1000)
n_lhs_sampled=$(expr $n_thousands \* $multiplier)

for i in $(seq 0 $(expr $n_thousands - 1));
do
    python egfr_microsim/model/sample_param_sets.py -e $exp_name --parameter_set_number 1000 -i $i --presampled_lhs_file "lhs_sample_${n_lhs_sampled}k.feather" --params_base_path $param_file 
done