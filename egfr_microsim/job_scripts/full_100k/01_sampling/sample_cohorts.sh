#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH -p sherrir

exp_name=$1
egfr_microsim_path=$2
source $HOME/.bashrc
conda activate egfr_microsim
cd $egfr_microsim_path

json_file='egfr_microsim/job_scripts/'$exp_name'/exp_params.json' 
## cohort_number
M=$(jq -r '.cohort_number' $json_file) 
## cohort_size
N=$(jq -r '.cohort_size' $json_file)  
## parameter_set_number
R=$(jq -r '.param_sets' $json_file) 

exp_name=$(jq -r '.exp_name' $json_file)


echo "initiated"
echo "experiment $exp_name"
echo "param_set_number $R"

python egfr_microsim/model/sample_cohort.py -e $exp_name --cohort_size $N --cohort_number $M