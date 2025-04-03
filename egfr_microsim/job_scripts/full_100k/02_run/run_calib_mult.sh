#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=400M
#SBATCH --mail-type=ALL

# array size: M*R/1000

exp_name=$1
egfr_microsim_path=$2
n_repeat=$3
source $HOME/.bashrc
conda activate egfr_microsim
cd $egfr_microsim_path

json_file='egfr_microsim/job_scripts/'$exp_name'/exp_params.json' 
param_file='job_scripts/'$exp_name'/params_base.json' 

## cohort_number
M=$(jq -r '.cohort_number' $json_file) 
## cohort_size
N=$(jq -r '.cohort_size' $json_file)  
## parameter_set_number
R=$(jq -r '.param_sets' $json_file)

equation=$(jq -r '.interv_under_eq' $json_file)

exp_name=$(jq -r '.exp_name' $json_file)

i=1
i1=$SLURM_ARRAY_TASK_ID
while [ $i -le $n_repeat ]
do
echo "initiated experiment $exp_name, task id $i1"
echo "cohort id" $(expr $i1 % $M)
echo "param id" $(expr $i1 / $M)
python egfr_microsim/calibration/calibration_batch.py --experiment_string $exp_name --cohort_id $(expr $i1 % $M) --parameter_set_id $(expr $i1 / $M) --parameter_set_number 1000 --interv_under_eq $equation --cohort_size $N --params_base_path $param_file
#echo "Number: $i"
i1=$((i1 + 1000))
i=$((i+1))
done 

# calculate failed jobs at the end of each batch script
# by the time the last one completes, failed jobs 
python egfr_microsim/calibration/failed_jobs.py -e $exp_name -o "missing_experiments" -r $R
