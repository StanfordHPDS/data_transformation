## Introduction

This repository contains code for the paper titled "A microsimulation-based framework for mitigating societal bias in primary care data" by Agata Foryciarz, Fernando Alarid-Escudero, Gabriela Basel, Marika Cusick, Robert L. Phillips, Andrew Bazemore, Alyce S. Adams and Sherri Rose.

The code generates eGFR trajectories across a person's lifetime, conditional on comorbidities, age, CKD stage and sex.

## Environment setup

Using this repo requires a local [anaconda installation](https://www.anaconda.com/docs/getting-started/anaconda/install).

```         
git clone https://github.com/StanfordHPDS/data_transformation
cd data_transformation
conda update -n base -c defaults conda

# set up and activate a conda environment
conda env create --name data_transformation -f environment.yml

conda activate data_transformation

# install this repo as a package
pip install -e .
```

## Cohort extraction

All SQL scripts necessary for extracting the cohort from an OMOP BigQuery dataset are provided under `egfr_microsim/cohort`. Running `python egfr_microsim/cohort/run_all.py` generates all calibration targets and cohort summary statistics necessary for running the model and reported in the manuscript. Cell counts below 11 are suppressed and converted to 11, so all summary data saved to CSV files by `run_all.py` can be exported to a non-secure environment.

Results of cohort extraction were included under

```
cd data_transformation

egfr_microsim/cohort/calibration_targets
egfr_microsim/cohort/data
```

## Experiment

The entire experiment can be reproduced by running 
```
egfr_microsim/job_scripts/full_100k.sh
```

Given the size of the experiment, it is best ran in a distributed way:
```
egfr_microsim/job_scripts/full_100k/run_experiment.sh
egfr_microsim/job_scripts/full_100k/run_analyze.sh
egfr_microsim/job_scripts/full_100k/gen_figures.sh
```

Experiment results will be saved in 

```
egfr_microsim/exp_results/full_100k
```