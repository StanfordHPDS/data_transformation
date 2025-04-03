import os
from google.cloud import bigquery
import time
import json

with open('gcloud_paths.json', 'r') as file:
    dataset_paths = json.load(file)

os.environ['GCLOUD_PROJECT'] = dataset_paths["gcloud_project"]
client = bigquery.Client()

save_dir = "data"
os.makedirs(save_dir, exist_ok=True)

save_dir_targets = "calibration_targets"
os.makedirs(save_dir_targets, exist_ok=True)

queries = [
    'in_cohort.sql',
    'race_eth.sql', # requires in_cohort to be complete
    'creatinine.sql',
    'acute_diag.sql',
    'eGFR_generated.sql', # requires creatinine to be complete
    'eGFR_gen_incl_noacute.sql', # requires eGFR_generated and acute_diag to be complete
    'in_cohort_eGFR.sql',
    'diabetes.sql',
    'hypertension.sql',
    'CKD_diag.sql',
    'sdvi.sql',
    'cohort_2017_2018.sql',
]

target_queries = [
    'calib_targ_sex.sql',
    'calib_targ_HT.sql',
    'calib_targ_DM.sql',
    'calib_targ_HT_combined.sql',
    'calib_targ_DM_combined.sql',
    'calib_targ_ICE.sql',
    'calib_targ_SDI.sql'
]

summary_queries = [
    'calib_targ_sex.sql',
    'calib_targ_HT.sql',
    'calib_targ_DM.sql',
    'calib_targ_HT_combined.sql',
    'calib_targ_DM_combined.sql',
    'calib_targ_ICE.sql',
    'calib_targ_SDI.sql',
    'cohort_summary.sql',
    'extraction_flowchart.sql'
]

# generate cohort
for el in queries:
    with open(os.path.join(('sql_scripts'), el), 'r') as file:
        query = file.read().format_map(dataset_paths)
        client.query(query)
        # allow job to complete before the next job is ran
        time.sleep(10)

# generate data summaries
for el in summary_queries:
    with open(os.path.join(('summary_scripts'), el), 'r') as file:
        query = file.read().format_map(dataset_paths)
        client.query(query)
        time.sleep(10)

read_script_sql = """
    SELECT * FROM {table_location}
"""

# save data summaries
for el in summary_queries:
    table_name = el.split('.')[0]
    table_location = '.'.join((dataset_paths["destination_dataset"], table_name))
    query = read_script_sql.format_map({"table_location": table_location})
    query_job = client.query(query)
    df = query_job.to_dataframe()
    if el in target_queries:
        df.to_csv(os.path.join(save_dir_targets, "".join((table_name, ".csv"))))
    else:
        df.to_csv(os.path.join(save_dir, "".join((table_name, ".csv"))))