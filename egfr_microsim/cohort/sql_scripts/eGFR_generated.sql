DECLARE filter_date DATE;
SET filter_date = DATE("2017-01-01");

-- define a javascript function to calculate eGFR
-- scraped from https://www.kidney.org/professionals/kdoqi/gfr_calculator

CREATE TEMP FUNCTION egfr(scr FLOAT64, age FLOAT64, sex STRING)
RETURNS FLOAT64
LANGUAGE js
AS r'''
if ((sex != 'F') & (sex != 'M')) {{
  return null;
}}
var kappa = (sex == "F" ? 0.7 : 0.9);
var alpha = (sex == "F" ? -0.241 : -0.302);
var result_1 = 
  142 *
  Math.pow(Math.min(scr / kappa, 1), alpha) *
  Math.pow(Math.max(scr / kappa, 1), -1.200) *
  Math.pow(0.9938, age) *
  (sex == "F" ? 1.012 : 1)
  return result_1.toFixed(3);
''';

-- define a javascript function to calculate G-stage from eGFR
-- scraped from https://www.kidney.org/professionals/kdoqi/gfr_calculator
CREATE TEMP FUNCTION ckd_stage(gfr FLOAT64)
RETURNS STRING
LANGUAGE js
AS r'''
var result_5 = "";
if(gfr >= 90) result_5 = "G1";
else if (gfr >= 60) result_5 = "G2";
else if (gfr >= 45) result_5 = "G3a";
else if (gfr >= 30) result_5 = "G3b";
else if (gfr >= 15) result_5 = "G4";
else result_5 = "G5";
return result_5;
''';


-- generate a demographics table

CREATE OR REPLACE TABLE {destination_dataset}.eGFR_generated AS (

-- select patients with at least one creatinine record since filter_date
WITH vars AS (
    SELECT DISTINCT(person_id) as patients_seen_since_2017
    FROM {destination_dataset}.creatinine
   -- WHERE measurement_date >= filter_date
),

-- create a new table that will include creatinine calculation results
-- only uncluding individuals with binary sex and non-missing date of birth
creatinine_demographics AS (
    SELECT 
      creat.person_id, 
      creat.measurement_date,
      creat.care_site_id,
      creat.value_as_number as creatinine,
      patients.gender_string as sex,
      patients.birth_datetime,
      -- TODO: now, we're using a float value for age - should we use int?
      DATE_DIFF(creat.measurement_date, patients.birth_datetime, DAY)/365.25 as age
    FROM {destination_dataset}.creatinine as creat 
    LEFT JOIN {destination_dataset}.in_cohort as patients
    ON creat.person_id = patients.person_id
    WHERE creat.person_id IN (SELECT patients_seen_since_2017 FROM vars)
),

-- first, just calculate eGFR
egfr_calc AS (
  SELECT
    *,
    egfr(creatinine, age, sex) as egfr,
  FROM creatinine_demographics
  ),

  -- then, find the date of the first observation for any given patient
  -- t0: first creatinine measurement recorded for patient
  -- t1: first creatinine measurement recorded for patient on or after filter_date
  -- n0: total number of creatinine measurements recorded for patient
  -- n1: total number of creatinine measurements recorded for patient on and after filter_date
  first_patient_obs AS (
  SELECT 
      person_id, 
      MIN(measurement_date) as t0,
      MIN(case when measurement_date >= filter_date then measurement_date end) AS t1,
      COUNT(*) as n0,
      COUNT(case when measurement_date >= filter_date then 1 end) AS n1,
      DATE_DIFF(MAX(measurement_date), MIN(measurement_date), DAY) as obs_period_days,
      DATE_DIFF(MAX(measurement_date), 
               MIN(case when measurement_date >= filter_date then measurement_date end), DAY)
               as obs_period_days_since_t1
  FROM egfr_calc
  GROUP BY person_id
  )

  -- add staging and days since first observation
  SELECT
    a.person_id,
    a.measurement_date,
    a.care_site_id, 
    a.creatinine,
    a.sex,
    a.age,
    b.t0,
    b.t1,
    DATE_DIFF(a.measurement_date,b.t0, DAY)+1 as days_since_t0,
    DATE_DIFF(a.measurement_date,b.t1, DAY)+1 as days_since_t0_5yr,
    ckd_stage(a.egfr) as g_stage,
    b.n0 as n_total_obs,
    b.n1 as n_total_obs_5yr,
    obs_period_days,
    obs_period_days_since_t1,
    egfr

  FROM egfr_calc AS a
  JOIN first_patient_obs as b
  ON a.person_id = b.person_id  
)

