CREATE OR REPLACE TABLE {destination_dataset}.cohort_2017_2018 AS (

WITH distinct_egfr AS (SELECT
  DISTINCT *
FROM {destination_dataset}.eGFR_gen_incl_noacute),

unique_eGFR AS (
  SELECT 
    person_id,
    measurement_date,
    ANY_VALUE(g_stage) as g_stage,
    ANY_VALUE(sex) AS sex,
    ANY_VALUE(age) AS age,
    ANY_VALUE(egfr) AS egfr,
    ANY_VALUE(obs_period_days_since_t1) AS days_to_last_egfr,
    ANY_VALUE(n_total_obs_5yr) AS n_total_obs_5yr
  FROM distinct_egfr
  WHERE measurement_date < DATE('2018-01-01')
  AND measurement_date > DATE('2017-01-01')
  GROUP BY person_id, measurement_date
),

unique_first_eGFR AS (
  SELECT 
    person_id,
    ANY_VALUE(g_stage) as g_stage,
    ANY_VALUE(sex) AS sex,
    ANY_VALUE(age) AS age,
    ANY_VALUE(egfr) AS egfr,
    ANY_VALUE(measurement_date) AS egfr_date,
    ANY_VALUE(days_to_last_egfr) AS days_to_last_egfr,
    ANY_VALUE(n_total_obs_5yr) AS n_total_obs_5yr
  FROM (
    SELECT
      *,
      FIRST_VALUE(measurement_date) OVER (PARTITION BY person_id ORDER BY measurement_date ASC) AS first_measurement_time
  FROM unique_eGFR)
  WHERE measurement_date=first_measurement_time
  GROUP BY person_id
),

first_diabetes AS (
 SELECT
   person_id,
  --  1 AS ever_diabetes,
   IF(MIN(start_date) < DATE('2018-01-01'), 1, 0) AS diabetes,
   MIN(start_date) AS diabetes_start_date
 FROM {destination_dataset}.diabetes
 GROUP BY person_id
),

first_CKD AS (
 SELECT
   person_id,
  -- 1 AS ever_CKD_code,
   IF(MIN(condition_start_datetime) < DATE('2018-01-01'), 1, 0) AS CKD_code,
   MIN(condition_start_datetime) AS CKD_start_date
 FROM {destination_dataset}.CKD_diag
 GROUP BY person_id
),

first_hypertension AS (
 SELECT
   person_id,
  --  1 AS ever_hypertension,
   IF(MIN(start_date) < DATE('2018-01-01'), 1, 0) AS hypertension,
   MIN(start_date) AS hypertension_start_date
 FROM {destination_dataset}.hypertension
 GROUP BY person_id
),

sdi AS (
  SELECT
    person_id,
    ROUND(SDI_CT,2) AS SDI,
    ROUND(ICE_INC_WNH_CT,2) AS ICE,
    NTILE(3) OVER (ORDER BY SDI_CT DESC) AS SDI_q3,
    NTILE(3) OVER (ORDER BY ICE_INC_WNH_CT ASC) AS ICE_q3
  FROM {destination_dataset}.sdvi
  WHERE SDI_CT IS NOT NULL
),

race_eth AS (
  SELECT
    person_id,
    CASE
        WHEN race_5_name = "No matching concept" THEN "Unknown"
        ELSE race_5_name
    END AS race,
    CASE
        WHEN ethnicity_concept_name = "No matching concept" THEN "Unknown"
        ELSE ethnicity_concept_name
    END AS ethnicity,
    race_omb_name AS race_omob
  FROM {destination_dataset}.race_eth
  -- WHERE SDI_CT IS NOT NULL
)

SELECT
  b.person_id,
  b.sex,
  ROUND(b.age, 2) AS age,
  ROUND(b.egfr, 2) AS egfr,
  b.g_stage,
  b.egfr_date,
  b.days_to_last_egfr,
  b.n_total_obs_5yr,
  c.diabetes,
  c.diabetes_start_date,
  d.CKD_code,
  d.CKD_start_date,
  e.hypertension,
  e.hypertension_start_date,
  f.SDI,
  f.SDI_q3,
  f.ICE,
  f.ICE_q3,
  g.race,
  g.ethnicity,
  g.race_omob
  -- b.first_egfr_date,
  -- b.n_obs,
  -- b.obs_period

FROM unique_first_eGFR AS b
LEFT JOIN first_diabetes AS c USING(person_id)
LEFT JOIN first_CKD AS d USING(person_id)
LEFT JOIN first_hypertension AS e USING(person_id)
LEFT JOIN sdi AS f USING(person_id)
LEFT JOIN race_eth AS g USING(person_id)
)