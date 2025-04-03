WITH table1 AS (
  SELECT 
    person_id,
    egfr,
    g_stage,
    IFNULL(diabetes, 0) AS diabetes,
    IFNULL(hypertension, 0) AS hypertension,
    sex,
    age,
  FROM {destination_dataset}.cohort_2017_2018
  WHERE age < 35 AND age >= 25
)

SELECT 
  "all" AS stratum,
  0 AS value,
  ROUND(AVG(egfr),2) AS egfr,
  ROUND(STDDEV(egfr),2) AS egfr_std,
  ROUND(AVG(age),2) AS mean_age,
  COUNT(*) AS n
FROM table1

UNION ALL 
    
SELECT 
  "diabetes" AS stratum,
CAST(diabetes AS STRING) AS value,
  ROUND(AVG(egfr),2) AS egfr,
  ROUND(STDDEV(egfr),2) AS egfr_std,
  ROUND(AVG(age),2) AS mean_age,
  COUNT(*) AS n
FROM table1
GROUP BY diabetes


UNION ALL 

SELECT 
  "hypertension" AS stratum,
  CAST(hypertension AS STRING) AS value,
  ROUND(AVG(egfr),2) AS egfr,
  ROUND(STDDEV(egfr),2) AS egfr_std,
  ROUND(AVG(age),2) AS mean_age,
  COUNT(*) AS n
FROM table1
GROUP BY hypertension


UNION ALL 

SELECT 
  "sex" AS stratum,
  sex AS value,
  ROUND(AVG(egfr),2) AS egfr,
  ROUND(STDDEV(egfr),2) AS egfr_std,
  ROUND(AVG(age),2) AS mean_age,
  COUNT(*) AS n
FROM table1
GROUP BY sex
ORDER BY stratum, value