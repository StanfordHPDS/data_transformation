CREATE OR REPLACE TABLE {destination_dataset}.extraction_flowchart AS (
    
SELECT 
"full_cohort" as variable,
COUNT(DISTINCT(person_id)) AS n
 from {source_dataset}.person

 UNION ALL

 SELECT 
"in_cohort" as variable,
COUNT(DISTINCT(person_id)) AS n
 from {destination_dataset}.in_cohort

 UNION ALL

 SELECT 
"creatinine" as variable,
COUNT(DISTINCT(person_id)) AS n
 from {destination_dataset}.creatinine

 UNION ALL

 SELECT 
"eGFR_generated" as variable,
COUNT(DISTINCT(person_id)) AS n
 from {destination_dataset}.eGFR_generated

UNION ALL

 SELECT 
"in_cohort_eGFR" as variable,
COUNT(DISTINCT(person_id)) AS n
 from {destination_dataset}.in_cohort_eGFR

UNION ALL

 SELECT 
"cohort_2017_2018" as variable,
COUNT(DISTINCT(person_id)) AS n
 from {destination_dataset}.cohort_2017_2018

ORDER BY n DESC
)