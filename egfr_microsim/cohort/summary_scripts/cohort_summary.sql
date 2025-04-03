CREATE OR REPLACE TABLE {destination_dataset}.cohort_summary AS (
    WITH demographics AS (
        SELECT 
            ROUND(AVG(age),2) as age_mean,
            ROUND(STDDEV(age), 2) as age_std,
            ROUND(AVG(egfr),2) as egfr_mean,
            COUNTIF(sex='M') AS male,
            COUNTIF(diabetes=1) AS diabetes,
            COUNTIF(hypertension=1) AS hypertension,
            COUNTIF(CKD_code=1) AS CKD_code,
            COUNT(*) AS total_n
        FROM 
    {destination_dataset}.cohort_2017_2018
    )

    SELECT 
    "total_n" as measure,
    "demographics" as type,
    total_n as value
    FROM demographics

    UNION ALL

    SELECT 
    "age_mean" as measure,
    "demographics" as type,
    age_mean as value
    FROM demographics

    UNION ALL

    SELECT 
    "age_std" as measure,
    "demographics" as type,
    age_std as value
    FROM demographics

    UNION ALL

    SELECT 
    "egfr_mean" as measure,
    "labs" as type,
    egfr_mean as value
    FROM demographics

    UNION ALL

    SELECT 
    "male" as measure,
    "demographics" as type,
    male as value
    FROM demographics

    UNION ALL

    SELECT 
    "diabetes" as measure,
    "diagnoses" as type,
    diabetes as value
    FROM demographics

    UNION ALL

    SELECT 
    "hypertension" as measure,
    "diagnoses" as type,
    hypertension as value
    FROM demographics

    UNION ALL

    SELECT 
    "CKD_code" as measure,
    "diagnoses" as type,
    CKD_code as value
    FROM demographics

    UNION ALL 

    SELECT 
    race,
    "race" as type,
    IF(COUNT(*) < 11, 11, COUNT(*)) AS n
    FROM 
    {destination_dataset}.cohort_2017_2018
    GROUP BY race

    UNION ALL
    SELECT 
    ethnicity,
    "ethnicity" as type,
    IF(COUNT(*) < 11, 11, COUNT(*)) AS n
    FROM 
    {destination_dataset}.cohort_2017_2018
    GROUP BY ethnicity

    UNION ALL
    SELECT 
    g_stage,
    "g_stage" as type,
    IF(COUNT(*) < 11, 11, COUNT(*)) AS n
    FROM 
    {destination_dataset}.cohort_2017_2018
    GROUP BY g_stage
)
