CREATE OR REPLACE TABLE {destination_dataset}.calib_targ_HT AS (

    WITH filled_with_zeros AS (
      SELECT 
        person_id,
        g_stage,
        IFNULL(hypertension, 0) AS hypertension,
        CASE 
            WHEN (age < 30) THEN 20
            WHEN (age >= 30 AND age < 40) THEN 30
            WHEN (age >= 40 AND age < 50) THEN 40
            WHEN (age >= 50 AND age < 60) THEN 50
            WHEN (age >= 60 AND age < 70) THEN 60
            WHEN (age >= 70 AND age < 80) THEN 70
            WHEN (age >= 80 AND age < 90) THEN 80
            ELSE 90
        END AS age_group,
      FROM {destination_dataset}.cohort_2017_2018
    ),
    
    full_table AS (
      SELECT
        age_group,
        g_stage,
        hypertension,
        IF(COUNT(*) < 11, 11, COUNT(*)) AS n,
        FROM filled_with_zeros
        GROUP BY age_group, g_stage, hypertension
      ),
    
    category_sums AS (
        SELECT
            age_group,
            hypertension,
            COUNT(*) AS group_sum,
        FROM filled_with_zeros
        GROUP BY age_group, hypertension
    )
        
    SELECT 
        *,
        ROUND(a.n/b.group_sum, 5) AS frac
    FROM full_table AS a
    JOIN category_sums AS b
    USING(age_group, hypertension)
    ORDER BY age_group, hypertension, g_stage
    
)
