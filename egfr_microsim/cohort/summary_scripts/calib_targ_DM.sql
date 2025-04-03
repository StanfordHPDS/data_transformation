CREATE OR REPLACE TABLE {destination_dataset}.calib_targ_DM AS (

    WITH filled_with_zeros AS (
      SELECT 
        person_id,
        g_stage,
        IFNULL(diabetes, 0) AS diabetes,
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
        diabetes,
        IF(COUNT(*) < 11, 11, COUNT(*)) AS n,
        FROM filled_with_zeros
        GROUP BY age_group, g_stage, diabetes
      ),

    category_sums AS (
        SELECT
            age_group,
            diabetes,
            COUNT(*) AS group_sum,
        FROM filled_with_zeros
        GROUP BY age_group, diabetes
    )
        
    SELECT 
        *,
        ROUND(a.n/b.group_sum, 5) AS frac
    FROM full_table AS a
    JOIN category_sums AS b
    USING(age_group, diabetes)
    ORDER BY age_group, diabetes, g_stage
)
