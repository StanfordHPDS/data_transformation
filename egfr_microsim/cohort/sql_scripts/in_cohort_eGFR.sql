DECLARE filter_date DATE;
SET filter_date = DATE("2017-01-01");

CREATE OR REPLACE TABLE {destination_dataset}.in_cohort_eGFR AS (

WITH spanning_1_year AS (

SELECT
    person_id,
    DATE_DIFF(
        MAX(measurement_date), 
        MIN(case when measurement_date >= filter_date then measurement_date end),
        DAY
    ) as obs_period_days_since_t1
    FROM som-nero-phi-sherrir-afc.agataf_omop_0225.eGFR_gen_incl_noacute
    GROUP BY person_id
),
    
selected_patients AS (
        SELECT 
            DISTINCT(person_id)
        FROM spanning_1_year
    -- require creatinine observations spanning >= 1 year
    -- following inclusion date
        WHERE obs_period_days_since_t1/365.25 >= 1
    )

    SELECT 
        b.*
    FROM selected_patients AS a
    INNER JOIN som-nero-phi-drehkopf-afc.agataf_omop.in_cohort AS b
    ON a.person_id=b.person_id 
)

