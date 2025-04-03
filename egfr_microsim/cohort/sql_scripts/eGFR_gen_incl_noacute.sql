CREATE OR REPLACE TABLE {destination_dataset}.eGFR_gen_incl_noacute AS (

    WITH eGFR_gen_incl AS (
        SELECT 
          a.*
        FROM {destination_dataset}.eGFR_generated AS a
        INNER JOIN {destination_dataset}.in_cohort_eGFR AS b
        ON a.person_id=b.person_id 
    ),
    
    exclude_eGFR AS (
      SELECT 
        a.person_id, 
        a.measurement_date,
        b.condition_start_date
      FROM eGFR_gen_incl AS a
      INNER JOIN {destination_dataset}.acute_diag AS b
      USING(person_id)
      WHERE (a.measurement_date <=  DATE_ADD(b.condition_start_date, INTERVAL 30 DAY)
      AND a.measurement_date >=  DATE_SUB(b.condition_start_date, INTERVAL 30 DAY))
    )


SELECT
   DISTINCT *
FROM
   eGFR_gen_incl a
WHERE
 a.egfr > 5 AND
   NOT EXISTS (SELECT * FROM exclude_eGFR b
     WHERE
         a.person_id=b.person_id AND a.measurement_date=b.measurement_date)
)

 