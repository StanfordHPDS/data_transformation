DECLARE filter_date DATE;
DECLARE codes ARRAY <STRING>;

SET filter_date = DATE("2017-01-01");

SET codes = ['11041-1', '11042-9', '16188-5', '16189-3', '2160-0', '35203-9', '39955-0', '39956-8', '39957-6', '39958-4', '39959-2', '39960-0', '39961-8', '39962-6', '39963-4', '39964-2', '39965-9', '39966-7', '39967-5', '39968-3', '39969-1', '39970-9', '39971-7', '39972-5', '39973-3', '39974-1', '39975-8', '39976-6', '40248-7', '40249-5', '40250-3', '40251-1', '40252-9', '40253-7', '40254-5', '40255-2', '40256-0', '40257-8', '40258-6', '44784-7', '54052-6', '57811-2', '67764-1', '72271-0', '74256-9'];

CREATE OR REPLACE TABLE {destination_dataset}.creatinine AS (
    WITH omop_concepts AS (
        SELECT 
            concept_id, concept_name, vocabulary_id
        FROM
            {source_dataset}.concept
        WHERE 
            concept_code IN UNNEST(codes)
    ),

    cohort_measurements AS (
        SELECT
            a.person_id, 
            a.measurement_date,
            a.measurement_concept_id,
            a.value_as_number,
            b.care_site_id
        FROM {source_dataset}.measurement AS a
        INNER JOIN {destination_dataset}.in_cohort AS b
        ON a.person_id=b.person_id
    )
    
    SELECT 
        a.person_id, 
        a.measurement_date,
        a.measurement_concept_id,
        a.value_as_number,
        a.care_site_id,
        b.concept_name,
        b.vocabulary_id
    FROM cohort_measurements AS a
    INNER JOIN omop_concepts AS b
    ON a.measurement_concept_id=b.concept_id
    
    WHERE a.value_as_number < 73.8 
    AND a.value_as_number > 0
    AND a.value_as_number IS NOT NULL
    AND measurement_date >= filter_date

)


 