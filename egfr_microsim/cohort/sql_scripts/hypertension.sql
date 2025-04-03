DECLARE hyp_codes STRING;
SET hyp_codes = r'(^40[123]\.[0-9]*$)|(^40[123]$)';

CREATE OR REPLACE TABLE {destination_dataset}.hypertension AS (

  WITH concept_ids_list AS (
    SELECT
      DISTINCT(condition_concept_id),
    FROM {source_dataset}.condition_occurrence
    WHERE REGEXP_CONTAINS(condition_source_value, hyp_codes)
    AND condition_concept_id != 0
  ),
  
  omop_concepts AS (
    SELECT
      a.concept_id, 
      a.concept_code, 
      a.concept_name, 
      a.vocabulary_id
    FROM {source_dataset}.concept AS a
    INNER JOIN concept_ids_list AS b
    ON a.concept_id=b.condition_concept_id
  ),

  has_diag AS (
    SELECT
      a.*,
      b.concept_name,
      b.vocabulary_id
    FROM {source_dataset}.condition_occurrence AS a
    INNER JOIN omop_concepts AS b
    ON a.condition_concept_id=b.concept_id
  )


  SELECT
    a.condition_occurrence_id as occurrence_id,
    a.person_id,
    a.condition_concept_id as concept_id,
    a.concept_name,
    a.condition_source_value,
    'condition_occurence' as source_table,
    a.load_table_id,
    a.condition_start_date as start_date,
    a.condition_end_date as end_date,
    a.vocabulary_id

  FROM has_diag AS a
  INNER JOIN {destination_dataset}.in_cohort_eGFR as b
  ON a.person_id=b.person_id
)