DECLARE diab2_codes STRING;
DECLARE diab1_codes STRING;
DECLARE diab_other_codes STRING;

SET diab2_codes = r'(^250\.[0-9][02]$)|(^E11[\.][0-9]*$)';
SET diab1_codes = r'(^250\.[0-9][13]$)|(^E10[\.][0-9]*$)';
SET diab_other_codes = r'(^357\.23$)|(^62\.0[1-7]+$)|(^E08\.42$)|(^E13\.42$)';

CREATE OR REPLACE TABLE {destination_dataset}.diabetes AS (

  WITH concept_ids_list AS (
    SELECT
      DISTINCT(condition_concept_id),
      2 as diabetes_type
    FROM {source_dataset}.condition_occurrence
    WHERE REGEXP_CONTAINS(condition_source_value, diab2_codes)
    AND condition_concept_id != 0

    UNION ALL

    SELECT
      DISTINCT(condition_concept_id),
      1 as diabetes_type
    FROM {source_dataset}.condition_occurrence
    WHERE REGEXP_CONTAINS(condition_source_value, diab1_codes)
    AND condition_concept_id != 0

    UNION ALL
         
    SELECT
      DISTINCT(condition_concept_id),
      -- 0 as unknown type
      0 as diabetes_type
    FROM {source_dataset}.condition_occurrence
    WHERE REGEXP_CONTAINS(condition_source_value, diab_other_codes)
    AND condition_concept_id != 0
  ),
  
  omop_concepts AS (
    SELECT
      a.concept_id, 
      a.concept_code, 
      a.concept_name, 
      a.vocabulary_id, 
      b.diabetes_type
    FROM {source_dataset}.concept AS a
    INNER JOIN concept_ids_list AS b
    ON a.concept_id=b.condition_concept_id
    WHERE 
    (b.diabetes_type = 2 AND REGEXP_CONTAINS(a.concept_name, r'([Tt]ype 2)|([Tt]ype II)'))
    OR
    (b.diabetes_type = 1 AND REGEXP_CONTAINS(a.concept_name, r'([Tt]ype 1)|([Tt]ype I$)|([Tt]ype I )'))
    OR
    (b.diabetes_type = 0 AND REGEXP_CONTAINS(a.concept_name, r'[Dd]iabet'))    
  ),


  has_diag AS (
    SELECT
      a.*,
      b.concept_name,
      b.vocabulary_id,
      b.diabetes_type
    FROM {source_dataset}.condition_occurrence AS a
    INNER JOIN omop_concepts AS b
    ON a.condition_concept_id=b.concept_id
  )


  SELECT
    a.condition_occurrence_id as occurrence_id,
    a.person_id,
    a.condition_concept_id as concept_id,
    a.diabetes_type,
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

