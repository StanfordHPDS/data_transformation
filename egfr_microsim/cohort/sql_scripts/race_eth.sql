CREATE OR REPLACE TABLE {destination_dataset}.race_eth AS (

WITH mapped_categories_1997 AS (
	SELECT
		_source_.concept_id AS concept_id,
		_source_.concept_name AS concept_name,
		t1.concept_id_1 AS concept_id_1,
		t1.concept_id_2 AS concept_id_2,
		t1.relationship_id AS relationship_id
	FROM
		{source_dataset}.concept AS _source_
		-- Step 1 (join):
		INNER JOIN {source_dataset}.concept_relationship AS t1 ON _source_.concept_id = t1.concept_id_1
	WHERE
		-- Step 2 (rowFilter):
		(
			_source_.domain_id = 'Race'
			AND t1.relationship_id = 'Is a'
			AND _source_.vocabulary_id IN ('Race', 'AFC_race', 'AFC_ethnicity')
		)
),

category_labels_1997 AS (
SELECT
  -- Step 2 (rename):
  _source_.concept_name AS name_id_1,
  _source_.concept_id_1 AS concept_id_1,
  _source_.concept_id_2 AS concept_id_2,
  _source_.relationship_id AS relationship_id,
  t1.concept_name AS name_id_2
FROM
  mapped_categories_1997 AS _source_
  -- Step 1 (join):
  LEFT JOIN {source_dataset}.concept AS t1 ON _source_.concept_id_2 = t1.concept_id
),

in_cohort_race_ethnicity AS (
  SELECT 
    b.person_id, 
    b.race_concept_id, 
    b.ethnicity_concept_id,
    b.race_source_value
  FROM {destination_dataset}.in_cohort AS a
  LEFT JOIN {source_dataset}.person AS b USING(person_id)
),

get_concept_ids AS (
  SELECT
    a.person_id,
    a.race_concept_id,
    a.ethnicity_concept_id,
    b.concept_id,
    b.concept_name,
    a.race_source_value
  FROM in_cohort_race_ethnicity AS a
    LEFT JOIN {source_dataset}.concept AS b ON a.race_concept_id=b.concept_id OR a.ethnicity_concept_id=b.concept_id
),

divide_race_ethnicity AS (

SELECT 
    a.person_id,
    a.race_concept_id,
    a.ethnicity_concept_id,
    a.concept_name AS race_concept_name,
    b.concept_name AS ethnicity_concept_name,
    LOWER(a.race_source_value) AS race_source_value_lower
 FROM 
get_concept_ids AS a 
INNER JOIN get_concept_ids AS b USING(person_id)
WHERE a.concept_id=a.race_concept_id AND b.concept_id=b.ethnicity_concept_id
 )

	SELECT
		-- Step 2 (create variables):
		_source_.person_id AS person_id,
		_source_.race_concept_id AS race_concept_id,
		_source_.race_concept_name AS race_concept_name,
		_source_.ethnicity_concept_id AS ethnicity_concept_id,
		_source_.ethnicity_concept_name AS ethnicity_concept_name,
		t1.concept_id_1 AS concept_id_1,
        CASE
			WHEN _source_.race_concept_id IN (38003615, 38003616) THEN 'Middle Eastern or Northern African'
    		WHEN _source_.race_source_value_lower IN ('other', 'others', '2131-1', 'other race') THEN 'Other'
            WHEN t1.name_id_2 IS NULL THEN _source_.race_concept_name
            WHEN name_id_2 = 'No matching concept' THEN 'Unknown'
            ELSE name_id_2
        END as race_5_name,
        CASE
            WHEN _source_.ethnicity_concept_id = 38003563 THEN 'Hispanic or Latino'
			WHEN _source_.race_concept_id IN (38003615, 38003616) THEN 'Middle Eastern or Northern African'
    		WHEN _source_.race_source_value_lower IN ('other', 'others', '2131-1', 'other race') THEN 'Other'
            WHEN t1.name_id_2 IS NULL THEN _source_.race_concept_name
            WHEN name_id_2 = 'No matching concept' THEN 'Unknown'
            ELSE name_id_2
        END as race_omb_name
	FROM
		divide_race_ethnicity AS _source_
		-- Step 1 (join):
		LEFT JOIN category_labels_1997 AS t1 ON _source_.race_concept_id = t1.concept_id_1
) 
