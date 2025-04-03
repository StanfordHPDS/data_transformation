CREATE OR REPLACE TABLE {destination_dataset}.in_cohort AS (

    WITH people AS (
        SELECT 
            person.person_id,
            person.care_site_id,
            CASE person.gender_concept_id
                WHEN 8532 THEN 'F'
                WHEN 8507 THEN 'M'
                ELSE NULL
                END
                AS gender_string,
            person.birth_datetime,
            DATETIME_DIFF('2017-01-01', person.birth_datetime, YEAR) AS age_at_inclusion,
            person.race_concept_id,
            person.ethnicity_concept_id,
            person.location_id,
            person.provider_id,
            obs_period.observation_period_start_date,
            obs_period.observation_period_end_date,
            DATE_DIFF(
            obs_period.observation_period_end_date,
            obs_period.observation_period_start_date, DAY
            ) AS t_obs_day,
            DATE_DIFF(
            obs_period.observation_period_end_date,
            '2017-01-01', DAY
            ) AS t_obs_day_post_2017
        FROM {source_dataset}.person as person
        JOIN {source_dataset}.observation_period as obs_period
        ON person.person_id = obs_period.person_id
    )

    SELECT * FROM people
    WHERE age_at_inclusion >= 18
    AND t_obs_day >= 365.25 
    AND t_obs_day_post_2017 >= 365.25 
    AND gender_string in UNNEST(['M', 'F'])

)
