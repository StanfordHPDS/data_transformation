CREATE OR REPLACE TABLE som-nero-phi-sherrir-afc.agataf_omop_0225.sdvi AS (
    WITH person_location AS (
    SELECT b.person_id, b.location_id
    FROM agataf_omop_0225.in_cohort_eGFR AS a
    
    LEFT JOIN {source_dataset}.person AS b
    ON a.person_id=b.person_id
    ),
    
    person_ct AS (
    SELECT 
    a.person_id, 
    b.census_block_group, 
    -- the first 11 digits of the census block group are the census track IDs
    SUBSTR(b.census_block_group, 0, 11) AS CT_id
    FROM person_location AS a
    LEFT JOIN som-nero-phi-sherrir-afc.afc0524phiomop.location AS b
    ON a.location_id=b.location_id
    ),

    SDVIs AS (
    SELECT DISTINCT
    CONCAT(fips_state, fips_county, tract) as CT,
    SDI_CT, 
    ICE_INC_WNH_CT
    FROM {source_dataset_non_omop}.GeneratedPatientBaseline
    )
    
    SELECT 
        a.person_id,
        b.CT,
        b.SDI_CT,
        b.ICE_INC_WNH_CT
    FROM person_ct AS a
    LEFT JOIN SDVIs AS b
    ON a.CT_id=b.CT
)