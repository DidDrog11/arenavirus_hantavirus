# Project ArHa: Impute Non-Detections
# 03_08_impute_nondetections.R
# Purpose: This script creates "imputed non-detection" records for species
# that were present in a study but not detected at a specific site, based on
# a reported sampling effort.

combined_data <- read_rds(here("data", "data_cleaning", "03_07_output.rds"))

# Step 1: Identify records with reported trap effort and sampling dates
# These are the records where we have confidence in "non-detection"
records_with_effort <- combined_data$host %>%
  filter(!str_detect(trap_nights_status, "unspecified|not_reported") & !is.na(trap_nights_status)
         & !is.na(start_date) & !is.na(end_date))

# Step 2: Create a table to store the new non-detection records
all_imputed_records <- tibble()

# Create a master species lookup table from the cleaned host data
master_species_lookup <- combined_data$host %>%
  distinct(resolved_species, resolved_genus, resolved_family, resolved_order, resolved_class, gbif_id, taxon_rank) %>%
  mutate(extracted_name_raw = NA_character_,
         extracted_genus_raw = NA_character_)

# Step 3: Iterate through each study with reported effort
# Iterate through each study with reported effort
for (cur_study_id in unique(records_with_effort$study_id)) {
  
  # Get all records for the current study
  study_records <- records_with_effort %>%
    filter(study_id == cur_study_id)
  
  # Identify all species detected within this study
  detected_species_in_study <- study_records %>%
    drop_na(resolved_species) %>%
    distinct(resolved_species) %>%
    pull(resolved_species)
  
  # Identify all unique sampling sites within this study
  unique_sites <- study_records %>%
    distinct(start_date, end_date, NAME_0, locality_raw, verbatimLocality_raw)
  
  # Iterate through each unique site
  for (i in 1:nrow(unique_sites)) {
    site <- unique_sites[i,]
    
    # Get species detected at this specific site
    species_detected_at_site <- study_records %>%
      filter(
        start_date == site$start_date,
        end_date == site$end_date,
        NAME_0 == site$NAME_0,
        locality_raw == site$locality_raw,
        coalesce(verbatimLocality_raw, "NA") == coalesce(site$verbatimLocality_raw, "NA")
      ) %>%
      drop_na(resolved_species) %>%
      distinct(resolved_species) %>%
      pull(resolved_species)
    
    # Find species that were in the study but not at this site
    species_not_detected <- setdiff(detected_species_in_study, species_detected_at_site)
    
    # If there are non-detections to impute, create the records
    if (length(species_not_detected) > 0) {
      
      site_template <- study_records %>%
        filter(
          start_date == site$start_date &
            end_date == site$end_date &
            NAME_0 == site$NAME_0 &
            coalesce(locality_raw, "NA") == coalesce(site$locality_raw, "NA") &
            coalesce(verbatimLocality_raw, "NA") == coalesce(site$verbatimLocality_raw, "NA")
        ) %>%
        head(1) %>%
        select(
          study_id, country_raw,
          country_processed,
          NAME_0,
          iso3_code_processed,
          GID_0,
          country_match_status,
          dist_from_expected,
          decimalLatitude,
          decimalLongitude,
          coord_status,
          locality_raw,
          verbatimLocality_raw,
          coordinate_resolution_raw,
          coordinate_resolution_processed,
          adm_level_num,
          adm1_name,
          adm1_id,
          adm2_name,
          adm2_id,
          adm3_name,
          adm3_id,
          event_date_raw,
          temporal_resolution,
          start_date,
          end_date,
          trap_nights_status,
          trap_nights_clean,
          trap_effort_raw,
          trap_effort_resolution_clean,
          trap_effort_resolution_raw)
      
      new_records <- tibble(resolved_species = species_not_detected) %>%
        mutate(
          study_id = site_template$study_id,
          individual_count = 0,
          detection_status = "imputed_non_detection",
          country_raw = site_template$country_raw,
          country_processed = site_template$country_processed,
          NAME_0 = site_template$NAME_0,
          iso3_code_processed = site_template$iso3_code_processed,
          GID_0 = site_template$GID_0,
          country_match_status = site_template$country_match_status,
          dist_from_expected = site_template$dist_from_expected,
          decimalLatitude = site_template$decimalLatitude,
          decimalLongitude = site_template$decimalLongitude,
          coord_status = site_template$coord_status,
          locality_raw = site_template$locality_raw,
          verbatimLocality_raw = site_template$verbatimLocality_raw,
          coordinate_resolution_raw = site_template$coordinate_resolution_raw,
          coordinate_resolution_processed = site_template$coordinate_resolution_processed,
          adm_level_num = site_template$adm_level_num,
          adm1_name = site_template$adm1_name,
          adm1_id = site_template$adm1_id,
          adm2_name = site_template$adm2_name,
          adm2_id = site_template$adm2_id,
          adm3_name = site_template$adm3_name,
          adm3_id = site_template$adm3_id,
          event_date_raw = site_template$event_date_raw,
          temporal_resolution = site_template$temporal_resolution,
          start_date = site_template$start_date,
          end_date = site_template$end_date,
          trap_nights_status = site_template$trap_nights_status,
          trap_nights_clean = site_template$trap_nights_clean,
          trap_effort_raw = site_template$trap_effort_raw,
          trap_effort_resolution_clean = site_template$trap_effort_resolution_clean,
          trap_effort_resolution_raw = site_template$trap_effort_resolution_raw
        )
      
      all_imputed_records <- bind_rows(all_imputed_records, new_records)
    }
  }
}

# Step 4: Combine the original data with the new imputed records
imputed_non_detections_df <- bind_rows(all_imputed_records) %>%
  left_join(master_species_lookup, by = "resolved_species") %>%
  mutate(rodent_record_id = fct(paste0("imputed_", row_number())),
         detection_status = "imputed")

# Finalize the original data by adding a detection status flag
cleaned_host_final <- combined_data$host %>%
  mutate(detection_status = "reported") %>%
  # Combine with the imputed data
  bind_rows(imputed_non_detections_df) %>%
  mutate(detection_status = fct(detection_status, levels = c("reported", "imputed"))) %>%
  group_by(study_id, start_date, end_date, NAME_0, locality_raw, verbatimLocality_raw, decimalLatitude, decimalLongitude) %>%
  mutate(site_id = cur_group_id()) %>%
  ungroup() %>%
  arrange(study_id, site_id, resolved_species, resolved_genus, resolved_family) %>%
  mutate(rodent_record_id = fct_inorder(rodent_record_id)) %>%
  arrange(rodent_record_id) %>%
  select(-site_id) %>%
  relocate("detection_status", .before = individual_count) %>%
  relocate(contains("trap_"), .after = individual_count)
    
# --- Finalize and Save ---
combined_data$host <- cleaned_host_final
write_rds(combined_data, here::here("data", "data_cleaning", "03_08_output.rds"))
