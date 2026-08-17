# Project ArHa: Locality and Areal Taxonomy
# 03_05_clean_locality.R
# Purpose: This script standardizes locality names and enriches all records with
# a full administrative hierarchy (ADM1, ADM2, ADM3) based on their coordinates.

# Load the combined data from the previous script
combined_data <- read_rds(here("data", "data_cleaning", "03_04_output.rds"))

# Step 1: Harmonise Coordinate Resolution ---------------------------------

cleaning_coord_res <- combined_data$host %>%
  mutate(
    coordinate_resolution_clean = case_when(
      # High-precision / site-level data
      str_detect(tolower(coordinate_resolution), "site|trap|station|field") ~ "site",
      str_detect(tolower(coordinate_resolution), "verbatim|exact|given|stated|reported|madeira") ~ "site",
      str_detect(coordinate_resolution, "\\d+(km|m)") ~ "site",
      str_detect(tolower(coordinate_resolution), "average|estimated|locality|location") ~ "adm3",
      
      # Administrative units
      str_detect(tolower(coordinate_resolution), "municipality|parish|island|prefecture|subdistrict") ~ "adm3",
      str_detect(tolower(coordinate_resolution), "county|district|commune|department|banner") ~ "adm2",
      str_detect(tolower(coordinate_resolution), "state|province|region|krai|oblast|republic|vasterbotten") ~ "adm1",
      str_detect(tolower(coordinate_resolution), "country|continent") ~ "country",
      
      # Natural areas and landmarks
      str_detect(tolower(coordinate_resolution), "mountain|river|lake|wetland|coast|plateau|forest|scenic|valley") ~ "site",
      str_detect(tolower(coordinate_resolution), "national park|national reserve|nature reserve|national wildlife refuge|wildlife management area|wildlife sanctuary|reserve|preserve") ~ "adm2",
      
      # Human settlement types
      str_detect(tolower(coordinate_resolution), "village|settlement|fishing village|community|rural") ~ "village",
      str_detect(tolower(coordinate_resolution), "town|center") ~ "town",
      str_detect(tolower(coordinate_resolution), "city|metropolitan|alexandria|djibouti|urban|airport|neighborhood|neighbourhood|complex|port") ~ "city",
      
      # Study-specific or ambiguous terms
      str_detect(tolower(coordinate_resolution), "study|surveillance|4 corner|basin|area|domain") ~ "study_area",
      str_detect(tolower(coordinate_resolution), "average|center|location|regions|long term ecological monitoring") ~ "unspecified",
      
      # Default catch-all
      TRUE ~ "unspecified"
      ),
    coordinate_resolution_clean = fct(coordinate_resolution_clean, levels =  c("site", "village", "town", "city", "study_area", "adm3", "adm2", "adm1", "country", "unspecified"))) %>%
  rename(coordinate_resolution_raw = coordinate_resolution,
         locality_raw = locality,
         verbatimLocality_raw = verbatimLocality) %>%
  group_by(study_id, locality_raw, decimalLatitude, decimalLongitude, NAME_0, GID_0, coordinate_resolution_clean)

cleaning_coord_res$spatial_id = group_indices(cleaning_coord_res)

spatial_id <- cleaning_coord_res %>%
  distinct(study_id, locality_raw, decimalLatitude, decimalLongitude, NAME_0, GID_0, coordinate_resolution_clean, spatial_id)

# Step 2: Identify and Download GADM Files --------------------------------
# Create a list of unique countries and the max admin level needed for each
gadm_levels_needed <- cleaning_coord_res %>%
  ungroup() %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  mutate(
    adm_level_num = case_when(
      coordinate_resolution_clean == "adm1" ~ 1,
      coordinate_resolution_clean == "adm2" ~ 2,
      coordinate_resolution_clean == "adm3" ~ 3,
      coordinate_resolution_clean == "city" ~ 3,
      coordinate_resolution_clean == "town" ~ 3,
      coordinate_resolution_clean == "village" ~ 3,
      coordinate_resolution_clean == "site" ~ 3,
      coordinate_resolution_clean == "study_area" ~ 3,
      coordinate_resolution_clean == "country" ~ 0,
      TRUE ~ 0
    )
  ) %>%
  group_by(GID_0) %>%
  summarise(max_level = max(adm_level_num, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(GID_0))

# Download the necessary GADM files
gadm_files_path <- here("data")

# Create a list to store the downloaded SpatVectors
gadm_spat_list <- list()

pwalk(gadm_levels_needed, function(GID_0, max_level) {
  country_code <- GID_0
  current_level <- max_level
  
  # Step down the administrative levels until a file is found or we reach level 0
  while (current_level >= 0) {
    message(paste("Attempting to download GADM level", current_level, "for", country_code))
    
    # Try to download the file and store the result
    gadm_result <- try(
      gadm(country = country_code, level = current_level, path = gadm_files_path, version = "4.1"),
      silent = TRUE
    )
    
    # Check if the download was successful
    if (!inherits(gadm_result, "try-error") && !is.null(gadm_result)) {
      message(paste("Successfully downloaded GADM level", current_level, "for", country_code))
      gadm_spat_list[[country_code]] <<- gadm_result 
      break # Exit the loop if successful
    } else {
      message(paste0("GADM level ", current_level, " not found for ", country_code, ". Stepping down."))
      current_level <- current_level - 1
    }
  }
})

# GADM level obtained
gadm_levels_obtained <- tibble(
  GID_0 = names(gadm_spat_list),
  levels = map_int(gadm_spat_list, ~ {
    column_names <- names(.x)
    gid_levels <- str_extract(column_names, "GID_(\\d)")
    max_level <- suppressWarnings(max(as.numeric(str_extract(na.omit(gid_levels), "\\d")), na.rm = TRUE))
    if(is.infinite(max_level)) { # Handle countries with no ADM levels
      0L
    } else {
      max_level
    }
  })
)

# --- Perform spatial join using each country's single GADM file ---
spatial_id <- spatial_id %>%
  mutate(
    adm_level_num = case_when(
      coordinate_resolution_clean == "adm1" ~ 1,
      coordinate_resolution_clean == "adm2" ~ 2,
      coordinate_resolution_clean == "adm3" ~ 3,
      coordinate_resolution_clean == "city" ~ 3,
      coordinate_resolution_clean == "town" ~ 3,
      coordinate_resolution_clean == "village" ~ 3,
      coordinate_resolution_clean == "site" ~ 3,
      coordinate_resolution_clean == "study_area" ~ 3,
      coordinate_resolution_clean == "country" ~ 0,
      TRUE ~ 0
    )
  ) %>%
  left_join(gadm_levels_obtained, by = "GID_0") %>%
  mutate(
    adm_level_num = if_else(adm_level_num > levels, levels, adm_level_num)
  ) %>%
  select(-levels)

points_with_coords <- spatial_id %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

points_without_coords <- spatial_id %>%
  filter(is.na(decimalLatitude) | is.na(decimalLongitude))

# One extract call per country -- not per country/level combination
groups_to_join <- points_with_coords %>%
  filter(adm_level_num > 0) %>%
  group_by(GID_0) %>%
  group_split()

adm_join_results <- map(groups_to_join, function(group) {
  
  country_code <- unique(group$GID_0)
  gadm_spat <- gadm_spat_list[[country_code]]
  if (is.null(gadm_spat)) return(NULL)
  
  points_vect <- vect(group, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  
  adm_extracted <- terra::extract(gadm_spat, points_vect) %>%
    select(starts_with("GID_"), starts_with("NAME_"))
  
  # Ensure all three levels exist as columns even for countries whose file
  # doesn't reach that deep, so bind_cols stays consistent across countries
  for (col in c("GID_1", "NAME_1", "GID_2", "NAME_2", "GID_3", "NAME_3")) {
    if (!col %in% names(adm_extracted)) adm_extracted[[col]] <- NA_character_
  }
  
  group %>%
    bind_cols(adm_extracted %>% select(-any_of(c("GID_0", "NAME_0")))) %>%
    rename(adm1_name = GID_1, adm1_id = NAME_1,
           adm2_name = GID_2, adm2_id = NAME_2,
           adm3_name = GID_3, adm3_id = NAME_3) %>%
    select(any_of(c(
      "spatial_id", "study_id", "locality_raw",
      "decimalLatitude", "decimalLongitude",
      "NAME_0", "GID_0", "coordinate_resolution_clean", "adm_level_num",
      "adm1_name", "adm1_id", "adm2_name", "adm2_id", "adm3_name", "adm3_id"
    )))
})

adm_joined_combined <- do.call(bind_rows, adm_join_results) %>%
  distinct(spatial_id, .keep_all = TRUE)

host_data <- cleaning_coord_res %>%
  ungroup() %>%
  rename("spatial_id_host" = spatial_id) %>%
  left_join(adm_joined_combined %>%
              select("spatial_id_spatial" = spatial_id, adm_level_num, adm1_name, adm1_id, adm2_name, adm2_id, adm3_name, adm3_id), by = c("spatial_id_host" = "spatial_id_spatial")) %>%
  # Final selection of columns and renaming
  select(
    # Core IDs
    rodent_record_id, study_id,
    
    # Taxonomic Information
    resolved_species, resolved_genus, resolved_family, resolved_order, resolved_class,
    gbif_id, taxon_rank,
    extracted_name_raw, extracted_genus_raw,
    
    # Location Data
    country_raw, country_processed, NAME_0, iso3_code_processed, GID_0,
    country_match_status, dist_from_expected,
    decimalLatitude, decimalLongitude, coord_status,
    locality_raw, verbatimLocality_raw, coordinate_resolution_raw,
    coordinate_resolution_processed = coordinate_resolution_clean,
    adm_level_num,
    adm1_name, adm1_id,
    adm2_name, adm2_id,
    adm3_name, adm3_id,
    
    # Effort and Count
    individualCount, trapEffort, trapEffortResolution,
    
    # Raw data and other columns
    eventDate
  ) %>%
  arrange(rodent_record_id)

# --- Save Final Output ---

combined_data$host <- host_data

write_rds(combined_data, here("data", "data_cleaning", "03_05_output.rds"))
