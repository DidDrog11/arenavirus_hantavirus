# Project ArHa: Locality and Areal Taxonomy
# 0_05_clean_locality.R
# Purpose: This script standardizes locality names and enriches all records with
# a full administrative hierarchy (ADM1, ADM2, ADM3) based on their coordinates.
# It also performs a geocoding step for records with missing coordinates.

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

# --- Step 2: Dissolve polygons to create a full administrative hierarchy ---
#
# This step creates a nested list of SpatVector objects for each administrative level.
gadm_hierarchy_spat_list <- list()

# Iterate over each country in the GADM list
# This will download and add them to gadm_heirarchy_spat_list
# If they are already downloaded they will be read in
for (country_code in names(gadm_spat_list)) {
  
  max_level_obj <- gadm_spat_list[[country_code]]
  max_level_num <- max(as.numeric(str_extract(names(gadm_spat_list[[country_code]]), "\\d+")), na.rm = TRUE)
  
  if (max_level_num == 3) {
    # Keep ADM3
    message(paste("Processing ADM3 level for", country_code))
    gadm_adm3 <- max_level_obj
    gadm_hierarchy_spat_list[[country_code]][["adm3"]] <- gadm_adm3
    
    # Dissolve to ADM2
    message(paste("Downloading ADM2 for", country_code))
    gadm_adm2 <- gadm(country = country_code, level = "2", path = gadm_files_path, version = "4.1")
    gadm_hierarchy_spat_list[[country_code]][["adm2"]] <- gadm_adm2
    
    # Dissolve to ADM1
    message(paste("Downloading ADM1 for", country_code))
    gadm_adm1 <- gadm(country = country_code, level = "1", path = gadm_files_path, version = "4.1")
    gadm_hierarchy_spat_list[[country_code]][["adm1"]] <- gadm_adm1
    
    # Dissolve to ADM0
    message(paste("Downloading ADM0 for", country_code))
    gadm_adm0 <- gadm(country = country_code, level = "0", path = gadm_files_path, version = "4.1")
    gadm_hierarchy_spat_list[[country_code]][["adm0"]] <- gadm_adm0
    
  } else if (max_level_num == 2) {
    # Keep ADM2
    message(paste("Processing ADM2 level for", country_code))
    gadm_adm2 <-  max_level_obj
    gadm_hierarchy_spat_list[[country_code]][["adm2"]] <- gadm_adm2
    
    # Dissolve to ADM1
    message(paste("Downloading ADM1 for", country_code))
    gadm_adm1 <- gadm(country = country_code, level = "1", path = gadm_files_path, version = "4.1")
    gadm_hierarchy_spat_list[[country_code]][["adm1"]] <- gadm_adm1
    
    # Dissolve to ADM0
    message(paste("Downloading ADM0 for", country_code))
    gadm_adm0 <- gadm(country = country_code, level = "0", path = gadm_files_path, version = "4.1")
    gadm_hierarchy_spat_list[[country_code]][["adm0"]] <- gadm_adm0
    
  } else if (max_level_num == 1) {
    # Keep ADM1
    message(paste("Processing ADM1 for", country_code))
    gadm_adm1 <-  max_level_obj
    gadm_hierarchy_spat_list[[country_code]][["adm1"]] <- gadm_adm1
    
    # Dissolve to ADM0
    message(paste("Downloading ADM0 for", country_code))
    gadm_adm0 <- gadm(country = country_code, level = "0", path = gadm_files_path, version = "4.1")
    gadm_hierarchy_spat_list[[country_code]][["adm0"]] <- gadm_adm0
    
  } else if (max_level_num == 0) {
    # Keep ADM0
    message(paste("Processing ADM0 level for", country_code))
    gadm_adm0 <-  max_level_obj
    gadm_hierarchy_spat_list[[country_code]][["adm0"]] <- gadm_adm0
  }
}

adm2_list <- map(gadm_hierarchy_spat_list, ~ .x[["adm2"]])
adm2_list_clean <- compact(adm2_list)
gadm_adm2_combined <- vect(adm2_list_clean)
output_path <- here("data", "gadm", "gadm_adm2_combined.shp")
writeVector(gadm_adm2_combined, output_path, overwrite = TRUE)

# Prepare and Perform Spatial Join ----------------------------------------
# Add the adm_level_num column to the spatial_id table
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
      coordinate_resolution_clean == "study_area" ~ 1,
      coordinate_resolution_clean == "country" ~ 0,
      TRUE ~ 0
    )
  ) %>%
  left_join(gadm_levels_obtained, by = "GID_0") %>%
  # Use coalesce to find the correct administrative level
  mutate(
    adm_level_num = if_else(adm_level_num > levels, levels, adm_level_num) # Cap at the highest available level
  ) %>%
  select(-levels) # Remove the temporary column

# Separate records that have coordinates from those that don't
points_with_coords <- spatial_id %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude))

points_without_coords <- spatial_id %>%
  filter(is.na(decimalLatitude) | is.na(decimalLongitude))


# --- Step 3: Perform spatial join for each unique GADM level per country ---

# Create a list of data frames, one for each country/adm_level combination
groups_to_join <- points_with_coords %>%
  filter(adm_level_num > 0) %>%
  group_by(GID_0, adm_level_num) %>%
  group_split()

adm_join_results <- map(groups_to_join, function(group) {
  
  country_code <- unique(group$GID_0)
  adm_level <- unique(group$adm_level_num)
  adm_level_name <- paste0("adm", adm_level)
  
  # Check if the SpatVector object exists in the list
  if (!is.null(gadm_hierarchy_spat_list[[country_code]][[adm_level_name]])) {
    gadm_spat <- gadm_hierarchy_spat_list[[country_code]][[adm_level_name]]
    
    points_vect <- vect(group, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    
    adm_extracted <- terra::extract(gadm_spat, points_vect) %>%
      select(starts_with("GID_"), starts_with("NAME_"))
    
    results <- group %>%
      bind_cols(adm_extracted %>%
                  select(-GID_0)) %>%
      select(any_of(c(
        "spatial_id",
        "study_id",
        "locality_raw",
        "decimalLatitude", "decimalLongitude",
        "NAME_0", "GID_0", "coordinate_resolution_clean", "adm_level_num",
        "adm1_name" = "GID_1", "adm1_id" = "NAME_1",
        "adm2_name" = "GID_2", "adm2_id" = "NAME_2",
        "adm3_name" = "GID_3", "adm3_id" = "NAME_3"
      )))
    
    return(results)
  }
})

adm_joined_combined <- do.call(bind_rows, adm_join_results) %>%
  # Ensure the join key is unique
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
