# Project ArHa: Clean Rodent Spatial Data
# 03_04_clean_coords.R
# Purpose: This script performs initial cleaning and validation of rodent spatial data.
# It checks for illogical coordinates, standardizes country names, and performs
# a spatial join to validate that coordinates fall within the reported country.

combined_data <- read_rds(here("data", "data_cleaning", "03_03_output.rds"))

# --- Helper Functions ---

# Custom function to check if a point falls within one of several country boundaries
check_multiple_countries <- function(lat, lon, iso3_list, world_shapefile) {
  # If iso3_list is NA or empty, use coordinates to determine the country
  if (is.null(iso3_list) || all(is.na(iso3_list)) || length(iso3_list) == 0) {
    point <- vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
    extracted <- terra::extract(world_shapefile, point)
    return(list(coords_in_country = !is.na(extracted$GID_0), iso3_code = extracted$GID_0))
  }
  
  # Subset world_shapefile for the composite countries
  composite_boundary <- world_shapefile[world_shapefile$GID_0 %in% iso3_list, ]
  if (nrow(composite_boundary) == 0) return(list(coords_in_country = FALSE, iso3_code = NA_character_))
  
  # Create a single composite polygon
  composite_union <- aggregate(composite_boundary, dissolve = TRUE)
  
  # Create a point object
  point <- vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
  # Check if the point falls within the composite boundary
  in_country <- relate(point, composite_union, "intersects")
  
  return(list(coords_in_country = as.logical(in_country), iso3_code = NA_character_))
}

# Main function to clean and validate spatial data for a single version (v2 or v3)
clean_spatial_data <- function(raw_host_df, world_shapefile) {
  
  # Step 1: Initial logical checks on coordinates
  cleaned_df <- raw_host_df %>%
    mutate(
      decimalLatitude = as.numeric(decimalLatitude),
      decimalLongitude = as.numeric(decimalLongitude),
      coord_status = case_when(
        is.na(decimalLatitude) | is.na(decimalLongitude) ~ "missing",
        decimalLatitude < -90 | decimalLatitude > 90 | decimalLongitude < -180 | decimalLongitude > 180 ~ "illogical",
        TRUE ~ "valid"
      )
    )  %>% 
    # Clean and standardize country names
    mutate(
      country_raw = country,
      country_processed = case_when(
        str_detect(country_raw, "England|Scotland|Wales") ~ "United Kingdom",
        str_detect(country_raw, "California U.S.|Colorado") ~ "USA",
        str_detect(country_raw, "Finald") ~ "Finland",
        str_detect(country_raw, "Gemany") ~ "Germany",
        str_detect(country_raw, "Nambia") ~ "Namibia",
        str_detect(country_raw, "Phillipines") ~ "Philippines",
        str_detect(country_raw, "Russa|Siberia|Karelia") ~ "Russia",
        str_detect(country_raw, "Solvenia") ~ "Slovenia",
        str_detect(country_raw, "Uruaguay") ~ "Uruguay",
        str_detect(country_raw, "Venezuala|Venezeula") ~ "Venezuela",
        TRUE ~ country
      ),
      iso3_code_processed = if_else(
        str_detect(country_processed, ", |Multiple|Czechoslovakia"),
        NA_character_,
        countrycode(country_processed, "country.name", "iso3c")
      )
    )
  
  # Step 2: Create a unique spatial lookup table
  spatial_lookup <- cleaned_df %>%
    filter(coord_status == "valid") %>%
    distinct(country_processed, decimalLatitude, decimalLongitude) %>%
    # Add a new, unique spatial ID for each distinct point
    rowid_to_column("spatial_id")
  
  # Join the lookup table back to the original data
  cleaned_df <- cleaned_df %>%
    left_join(spatial_lookup, by = c("country_processed", "decimalLatitude", "decimalLongitude"))
  
  # Step 3: Perform spatial join on the unique subset
  
  # A. Separate into single and multiple country groups
  single_countries_subset <- cleaned_df %>%
    filter(!is.na(iso3_code_processed)) %>%
    distinct(spatial_id, .keep_all = TRUE)
  
  multiple_countries_subset <- cleaned_df %>%
    filter(is.na(iso3_code_processed)) %>%
    distinct(spatial_id, .keep_all = TRUE) %>%
    mutate(
      iso3_composite = str_split(country_processed, ", |Multiple|Czechoslovakia") %>%
        map(~ countrycode(.x, "country.name", "iso3c")) %>%
        map_chr(~ paste0(na.omit(.x), collapse = ", "))
    )
  
  # B. Validate single-country studies
  single_countries_vect <- vect(single_countries_subset, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
  country_mapping_single <- terra::extract(world_shapefile, single_countries_vect) %>%
    cbind(spatial_id = single_countries_vect$spatial_id)
  
  single_countries_checked <- single_countries_subset %>%
    left_join(as_tibble(country_mapping_single), by = c("spatial_id")) %>%
    mutate(
      country_match_status = fct(case_when(
        is.na(GID_0) ~ "point_not_in_matched_boundary",
        iso3_code_processed == GID_0 ~ "match",
        TRUE ~ "mismatch"
      ), levels = c("match", "mismatch", "point_not_in_matched_boundary")))
  
  # C. Add a distance-to-border check only for mismatched/no-match records
  mismatched_points <- single_countries_checked %>%
    filter(country_match_status %in% c("mismatch", "point_not_in_matched_boundary"))
  
  if (nrow(mismatched_points) > 0) {
    
    #Check if point is in a different country
    mismatched_points_vect <- mismatched_points %>%
      drop_na(decimalLongitude, decimalLatitude) %>%
      vect(geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    
    other_countries <- terra::extract(world_shapefile, mismatched_points_vect) %>%
      cbind(spatial_id = mismatched_points %>% drop_na(decimalLongitude, decimalLatitude) %>% pull(spatial_id))
    
    mismatched_points$dist_from_expected_border <- NA_real_
    
    for(i in 1:nrow(mismatched_points_vect)) {
      
      mismatched_points$dist_from_expected_border[i] = as.numeric(terra::distance(mismatched_points_vect[i, ], world_shapefile[world_shapefile$GID_0 == mismatched_points_vect[i, ]$iso3_code_processed], unit = "km"))
      
    }
    
    mismatched_points <- mismatched_points %>%
      mutate(dist_from_expected = case_when(dist_from_expected_border <= 50 ~ "less than 50km",
                                            dist_from_expected_border >50 ~ "greater than 50km",
                                            TRUE ~ "not calculated"),
             NAME_0 = case_when(dist_from_expected_border <= 50 ~ country_processed,
                                  dist_from_expected_border >50 & !is.na(NAME_0) ~ NAME_0,
                                  TRUE ~ country_processed),
             GID_0 = case_when(dist_from_expected_border <= 50 ~ iso3_code_processed,
                                    dist_from_expected_border >50 & !is.na(GID_0) ~ GID_0,
                                    TRUE ~ iso3_code_processed)) %>%
      select(rodent_record_id, study_id, spatial_id, GID_0, NAME_0, country_match_status, dist_from_expected)
    
    single_countries_checked <- single_countries_checked %>%
      filter(!spatial_id %in% mismatched_points$spatial_id) %>%
      bind_rows(mismatched_points) %>%
      arrange(spatial_id)
      
  } else {
    single_countries_checked <- single_countries_checked %>%
      mutate(dist_to_border_km = NA_real_, border_status = NA_character_)
  } 
  
  # D. Validate multi-country studies using the helper function
  multiple_countries_checked <- multiple_countries_subset %>%
    rowwise() %>%
    mutate(
      check_result = list(check_multiple_countries(decimalLatitude, decimalLongitude, str_split(iso3_composite, ", ") %>% unlist(), world_shapefile))
    ) %>%
    ungroup() %>%
    mutate(
      country_match_status = fct(map_chr(check_result, ~ if_else(.x$coords_in_country, "match", "mismatch")),
                                 levels = c("match", "mismatch", "unspecified_country")),
      country_match_status = if_else(is.na(country_match_status), "unspecified_country", country_match_status)
    ) %>%
    select(-check_result)
  
  # E. Combine the results of all checks
  unique_coords_final <- bind_rows(
    single_countries_checked %>% select(spatial_id, country_match_status, dist_from_expected, GID_0, NAME_0),
    multiple_countries_checked %>% mutate(dist_from_expected = NA_character_) %>% select(spatial_id, country_match_status, dist_from_expected)
  ) %>% 
    distinct(spatial_id, .keep_all = TRUE)
  
  # Step 5: Join the results back to the original full data frame
  final_df <- cleaned_df %>%
    left_join(unique_coords_final, by = c("spatial_id" = "spatial_id")) %>%
    select(
      # Core IDs
      rodent_record_id, study_id, 
      
      # Taxonomic Information
      resolved_species, resolved_genus, resolved_family, resolved_order, resolved_class,
      gbif_id, taxon_rank, 
      extracted_name_raw, extracted_genus_raw,
      
      # Location Data
      country_raw = country, country_processed, NAME_0,
      iso3_code_processed, GID_0,
      country_match_status, dist_from_expected,
      decimalLatitude, decimalLongitude, coord_status,
      locality, verbatimLocality, coordinate_resolution,
      
      # Effort and Count
      individualCount, trapEffort, trapEffortResolution,
      
      # Raw data and other columns
      eventDate,
      
      # Remove temporary ID
      -spatial_id
    ) %>%
    arrange(rodent_record_id)
  
  return(final_df)
}


# --- Apply Cleaning Function to Host Data ---
# Run the function on the combined data set
combined_data$host <- combined_data$host %>% clean_spatial_data(world_shapefile)

# --- Consolidate and Save Final Output ---
write_rds(combined_data, here("data", "data_cleaning", "03_04_output.rds"))

# --- Summary Report ---
#
# Generate a summary of spatial cleaning status
summary_spatial <- combined_data$host %>%
  group_by(coord_status) %>%
  summarise(
    n_records = n(),
    .groups = "drop"
  )

mismatched_records <- combined_data$host %>%
  filter(country_match_status %in% c("mismatch", "no_country_boundary")) %>%
  select(rodent_record_id, study_id, country_raw, country_match_status, dist_from_expected)

cat("\n--- Spatial Cleaning Summary Report ---\n")
cat(paste("Total records processed:", nrow(combined_data$host), "\n\n"))
cat("Coordinate Status:\n")
print(summary_spatial)
cat("\nRecords with mismatched coordinates (mismatch or no_country_boundary):\n")
print(mismatched_records)
