# Project ArHa: Clean Rodent Spatial Data
# 03_04_clean_coords.R
# Purpose: This script performs initial cleaning and validation of rodent spatial data.
# It checks for illogical coordinates, standardizes country names, and validates that
# coordinates fall within the reported country, using the full-resolution GADM 4.1
# outline of that country. Points outside it that were checked against the source
# paper and kept are listed in data/corrections/coordinate_reviews.csv.

combined_data <- read_rds(here("data", "data_cleaning", "03_03_output.rds"))

# --- Helper Functions ---

# Key for matching a record to coordinate_reviews.csv: study, stated country and coordinate to 5 dp
coordinate_key <- function(study_id, iso3c, lat, lon) {
  paste(study_id, iso3c, sprintf("%.5f", as.numeric(lat)), sprintf("%.5f", as.numeric(lon)))
}

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
  
  # B. Validate single-country studies against the full-resolution GADM 4.1 outline of the
  # stated country (the files 03_05 uses). The coarse world map put many coastal and island
  # points outside their country; it is now used only to name the country a point falls in.
  single_countries_checked <- single_countries_subset |>
    filter(!is.na(spatial_id)) |>
    mutate(inside_country = NA, dist_from_expected_border = NA_real_)

  for (iso in unique(single_countries_checked$iso3_code_processed)) {
    country_outline <- gadm(country = iso, level = 0, path = here("data"), version = "4.1")
    if (is.null(country_outline)) stop("No GADM 4.1 outline for ", iso)
    rows <- which(single_countries_checked$iso3_code_processed == iso)
    points <- vect(single_countries_checked[rows, ], geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    single_countries_checked$inside_country[rows] <- relate(points, country_outline, "intersects")[, 1]

    # C. Distance to the stated country's border, for points outside it
    for (k in which(!single_countries_checked$inside_country[rows])) {
      single_countries_checked$dist_from_expected_border[rows[k]] <- min(distance(points[k], country_outline, unit = "km"))
    }
  }

  # The country each outside point falls in on the world map (NA at sea)
  outside <- which(!single_countries_checked$inside_country)
  single_countries_checked$point_GID_0 <- NA_character_
  single_countries_checked$point_NAME_0 <- NA_character_
  if (length(outside) > 0) {
    in_world <- terra::extract(world_shapefile, vect(single_countries_checked[outside, ], geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326"))
    stopifnot(nrow(in_world) == length(outside))
    single_countries_checked$point_GID_0[outside] <- in_world$GID_0
    single_countries_checked$point_NAME_0[outside] <- in_world$NAME_0
  }

  single_countries_checked <- single_countries_checked |>
    mutate(dist_from_expected = case_when(dist_from_expected_border <= 50 ~ "less than 50km",
                                          dist_from_expected_border > 50 ~ "greater than 50km"))

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
    single_countries_checked |> select(spatial_id, inside_country, dist_from_expected_border, dist_from_expected, point_GID_0, point_NAME_0),
    multiple_countries_checked %>% mutate(dist_from_expected = NA_character_) %>% select(spatial_id, country_match_status, dist_from_expected)
  ) %>%
    distinct(spatial_id, .keep_all = TRUE)

  # Step 5: Join the results back to the original full data frame and apply the reviews.
  # A point outside its stated country is "reviewed" if it is a confirmed, accepted row of
  # coordinate_reviews.csv, otherwise "mismatch". Only unreviewed points more than 50 km out
  # take the country they fall in as GID_0/NAME_0, which sets their admin units in 03_05.
  final_df <- cleaned_df %>%
    left_join(unique_coords_final, by = c("spatial_id" = "spatial_id")) |>
    mutate(
      review_key = coordinate_key(study_id, iso3_code_processed, decimalLatitude, decimalLongitude),
      country_match_status = case_when(
        !is.na(country_match_status) ~ as.character(country_match_status), # multi-country records
        is.na(inside_country) ~ NA_character_,
        inside_country ~ "match",
        review_key %in% accepted_reviews$review_key ~ "reviewed",
        TRUE ~ "mismatch"),
      far_unreviewed = country_match_status == "mismatch" & dist_from_expected_border > 50 & !is.na(point_GID_0),
      GID_0 = case_when(is.na(iso3_code_processed) ~ NA_character_,
                        far_unreviewed ~ point_GID_0,
                        TRUE ~ iso3_code_processed),
      NAME_0 = case_when(is.na(iso3_code_processed) ~ NA_character_,
                         far_unreviewed ~ point_NAME_0,
                         TRUE ~ coalesce(world_shapefile$NAME_0[match(iso3_code_processed, world_shapefile$GID_0)], country_processed))) %>%
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


# --- Coordinate Reviews ---
# Points outside their stated country that were checked against the source paper. Rows with
# outcome "accepted" and confirmed = TRUE are marked as reviewed instead of being flagged again;
# "corrected" rows record points fixed in field_corrections.csv and are not used here.
coordinate_reviews <- read_csv(here("data", "corrections", "coordinate_reviews.csv"), col_types = cols(.default = col_character())) |>
  mutate(review_key = coordinate_key(study_id, iso3c, decimalLatitude, decimalLongitude))
accepted_reviews <- coordinate_reviews |> filter(outcome == "accepted", confirmed == "TRUE")

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

# Points outside their stated country that are not on the confirmed review list: check these
# against the source paper, then correct them or add them to coordinate_reviews.csv. Repeat
# extractions that 07 drops (not canonical in 02_01) are left out, as they are never released.
noncanonical_ids <- combined_data$extraction_status |> filter(!is_canonical) |> pull(study_id)
unreviewed_points <- combined_data$host |>
  filter(country_match_status == "mismatch", !is.na(iso3_code_processed), !study_id %in% noncanonical_ids) |>
  count(study_id, country_processed, decimalLatitude, decimalLongitude, dist_from_expected, name = "n_records")

# Unconfirmed review rows that match a point outside its country (dry run), and confirmed rows
# that no longer match one (the coordinate or country has since changed)
outside_keys <- combined_data$host |>
  filter(country_match_status %in% c("reviewed", "mismatch")) |>
  mutate(review_key = coordinate_key(study_id, iso3_code_processed, decimalLatitude, decimalLongitude)) |>
  pull(review_key) |> unique()
pending_reviews <- coordinate_reviews |> filter(outcome == "accepted", confirmed != "TRUE", review_key %in% outside_keys)
stale_reviews <- accepted_reviews |> filter(!review_key %in% outside_keys)

cat("\n--- Spatial Cleaning Summary Report ---\n")
cat(paste("Total records processed:", nrow(combined_data$host), "\n\n"))
cat("Coordinate Status:\n")
print(summary_spatial)
cat("\nCountry check (GADM 4.1):\n")
print(count(combined_data$host, country_match_status))
cat("\nPoints outside their stated country and not reviewed:\n")
print(unreviewed_points, n = Inf)
cat("\nDRY RUN: unconfirmed reviews matching a point outside its country:\n")
print(pending_reviews |> select(study_id, iso3c, decimalLatitude, decimalLongitude, reason), n = Inf)
if (nrow(stale_reviews) > 0) {
  warning(nrow(stale_reviews), " confirmed coordinate reviews no longer match a point outside its country: ",
          paste(stale_reviews$review_key, collapse = "; "), call. = FALSE)
}
