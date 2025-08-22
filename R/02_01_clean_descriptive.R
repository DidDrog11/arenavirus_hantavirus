# Project ArHa: Initial Cleaning of Descriptive Data
# 02_01_clean_descriptives.R
# Purpose: This script consolidates descriptive metadata from v2 and v3,
# and performs initial cleaning and standardization of key fields.
# Output: A single, consolidated data frame named `descriptives_cleaned`,
# ready for use in subsequent analysis and final consolidation.
# 
# Load the raw data from 00_load_data.R
# e.g., combined_data_v2, combined_data_v3
combined_data_v2 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
combined_data_v3 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))

# Descriptives v2 ---------------------------------------------------------

descriptives_v2 <- combined_data_v2$descriptive

descriptives_v2_cleaned <- descriptives_v2 %>%
  
  # Remove unnecessary columns first to simplify subsequent code
  select(-data, -data_extractor) %>%
  
  # Clean and flag identifiedBy (first author)
  rename(author_first_last_raw = identifiedBy) %>%
  mutate(
    author_status = if_else(is.na(author_first_last_raw) | str_trim(author_first_last_raw) == "", "missing", "present")
  ) %>%
  
  # Clean and parse publication_year
  rename(publication_year_raw = publication_year) %>%
  mutate(
    publication_year_clean = as.numeric(publication_year_raw),
    publication_year_status = if_else(
      is.na(publication_year_clean) | publication_year_clean < 1900,
      "invalid_year",
      "valid_year"
    )
  ) %>%
  
  # Clean and parse sampling_effort
  rename(sampling_effort_raw = sampling_effort) %>%
  mutate(
    sampling_effort_raw = str_trim(sampling_effort_raw), # Clean whitespace
    sampling_effort_status = case_when(
      is.na(sampling_effort_raw) | str_detect(tolower(sampling_effort_raw), "not reported|ifa|pcr|^-?$|n/a") ~ "not_reported",
      str_detect(sampling_effort_raw, "trapping nights|trap nights") ~ "trap_nights",
      str_detect(tolower(sampling_effort_raw), "session|sessions") ~ "sessions",
      str_detect(sampling_effort_raw, "month|months|week|weeks|year|years|annual|sample|season|\\d{4}") ~ "time_period",
      TRUE ~ "unspecified"
    ),
    trap_nights = if_else(
      sampling_effort_status == "trap_nights",
      as.numeric(str_remove_all(str_extract(sampling_effort_raw, "\\d+(,\\d+)*(?=\\s*trapping nights|\\s*trap nights)"), ",")),
      NA_real_
    ),
    n_sessions = if_else(
      sampling_effort_status == "sessions",
      as.numeric(str_extract(sampling_effort_raw, "\\d+")),
      NA_real_
    ),
    start_year = if_else(
      sampling_effort_status == "time_period",
      as.numeric(str_extract(sampling_effort_raw, "\\d{4}")),
      NA_real_
    ),
    end_year = if_else(
      sampling_effort_status == "time_period",
      as.numeric(str_extract(sampling_effort_raw, "\\d{4}$")),
      NA_real_
    ),
    sampling_effort_clean = if_else(
      sampling_effort_status == "time_period",
      sampling_effort_raw,
      NA_character_
    )
  ) %>%
  
  # Clean and standardize data_access
  rename(data_access_raw = data_access) %>%
  mutate(
    data_access_level = case_when(
      str_detect(tolower(data_access_raw), "summarised") ~ "summarised",
      str_detect(tolower(data_access_raw), "individual") ~ "individual",
      TRUE ~ "unspecified"
    )
  ) %>%
  
  # Clean and flag linked_manuscripts (as a DOI)
  rename(linked_manuscripts_raw = linked_manuscripts) %>%
  mutate(
    linked_manuscripts_status = case_when(
      is.na(linked_manuscripts_raw) | str_trim(linked_manuscripts_raw) == "" ~ "missing",
      str_detect(linked_manuscripts_raw, "doi\\.org") ~ "valid_doi",
      str_detect(linked_manuscripts_raw, "ft_") ~ "other_citations",
      TRUE ~ "unspecified"
    )
  ) %>%
  
  # Final selection and reordering of columns
  select(
    # Core IDs
    study_id, full_text_id,
    
    # Core Metadata
    dataset_name = datasetName,
    publication_year = publication_year_clean, publication_year_status,
    
    # Author Info
    author_first_last = author_first_last_raw, author_status,
    
    # Sampling Effort
    sampling_effort_raw, sampling_effort_status, trap_nights, n_sessions,
    sampling_effort_clean, start_year, end_year,
    
    # Data Access & Links
    data_access_level, data_access_raw,
    linked_manuscripts_raw, linked_manuscripts_status,
    
    # QC & Notes
    data_checker, notes
  )


# Descriptives v3 ---------------------------------------------------------

descriptives_v3 <- combined_data_v3$descriptive

# 1. Join with master citation list to get publication year, author, and title
descriptives_v3_cleaned <- descriptives_v3 %>%
  left_join(
    combined_data$citations$all_citations %>% dplyr::select(
      full_text_id, `Publication Year`, Author, Title),
    # It's good practice to specify relationship to silence the warning if duplicates are expected
    by = "full_text_id", relationship = "many-to-many"
  ) %>%
  
  # 2. Standardize columns from the join
  rename(
    publication_year_raw = `Publication Year`,
    author_first_last_raw = Author
  ) %>%
  
  # --- Add v2-specific placeholder columns for a clean merge later ---
  mutate(
    author_status = NA_character_,
    data_checker = NA_character_
  ) %>%
  
  # 3. Rename columns to a standardized format
  rename(
    sampling_effort_raw = sampling_effort,
    data_access_raw = data_access,
    data_resolution_raw = data_resolution,
    linked_manuscripts_raw = linked_manuscripts,
    dataset_name = datasetName
  ) %>%
  
  # --- Clean and flag publication_year ---
  mutate(
    publication_year_clean = as.numeric(publication_year_raw),
    publication_year_status = if_else(
      is.na(publication_year_clean) | publication_year_clean < 1900,
      "invalid_year",
      "valid_year"
    )
  ) %>%
  
  # --- Clean and flag author name ---
  mutate(
    # First, isolate the first author's name
    first_author_only = str_trim(str_extract(author_first_last_raw, "^[^;]+")),
    # Then, clean the format of just that name
    author_first_last = if_else(
      str_detect(first_author_only, "^[^,]+, "),
      str_replace(first_author_only, "^([^,]+),\\s*(.*)$", "\\2 \\1"),
      first_author_only
    ),
    author_status = if_else(is.na(first_author_only), "missing", "present")
  ) %>%
  
  # 4. Correctly place coalesce() inside a mutate()
  mutate(
    dataset_name_clean = coalesce(dataset_name, Title)
  ) %>%
  
  # Clean and parse sampling_effort
  mutate(
    sampling_effort_raw = str_trim(sampling_effort_raw),
    sampling_effort_status = case_when(
      is.na(sampling_effort_raw) | str_detect(tolower(sampling_effort_raw), "not reported|ifa|pcr|^-?$|n/a") ~ "not_reported",
      str_detect(tolower(sampling_effort_raw), "trapping nights|trap nights|trap-nights|trapnight|trap nigths|trap night") ~ "trap_nights",
      str_detect(tolower(sampling_effort_raw), "session|sessions") ~ "sessions",
      str_detect(tolower(sampling_effort_raw), "month|months|week|weeks|year|years|annual|sample|season|\\d{4}") ~ "time_period",
      TRUE ~ "unspecified"
    ),
    trap_nights = if_else(
      sampling_effort_status == "trap_nights",
      map_dbl(
        str_extract_all(sampling_effort_raw, "\\d+(,\\d+)*"),
        ~ suppressWarnings(max(as.numeric(str_remove_all(.x, ","))))
      ),
      NA_real_
    ),
    n_sessions = if_else(
      sampling_effort_status == "sessions",
      map_dbl(
        str_extract_all(sampling_effort_raw, "\\d+(,\\d+)*"),
        ~ suppressWarnings(max(as.numeric(str_remove_all(.x, ","))))
      ),
      NA_real_
    ),
    start_year = if_else(
      sampling_effort_status == "time_period",
      as.numeric(str_extract(sampling_effort_raw, "\\d{4}")),
      NA_real_
    ),
    end_year = if_else(
      sampling_effort_status == "time_period",
      as.numeric(str_extract(sampling_effort_raw, "\\d{4}$")),
      NA_real_
    ),
    sampling_effort_clean = if_else(
      sampling_effort_status == "time_period",
      sampling_effort_raw,
      NA_character_
    )
  ) %>%
  
  # Clean and standardize data_access
  mutate(
    data_access_level = case_when(
      str_detect(tolower(data_access_raw), "individual") ~ "individual",
      str_detect(tolower(data_access_raw), "summari|summariz|sumariz|summatiz|summmariz") ~ "summarised",
      TRUE ~ "unspecified"
    )
  ) %>%
  
  # Clean and standardize data_resolution
  mutate(
    data_resolution_clean = case_when(
      str_detect(tolower(data_resolution_raw), "individual|species") ~ "individual",
      str_detect(tolower(data_resolution_raw), "site-session") ~ "site-session",
      str_detect(tolower(data_resolution_raw), "site-season") ~ "site-season",
      str_detect(tolower(data_resolution_raw), "session") ~ "session",
      str_detect(tolower(data_resolution_raw), "site|stite") ~ "site",
      str_detect(tolower(data_resolution_raw), "study") ~ "study",
      TRUE ~ "unspecified"
    )
  ) %>%
  
  # Clean and flag linked_manuscripts
  mutate(
    linked_manuscripts_status = case_when(
      is.na(linked_manuscripts_raw) | str_trim(linked_manuscripts_raw) == "" ~ NA_character_,
      str_detect(linked_manuscripts_raw, "doi\\.org") ~ "valid_doi",
      str_detect(linked_manuscripts_raw, "ft_") ~ "other_citations",
      TRUE ~ NA_character_
    )
  ) %>%
  
  mutate(
    author_status = forcats::fct_relevel(author_status, "present", "missing"),
    sampling_effort_status = forcats::fct_relevel(sampling_effort_status, "trap_nights", "sessions", "time_period", "unspecified", "not_reported"),
    data_access_level = forcats::fct_relevel(data_access_level, "individual", "summarised", "unspecified"),
    data_resolution_clean = forcats::fct_relevel(data_resolution_clean, "individual", "site-session", "site-season", "site", "study", "session", "unspecified")
  ) %>%
  
  # Final selection and reordering of columns
  select(
    # Core IDs
    study_id, full_text_id,
    
    # Core Metadata
    dataset_name = dataset_name_clean,
    publication_year = publication_year_clean, publication_year_status,
    
    # Author Info
    author_first_last, author_status,
    
    # Sampling Effort
    sampling_effort_raw, sampling_effort_status, trap_nights, n_sessions,
    sampling_effort_clean, start_year, end_year,
    
    # Data Access & Links
    data_access_level, data_access_raw,
    data_resolution_clean, data_resolution_raw,
    linked_manuscripts_raw, linked_manuscripts_status,
    
    # QC & Notes
    data_checker, notes,
    
    # Remove temporary columns
    -Title, -author_first_last_raw, -first_author_only
  )


# Combined Descriptives ---------------------------------------------------

descriptives_v2_cleaned <- descriptives_v2_cleaned %>%
  mutate(
    # The v2 script should also have these mutate calls, but this ensures safety
    author_status = forcats::fct_relevel(author_status, "present", "missing"),
    sampling_effort_status = forcats::fct_relevel(sampling_effort_status, "trap_nights", "sessions", "time_period", "unspecified", "not_reported"),
    data_access_level = forcats::fct_relevel(data_access_level, "individual", "summarised", "unspecified")
  )

# --- Combine the cleaned data frames ---
descriptives_final <- bind_rows(descriptives_v3_cleaned, descriptives_v2_cleaned) %>%
  distinct(full_text_id, .keep_all = TRUE) %>%
  mutate(
    full_text_id = factor(
      full_text_id, 
      levels = levels(combined_data$citations$all_citations$full_text_id)
    )
  ) %>%
  arrange(full_text_id, study_id)

# Finalize the combined list of dataframes
combined_data <- read_rds(here("data", "data_cleaning", "01_01_output.rds"))
combined_data$descriptives_cleaned <- descriptives_final

write_rds(combined_data, here("data", "data_cleaning", "02_01_output.rds"))
