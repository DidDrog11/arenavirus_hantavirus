# Project ArHa: Impute Trapping Effort
#
# Purpose: This script assesses the feasibility of imputing `trap_nights_clean`
# and performed statistical imputation.

# Will review the need for this

# Load the cleaned data from the previous script
combined_data <- read_rds(here::here("data", "data_cleaning", "03_08_output.rds"))
cleaned_host <- combined_data$host

# --- Step 1: Create Descriptive Tables to Assess Data Sparsity ---

# Summary of records with reported trap_nights (the potential training data)
training_data_summary <- cleaned_host %>%
  filter(detection_status == "reported") %>%
  filter(!is.na(trap_nights_clean)) %>%
  group_by(NAME_0, resolved_species, resolved_genus, trap_effort_resolution_clean, coordinate_resolution_processed, temporal_resolution) %>%
  summarise(
    n_records = n(),
    min_trap_nights = min(trap_nights_clean, na.rm = TRUE),
    max_trap_nights = max(trap_nights_clean, na.rm = TRUE),
    .groups = "drop"
  )

# Species with no trapping effort data
species_sampling_effort <- cleaned_host %>%
  filter(detection_status == "reported") %>%
  filter(is.na(trap_nights_clean)) %>%
  group_by(resolved_species, resolved_genus) %>%
  summarise(n_records = n(),
            .groups = "drop") %>%
  filter(resolved_species %in% training_data_summary$resolved_species &
           resolved_genus %in% training_data_summary$resolved_genus)

species_without_sampling_effort <- cleaned_host %>%
  filter(detection_status == "reported") %>%
  filter(is.na(trap_nights_clean)) %>%
  group_by(resolved_species, resolved_genus) %>%
  summarise(n_records = n(),
            .groups = "drop") %>%
  filter(!resolved_species %in% training_data_summary$resolved_species |
           !resolved_genus %in% training_data_summary$resolved_genus)

# Summary of records with missing trap_nights (the data to be imputed)
imputation_data_summary <- cleaned_host %>%
  filter(detection_status == "reported") %>%
  filter(is.na(trap_nights_clean)) %>%
  group_by(NAME_0, resolved_species) %>%
  summarise(
    n_records = n(),
    .groups = "drop"
  )
