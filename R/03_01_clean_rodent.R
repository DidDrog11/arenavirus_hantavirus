# Project ArHa: Initial Rodent Data Cleaning
# 03_01_clean_rodent.R
# Purpose: This script performs the initial data preparation and quality control
# for the rodent data. It loads the raw data, checks for unique record IDs,
# and captures the raw number of records for subsequent validation steps.

# Load the raw data from 00_load_data.R
# e.g., combined_data_v2, combined_data_v3
combined_data_v2 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
combined_data_v3 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))

# --- Rodent Data Preparation and QC ---

# Get the raw lengths of the rodent data frames for validation later
v2_raw_length <- nrow(combined_data_v2$rodent)
v3_raw_length <- nrow(combined_data_v3$rodent)

message(paste("Raw v2 rodent records:", v2_raw_length))
message(paste("Raw v3 rodent records:", v3_raw_length))


# Check for unique record IDs across both v2 and v3
rodent_uid <- bind_rows(
  combined_data_v2$rodent %>% mutate(source = "v2"),
  combined_data_v3$rodent %>% mutate(source = "v3")
) %>%
  select(rodent_id = rodent_record_id, source) %>%
  group_by(rodent_id) %>%
  mutate(n = n()) %>%
  ungroup()

message("\nRodent record ID uniqueness check:")
print(table(id_count = rodent_uid$n))

# Check that individualCounts are whole numebrs
rodents_detected <- bind_rows(
  combined_data_v2$rodent %>% mutate(source = "v2"),
  combined_data_v3$rodent %>% mutate(source = "v3")
) %>%
  group_by(rodent_record_id) %>%
  mutate(n_check = as.integer(individualCount) == as.numeric(individualCount))

# --- Finalize and Save ---

# Create a list of key outputs to save for the next scripts
rodent_qc_data <- list(
  rodent_uid_check = rodent_uid,
  v2_raw_length = v2_raw_length,
  v3_raw_length = v3_raw_length
)

rodent_initial_qc <- write_rds(rodent_qc_data, here("data", "clean_data", "rodent_initial_qc.rds"))
