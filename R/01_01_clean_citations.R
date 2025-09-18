# Project ArHa: Initial Citation Cleaning
# 01_01_clean_citations.R
# This script loads the raw data from disk and performs the first stage of cleaning:
# consolidating and processing citation and descriptive metadata.

# Load the raw data from 00_load_data.R
# e.g., combined_data_v2, combined_data_v3
combined_data_v2 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
combined_data_v3 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))

# Cleaned data will be consolidated as combined_data
combined_data <- list()

# --- Initial Citation Cleaning and Consolidation ---

# 1. Create a master list of all citations (all_citations).
# This combines citation info from v2 (which has no duplicates) and v3 (which has a separate sheet).
# This table will serve as the master lookup for all other records.
all_citations <- combined_data_v3$inclusion_full_text %>%
  full_join(bind_rows(
    combined_data_v2$descriptive %>%
      dplyr::select(study_id, full_text_id),
    combined_data_v3$descriptive %>%
      dplyr::select(study_id, full_text_id)
  ), by = "full_text_id") %>%
  select(-any_of(c("study_id.x", "study_id.y"))) %>%
  mutate(full_text_id = forcats::fct_inorder(full_text_id)) %>%
  relocate(study_id, .before = 1)


# 2. Consolidate descriptive metadata from v2 and v3.
# This prioritizes the newer v3 data and combines it with non-overlapping v2 data.
descriptives_consolidated <- combined_data_v3$descriptive %>%
  bind_rows(combined_data_v2$descriptive %>%
              filter(!full_text_id %in% combined_data_v3$descriptive)) %>%
  mutate(is_manual = str_detect(full_text_id, "ft_m_"),
         full_text_id_number = as.numeric(str_extract(full_text_id, "\\d+"))) %>%
  arrange(is_manual, full_text_id_number) %>%
  mutate(full_text_id = fct_inorder(full_text_id)) %>%
  select(-is_manual, -full_text_id_number, -study_id)

# 3. Create the 'extractions' table by joining citations with descriptive metadata.
# This table contains all successfully extracted and consolidated records.
extractions <- all_citations %>%
  distinct(full_text_id, study_id, processed, decision, reason) %>%
  drop_na(study_id) %>%
  group_by(full_text_id) %>%
  mutate(n_records = n()) %>%
  right_join(descriptives_consolidated,
             by = c("full_text_id")) %>%
  distinct(full_text_id, study_id, .keep_all = TRUE) %>%
  group_by(full_text_id) %>%
  mutate(n_records = n()) %>%
  relocate(study_id, .before = 1) %>%
  ungroup() %>%
  mutate(is_duplicated = n_records > 1)

# 4. Categorize the extracted data for quality control and reporting.
# These tables help assess the completeness of the review and flag issues.
multiple_extractions <- extractions %>%
  filter(n_records > 1)

single_extractions <- extractions %>%
  filter(n_records == 1)

no_extractions <- all_citations %>%
  filter(is.na(study_id)) %>%
  filter(!str_detect(decision, "Exclude"))

excluded <- all_citations %>%
  filter(is.na(study_id)) %>%
  filter(str_detect(decision, "Exclude")) %>%
  distinct(full_text_id, .keep_all = TRUE) %>%
  select(full_text_id, decision, reason) %>%
  left_join(all_citations, by = c("full_text_id", "decision", "reason"))

# 5. Finalize the consolidated data list.
# This list contains all the processed tables from this script.
combined_data$citations <- list(
  all_citations = all_citations,
  extractions = extractions,
  multiple_extractions = multiple_extractions,
  single_extractions = single_extractions,
  no_extractions = no_extractions,
  excluded = excluded
)

write_rds(combined_data, here("data", "data_cleaning", "01_01_output.rds"))
