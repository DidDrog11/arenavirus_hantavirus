# 05_01_clean_sequences.R
# Load the raw data from 00_load_data.R
# e.g., combined_data_v2, combined_data_v3
combined_data_v2 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
combined_data_v3 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))

sequence_data_v2 <- combined_data_v2$sequences
sequence_data_v3 <- combined_data_v3$sequences

# Load the cleaned data from the previous step (04_01_clean_pathogen.R)
combined_data <- read_rds(here("data", "data_cleaning", "04_01_output.rds"))
host_data <- combined_data$host
pathogen_data <- combined_data$pathogen

# --- Step 1: Consolidate v2 and v3 Sequence Data ---

# Standardise column names and types before binding.
sequence_v2_std <- sequence_data_v2 %>%
  drop_na(sequence_record_id) %>%
  mutate(associatedTaxa = coalesce(associatedTaxa, host_genus)) %>%
  select(
    sequence_record_id,
    associated_pathogen_record_id,
    associated_rodent_record_id,
    study_id,
    host_species_raw = associatedTaxa,
    pathogen_species_raw = scientificName,
    sequenceType,
    accession_number,
    method,
    note
  )

sequence_v3_std <- sequence_data_v3 %>%
  drop_na(sequence_record_id) %>%
  select(
    sequence_record_id,
    associated_pathogen_record_id,
    associated_rodent_record_id,
    study_id,
    host_species_raw = associatedTaxa,
    pathogen_species_raw = scientificName,
    sequenceType,
    accession_number,
    method,
    note
  )

sequence_combined <- bind_rows(sequence_v2_std, sequence_v3_std) %>%
  # Records whose host or pathogen key was constructed from a missing ID
  # arrive as the literal string "<extractor>_NA" (e.g. a Host-type row with
  # no real pathogen link still carries associated_pathogen_record_id =
  # "david_NA"). Recode to true
  # NA so the missingness is explicit downstream.
  mutate(
    associated_rodent_record_id = if_else(
      str_detect(associated_rodent_record_id, "_NA$"),
      NA_character_,
      associated_rodent_record_id
    ),
    associated_pathogen_record_id = if_else(
      str_detect(associated_pathogen_record_id, "_NA$"),
      NA_character_,
      associated_pathogen_record_id
    )
  )

# --- Step 2: QC - Ensure Unique Identifiers and Clean Key Columns ---

# Check for duplicate sequence_record_id values
duplicate_ids <- sequence_combined %>%
  group_by(sequence_record_id) %>%
  filter(n() > 1)

if (nrow(duplicate_ids) > 0) {
  message(paste("\nERROR:", nrow(duplicate_ids), "duplicate sequence_record_ids found. Manual review required."))
  print(duplicate_ids)
  # Consider adding stop("Duplicate IDs found.") here for a strict pipeline
} else {
  message("\nAll sequence_record_ids are unique.")
}

# Clean the accession_number column
sequence_clean <- sequence_combined %>%
  mutate(accession_number = str_trim(accession_number))

# --- Step 3: QC - Check for Orphaned Sequence Records ---

orphaned_by_pathogen <- sequence_clean %>%
  filter(sequenceType == "Pathogen") %>%
  anti_join(pathogen_data, by = c("associated_pathogen_record_id" = "pathogen_record_id"))

if (nrow(orphaned_by_pathogen) > 0) {
  message(paste(nrow(orphaned_by_pathogen), "sequence records have no corresponding pathogen record."))
}

orphaned_by_host <- sequence_clean %>%
  filter(sequenceType == "Host") %>%
  anti_join(host_data, by = c("associated_rodent_record_id" = "rodent_record_id"))

if (nrow(orphaned_by_host) > 0) {
  message(paste(nrow(orphaned_by_host), "sequence records have no corresponding host record."))
}

# --- Step 4: Join Canonical Names, Drop Orphans ---

sequence_flagged <- sequence_clean %>%
  filter(!sequence_record_id %in% orphaned_by_pathogen$sequence_record_id &
           !sequence_record_id %in% orphaned_by_host$sequence_record_id) %>%
  left_join(host_data %>%
              select(rodent_record_id, host_species_clean = resolved_species),
            by = c("associated_rodent_record_id" = "rodent_record_id")
  ) %>%
  left_join(pathogen_data %>%
              select(pathogen_record_id, pathogen_name_clean),
            by = c("associated_pathogen_record_id" = "pathogen_record_id")
  )

# --- Step 5: Finalize and Save ---

sequence_final <- sequence_flagged %>%
  select(sequence_record_id,
         study_id,
         associated_rodent_record_id,
         associated_pathogen_record_id,
         accession_number,
         host_species_clean,
         pathogen_species_clean = pathogen_name_clean,
         sequenceType,
         method,
         note)

combined_data$sequence <- sequence_final

write_rds(combined_data, here("data", "data_cleaning", "05_01_output.rds"))
