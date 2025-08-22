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

# Standardize column names and types before binding.
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
  # v3 names are already quite close to our target standard
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

sequence_combined <- bind_rows(sequence_v2_std, sequence_v3_std)

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

# Check for sequence records that don't link to a valid pathogen record
orphaned_by_pathogen <- sequence_clean %>%
  filter(sequenceType == "Pathogen") %>%
  anti_join(pathogen_data, by = c("associated_pathogen_record_id" = "pathogen_record_id"))

if (nrow(orphaned_by_pathogen) > 0) {
  message(paste(nrow(orphaned_by_pathogen), "sequence records have no corresponding pathogen record."))
  orphaned_by_pathogen = orphaned_by_pathogen
}

# Check for sequence records that don't link to a valid host record
orphaned_by_host <- sequence_clean %>%
  filter(sequenceType == "Host") %>%
  anti_join(host_data, by = c("associated_rodent_record_id" = "rodent_record_id"))

if (nrow(orphaned_by_host) > 0) {
  message(paste(nrow(orphaned_by_host), "sequence records have no corresponding host record."))
  orphaned_by_host = orphaned_by_host
}

# --- Step 4: Flag Discrepancies with Canonical Data ---

# Join with canonical host and pathogen data to compare names
sequence_flagged <- sequence_clean %>%
  filter(!sequence_record_id %in% orphaned_by_pathogen$sequence_record_id |
           !sequence_record_id %in% orphaned_by_host$sequence_record_id) %>%
  # Join canonical host species name
  left_join(host_data %>%
    select(rodent_record_id, host_species_clean = resolved_species),
    by = c("associated_rodent_record_id" = "rodent_record_id")
  ) %>%
  # Join canonical pathogen species name
  left_join(pathogen_data %>%
    select(pathogen_record_id, pathogen_species_clean = pathogen_name_clean),
    by = c("associated_pathogen_record_id" = "pathogen_record_id")
  ) %>%
  # Create flagging columns
  mutate(
    flag_host_mismatch =  case_when(host_species_raw == host_species_clean ~ FALSE,
                                    is.na(host_species_raw) ~ NA,
                                    TRUE ~ TRUE),
    flag_pathogen_mismatch = case_when(pathogen_species_raw == pathogen_species_clean ~ FALSE,
                                       is.na(pathogen_species_raw) ~ NA,
                                       TRUE ~ TRUE)
  )

# --- Step 5: Finalize and Save ---
sequence_final <- sequence_flagged %>%
  select(sequence_record_id,
         study_id,
         associated_rodent_record_id,
         associated_pathogen_record_id,
         accession_number,
         host_species_raw,
         host_species_clean,
         flag_host_mismatch,
         pathogen_species_raw,
         pathogen_species_clean,
         flag_pathogen_mismatch,
         sequenceType,
         method,
         note)

# Add the cleaned sequence data frame to the main data object
combined_data$sequence <- sequence_final

# Save the updated combined_data object
write_rds(combined_data, here("data", "data_cleaning", "05_01_output.rds"))
