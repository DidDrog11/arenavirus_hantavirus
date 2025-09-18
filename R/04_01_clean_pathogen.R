# Project ArHa: Clean Pathogen Data
# 04_01_clean_pathogen.R
# Purpose: This script consolidates raw pathogen data and standardizes
# pathogen names, assay types, and data counts. It maintains a relational
# link to the host data.

# Load the raw data from 00_load_data.R
# e.g., combined_data_v2, combined_data_v3
combined_data_v2 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
combined_data_v3 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))

pathogen_data_v2 <- combined_data_v2$pathogen
pathogen_data_v3 <- combined_data_v3$pathogen

# Load the currently cleaned data from 03_07
combined_data <- read_rds(here("data", "data_cleaning", "03_07_output.rds"))
host_data <- combined_data$host

# --- Step 1: Consolidate v2 and v3 pathogen data and clean names ---

pathogen_data_v2_clean <- pathogen_data_v2 %>%
  # Rename columns to a consistent format
  rename(
    pathogen_name_raw = scientificName,
    assay_raw = identificationRemarks,
    number_tested = occurrenceRemarks,
    number_positive = organismQuantity
  ) %>%
  # Select only the columns needed
  select(
    pathogen_record_id, associated_rodent_record_id, study_id,
    pathogen_name_raw, family, 
    assay_raw,
    number_tested, number_positive, number_negative, number_inconclusive,
    note
  )

pathogen_data_v3_clean <- pathogen_data_v3 %>%
  # Rename columns to a consistent format
  rename(
    pathogen_name_raw = scientificName,
    assay_raw = assay,
    number_tested = tested,
    number_positive = positive,
    number_negative = negative
  ) %>%
  # Select only the columns needed
  select(
    pathogen_record_id, associated_rodent_record_id, study_id,
    pathogen_name_raw, family,
    assay_raw,
    number_tested, number_positive, number_negative, number_inconclusive,
    note
  )

# Combine the two cleaned data frames
pathogen_combined <- bind_rows(pathogen_data_v2_clean, pathogen_data_v3_clean)

# --- QC Step: Check for orphaned pathogen records ---

# Find associated_rodent_record_ids that do not exist in the host table
orphaned_pathogens <- pathogen_combined %>%
  anti_join(host_data, by = c("associated_rodent_record_id" = "rodent_record_id"))

if (nrow(orphaned_pathogens) > 0) {
  message(paste("\nWARNING:", nrow(orphaned_pathogens), "pathogen records have no corresponding host record."))
  print(orphaned_pathogens)
} else {
  message("\nAll pathogen records successfully matched to a host record.")
}

# --- Step 2: Create a manual virus matching dictionary ---
virus_mapping <- read_csv(here("data", "matching", "virus_names_matching.csv"))

viruses_missing <- pathogen_combined %>%
  drop_na(pathogen_name_raw) %>%
  distinct(pathogen_name_raw) %>%
  filter(!pathogen_name_raw %in% virus_mapping$virus)

virus_matched_long <- virus_mapping %>%
  pivot_longer(cols = starts_with("virus_"),
               names_to = "virus_number",
               values_to = "matched_virus",
               values_drop_na = TRUE) %>%
  select(-virus_number) %>%
  distinct(virus, clean_name, matched_virus, taxonomic_level) %>%
  group_by(virus) %>%
  summarise(
    clean_name = str_c(unique(clean_name), collapse = ", "),
    matched_virus = str_c(unique(matched_virus), collapse = ", "),
    taxonomic_level = str_c(unique(taxonomic_level), collapse = ", "),
    .groups = "drop"
  ) %>%
  rename("virus_clean" = clean_name) %>%
  mutate(family_clean = case_when(
    str_detect(virus_clean, regex("Catarina|Kodoko|Middle Pease River|Patawa|Pinhal|Pirital", ignore_case = TRUE)) ~ "Arenaviridae",
    str_detect(virus_clean, regex("Mammarenvirus", ignore_case = TRUE)) ~ "Arenaviridae",
    str_detect(virus_clean, regex("Altai|Ape Aime|Araucaria|Azagny|Boginia|Calabazo|Dahonggou|Gou|Thottimvirus|Jabora|Juquitiba|Lohja|Maripa|Mayotte|Mobatvirus|Oran", ignore_case = TRUE)) ~ "Hantaviridae",
    str_detect(virus_clean, regex("Hantavirus", ignore_case = TRUE)) ~ "Hantaviridae",
    str_detect(virus_clean, regex("orthohantavirus", ignore_case = TRUE)) ~ "Hantaviridae",
    str_detect(virus_clean, regex("Cowpox", ignore_case = TRUE)) ~ "Poxviridae",
    str_detect(virus_clean, regex("Malacky virus", ignore_case = TRUE)) ~ "Peribunyaviridae",
    str_detect(virus_clean, regex("Rocahepevirus", ignore_case = TRUE)) ~ "Hepaviridae",
    str_detect(virus_clean, regex("arenavirus", ignore_case = TRUE)) ~ "Arenaviridae",
    TRUE ~ "Check"
  ))

virus_long_unnested <- virus_matched_long %>%
  mutate(matched_virus = str_split(matched_virus, ", ")) %>%
  unnest(matched_virus)

unique_viruses <- virus_long_unnested %>%
  drop_na(matched_virus) %>%
  distinct(matched_virus) %>%
  pull(matched_virus)

# --- Step 3: Match to ICTV via NCBI ---
ictv_ids_safe <- safely(get_ids)
ictv_hierarchy_safe <- safely(classification)

ictv_ids <- ictv_ids_safe(unique_viruses, db = "ncbi")

if (is.null(ictv_ids$error)) {
  # Get the full taxonomic hierarchy for each ID
  ictv_hierarchy <- ictv_hierarchy_safe(ictv_ids$result, db = "ncbi")
  
  # Check for errors in the hierarchy call
  if (is.null(ictv_hierarchy$error)) {
    ictv_hierarchy_df <- do.call(rbind, ictv_hierarchy) %>%
      select(query, rank, name) %>%
      pivot_wider(names_from = rank, values_from = name, id_cols = query, values_fn = list) %>%
      select(query, species, genus, subfamily, family) %>%
      mutate(across(c(family, genus, any_of("species"), any_of("subfamily")), ~ map_chr(., ~ if(is.null(.x)) NA_character_ else .x[[1]]))) %>%
      mutate(query = str_remove(query, "ncbi.")) %>%
      # Select the hierarchy columns we care about
      rename(ncbi_species = species) %>%
      select(query, family, any_of("subfamily"), genus, ncbi_species)
    
    ictv_lookup <- tibble(
      pathogen_name_clean = unique_viruses,
      ncbi_id = as.character(ictv_ids$result$ncbi)
    ) %>%
      left_join(ictv_hierarchy_df, by = c("ncbi_id" = "query")) %>%
      select(pathogen_name_clean, ncbi_id, family, subfamily, genus, ncbi_species)
    
    message("\nVirus taxonomy matching complete. Joining to main data frame.")
  }
}

pathogen_unnested_with_ictv <- virus_long_unnested %>%
  left_join(ictv_lookup, by = c("matched_virus" = "pathogen_name_clean"))

nested_matches <- pathogen_unnested_with_ictv %>%
  # Group by the columns that define a single, broad pathogen name
  group_by(virus) %>%
  # Nest the detailed match information into a single list-column
  nest(nested_matches = c(virus_clean, matched_virus, ncbi_species, ncbi_id, genus, subfamily, family)) %>%
  ungroup()

# --- Step 4: Join to pathogen data ---

pathogen_final <- pathogen_combined %>%
  left_join(nested_matches, by = c("pathogen_name_raw" = "virus")) %>%
  mutate(
    # Use map_chr to iterate through the list-column
    family_clean = map_chr(nested_matches, ~ {
      # Handle NULL cases
      if (is.null(.x) || nrow(.x) == 0) {
        return(NA_character_)
      }
      # Get the unique genus name from the tibble
      family_name <- unique(.x$family)
      # Coalesce to handle potential NA values
      coalesce(family_name, NA_character_)
    }),
    family_clean = case_when(is.na(family_clean) & str_detect(family, "hanta|Hanta") ~ "Hantaviridae",
                             is.na(family_clean) & str_detect(family, "arena|Arena") ~ "Arenaviridae",
                             TRUE ~ family_clean),
    taxonomic_level = map_chr(nested_matches, ~ {
      # Handle NULL cases
      if (is.null(.x) || nrow(.x) == 0) {
        return(NA_character_)
      }
      # If a single species has been matched it is species, otherwise set as genus or family
      if(length(.x$ncbi_species) == 1) "species" else if(length(unique(.x$genus)) == 1) "genus" else "family"
    }),
    taxonomic_level = case_when(is.na(taxonomic_level) & !is.na(family_clean) ~ "family",
                                TRUE ~ taxonomic_level),
    pathogen_name_clean = map_chr(nested_matches, ~ {
      # Handle NULL cases
      if (is.null(.x) || nrow(.x) == 0) {
        return(NA_character_)
      }
      # Get the unique genus name from the tibble
      if(length(.x$matched_virus) == 1) unique(.x$matched_virus) else str_c(.x %>%
                                                                            distinct(matched_virus) %>%
                                                                            pull(matched_virus), collapse = ", ")
    }),
    ncbi_species_name = map_chr(nested_matches, ~ {
      # Handle NULL cases
      if (is.null(.x) || nrow(.x) == 0) {
        return(NA_character_)
      }
      # Get the unique genus name from the tibble
      if(length(.x$ncbi_species) == 1) unique(.x$ncbi_species) else str_c(.x %>%
                                                                            distinct(ncbi_species) %>%
                                                                            pull(ncbi_species), collapse = ", ")
    }),
    ncbi_id = map_chr(nested_matches, ~ {
      # Handle NULL cases
      if (is.null(.x) || nrow(.x) == 0) {
        return(NA_character_)
      }
      # Get the unique genus name from the tibble
      if(length(.x$ncbi_id) == 1) unique(.x$ncbi_id) else str_c(.x %>%
                                                                  distinct(ncbi_id) %>%
                                                                  pull(ncbi_id), collapse = ", ")
    }),
    family_clean = fct(family_clean, levels = c("Arenaviridae", "Hantaviridae", "Peribunyaviridae", "Hepeviridae", "Poxviridae")),
    taxonomic_level = fct(taxonomic_level, levels = c("species", "genus", "family"))
    
  ) %>%
  select(pathogen_record_id, associated_rodent_record_id, study_id, ncbi_species_name, ncbi_id, pathogen_name_clean, pathogen_name_raw, family_clean, family, taxonomic_level, pathogen_taxonomy = nested_matches, assay_raw, number_tested, number_positive, number_inconclusive, note)

# --- Step 5: Clean assay type ---

pathogen_assay <- pathogen_final %>%
  mutate(
    assay_clean = fct(case_when(
      # High-priority: Culture (most definitive evidence)
      str_detect(tolower(assay_raw), "culture|viral culture") ~ "Virus Culture",
      # High-priority: PCR and sequencing (direct evidence)
      str_detect(tolower(assay_raw), "pcr|rt-pcr") ~ "PCR",
      str_detect(tolower(assay_raw), "sequencing|illumina") ~ "Sequencing",
      # Mid-priority: Serology and related methods (indirect evidence)
      str_detect(tolower(assay_raw), "serology|elisa|elsia|antibody|antigen|ifa|eia|immuno|frnt|hdp|igg|complement") ~ "Serology",
      str_detect(tolower(assay_raw), "western blot|immunoblot|ib") ~ "Western Blot",
      # Low-priority: Other and missing
      str_detect(tolower(assay_raw), "rapid field test") ~ "Other",
      TRUE ~ "Missing"
    ), levels = c("Virus Culture", "PCR", "Sequencing", "Serology", "Western Blot", "Other", "Missing")
    )) %>%
  relocate(assay_clean, .before = assay_raw)

# --- Step 6: Match to ICTV via NCBI ---

pathogen_n <- pathogen_assay %>%
  left_join(host_data %>%
              filter(rodent_record_id %in% pathogen_assay$associated_rodent_record_id) %>%
              select(rodent_record_id, individual_count),
            by = c("associated_rodent_record_id" = "rodent_record_id")) %>%
  mutate(
    non_integer = (number_tested %% 1 != 0 & !is.na(number_tested)) |
      (number_positive %% 1 != 0 & !is.na(number_positive)),
    tested_detected = number_tested > individual_count,
    positive_tested = number_positive > number_tested,
    number_negative = case_when(number_tested >= number_positive ~ number_tested - number_positive,
                                TRUE ~ NA_real_)
  ) %>%
  relocate(number_negative, .after = number_positive) %>%
  select(-individual_count)

# --- Step 7: Add to combined data and save ---

combined_data$pathogen <- pathogen_n

write_rds(combined_data, here::here("data", "data_cleaning", "04_01_output.rds"))
