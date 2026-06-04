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
virus_mapping <- read_delim(here("data", "matching", "virus_names_matching.csv"), delim = ";")

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


# Match using Python script -----------------------------------------------
# 1. Setup Paths
dir.create(here("data", "temp"), showWarnings = FALSE)
input_tsv  <- here("data", "temp", "taxonomy_input.tsv")
output_tsv <- here("data", "temp", "taxonomy_output.tsv")
script_path <- here("Python", "Taxonomy", "TaxonomyReconciler.py")
lookup_arena <- here("data", "matching", "arenaviridae_ncbi_taxonomy_lookup.tsv")
lookup_hanta <- here("data", "matching", "hantaviridae_ncbi_taxonomy_lookup.tsv")

input_base <- virus_long_unnested |>
  distinct(matched_virus, family_clean) |>
  mutate(Organism_Name = str_squish(matched_virus),
         Family = family_clean,
         input_variant = "Original")

input_cleaned <- input_base |>
  mutate(
    Organism_Name_Clean = Organism_Name,
    # Map Common Names -> Scientific Binomials (ICTV format)
    Organism_Name_Clean = case_when(
      # Hantaviruses
      str_detect(Organism_Name_Clean, "Sin Nombre") ~ "Orthohantavirus sinnombreense",
      str_detect(Organism_Name_Clean, "Andes") ~ "Orthohantavirus andesense",
      str_detect(Organism_Name_Clean, "Hantaan") ~ "Orthohantavirus hantanense",
      str_detect(Organism_Name_Clean, "Dobrava") ~ "Orthohantavirus dobravaense",
      str_detect(Organism_Name_Clean, "Puumala") ~ "Orthohantavirus puumalaense",
      str_detect(Organism_Name_Clean, "Seoul") ~ "Orthohantavirus seoulense",
      str_detect(Organism_Name_Clean, "Tula") ~ "Orthohantavirus tulaense",
      str_detect(Organism_Name_Clean, "Black Creek") ~ "Orthohantavirus blackcreekense",
      str_detect(Organism_Name_Clean, "Bayou") ~ "Orthohantavirus bayoui",
      str_detect(Organism_Name_Clean, "Choclo") ~ "Orthohantavirus chocloense",
      # Arenaviruses
      str_detect(Organism_Name_Clean, "Lassa") ~ "Mammarenavirus lassaense",
      str_detect(Organism_Name_Clean, "Guanarito") ~ "Mammarenavirus guanaritoense",
      str_detect(Organism_Name_Clean, "Junin") ~ "Mammarenavirus juninense",
      str_detect(Organism_Name_Clean, "Machupo") ~ "Mammarenavirus machupoense",
      str_detect(Organism_Name_Clean, "Chapare") ~ "Mammarenavirus chapareense",
      str_detect(Organism_Name_Clean, "Gairo") ~ "Mammarenavirus gairoense",
      str_detect(Organism_Name_Clean, "Whitewater") ~ "Mammarenavirus whitewaterense",
      str_detect(Organism_Name_Clean, "Lymphocytic|choriomeningitidis") ~ "Mammarenavirus choriomeningitidis",
      # Fallback: Try generic Regex for others
      TRUE ~ Organism_Name_Clean
    ),
    # Generic Cleanup if no specific rule matched
    Organism_Name_Clean = if_else(Organism_Name_Clean == Organism_Name, 
                                  Organism_Name_Clean |> 
                                    str_replace(" virus", "") |> # Remove 'virus' to try and match the genus/species root
                                    str_to_title(),
                                  Organism_Name_Clean),
    input_variant = "Scientific_Fallback") |>
  select(matched_virus, Organism_Name = Organism_Name_Clean, Family, input_variant) |>
  distinct()

python_input_combined <- bind_rows(input_base, input_cleaned) |>
  select(matched_virus, Organism_Name, Family, input_variant)

write_tsv(python_input_combined |> select(Organism_Name, Family), input_tsv)
message(paste("Sending", nrow(python_input_combined), "variants to Python."))

cmd_args <- c(paste("python", shQuote(script_path)),
              paste("--input", shQuote(input_tsv)),
              paste("--output", shQuote(output_tsv)),
              paste0("--lookup Arenaviridae=", shQuote(lookup_arena)))

if (file.exists(lookup_hanta)) {
  cmd_args <- c(cmd_args, paste0("--lookup Hantaviridae=", shQuote(lookup_hanta)))
}

exit_code <- system(paste(cmd_args, collapse = " "))

python_results <- read_tsv(output_tsv, show_col_types = FALSE)

# --- Step 3: Match to ICTV via NCBI ---
ictv_ids_safe <- safely(get_ids)
ictv_hierarchy_safe <- safely(classification)

names_to_query <- python_results |>
  filter(!is.na(assigned_name)) |>
  filter(!str_detect(assigned_name, "Unclassified")) |>
  filter(!str_detect(assigned_name, "family_not_supported")) |>
  distinct(assigned_name) |>
  pull(assigned_name)

ictv_ids <- ictv_ids_safe(names_to_query, db = "ncbi")

found_ids <- ictv_ids$result$ncbi
valid_indices <- !is.na(found_ids)
ids_to_classify <- found_ids[valid_indices]

ictv_hierarchy <- ictv_hierarchy_safe(ids_to_classify, db = "ncbi")

hierarchy_df <- do.call(rbind, ictv_hierarchy) |>
  select(query, rank, name) |>
  pivot_wider(names_from = rank, values_from = name, id_cols = query, values_fn = list) |>
  mutate(across(everything(), ~ map_chr(., ~ if(is.null(.x)) NA_character_ else .x[[1]]))) |>
  rename(ncbi_id = query, ncbi_species = species, ncbi_genus = genus, ncbi_family = family) |>
  mutate(ncbi_id = str_remove(ncbi_id, "result.")) |>
  select(ncbi_id, ncbi_family, any_of("subfamily"), ncbi_genus, ncbi_species)

id_map <- tibble(assigned_name = names_to_query[valid_indices],
                 ncbi_id = as.character(ids_to_classify)) |>
  left_join(hierarchy_df, by = "ncbi_id")

full_taxonomy_map <- python_results |>
  select(Organism_Name, assigned_name, assigned_family, assignment_source) |>
  left_join(id_map, by = "assigned_name") |>
  mutate(ncbi_family = coalesce(ncbi_family, assigned_family))

virus_mapping_final <- virus_long_unnested |>
  mutate(Organism_Name = str_squish(matched_virus)) |> 
  left_join(full_taxonomy_map, by = "Organism_Name")

virus_mapping_flat <- virus_mapping_final |>
  group_by(virus) |>
  summarise(
    # Combined Names (for reference)
    matched_virus_list = paste(unique(matched_virus), collapse = " | "),
    pathogen_name_clean_list = paste(unique(assigned_name), collapse = " | "),
    
    # Consensus Taxonomy
    # If all matches belong to the same Family, keep it. Else NA.
    family_clean = if(n_distinct(ncbi_family, na.rm = TRUE) == 1) first(na.omit(ncbi_family)) else NA_character_,
    
    # If all matches belong to the same Genus, keep it. Else NA.
    ncbi_genus = if(n_distinct(ncbi_genus, na.rm = TRUE) == 1) first(na.omit(ncbi_genus)) else NA_character_,
    
    # Species is ONLY kept if all inputs map to the exact same species ID
    ncbi_species_name = if(n_distinct(ncbi_species, na.rm = TRUE) == 1) first(na.omit(ncbi_species)) else NA_character_,
    
    # IDs
    # If ambiguous, we concatenate them so we don't lose the trace
    ncbi_id = paste(unique(na.omit(ncbi_id)), collapse = ";"),
    
    # Recalculate Taxonomic Level based on the Consensus above
    taxonomic_level = case_when(!is.na(ncbi_species_name) ~ "species",
                                !is.na(ncbi_genus) ~ "genus",
                                !is.na(family_clean) ~ "family",
                                TRUE ~ "unclassified"),
    
    # E. Nest the details
    taxonomy_details = list(pick(everything())),
    .groups = "drop")

pathogen_final <- pathogen_combined |>
  left_join(virus_mapping_flat, by = c("pathogen_name_raw" = "virus")) |>
  mutate(family_clean = case_when(is.na(family_clean) & str_detect(pathogen_name_raw, regex("hanta", ignore_case=T)) ~ "Hantaviridae",
                                  is.na(family_clean) & str_detect(pathogen_name_raw, regex("arena", ignore_case=T)) ~ "Arenaviridae",
                                  TRUE ~ family_clean),
         taxonomic_level = case_when(!is.na(ncbi_species_name) ~ "species",
                                     !is.na(ncbi_genus) ~ "genus",
                                     !is.na(family_clean) ~ "family",
                                     TRUE ~ "unclassified"),
    # Factor levels
    family_clean = fct(family_clean, levels = c("Arenaviridae", "Hantaviridae", "Peribunyaviridae", "Hepeviridae", "Poxviridae")),
    taxonomic_level = fct(taxonomic_level, levels = c("species", "genus", "family", "unclassified"))) |>
  select(pathogen_record_id, associated_rodent_record_id, study_id, 
         pathogen_name_raw, 
         pathogen_name_clean = pathogen_name_clean_list,  
         ncbi_species_name, ncbi_genus,  ncbi_id, family_clean, family, taxonomic_level, taxonomy_details, # The list column
         assay_raw, number_tested, number_positive, number_inconclusive, note)

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

combined_data$pathogen <- pathogen_n

write_rds(combined_data, here::here("data", "data_cleaning", "04_01_output.rds"))
