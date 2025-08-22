# Project ArHa: Clean Rodent Taxonomic Data
# 03_03_clean_rodent_taxa.R
# Purpose: This script takes the raw rodent data and joins it with pre-built
# taxonomic lookup tables (from manual cleaning and GBIF).

combined_data_v2 <- read_rds(here::here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
combined_data_v3 <- read_rds(here::here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))

# Load the taxonomic lookup tables
rodent_names <- read_rds(here("data", "matching", "rodent_names_manual.rds"))
species_names <- read_rds(here("data", "matching", "gbif_species_names.rds"))
genus_names <- read_rds(here("data", "matching", "gbif_genus_names.rds"))
higher_taxa_table <- read_rds(here("data", "matching", "gbif_higher_taxa.rds"))

# Load initial QC data for validation checks
rodent_initial_qc <- read_rds(here("data", "clean_data", "rodent_initial_qc.rds"))
v2_raw_length <- rodent_initial_qc$v2_raw_length
v3_raw_length <- rodent_initial_qc$v3_raw_length

# Join Classifications to the Raw Data ------------------------------------
clean_and_join_rodent_data <- function(raw_host_df) {
  
  # Step 1: Join with manual cleaning table and extract cleaned genus
  cleaned_df <- raw_host_df %>%
    mutate(
      # Coalesce scientificName and genus, handling the case where genus might not exist
      scientificName_coalesced = if("genus" %in% names(.)) {
        coalesce(scientificName, genus)
      } else {
        scientificName
      }
    ) %>%
    mutate(extracted_name = str_to_sentence(str_squish(scientificName_coalesced))) %>%
    left_join(rodent_names, by = c("extracted_name" = "rodent_names")) %>%
    mutate(cleaned_rodent_names = coalesce(cleaned_rodent_names, extracted_name)) %>%
    mutate(genus_clean = str_split(cleaned_rodent_names, " ", simplify = TRUE)[, 1])
  
  # Step 2: Perform all joins and consolidate in a single pipeline
  final_df <- cleaned_df %>%
    # Join species data
    left_join(species_names %>%
                select(rodent_names, resolved_species = resolved_name, resolved_genus = genus, resolved_family = family, resolved_order = order, resolved_class = class, gbif_id),
              by = c("extracted_name" = "rodent_names")) %>%
    rename(gbif_id_species = gbif_id,
           resolved_species_species = resolved_species,
           resolved_genus_species = resolved_genus,
           resolved_family_species = resolved_family,
           resolved_order_species = resolved_order,
           resolved_class_species = resolved_class) %>%
    
    # Join genus data
    left_join(genus_names %>%
                select(rodent_genus,  resolved_genus = genus, resolved_family = family, resolved_order = order, resolved_class = class, gbif_id),
              by = c("genus_clean" = "rodent_genus")) %>%
    rename(gbif_id_genus = gbif_id,
           resolved_genus_genus = resolved_genus,
           resolved_family_genus = resolved_family,
           resolved_order_genus = resolved_order,
           resolved_class_genus = resolved_class) %>%
    
    # Join higher taxa data
    left_join(higher_taxa_table %>%
                select(extracted_name, resolved_genus = genus, resolved_family = family, resolved_order = order, resolved_class = class, gbif_id), by = c("cleaned_rodent_names" = "extracted_name")) %>%
    rename(gbif_id_higher = gbif_id,
           resolved_class_higher = resolved_class,
           resolved_order_higher = resolved_order,
           resolved_family_higher = resolved_family,
           resolved_genus_higher = resolved_genus) %>%
    
    # Step 3: Consolidate the GBIF-resolved taxonomic information
    mutate(
      gbif_id = coalesce(gbif_id_species, gbif_id_genus, gbif_id_higher),
      resolved_species = resolved_species_species,
      resolved_genus = coalesce(resolved_genus_species, resolved_genus_genus, resolved_genus_higher),
      resolved_family = coalesce(resolved_family_species, resolved_family_genus, resolved_family_higher),
      resolved_order = coalesce(resolved_order_species, resolved_order_genus, resolved_order_higher),
      resolved_class = coalesce(resolved_class_species, resolved_class_genus, resolved_class_higher),
      taxon_rank = fct(case_when(
        !is.na(resolved_species) ~ "species",
        !is.na(resolved_genus) ~ "genus",
        !is.na(resolved_family) ~ "family",
        !is.na(resolved_order) ~ "order",
        !is.na(resolved_class) ~ "class",
        TRUE ~ "unresolved"),
        levels = c("species", "genus", "family", "order", "class", "unresolved"))
    ) %>%
    
    # Final selection and reordering of columns
    select(rodent_record_id, study_id, eventDate, resolved_species, resolved_genus, resolved_family, resolved_order, resolved_class, gbif_id, taxon_rank, extracted_name_raw = extracted_name, extracted_genus_raw = genus_clean,
           locality, country, verbatimLocality, coordinate_resolution, decimalLatitude, decimalLongitude, individualCount, trapEffort, trapEffortResolution)
  
  return(final_df)
}

# Apply the function to v2 and v3 data
v2_taxa_cleaned <- clean_and_join_rodent_data(combined_data_v2$rodent)
v3_taxa_cleaned <- clean_and_join_rodent_data(combined_data_v3$rodent)

# QC
rodent_qc_data$v2_raw_length == nrow(v2_taxa_cleaned)
rodent_qc_data$v3_raw_length == nrow(v3_taxa_cleaned)

# --- Combine Species and Genus Level Matching ---
rodent_taxa_cleaned <- bind_rows(v2_taxa_cleaned, v3_taxa_cleaned) %>%
  arrange(rodent_record_id)

combined_data <- read_rds(here("data", "data_cleaning", "02_01_output.rds"))
combined_data$host <- rodent_taxa_cleaned

write_rds(combined_data, here("data", "data_cleaning", "03_03_output.rds"))
