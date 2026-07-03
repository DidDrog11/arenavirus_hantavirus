# Project ArHa: 07_01_finalize_database.R
#
# Purpose: This is the final script in the cleaning pipeline. It takes the fully
# processed data objects and applies a final layer of polish to ensure all
# tables have consistent, logical, and clean column names and a thematic column order.
# The output of this script is the final, distributable database object.
# 

combined_data <- read_rds(here("data", "data_cleaning", "05_02_output.rds"))

# --- 1. Clean the Citations Table ---
message("Cleaning citations table...")
citations_final <- combined_data$citations$all_citations %>%
  # Standardize names to snake_case
  clean_names() %>%
  # Rename for clarity and consistency
  rename(
    zotero_key = key
  ) %>%
  # Reorder columns thematically
  select(
    # Core IDs
    study_id,
    full_text_id,
    zotero_key,
    doi,
    # Citation Info
    author,
    publication_year,
    title,
    publication_title,
    volume,
    issue,
    # Extraction Metadata
    decision,
    reason,
    processed,
    extractor,
    # Other
    rights,
    extra
  )

# --- 2. Clean the Descriptives Table ---
message("Cleaning descriptives table...")
descriptives_final <- combined_data$descriptives_cleaned %>%
  clean_names() %>%
  rename(
    publication_title = dataset_name,
    author_first = author_first_last
  ) %>%
  select(
    # Core IDs
    study_id,
    full_text_id,
    # Study Info
    publication_title,
    publication_year,
    publication_year_status,
    author_first,
    author_status,
    # Sampling Info
    sampling_effort_raw,
    sampling_effort_clean,
    sampling_effort_status,
    trap_nights,
    n_sessions,
    # Data Info
    data_access_level,
    data_resolution_clean,
    data_access_raw,
    data_resolution_raw,
    # Provenance
    notes,
    data_checker
  )

# --- 3. Clean the Host Table ---
message("Cleaning host table...")
host_final <- combined_data$host %>%
  # Standardize names first
  clean_names() %>%
  # Rename for clarity and consistency
  rename(
    host_record_id = rodent_record_id,
    host_species = resolved_species,
    host_genus = resolved_genus,
    host_family = resolved_family,
    host_order = resolved_order,
    host_class = resolved_class,
    host_species_original = extracted_name_raw,
    host_genus_original = extracted_genus_raw,
    number_of_hosts = individual_count,
    country = country_processed,
    iso3c = iso3_code_processed,
    latitude = decimal_latitude,
    longitude = decimal_longitude,
    gadm_country = name_0,
    gadm_adm1 = adm1_name,
    gadm_adm2 = adm2_name,
    gadm_adm3 = adm3_name
  ) %>%
  # Reorder columns thematically
  select(
    # Record IDs
    host_record_id,
    study_id,
    # Taxonomy (Cleaned)
    host_species, host_genus, host_family, host_order, host_class,
    taxon_rank, gbif_id,
    # Geography (Cleaned)
    country, iso3c, latitude, longitude, coord_status,
    dist_from_expected, coordinate_resolution_processed,
    gadm_country, gadm_adm1, gadm_adm2, gadm_adm3,
    # Temporal (Cleaned)
    start_date, end_date, temporal_resolution,
    # Sampling Effort
    number_of_hosts,
    trap_nights_clean, trap_nights_status,
    # Raw/Source Data for Provenance
    host_species_original, host_genus_original,
    country_raw, locality_raw, verbatim_locality_raw,
    event_date_raw
  )

# --- 4. Clean the Pathogen Table ---
message("Cleaning pathogen table...")
pathogen_final <- combined_data$pathogen |>
  janitor::clean_names() |>
  rename(
    host_record_id = associated_rodent_record_id,
    pathogen_species_ncbi = ncbi_species_name,
    pathogen_genus_ncbi = ncbi_genus,
    pathogen_species_cleaned = pathogen_name_clean,
    pathogen_species_original = pathogen_name_raw,
    pathogen_family = family_clean,
    pathogen_family_original = family,
    assay = assay_clean,
  ) |>
  select(
    # Record IDs
    pathogen_record_id,
    host_record_id,
    study_id,
    # Taxonomy
    pathogen_species_cleaned,   # The Python Standardized Name
    pathogen_species_ncbi,      # The Official NCBI Species Node (NA for ambiguous)
    pathogen_genus_ncbi,        # The Official NCBI Genus
    pathogen_family,
    taxonomic_level,
    ncbi_id,
    taxonomy_details,           # List column with the specific details
    
    # Assay Info
    assay,
    
    # Test Counts
    number_tested,
    number_positive,
    number_negative,
    number_inconclusive,
    
    # QC Flags
    non_integer,
    tested_detected,
    positive_tested,
    
    # Raw/Source Data for Provenance
    pathogen_species_original,
    pathogen_family_original,
    assay_raw,
    note
  )

# --- 5. Clean the Sequence Table ---
message("Cleaning sequence table...")
sequence_final <- combined_data$sequence %>%
  clean_names() %>%
  # Rename for clarity, consistency, and to remove redundant 'ncbi_' prefix
  rename(
    host_record_id = associated_rodent_record_id,
    pathogen_record_id = associated_pathogen_record_id,
    resolution_status = ncbi_resolved_status,
    accession_primary = ncbi_accession_primary,
    seq_length = ncbi_seq_length,
    strandedness = ncbi_strandedness,
    molecule_type = ncbi_molecule_type,
    definition = ncbi_definition,
    host = ncbi_host,
    pathogen = ncbi_pathogen,
    create_date = ncbi_create_date,
    update_date = ncbi_update_date,
    isolate = ncbi_isolate,
    country_raw = ncbi_country_raw,
    country = ncbi_country_clean,
    iso3c = ncbi_iso3,
    locality = ncbi_locality,
    collection_date = ncbi_collection_date,
    gene_name = ncbi_gene_name,
    protein_name = ncbi_protein_name,
    host_species_original = host_species_raw,
    pathogen_species_original = pathogen_species_raw,
    query_accession = accession_number
  ) %>%
  select(
    # Record IDs
    sequence_record_id,
    host_record_id,
    pathogen_record_id,
    study_id,
    sequence_type,
    # Query & Resolution Info
    query_accession,
    resolution_status,
    # NCBI Record Info
    accession_primary,
    definition,
    seq_length,
    molecule_type,
    strandedness,
    # Taxonomy
    host,
    pathogen,
    isolate,
    gene_name,
    protein_name,
    # Geography
    country,
    iso3c,
    locality,
    country_raw,
    # Temporal
    collection_date,
    create_date,
    update_date,
    # Original raw data for provenance
    host_species_original, host_species_clean,
    pathogen_species_original, pathogen_species_clean
  )

# --- 6. Assemble Final Database Object and Save ---
message("Assembling final database object...")
# Create the final list object with consistent names
arha_db <- list(
  citations = citations_final,
  descriptives = descriptives_final,
  host = host_final,
  pathogen = pathogen_final,
  sequence = sequence_final
)

# Define a file name with the current date for versioning
final_db_path <- here("data", "database", paste0("Project_ArHa_database_", Sys.Date(), ".rds"))

# Create the directory if it doesn't exist
if (!dir.exists(here("data", "database"))) {
  dir.create(here("data", "database"))
}

# Save the final object
write_rds(arha_db, final_db_path)

# --- 7. Data Release: Darwin Core Mapping and Parquet Export ---
message("Packaging for data release...")

pkgs <- c("arrow", "fs", "jsonlite")
pacman::p_load(pkgs, character.only = T)

# Define target directories
release_path <- here("data", "data_release")
dir_create(path(release_path, "parquet"))
dir_create(path(release_path, "csv"))

event_map <- tibble(old_study_id = as.character(unique(c(descriptives_final$study_id, host_final$study_id, pathogen_final$study_id)))) |>
  drop_na() |>
  mutate(new_event_id = paste0("event_", sprintf("%04d", row_number())))

host_map <- tibble(old_host_id = as.character(unique(host_final$host_record_id))) |>
  drop_na() |>
  mutate(new_host_id = paste0("host_", sprintf("%05d", row_number())))

pathogen_map <- tibble(old_pathogen_id = as.character(unique(pathogen_final$pathogen_record_id))) |>
  drop_na() |>
  mutate(new_pathogen_id = paste0("path_", sprintf("%05d", row_number())))

# Apply via Relational Joins
citations_mapped <- citations_final |>
  mutate(study_id = as.character(study_id)) |>
  left_join(event_map, by = c("study_id" = "old_study_id")) |>
  mutate(study_id = new_event_id) |>
  select(-new_event_id)

descriptives_mapped <- descriptives_final |>
  mutate(study_id = as.character(study_id)) |>
  left_join(event_map, by = c("study_id" = "old_study_id")) |>
  mutate(study_id = new_event_id) |>
  select(-new_event_id)

host_mapped <- host_final |>
  mutate(study_id = as.character(study_id), 
         host_record_id = as.character(host_record_id)) |>
  left_join(event_map, by = c("study_id" = "old_study_id")) |>
  left_join(host_map, by = c("host_record_id" = "old_host_id")) |>
  mutate(study_id = coalesce(new_event_id, study_id), 
         host_record_id = coalesce(new_host_id, host_record_id)) |>
  select(-new_event_id, -new_host_id)

pathogen_mapped <- pathogen_final |>
  mutate(study_id = as.character(study_id),
         host_record_id = as.character(host_record_id),
         pathogen_record_id = as.character(pathogen_record_id)) |>
  left_join(event_map, by = c("study_id" = "old_study_id")) |>
  left_join(host_map, by = c("host_record_id" = "old_host_id")) |>
  left_join(pathogen_map, by = c("pathogen_record_id" = "old_pathogen_id")) |>
  mutate(study_id = coalesce(new_event_id, study_id),
         host_record_id = coalesce(new_host_id, host_record_id),
         pathogen_record_id = coalesce(new_pathogen_id, pathogen_record_id)) |>
  select(-new_event_id, -new_host_id, -new_pathogen_id)

sequence_mapped <- sequence_final |>
  mutate(study_id = as.character(study_id),
         host_record_id = as.character(host_record_id),
         pathogen_record_id = as.character(pathogen_record_id)) |>
  left_join(event_map, by = c("study_id" = "old_study_id")) |>
  left_join(host_map, by = c("host_record_id" = "old_host_id")) |>
  left_join(pathogen_map, by = c("pathogen_record_id" = "old_pathogen_id")) |>
  mutate(study_id = coalesce(new_event_id, study_id),
         host_record_id = coalesce(new_host_id, host_record_id),
         pathogen_record_id = coalesce(new_pathogen_id, pathogen_record_id)) |>
  select(-new_event_id, -new_host_id, -new_pathogen_id)

citations_tmp <- citations_final

sampling_events <- descriptives_mapped |>
  left_join( citations_tmp |> select( full_text_id, author, journal = publication_title, doi ), by = "full_text_id" ) |>
  distinct() |>
  mutate(associatedReferences = case_when(!is.na( doi ) & !is.na( publication_year ) & !is.na( journal ) & !is.na( author ) ~ paste0( author, " (", publication_year, "). ", publication_title, ". ", journal, ". DOI: ", str_trim( doi ) ),
                                          is.na( doi ) & !is.na( publication_year ) & !is.na( journal ) & !is.na( author ) ~ paste0( author, " (", publication_year, "). ", publication_title, ". ", journal, "." ),
                                          is.na( journal ) & !is.na( author ) & !is.na( publication_year ) & !is.na( doi ) ~ paste0( author, " (", publication_year, "). ", publication_title, ". DOI: ", str_trim( doi ) ),
                                          is.na( journal ) & !is.na( author ) & !is.na( publication_year ) & is.na( doi ) ~ paste0( author, " (", publication_year, "). ", publication_title, "."),
                                          !is.na( author ) & is.na( publication_year ) & !is.na( journal ) & !is.na( doi ) ~ paste0( author, publication_title, ". ", journal, ". DOI: ", str_trim( doi ) ),
                                          !is.na( author ) & is.na( publication_year ) & !is.na( doi ) ~ paste0( author, publication_title, ". DOI: ", str_trim( doi ) ),
                                          !is.na( author ) & is.na( publication_year ) & is.na( doi ) ~ paste0( author, publication_title, "."),
                                          !is.na( author ) ~ paste0( author, " (", publication_year, "). ", publication_title ),
                                          TRUE ~ publication_title)) |>
  select(eventID = study_id,
         associatedReferences,
         year = publication_year,
         samplingProtocol = sampling_effort_status,
         sampleSizeValue = trap_nights,
         samplingEffort = n_sessions,
         dataAccessLevel = data_access_level,
         dataResolution = data_resolution_clean)

write_parquet(sampling_events, here("arha_app", "data", "parquet", "sampling_events.parquet"))
write_csv(sampling_events, here("data", "data_release", "csv", "sampling_events.csv"))

host_occurrences <- host_mapped |>
  mutate(
    host_genus = case_when(
      !is.na(host_genus) ~ host_genus,
      stringr::str_starts(host_species_original, "N ")  ~ "Neotoma",
      stringr::str_starts(host_species_original, "C ")  ~ "Chaetodipus",
      stringr::str_starts(host_species_original, "D ")  ~ "Dipodomys",
      stringr::str_starts(host_species_original, "O ")  ~ "Onychomys",
      stringr::str_starts(host_species_original, "P ")  ~ "Peromyscus",
      host_species_original == "\"Social rat\""          ~ "Verbatim Extraction - Social Rat",
      TRUE                                              ~ host_genus
    ),
    taxon_rank = case_when(
      !is.na(taxon_rank) ~ taxon_rank,
      host_genus == "Verbatim Extraction - Social Rat"  ~ "unclassified",
      !is.na(host_genus)                                ~ "species",
      TRUE                                              ~ taxon_rank
    )
  ) |>
  mutate(
    eventDate = case_when(
      !is.na(start_date) & !is.na(end_date) & start_date != end_date ~ paste0(format(start_date, "%Y-%m-%d"), "/", format(end_date, "%Y-%m-%d")),
      !is.na(start_date) ~ format(start_date, "%Y-%m-%d"),
      TRUE ~ NA_character_
    )
  ) |>
  rename(
    occurrenceID = host_record_id,
    eventID = study_id,
    individualCount = number_of_hosts,
    scientificName = host_species,
    genus = host_genus,
    family = host_family,
    order = host_order,
    class = host_class,
    taxonRank = taxon_rank,
    taxonID = gbif_id,
    countryCode = iso3c,
    decimalLatitude = latitude,
    decimalLongitude = longitude,
    stateProvince = gadm_adm1,
    county = gadm_adm2,
    municipality = gadm_adm3,
    verbatimIdentification = host_species_original,
    verbatimLocality = verbatim_locality_raw,
    verbatimEventDate = event_date_raw
  ) |>
  select(
    occurrenceID, eventID, scientificName, genus, family, order, class, taxonRank, taxonID,
    individualCount, eventDate, country, countryCode, stateProvince, county, municipality,
    decimalLatitude, decimalLongitude, verbatimIdentification, verbatimLocality, verbatimEventDate,
    coord_status, dist_from_expected, coordinate_resolution_processed,
    temporal_resolution, trap_nights_clean, trap_nights_status
  )

write_parquet(host_occurrences, here("arha_app", "data", "parquet", "host_occurrences.parquet"))
write_csv(host_occurrences, here("data", "data_release", "csv", "host_occurrences.csv"))

# Pathogen to DwC Extended MeasurementOrFact
pathogen_mof <- pathogen_mapped |>
  # Harmonise missing classifications using historical family entries
  mutate(
    pathogen_family_fixed = case_when(
      !is.na(pathogen_family) ~ pathogen_family,
      pathogen_family_original %in% c("Arenaviridae", "Mammarenaviridae") ~ "Arenaviridae",
      pathogen_family_original == "Hantaviridae"                           ~ "Hantaviridae",
      TRUE                                                                 ~ pathogen_family
    ),
    pathogen_genus_fixed = case_when(
      !is.na(pathogen_genus_ncbi) ~ pathogen_genus_ncbi,
      pathogen_family_original == "Mammarenavirus"                         ~ "Mammarenavirus",
      TRUE                                                                 ~ pathogen_genus_ncbi
    ),
    taxonomic_level_fixed = case_when(
      !is.na(taxonomic_level) ~ taxonomic_level,
      pathogen_family_original %in% c("Arenaviridae", "Hantaviridae", "Mammarenaviridae") ~ "family",
      pathogen_family_original == "Mammarenavirus"                                        ~ "genus",
      pathogen_family_original == "Bacterial"                                             ~ "unclassified",
      TRUE                                                                                ~ taxonomic_level
    )
  ) |>
  # Construct final release fields
  mutate(
    taxonomy_details_json = sapply(taxonomy_details, function(x) {
      if(is.null(x) || length(x) == 0) return(NA_character_)
      jsonlite::toJSON(x, auto_unbox = TRUE) 
    }),
    measurementType = paste(pathogen_species_cleaned, "detection"),
    measurementRemarks = paste0(
      "Number tested: ", number_tested,
      " | Verbatim pathogen: ", pathogen_species_original, 
      " | Verbatim assay: ", assay_raw,
      " | Notes: ", note
    )
  ) |>
  rename(
    measurementID = pathogen_record_id,
    occurrenceID = host_record_id,
    eventID = study_id,
    measurementMethod = assay,
    measurementValue = number_positive,
    scientificName_pathogen = pathogen_species_ncbi,
    genus_pathogen = pathogen_genus_fixed,
    family_pathogen = pathogen_family_fixed,
    taxonRank_pathogen = taxonomic_level_fixed,
    taxonID_pathogen = ncbi_id,
    dynamicProperties = taxonomy_details_json
  ) |>
  select(
    measurementID, occurrenceID, eventID, measurementType, measurementValue,
    measurementMethod, measurementRemarks, scientificName_pathogen, genus_pathogen,
    family_pathogen, taxonRank_pathogen, taxonID_pathogen, dynamicProperties,
    number_tested, number_negative, number_inconclusive, tested_detected,
    positive_tested, pathogen_species_cleaned, pathogen_family_original
  )

write_parquet(pathogen_mof, here("arha_app", "data", "parquet", "pathogen_mof.parquet"))
write_csv(pathogen_mof, here("data", "data_release", "csv", "pathogen_mof.csv"))

# Create a clean, padded ID for the sequences
sequence_clean <- sequence_mapped |>
  mutate(clean_seq_id = paste0("seq_", sprintf("%05d", row_number())),
         # Convert accession to a resolvable URI
         ncbi_uri = paste0("https://www.ncbi.nlm.nih.gov/nuccore/", accession_primary))

# Build Host -> Sequence relationships
seq_rel_host <- sequence_clean |>
  filter(!is.na(host_record_id)) |>
  select(resourceRelationshipID = clean_seq_id,
         resourceID = host_record_id,
         relatedResourceID = ncbi_uri,
         sequence_type,
         gene_name,
         query_accession,
         accession_primary) |>
  mutate(relationshipOfResource = "host has sequence")

# Build Pathogen -> Sequence relationships
seq_rel_pathogen <- sequence_clean |>
  filter(!is.na(pathogen_record_id)) |>
  select(resourceRelationshipID = clean_seq_id,
         resourceID = pathogen_record_id,
         relatedResourceID = ncbi_uri,
         sequence_type,
         gene_name,
         query_accession,
         accession_primary) |>
  mutate(relationshipOfResource = "pathogen has sequence")

# Bind them to single DwC extension table
resource_relationships <- bind_rows(seq_rel_host, seq_rel_pathogen) |>
  select(resourceRelationshipID,
         resourceID,
         relatedResourceID,
         relationshipOfResource,
         # Internal App fields
         sequence_type,
         gene_name,
         query_accession,
         accession_primary)

# Save the Resource Relationship Extension
write_parquet(resource_relationships, here("arha_app", "data", "parquet", "resource_relationships.parquet"))
write_csv(resource_relationships, here("data", "data_release", "csv", "resource_relationships.csv"))
