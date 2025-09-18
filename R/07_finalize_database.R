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
pathogen_final <- combined_data$pathogen %>%
  clean_names() %>%
  rename(
    host_record_id = associated_rodent_record_id,
    pathogen_species_ncbi = ncbi_species_name,
    pathogen_species_cleaned = pathogen_name_clean,
    pathogen_species_original = pathogen_name_raw,
    pathogen_family = family_clean,
    pathogen_family_original = family,
    assay = assay_clean,
    number_positive = number_positive
  ) %>%
  select(
    # Record IDs
    pathogen_record_id,
    host_record_id,
    study_id,
    # Taxonomy (Cleaned)
    pathogen_species_cleaned,
    pathogen_family,
    taxonomic_level,
    pathogen_species_ncbi,
    ncbi_id,
    # Assay Info
    assay,
    # Test Counts
    number_tested,
    number_positive,
    number_negative,
    number_inconclusive,
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
