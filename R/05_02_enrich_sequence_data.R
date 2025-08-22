# Project ArHa: Enrich Sequence Data via NCBI API
# 05_02_enrich_sequence_data.R
# Purpose: This script uses accession numbers from the cleaned sequence data
# to query the NCBI Nucleotide database. It fetches metadata for each sequence
# and joins this enriched information back into the main data frame.
# This process is checkpointed to avoid redundant API calls.

# --- Source the API function ---
# This makes our custom function available in this script.
source(here("R", "fetch_accession_data.R"))

# --- Load Data ---
combined_data <- read_rds(here("data", "data_cleaning", "05_01_output.rds"))
sequence_data <- combined_data$sequence

# --- Step 1: Prepare Accession Numbers for Query ---

# As requested, we only query the unique, non-missing accession numbers to be efficient.
unique_accessions <- sequence_data %>%
  filter(!is.na(accession_number) & accession_number != "") %>%
  distinct(accession_number) %>%
  pull(accession_number)

message(paste("Found", length(unique_accessions), "unique accession numbers to query."))

# --- Step 2: Fetch or Load NCBI Metadata (Checkpointing) ---

# Define a path for our checkpoint file
checkpoint_path <- here("data", "data_cleaning", "ncbi_metadata_cache.rds")

if (file.exists(checkpoint_path)) {
  # If the cache file exists, load the data from it instead of calling the API
  message("Loading NCBI metadata from local cache file...")
  ncbi_metadata <- read_rds(checkpoint_path)
} else {
  # If the cache file does not exist, run the API query
  message("Local cache not found. Fetching data from NCBI...")
  
  # Set in Renviron
  api_key <- Sys.getenv("ENTREZ_KEY")
  
  set_entrez_key(api_key)
  
  ncbi_metadata <- fetch_accession_data(unique_accessions)
  
  # Save the results to the cache file for future runs
  write_rds(ncbi_metadata, file = checkpoint_path)
  message("NCBI metadata saved to local cache: ", checkpoint_path)
}

# --- Step 3: Join Enriched Data Back to Sequence Table ---

if (nrow(ncbi_metadata) > 0) {
  sequence_enriched <- sequence_data %>%
    left_join(ncbi_metadata %>%
                distinct(), by = c("accession_number" = "query_accession"))
  
  message("Successfully joined NCBI metadata to the sequence table.")
} else {
  warning("NCBI metadata is empty. No new data was joined.")
  sequence_enriched <- sequence_data # Proceed with original data
}

# --- Step 4: Download Sequences to FASTA File ---

# Source the new function
source(here("R", "download_sequences.R"))

# Get the list of accession numbers that were successfully found
resolved_accessions <- sequence_enriched %>%
  filter(ncbi_resolved_status == "resolved") %>%
  pull(ncbi_accession_version)

# Define the output path
fasta_output_path <- here("data", "data_cleaning", "project_arha_sequences.fasta")

# Call the function to download the sequences
download_sequences_as_fasta(accession_vec = resolved_accessions, 
                            output_file = fasta_output_path)

# --- Step 5: Finalize and Save ---

# Clean the sequence data frame in our main data object
combined_data$sequence <- sequence_enriched %>%
  mutate(ncbi_create_date = parse_date_time(ncbi_create_date, orders = c("%d-%b-%y")),
         ncbi_update_date = parse_date_time(ncbi_update_date, orders = c("%d-%b-%y")),
         ncbi_collection_date_clean = parse_date_time(ncbi_collection_date, orders = c("%d-%b-%y", "%y", "%d-%m-%y", "%y-%m-%d", "%b-%y")),
         ncbi_country_clean = case_when(
           str_detect(ncbi_country_raw, regex("USA|United States", ignore_case = TRUE)) ~ "United States",
           str_detect(ncbi_country_raw, regex("UK|United Kingdom", ignore_case = TRUE)) ~ "United Kingdom",
           str_detect(ncbi_country_raw, regex("Viet Nam", ignore_case = TRUE)) ~ "Vietnam",
           str_detect(ncbi_country_raw, regex("South Korea|Korea", ignore_case = TRUE)) ~ "South Korea",
           str_detect(ncbi_country_raw, regex("Cote d'Ivoire", ignore_case = TRUE)) ~ "Cote d'Ivoire",
           # For others, extract the first word/phrase before a colon or comma
           TRUE ~ str_extract(ncbi_country_raw, "^[^:,]+")
         ),
         ncbi_iso3 = countrycode(sourcevar = ncbi_country_clean, 
                             origin = 'country.name', 
                             destination = 'iso3c'),
         ncbi_locality = str_remove(ncbi_country_raw, regex(paste0("^", ncbi_country_clean), ignore_case = TRUE)) %>%
           str_remove("^[:punct:]*\\s*") %>%
           # If the removal results in an empty string, turn it into NA
           if_else(. == "", NA_character_, .),
         sequence_record_id = fct_inorder(sequence_record_id),
         ncbi_host = case_when(sequenceType == "Host" ~ ncbi_organism,
                               sequenceType == "Pathogen" ~ ncbi_host,
                               TRUE ~ ncbi_host),
         ncbi_pathogen = case_when(sequenceType == "Pathogen" ~ ncbi_organism,
                                   sequenceType == "Host" ~ NA,
                                   TRUE ~ ncbi_organism),
         pathogen_species_clean = case_when(sequenceType == "Pathogen" ~ pathogen_species_clean,
                                   sequenceType == "Host" ~ NA,
                                   TRUE ~ pathogen_species_clean)
  ) %>%
  select(sequence_record_id, associated_rodent_record_id, associated_pathogen_record_id, study_id, accession_number, ncbi_resolved_status, sequence_type = sequenceType,
         host_species_raw, host_species_clean, flag_host_mismatch, ncbi_host,
         pathogen_species_raw, pathogen_species_clean, flag_pathogen_mismatch, ncbi_pathogen,
         ncbi_accession_primary, ncbi_accession_version, ncbi_seq_length, ncbi_strandedness, ncbi_molecule_type,
         ncbi_definition, ncbi_isolate, ncbi_gene_name, ncbi_protein_name,
         ncbi_country_raw, ncbi_country_clean, ncbi_iso3, ncbi_locality,
         ncbi_collection_date = ncbi_collection_date_clean, ncbi_create_date, ncbi_update_date
         )

# Save the updated combined_data object
write_rds(combined_data, here("data", "data_cleaning", "05_02_output.rds"))
