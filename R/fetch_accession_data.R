# R/fetch_accession_data.R

#' Fetch and Parse Metadata from NCBI Nucleotide Database (Ultra-Robust Version)
#'
#' This function takes a vector of NCBI accession numbers and queries the Entrez
#' API. It is designed for very large queries and ensures that all input
#' accession numbers are represented in the output, flagging those not found
#' and retaining the original query accession for provenance.
#'
#' @param accession_vec A character vector of valid NCBI accession numbers.
#' @return A tibble with one row for each unique input accession. Columns contain
#'   the original query, parsed metadata or NA for unfound records, plus a status flag.

fetch_accession_data <- function(accession_vec) {
  
  if (!requireNamespace("rentrez", quietly = TRUE) ||
      !requireNamespace("xml2", quietly = TRUE) ||
      !requireNamespace("purrr", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("stringr", quietly = TRUE)) {
    stop("This function requires 'rentrez', 'xml2', 'purrr', 'dplyr', and 'stringr'.")
  }
  
  scaffold <- tibble::tibble(input_accession = unique(accession_vec))
  
  parse_records <- function(entrez_records_xml) {
    if (is.null(entrez_records_xml)) return(tibble())
    all_records <- xml2::xml_find_all(entrez_records_xml, ".//GBSeq")
    list_of_single_row_tibbles <- purrr::map(all_records, function(rec) {
      extract_text <- function(xpath) {
        node <- xml2::xml_find_first(rec, xpath)
        if (is.na(node)) NA_character_ else xml2::xml_text(node)
      }
      extract_qualifiers_collapsed <- function(feature_key, qual_name) {
        xpath <- paste0(".//GBFeature[GBFeature_key='", feature_key, "']/GBFeature_quals/GBQualifier[GBQualifier_name='", qual_name, "']/GBQualifier_value")
        nodes <- xml2::xml_find_all(rec, xpath)
        if (length(nodes) == 0) return(NA_character_) else return(paste(xml2::xml_text(nodes), collapse = "; "))
      }
      tibble::tibble(
        ncbi_accession_version = extract_text(".//GBSeq_accession-version"),
        ncbi_accession_primary = extract_text(".//GBSeq_primary-accession"),
        ncbi_seq_length        = as.integer(extract_text(".//GBSeq_length")),
        ncbi_strandedness      = extract_text(".//GBSeq_strandedness"),
        ncbi_molecule_type     = extract_text(".//GBSeq_moltype"),
        ncbi_definition        = extract_text(".//GBSeq_definition"),
        ncbi_organism          = extract_text(".//GBSeq_organism"),
        ncbi_create_date       = extract_text(".//GBSeq_create-date"),
        ncbi_update_date       = extract_text(".//GBSeq_update-date"),
        ncbi_isolate           = extract_qualifiers_collapsed("source", "isolate"),
        ncbi_host              = extract_qualifiers_collapsed("source", "host"),
        ncbi_country_raw       = {
          country_std <- extract_qualifiers_collapsed("source", "country")
          country_geo <- extract_qualifiers_collapsed("source", "geo_loc_name")
          # Combine non-NA results, separated by a semicolon
          combined <- paste(stats::na.omit(c(country_std, country_geo)), collapse = "; ")
          # Return NA if the result is an empty string, otherwise return the combined string
          if (nzchar(combined)) combined else NA_character_
        },
        ncbi_collection_date   = extract_qualifiers_collapsed("source", "collection_date"),
        ncbi_gene_name         = extract_qualifiers_collapsed("gene", "gene"),
        ncbi_protein_name      = extract_qualifiers_collapsed("CDS", "product")
      )
    })
    dplyr::bind_rows(list_of_single_row_tibbles)
  }
  
  fetch_from_history <- function(history_object, n_records) {
    # ... (This helper function remains unchanged) ...
    if (n_records == 0) return(tibble())
    fetch_batch_size <- 500
    start_indices <- seq(0, n_records - 1, by = fetch_batch_size)
    purrr::map(start_indices, function(start_index) {
      message(paste("...fetching batch starting at record", start_index + 1, "of", n_records))
      batch_xml_text <- rentrez::entrez_fetch(db = "nuccore", web_history = history_object, rettype = "xml", retmax = fetch_batch_size, retstart = start_index, parsed = FALSE)
      if(is.null(batch_xml_text) || nchar(batch_xml_text) == 0) return(tibble())
      batch_xml <- xml2::read_xml(batch_xml_text)
      parse_records(batch_xml)
    }) %>% dplyr::bind_rows()
  }
  
  # --- Main Logic ---
  # ... (The main logic remains unchanged) ...
  n_ids <- nrow(scaffold)
  message(paste("Querying NCBI for", n_ids, "unique accession numbers..."))
  post_batch_size <- 200
  id_chunks <- split(scaffold$input_accession, ceiling(seq_along(scaffold$input_accession) / post_batch_size))
  message(paste("Splitting into", length(id_chunks), "chunks for POSTing."))
  list_of_histories <- purrr::map(id_chunks, function(chunk) { rentrez::entrez_post(db = "nuccore", id = chunk) })
  list_of_tibbles <- list()
  for (i in seq_along(list_of_histories)) {
    history_obj  <- list_of_histories[[i]]
    chunk_of_ids <- id_chunks[[i]]
    message(paste0("\n--- Processing POST chunk ", i, " of ", length(list_of_histories), " (", length(chunk_of_ids), " IDs) ---"))
    list_of_tibbles[[i]] <- fetch_from_history(history_object = history_obj, n_records = length(chunk_of_ids))
  }
  resolved_metadata <- dplyr::bind_rows(list_of_tibbles)
  
  # --- Final Join Logic ---
  # ... (This logic remains unchanged) ...
  resolved_metadata_with_key <- resolved_metadata %>%
    dplyr::mutate(join_key = stringr::str_remove(ncbi_accession_version, "\\.\\d+$"))
  duplicates <- resolved_metadata_with_key %>%
    dplyr::group_by(join_key) %>% dplyr::filter(dplyr::n() > 1)
  if(nrow(duplicates) > 0){
    warning("Multiple versions of the same primary accession were found in the NCBI results for: ", 
            paste(unique(duplicates$join_key), collapse=", "))
  }
  final_output <- scaffold %>%
    dplyr::left_join(resolved_metadata_with_key, by = c("input_accession" = "join_key")) %>%
    dplyr::rename(query_accession = input_accession) %>%
    dplyr::mutate(
      ncbi_resolved_status = dplyr::if_else(is.na(ncbi_seq_length), "unresolved", "resolved")
    ) %>%
    dplyr::select(
      query_accession,
      ncbi_resolved_status,
      ncbi_accession_version,
      ncbi_accession_primary,
      dplyr::everything()
    )
  
  n_resolved <- sum(final_output$ncbi_resolved_status == "resolved")
  message(paste("\nSuccessfully resolved and parsed metadata for", n_resolved, "of", n_ids, "requested records."))
  
  return(final_output)
}
