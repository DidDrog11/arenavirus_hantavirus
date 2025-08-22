# R/download_sequences.R

#' Download NCBI Nucleotide Sequences as a FASTA file
#'
#' Takes a vector of NCBI accession numbers and downloads the corresponding
#' sequences in FASTA format, saving them to a single file. Uses a robust,
#' doubly-batched method to handle very large numbers of requests efficiently
#' and avoid API errors.
#'
#' @param accession_vec A character vector of valid NCBI accession.version numbers.
#' @param output_file The file path where the FASTA file should be saved.
#' @return Invisibly returns the file path. The main result is the created file.

download_sequences_as_fasta <- function(accession_vec, output_file) {
  
  if (!requireNamespace("rentrez", quietly = TRUE) ||
      !requireNamespace("purrr", quietly = TRUE)) {
    stop("This function requires the 'rentrez' and 'purrr' packages.")
  }
  
  # Ensure the output file is empty before we start appending
  if (file.exists(output_file)) {
    message("Output file exists. Deleting to start fresh.")
    file.remove(output_file)
  }
  
  # --- Helper to fetch from a single history object and append to file ---
  fetch_and_write_fasta <- function(history_object, n_records) {
    if (n_records == 0) return(NULL)
    
    fetch_batch_size <- 500
    start_indices <- seq(0, n_records - 1, by = fetch_batch_size)
    
    # Use purrr::walk for side-effects (writing to a file)
    purrr::walk(start_indices, function(start_index) {
      message(paste("...fetching batch starting at record", start_index + 1, "of", n_records))
      
      fasta_text <- rentrez::entrez_fetch(db = "nuccore", 
                                          web_history = history_object, 
                                          rettype = "fasta",
                                          retmax = fetch_batch_size,
                                          retstart = start_index)
      
      # Append the downloaded text to the output file
      if (!is.null(fasta_text) && nchar(fasta_text) > 0) {
        write(fasta_text, file = output_file, append = TRUE)
      }
    })
  }
  
  # --- Main Logic ---
  n_ids <- length(accession_vec)
  message(paste("Downloading", n_ids, "sequences..."))
  
  post_batch_size <- 200
  id_chunks <- split(accession_vec, ceiling(seq_along(accession_vec) / post_batch_size))
  message(paste("Splitting into", length(id_chunks), "chunks for POSTing to NCBI History."))
  
  # POST each chunk to the history server
  list_of_histories <- purrr::map(id_chunks, function(chunk) {
    rentrez::entrez_post(db = "nuccore", id = chunk)
  })
  
  # Loop through each history object and process it
  for (i in seq_along(list_of_histories)) {
    history_obj  <- list_of_histories[[i]]
    chunk_of_ids <- id_chunks[[i]]
    
    message(paste0("\n--- Processing POST chunk ", i, " of ", length(list_of_histories), " (", length(chunk_of_ids), " IDs) ---"))
    
    fetch_and_write_fasta(history_object = history_obj, n_records = length(chunk_of_ids))
  }
  
  message("\nFASTA file created successfully at: ", output_file)
  invisible(output_file)
}