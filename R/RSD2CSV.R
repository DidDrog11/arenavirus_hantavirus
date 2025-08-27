#Code for extracting the databases in ArHa as CSV that can be used outside of R
extract_db <- function(pathToInput, databaseToExtract, outputPath) {
  library(readr)
  library(dplyr)
  
  # Load the .rds object
  df <- readRDS(pathToInput)
  
  # Extract the specified element (assuming it’s a data frame inside the list)
  if (!databaseToExtract %in% names(df)) {
    stop(paste("Database", databaseToExtract, "not found in the input file"))
  }
  
  dfProduct <- df[[databaseToExtract]]
  
  # Write to CSV in UTF-8 encoding
  write_csv(dfProduct, file = outputPath)
}

