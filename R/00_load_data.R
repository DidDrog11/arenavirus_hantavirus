# 00_load_data.R
# This script loads the packages and raw data
# Set should_update_raw_data to TRUE to force a download of the current data

if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "countrycode",
  "geodata",
  "googledrive",
  "googlesheets4",
  "here",
  "janitor",
  "readr",
  "readxl",
  "rentrez",
  "stringdist",
  "terra",
  "taxize",
  "tidyverse",
  "tidyterra"
)

pacman::p_load(pkgs, character.only = T)

# Where the extraction data come from. "google" downloads the Google Sheets and saves a
# snapshot of every sheet as CSV in data/raw_data/, with the column types in
# <analysis_date>_column_types.csv. "csv" compiles from that snapshot instead, so the
# build can be repeated without access to the Sheets (the snapshot is archived with
# each release).
raw_data_source <- "csv"

# Searches of PubMed and Web of Science were run on 2023-10-16 (published protocol):
# 4,020 records, plus 11 from other sources, gave 2,755 after de-duplication
search_date <- "2023-10-16"
analysis_date <- "2026-10-20"

# --- Configuration Section ---

# Sheets to download for v2 data
v2_sheets <- c("inclusion_full_text", "descriptive", "rodent", 
               "pathogen", "sequences", "known_zoonoses")

v2_url <- "https://docs.google.com/spreadsheets/d/14eW_YwSP6EWWuDrnsvDX-vwTi-7KnVyRdYL8FDqYivk/edit?usp=sharing"

v2_col_types <- c("cccccciccccccccccccccc", "cccciccccccc", "ccccccccccdddcc", 
                  "ccccccccddccccddddc", "ccccccccccccc", "ccclccc")

names(v2_col_types) <- v2_sheets

# Sheets to download for v3 data
v3_sheets <- c("descriptive", "rodent", "pathogen", "sequences")

v3_urls <- list(david = "https://docs.google.com/spreadsheets/d/14Ghz07XOlZoaie8990mvsPbHfsHplxtfwscgf2-5KDk/edit?usp=sharing",
                ana = "https://docs.google.com/spreadsheets/d/10LsT76WZF4c1-LyvTKwGUsG5KH5VXIhAFcaALK6x1Cw/edit?usp=sharing",
                grant = "https://docs.google.com/spreadsheets/d/16pA8h33bwSMHsif6DaxfyeGKUUf0Ol9006oniEb2Gqg/edit?usp=sharing",
                harry = "https://docs.google.com/spreadsheets/d/1lKZitcCCP1iKGQdT8shZBewnkJoaal5flJOevOIFHzQ/edit?usp=sharing",
                ricardo = "https://docs.google.com/spreadsheets/d/1PVa2uF2-hG13XxGWOA_RAW0fZQ_rJKviX75McoNgsjQ/edit?usp=sharing")

v3_full_text_studies <- "https://docs.google.com/spreadsheets/d/1rh6H_N_f3NgLZoPRIjMKgkrYCTgaznReub_EhBctHxI/edit?usp=sharing"

v3_col_types <- c(descriptive = "c", rodent = "ccccccccdddcc", 
                  pathogen = "cccccccddddc", sequences = "cccccccccccc")

# --- Helper Functions ---

# Function to read and process a v2 sheet
process_v2_sheet <- function(sheet_name, col_type_string) {
  df <- read_sheet(v2_url, sheet = sheet_name, col_types = col_type_string)
  
  if (sheet_name == "inclusion_full_text") {
    df <- df %>% filter(processed == "y")
  }
  
  if (sheet_name == "descriptive") {
    df <- df %>%
      drop_na(study_id) %>%
      mutate(study_id = fct_inorder(paste0("david_", str_extract(study_id, "\\d+"))))
  }
  
  if (sheet_name == "rodent") {
    df <- df %>%
      mutate(
        study_id = fct_inorder(paste0("david_", str_extract(study_id, "\\d+"))),
        rodent_record_id = fct_inorder(paste0("david_", str_extract(rodent_record_id, "\\d+")))
      )
  }
  
  if (sheet_name == "pathogen") {
    df <- df %>%
      mutate(
        study_id = fct_inorder(paste0("david_", str_extract(study_id, "\\d+"))),
        associated_rodent_record_id = fct_inorder(paste0("david_", str_extract(associated_rodent_record_id, "\\d+"))),
        pathogen_record_id = fct_inorder(paste0("david_", str_extract(pathogen_record_id, "\\d+")))
      )
  }
  
  if (sheet_name == "sequences") {
    df <- df %>%
      mutate(
        study_id = fct_inorder(paste0("david_", str_extract(study_id, "\\d+"))),
        associated_rodent_record_id = fct_inorder(if_else(!is.na(associated_rodent_record_id), paste0("david_", str_extract(associated_rodent_record_id, "\\d+")), NA_character_)),
        associated_pathogen_record_id = fct_inorder(if_else(!is.na(associated_pathogen_record_id), paste0("david_", str_extract(associated_pathogen_record_id, "\\d+")), NA_character_))
      )
  }
  
  if (sheet_name == "known_zoonoses") {
    # No mutations needed for this sheet
  }

  return(df)
}

# Function to read and process v3 sheets for all extractors
process_v3_sheet <- function(sheet_name) {
  
  # Read and process for each extractor, then bind rows
  lapply(names(v3_urls), function(name) {
    url <- v3_urls[[name]]
    
    # Read sheet with specific range if needed
    if (sheet_name %in% c("rodent")) {
      df <- read_sheet(url, sheet = sheet_name, range = paste0(sheet_name, "!A:M"), col_types = v3_col_types[[sheet_name]])
    } else {
      if (sheet_name %in% c("pathogen", "sequences")) {
        df <- read_sheet(url, sheet = sheet_name, range = paste0(sheet_name, "!A:L"), col_types = v3_col_types[[sheet_name]])
      } else {
      df <- read_sheet(url, sheet = sheet_name, col_types = v3_col_types[[sheet_name]])
      }
    }
    
    # --- Sheet-Specific Transformations using if/else if ---
    if (sheet_name == "descriptive") {
      df <- df %>%
        drop_na(study_id) %>%
        mutate(study_id = fct_inorder(paste0(name, "_", str_extract(study_id, "\\d+")))) %>%
        select(any_of(c("study_id", "full_text_id", "datasetName", "sampling_effort", "data_access", "data_resolution", "linked_manuscripts", "notes")))
      
    } else if (sheet_name == "rodent") {
      df <- df %>%
        drop_na(study_id) %>%
        mutate(
          study_id = fct_inorder(paste0(name, "_", str_extract(study_id, "\\d+"))),
          rodent_record_id = fct_inorder(paste0(name, "_", str_extract(rodent_record_id, "\\d+")))
        )
      
    } else if (sheet_name == "pathogen") {
      df <- df %>%
        drop_na(study_id, pathogen_record_id) %>%
        mutate(
          study_id = fct_inorder(paste0(name, "_", str_extract(study_id, "\\d+"))),
          associated_rodent_record_id = fct_inorder(paste0(name, "_", str_extract(associated_rodent_record_id, "\\d+"))),
          pathogen_record_id = fct_inorder(paste0(name, "_", str_extract(pathogen_record_id, "\\d+")))
        )
      
    } else if (sheet_name == "sequences") {
      df <- df %>%
        drop_na(study_id) %>%
        mutate(
          study_id = fct_inorder(paste0(name, "_", str_extract(study_id, "\\d+"))),
          associated_rodent_record_id = fct_inorder(if_else(!is.na(associated_rodent_record_id), paste0(name, "_", str_extract(associated_rodent_record_id, "\\d+")), NA_character_)),
          associated_pathogen_record_id = fct_inorder(if_else(!is.na(associated_pathogen_record_id), paste0(name, "_", str_extract(associated_pathogen_record_id, "\\d+")), NA_character_)),
          sequence_record_id = fct_inorder(if_else(!is.na(sequence_record_id), paste0(name, "_", str_extract(sequence_record_id, "\\d+")), NA_character_))
        )
    }
    
    # Return the processed data frame
    return(df)
  }) %>% bind_rows()
}


# CSV snapshot of the downloaded sheets: one file per sheet, and a table of column types
# so that the "csv" route rebuilds exactly the objects the download produced.
# - Missing values are written as empty cells, because some cells hold the text "NA"
#   (the Sheets give no empty strings, so nothing else is written as an empty cell).
# - Decimal numbers are read as text and converted with as.numeric(), which reads
#   them exactly; readr's parser can be off in the last digit.
# - ID columns are factors in order of first appearance, which fct_inorder() restores.
snapshot_file <- function(version, sheet) {
  here("data", "raw_data", paste0(analysis_date, "_", version, "_", sheet, ".csv"))
}
types_file <- here("data", "raw_data", paste0(analysis_date, "_column_types.csv"))
type_codes <- c(character = "c", factor = "c", integer = "i", numeric = "c", logical = "l")

write_snapshot <- function(sheets, version) {
  for (s in names(sheets)) write_csv(sheets[[s]], snapshot_file(version, s), na = "")
  tibble(version = version,
         sheet   = rep(names(sheets), map_int(sheets, ncol)),
         column  = unlist(map(sheets, names), use.names = FALSE),
         type    = unlist(map(sheets, \(df) map_chr(df, \(col) class(col)[1])), use.names = FALSE))
}

read_snapshot <- function(version, column_types) {
  types <- column_types |> filter(version == !!version)
  map(set_names(unique(types$sheet)), function(s) {
    sheet_types <- types |> filter(sheet == s)
    stopifnot(all(sheet_types$type %in% names(type_codes)))
    read_csv(snapshot_file(version, s), col_types = paste(type_codes[sheet_types$type], collapse = ""), na = "") |>
      mutate(across(all_of(sheet_types$column[sheet_types$type == "numeric"]), as.numeric),
             across(all_of(sheet_types$column[sheet_types$type == "factor"]), fct_inorder))
  })
}

# Download data -----------------------------------------------------------

if (raw_data_source == "google") {

  # Process v2 sheets
  combined_data_v2 <- map2(v2_sheets, v2_col_types, ~ process_v2_sheet(.x, .y)) %>%
    set_names(v2_sheets)

  # Process v3 sheets
  combined_data_v3 <- map(v3_sheets, ~ process_v3_sheet(.x)) %>%
    set_names(v3_sheets)

  # Handle the separate v3 inclusion sheet. Read as text: guessed types turned Issue and
  # Volume into mixed lists and some issue numbers into dates.
  data_v3_full_text_studies <- read_sheet(v3_full_text_studies, col_types = "c") %>%
    select(full_text_id, processed, extractor, decision, reason, Key, `Publication Year`, Author, Title, `Publication Title`, DOI, Issue, Volume, Rights, Extra) %>%
    mutate(`Publication Year` = as.numeric(`Publication Year`))

  # Add the v3 inclusion sheet to the list
  combined_data_v3[["inclusion_full_text"]] <- data_v3_full_text_studies

  # Save the download as RDS and as the CSV snapshot
  write_rds(combined_data_v2, here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
  write_rds(combined_data_v3, here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))
  column_types <- bind_rows(write_snapshot(combined_data_v2, "v2"), write_snapshot(combined_data_v3, "v3"))
  write_csv(column_types, types_file)

} else if (raw_data_source == "csv") {

  column_types <- read_csv(types_file, col_types = "cccc")
  combined_data_v2 <- read_snapshot("v2", column_types)
  combined_data_v3 <- read_snapshot("v3", column_types)

} else {
  stop("raw_data_source must be \"google\" or \"csv\"")
}

# Download country shapefiles ---------------------------------------------

world_shapefile <- world(path = here("data"))

# Other settings ----------------------------------------------------------
project_crs <- "EPSG:4326"

# Apply tracked corrections ------------------------------------------------
# Corrections to extracted records live in data/corrections/ and are applied
# to a copy of the download; the cleaning scripts read the corrected copy.
source(here("R", "00_01_apply_corrections.R"))
