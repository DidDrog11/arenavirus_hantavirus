if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "googledrive",
  "googlesheets4",
  "here",
  "readxl",
  "stringdist",
  "taxize",
  "tidyverse"
)

pacman::p_load(pkgs, character.only = T)

# Searches were run manually on 2023-08-24 on PubMed, WoS and Google Scholar
# 2,448 citations were identified after de-duplication
search_date <- "2023-08-24"
analysis_date <- "2024-08-14"

if(length(list.files(here("data", "raw_data"), pattern = paste0(analysis_date, "_v2_data.rds"))) == 0) {
  
  data_v2 <- drive_get("https://docs.google.com/spreadsheets/d/14eW_YwSP6EWWuDrnsvDX-vwTi-7KnVyRdYL8FDqYivk/edit?usp=sharing")
  
  # Change this to use the MAIN inclusion full text document
  data_v2_full_text_studies <- data_v2 %>%
    read_sheet(sheet = "inclusion_full_text",
               col_types = "cccccciccccccccccccccc")
  write_csv(data_v2_full_text_studies, here("data", "raw_data", paste0(analysis_date, "_v2_inclusion_full_text.csv")))
  
  data_v2_descriptive <- data_v2 %>%
    read_sheet(sheet = "descriptive",
               col_types = "cccciccccccc")
  write_csv(data_v2_descriptive, here("data", "raw_data", paste0(analysis_date, "_v2_descriptive.csv")))
  data_v2_rodent <- data_v2 %>%
    read_sheet(sheet = "rodent",
               col_types = "ccccccccccdddcc")
  write_csv(data_v2_rodent, here("data", "raw_data", paste0(analysis_date, "_v2_rodent.csv")))
  data_v2_pathogen <- data_v2 %>%
    read_sheet(sheet = "pathogen",
               col_types = "ccccccccddccccddddc")
  write_csv(data_v2_pathogen, here("data", "raw_data", paste0(analysis_date, "_v2_pathogen.csv")))
  data_v2_sequences <- data_v2 %>%
    read_sheet(sheet = "sequences",
               col_types = "ccccccccccccc")
  write_csv(data_v2_sequences, here("data", "raw_data", paste0(analysis_date, "_v2_sequences.csv")))
  data_v2_known_zoonoses <- data_v2 %>%
    read_sheet(sheet = "known_zoonoses",
               col_types = "ccclccc")
  write_csv(data_v2_known_zoonoses, here("data", "raw_data", paste0(analysis_date, "_v2_known_zoonoses.csv")))
  
  combined_data_v2 <- list(citations = data_v2_full_text_studies,
                           studies = data_v2_descriptive,
                           host = data_v2_rodent,
                           pathogen = data_v2_pathogen,
                           sequence_data = data_v2_sequences,
                           zoonoses = data_v2_known_zoonoses)
  
  data_v3 <- drive_get("https://docs.google.com/spreadsheets/d/14Ghz07XOlZoaie8990mvsPbHfsHplxtfwscgf2-5KDk/edit?usp=sharing")
  data_v3_descriptive <- data_v3 %>%
    read_sheet(sheet = "descriptive",
               col_types = "cccccccc") %>%
    drop_na(study_id)
  write_csv(data_v3_descriptive, here("data", "raw_data", paste0(analysis_date, "_v3_descriptive.csv")))
  data_v3_rodent <- data_v3 %>%
    read_sheet(sheet = "rodent",
               col_types = "ccccccccdddcc")
  write_csv(data_v3_rodent, here("data", "raw_data", paste0(analysis_date, "_v3_rodent.csv")))
  data_v3_pathogen <- data_v3 %>%
    read_sheet(sheet = "pathogen",
               range = "A:L",
               col_types = "cccccccddddc") %>%
    drop_na(pathogen_record_id)
  write_csv(data_v3_pathogen, here("data", "raw_data", paste0(analysis_date, "_v3_pathogen.csv")))
  data_v3_sequences <- data_v3 %>%
    read_sheet(sheet = "sequences",
               range = "A:L",
               col_types = "cccccccccccc") %>%
    drop_na(sequence_record_id)
  write_csv(data_v3_sequences, here("data", "raw_data", paste0(analysis_date, "_v3_sequences.csv")))
  
  combined_data_v3 <- list(studies = data_v3_descriptive,
                           host = data_v3_rodent,
                           pathogen = data_v3_pathogen,
                           sequence_data = data_v3_sequences)
  
  write_rds(combined_data_v2, here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
  write_rds(combined_data_v3, here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))
  
} else {
  
  combined_data_v2 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
  combined_data_v3 <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))
  
}
