if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "countrycode",
  "geodata",
  "googledrive",
  "googlesheets4",
  "here",
  "readxl",
  "stringdist",
  "terra",
  "taxize",
  "tidyverse",
  "tidyterra"
)

pacman::p_load(pkgs, character.only = T)

# Searches were run manually on 2023-08-24 on PubMed, WoS and Google Scholar
# 2,448 citations were identified after de-duplication
search_date <- "2023-08-24"
analysis_date <- "2025-01-20"

if(length(list.files(here("data", "raw_data"), pattern = paste0(analysis_date, "_v2_data.rds"))) == 0) {
  
  data_v2 <- drive_get("https://docs.google.com/spreadsheets/d/14eW_YwSP6EWWuDrnsvDX-vwTi-7KnVyRdYL8FDqYivk/edit?usp=sharing")
  
  # Change this to use the MAIN inclusion full text document
  data_v2_full_text_studies <- data_v2 %>%
    read_sheet(sheet = "inclusion_full_text",
               col_types = "cccccciccccccccccccccc") %>%
    mutate(study_id = fct_inorder(paste0("david_", str_extract(record_number, "\\d+"))))
  write_csv(data_v2_full_text_studies, here("data", "raw_data", paste0(analysis_date, "_v2_inclusion_full_text.csv")))
  
  data_v2_descriptive <- data_v2 %>%
    read_sheet(sheet = "descriptive",
               col_types = "cccciccccccc") %>%
    mutate(study_id = fct_inorder(paste0("david_", str_extract(study_id, "\\d+"))))
  write_csv(data_v2_descriptive, here("data", "raw_data", paste0(analysis_date, "_v2_descriptive.csv")))
  
  data_v2_rodent <- data_v2 %>%
    read_sheet(sheet = "rodent",
               col_types = "ccccccccccdddcc") %>%
    mutate(study_id = fct_inorder(paste0("david_", str_extract(study_id, "\\d+"))),
           rodent_record_id = fct_inorder(paste0("david_", str_extract(rodent_record_id, "\\d+"))))
  write_csv(data_v2_rodent, here("data", "raw_data", paste0(analysis_date, "_v2_rodent.csv")))
  
  data_v2_pathogen <- data_v2 %>%
    read_sheet(sheet = "pathogen",
               col_types = "ccccccccddccccddddc") %>%
    mutate(study_id = fct_inorder(paste0("david_", str_extract(study_id, "\\d+"))),
           associated_rodent_record_id = fct_inorder(paste0("david_", str_extract(associated_rodent_record_id, "\\d+"))),
           pathogen_record_id = fct_inorder(paste0("david_", str_extract(pathogen_record_id, "\\d+"))))
  write_csv(data_v2_pathogen, here("data", "raw_data", paste0(analysis_date, "_v2_pathogen.csv")))
  
  data_v2_sequences <- data_v2 %>%
    read_sheet(sheet = "sequences",
               col_types = "ccccccccccccc") %>%
    mutate(study_id = fct_inorder(paste0("david_", str_extract(study_id, "\\d+"))),
           associated_rodent_record_id = fct_inorder(if_else(!is.na(associated_rodent_record_id), paste0("david_", str_extract(associated_rodent_record_id, "\\d+")), NA)),
           associated_pathogen_record_id = fct_inorder(if_else(!is.na(associated_pathogen_record_id), paste0("david_", str_extract(associated_pathogen_record_id, "\\d+")), NA)))
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
  
  data_v3 <- list()
  
  data_v3$david <- drive_get("https://docs.google.com/spreadsheets/d/14Ghz07XOlZoaie8990mvsPbHfsHplxtfwscgf2-5KDk/edit?usp=sharing")
  data_v3$ana <- drive_get("https://docs.google.com/spreadsheets/d/10LsT76WZF4c1-LyvTKwGUsG5KH5VXIhAFcaALK6x1Cw/edit?usp=sharing")
  data_v3$grant <- drive_get("https://docs.google.com/spreadsheets/d/16pA8h33bwSMHsif6DaxfyeGKUUf0Ol9006oniEb2Gqg/edit?usp=sharing")
  data_v3$harry <- drive_get("https://docs.google.com/spreadsheets/d/1lKZitcCCP1iKGQdT8shZBewnkJoaal5flJOevOIFHzQ/edit?usp=sharing")
  data_v3$ricardo <- drive_get("https://docs.google.com/spreadsheets/d/1PVa2uF2-hG13XxGWOA_RAW0fZQ_rJKviX75McoNgsjQ/edit?usp=sharing")
  #data_v3$sadie <- drive_get("https://docs.google.com/spreadsheets/d/1vZuLT6BKBOqgJ1tg9cz2o7dHzadd9mRortJNxENmI8c/edit?usp=sharing")
  #data_v3$steph <- drive_get("https://docs.google.com/spreadsheets/d/1MQgE6bZpiG4iJhkWhphSGhsYd7t8Kko66Hmo79WKbO0/edit?usp=sharing")
  
  data_v3_descriptive <- lapply(names(data_v3), function(name) {
    read_sheet(data_v3[[name]],
               sheet = "descriptive",
               col_types = "c") %>%
      drop_na(study_id) %>%
      mutate(study_id = fct_inorder(paste0(name, "_", str_extract(study_id, "\\d+")))) %>% # Extract numbers and prefix with name
      select(any_of(c("study_id", "full_text_id", "datasetName", "sampling_effort", "data_access", "data_resolution", "linked_manuscripts", "notes")))
  }) %>%
    bind_rows()
  
  write_csv(data_v3_descriptive, here("data", "raw_data", paste0(analysis_date, "_v3_descriptive.csv")))
  
  data_v3_rodent <- lapply(names(data_v3), function(name) {
    read_sheet(data_v3[[name]],
               range = "rodent!A:M",
               col_types = "ccccccccdddcc") %>%
      drop_na(study_id) %>%
      mutate(study_id = fct_inorder(paste0(name, "_", str_extract(study_id, "\\d+"))),
             rodent_record_id = fct_inorder(paste0(name, "_", str_extract(rodent_record_id, "\\d+"))))}) %>%
    bind_rows()
    
  write_csv(data_v3_rodent, here("data", "raw_data", paste0(analysis_date, "_v3_rodent.csv")))
  
  data_v3_pathogen <- lapply(names(data_v3), function(name) {
    read_sheet(data_v3[[name]],
               range = "pathogen!A:L",
               col_types = "cccccccddddc") %>%
      drop_na(study_id, pathogen_record_id) %>%
      mutate(study_id = fct_inorder(paste0(name, "_", str_extract(study_id, "\\d+"))),
             associated_rodent_record_id = fct_inorder(paste0(name, "_", str_extract(associated_rodent_record_id, "\\d+"))),
             pathogen_record_id = fct_inorder(paste0(name, "_", str_extract(pathogen_record_id, "\\d+"))),
             )}) %>%
    bind_rows()
  
  write_csv(data_v3_pathogen, here("data", "raw_data", paste0(analysis_date, "_v3_pathogen.csv")))
  
  data_v3_sequences <- lapply(names(data_v3), function(name) {
    read_sheet(data_v3[[name]],
               range = "sequences!A:L",
               col_types = "cccccccccccc") %>%
      drop_na(study_id) %>%
      mutate(study_id = fct_inorder(paste0(name, "_", str_extract(study_id, "\\d+"))),
             associated_rodent_record_id = fct_inorder(if_else(!is.na(associated_rodent_record_id), paste0(name, "_", str_extract(associated_rodent_record_id, "\\d+")), NA)),
             associated_pathogen_record_id = fct_inorder(if_else(!is.na(associated_pathogen_record_id), paste0(name, "_", str_extract(associated_pathogen_record_id, "\\d+")), NA))
      )}) %>%
    bind_rows()
  
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


# Download country shapefiles ---------------------------------------------

world_shapefile <- world(path = here("data"))
