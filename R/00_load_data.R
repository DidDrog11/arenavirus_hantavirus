if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "googledrive",
  "here",
  "readxl",
  "stringdist",
  "tidyverse"
)

pacman::p_load(pkgs, character.only = T)

# Searches were run manually on 2023-08-24 on PubMed, WoS and Google Scholar
# 2,448 citations were identified after de-duplication
search_date <- "2023-08-24"
analysis_date <- "2024-03-21"

if(length(list.files(here("data", "raw_data"), pattern = paste0(analysis_date, "_data.rds"))) == 0) {
  
  drive_download(file = "https://docs.google.com/spreadsheets/d/14eW_YwSP6EWWuDrnsvDX-vwTi-7KnVyRdYL8FDqYivk/edit?usp=sharing",
                 path = here("data", "raw_data", paste0(analysis_date, "_data.xlsx")),
                 overwrite = TRUE)
  
  # Studies that were assessed for inclusion on review of full text
  full_text_studies <- read_xlsx(path = here("data", "raw_data", paste0(analysis_date, "_data.xlsx")),
                                 sheet = "inclusion_full_text")
  
  # Studies included after review of full text
  included_studies <- read_xlsx(path = here("data", "raw_data", paste0(analysis_date, "_data.xlsx")),
                                sheet = "descriptive")
  
  # Rodent data extracted from included studies
  rodent_data <-  read_xlsx(path = here("data", "raw_data", paste0(analysis_date, "_data.xlsx")),
                            sheet = "rodent")
  
  # Pathogen data extracted from included studies
  pathogen_data <- read_xlsx(path = here("data", "raw_data", paste0(analysis_date, "_data.xlsx")),
                             sheet = "pathogen")
  
  # Pathogen and host sequence data extracted from included studies
  sequence_data <- read_xlsx(path = here("data", "raw_data", paste0(analysis_date, "_data.xlsx")),
                             sheet = "sequences")
  
  # The status of a pathogen as a zoonosis
  zoonosis_status <- read_xlsx(path = here("data", "raw_data", paste0(analysis_date, "_data.xlsx")),
                                            sheet = "known_zoonoses")
  
  combined_data <- list(citations = full_text_studies,
                        studies = included_studies,
                        host = rodent_data,
                        pathogen = pathogen_data,
                        sequence_data = sequence_data,
                        zoonoses = zoonosis_status)
  
  write_rds(combined_data, here("data", "raw_data", paste0(analysis_date, "_data.rds")))
  
} else {
  
  combined_data <- read_rds(here("data", "raw_data", paste0(analysis_date, "_data.rds")))
  
}
