if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "googledrive",
  "here",
  "tidyverse",
  "readxl"
)

pacman::p_load(pkgs, character.only = T)

search_date <- "2023-01-06"

if(length(list.files(here("data", "raw_data"), pattern = paste0(search_date, "_data.rds"))) == 0) {
  
  drive_download(file = "https://docs.google.com/spreadsheets/d/1BmohjsUC9rtwyULUPRrhaCjRrcOMNv5f9Zr1Qda0330/edit?usp=sharing",
                 path = here("data", "raw_data", paste0(search_date, "_data.xlsx")),
                 overwrite = TRUE)
  
  included_studies <- read_xlsx(path = here("data", "raw_data", paste0(search_date, "_data.xlsx")),
                                sheet = "descriptive")
  
  rodent_data <-  read_xlsx(path = here("data", "raw_data", paste0(search_date, "_data.xlsx")),
                            sheet = "rodent")
  
  pathogen_data <- read_xlsx(path = here("data", "raw_data", paste0(search_date, "_data.xlsx")),
                             sheet = "pathogen")
  
  sequence_data <- read_xlsx(path = here("data", "raw_data", paste0(search_date, "_data.xlsx")),
                             sheet = "pathogen_sequences")
  
  combined_data <- list(studies = included_studies,
                        rodent = rodent_data,
                        pathogen = pathogen_data,
                        sequence = sequence_data)
  
  write_rds(combined_data, here("data", "raw_data", paste0(search_date, "_data.rds")))
  
} else {
  
  combined_data <- read_rds(here("data", "raw_data", paste0(search_date, "_data.rds")))
  
}
