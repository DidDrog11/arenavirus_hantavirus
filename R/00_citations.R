if (!require("pacman")) install.packages("pacman")

pkgs <- c(
  "easyPubMed",
  "here",
  "methods",
  "tidyverse",
  "XML"
)

pacman::p_load(pkgs, character.only = T)

search_date <- "2023-01-06"

if(length(list.files(here("data", "search"), pattern = search_date)) >= 1) {
  
  articles <- read_rds(here("data", "search", paste0("citations_", search_date, ".rds")))
  
} else {
  
  pubmed_search <- get_pubmed_ids("rodent* AND (hantavir* OR arenavir*)")
  
  batch_pubmed_download("rodent* AND (hantavir* OR arenavir*)", dest_dir = here("data", "search"), dest_file_prefix = paste0("citations_", search_date), format = "xml", batch_size = 3000)
  
  pubmed_records <- fetch_pubmed_data(pubmed_search, format = "xml")
  articles_list <- articles_to_list(pubmed_data = here("data", "search", paste0("citations_", search_date, "01.txt")))
  articles_df <- lapply(articles_list, article_to_df, autofill = TRUE, max_chars = 1500) %>%
    bind_rows()
  
  deduplicated <- articles_df %>%
    group_by(pmid) %>%
    slice(1)
  
  write_rds(deduplicated, here("data", "search", paste0("citations_", search_date, ".rds")))
  write_csv(deduplicated, here("data", "search", paste0("citations_", search_date, ".csv")))
  
}

