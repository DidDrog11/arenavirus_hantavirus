if (!require("pacman")) install.packages("pacman")

pkgs <- c("here", "tidyverse", "DiagrammeR", "DiagrammeRsvg", "rsvg", "arrow")

pacman::p_load(pkgs, character.only = T)

analysis_date <- "2026-10-20" # as in 00_load_data.R

# 1. Search and title/abstract screening, from the published protocol (Simons et al. 2025,
# PRISMA figure): PubMed and Web of Science searched on 2023-10-16
n_database    <- 4020 # records identified through database searching
n_other       <- 11   # records identified from other sources at that time
n_deduplicated <- 2755 # records after de-duplication (all screened)
n_excluded_ta <- 1838 # records excluded at title and abstract
n_from_screening <- n_deduplicated - n_excluded_ta

# 2. Full-text stage, from the inclusion sheet (one row per report sought) and the extractions.
# Papers added by hand from reference lists have IDs starting "ft_m_"; those beyond the 11 in
# the protocol were added after screening.
inclusion <- read_rds(here("data", "raw_data", paste0(analysis_date, "_v3_data_corrected.rds")))$inclusion_full_text
extraction_status <- read_rds(here("data", "data_cleaning", "02_01_output.rds"))$extraction_status
sampling_events <- read_parquet(here("data", "data_release", "parquet", "sampling_events.parquet"))

# A report marked processed but never extracted was assessed and excluded, with no reason recorded;
# one never processed was not assessed
inclusion <- inclusion |>
  mutate(status = case_when(
    full_text_id %in% extraction_status$full_text_id ~ "extracted",
    str_detect(coalesce(decision, ""), regex("^Exclude", ignore_case = TRUE)) ~ "excluded",
    tolower(coalesce(processed, "")) == "y" ~ "excluded",
    TRUE ~ "not assessed"))

n_sought        <- nrow(inclusion)
n_added_later   <- n_sought - n_from_screening
stopifnot(n_added_later == sum(str_detect(inclusion$full_text_id, "^ft_m_"), na.rm = TRUE) - n_other)
n_not_assessed  <- sum(inclusion$status == "not assessed")
n_excluded_ft   <- sum(inclusion$status == "excluded")
n_extracted     <- sum(inclusion$status == "extracted")
# Whole studies removed in 05_03 as repeats of data in another included paper, or not primary data
n_in_release    <- nrow(sampling_events)
n_removed_after <- n_extracted - n_in_release

# 3. Full-text exclusion reasons, cleaned and standardised
excluded <- inclusion |>
  filter(status == "excluded") |>
  mutate(
    reason_clean = case_when(
      # 1. Full text / Abstract only
      str_detect(tolower(reason), "unable|abstract|conference") ~ "Full text unavailable / Abstract only",

      # 2. Experimental / Captive
      str_detect(tolower(reason), "experimental|colony|captive|laboratory|not wild") ~ "Experimental or captive study",

      # 3. Duplicate / Data previously published
      str_detect(tolower(reason), "elsewhere|duplicate|reanalysis|linked") ~ "Data previously published",

      # 4. Unextractable format
      str_detect(tolower(reason), "graphical|extract|individual level|species level|underlying data") ~ "Unextractable data format",

      # 5. No target pathogen (Prioritised over 'no host' to catch "no rodent associated pathogen")
      str_detect(tolower(reason), "pathogen|hanta|aren|tick") ~ "No target pathogen",

      # 6. No target host / No field sampling
      str_detect(tolower(reason), "human|owl|rodent|capture|sampling") ~ "No target host / No field sampling",

      # 7. No primary data
      str_detect(tolower(reason), "review|meta-analysis|letter|nomenclature|modelling|no data|testing") ~ "No primary data (e.g., review, modelling)",

      # Reason missing or not matched
      TRUE ~ "Other / not recorded"
    )
  )

exclusion_reasons <- excluded |>
  count(reason_clean) |>
  arrange(desc(n)) |>
  mutate(label = paste0(reason_clean, " (n = ", n, ")")) |>
  pull(label) |>
  paste(collapse = "\n")

fmt <- function(n) format(n, big.mark = ",")

# 4. Construct the Graphviz DOT syntax
prisma_dot <- paste0('
  digraph PRISMA {
    graph [layout = dot, rankdir = TB]
    node [shape = box, style = rounded, fontname = Arial, fontsize = 12, margin = "0.3,0.15"]
    edge [fontname = Arial]

    database [label = "Records identified from databases\n(n = ', fmt(n_database), ')"]
    other [label = "Records identified from other sources\n(n = ', fmt(n_other), ')"]
    deduplicated [label = "Records after de-duplication\n(n = ', fmt(n_deduplicated), ')"]
    screened [label = "Records screened\n(n = ', fmt(n_deduplicated), ')"]
    excluded_ta [label = "Records excluded at title/abstract\n(n = ', fmt(n_excluded_ta), ')"]
    added_later [label = "Reports added later\nfrom reference lists\n(n = ', fmt(n_added_later), ')"]
    sought [label = "Reports sought for full-text assessment\n(n = ', fmt(n_sought), ')"]
    not_assessed [label = "Reports not assessed\n(n = ', fmt(n_not_assessed), ')"]
    assessed [label = "Reports assessed for eligibility\n(n = ', fmt(n_sought - n_not_assessed), ')"]
    excluded_ft [label = "Reports excluded (n = ', fmt(n_excluded_ft), '):\n', exclusion_reasons, '"]
    extracted [label = "Reports from which data were extracted\n(n = ', fmt(n_extracted), ')"]
    removed_after [label = "Removed after extraction: data repeated\nfrom another included report,\nor not primary data (n = ', fmt(n_removed_after), ')"]
    included [label = "Studies in ArHa\n(n = ', fmt(n_in_release), ')"]

    # Vertical flow
    database -> deduplicated
    other -> deduplicated
    deduplicated -> screened
    screened -> sought [weight = 10]
    added_later -> sought
    sought -> assessed [weight = 10]
    assessed -> extracted [weight = 10]
    extracted -> included [weight = 10]

    # Horizontal exclusions
    screened -> excluded_ta
    sought -> not_assessed
    assessed -> excluded_ft
    extracted -> removed_after

    # Force horizontal alignment
    {rank = same; database; other}
    {rank = same; screened; excluded_ta}
    {rank = same; added_later; sought; not_assessed}
    {rank = same; assessed; excluded_ft}
    {rank = same; extracted; removed_after}
  }
')

# 5. Render and export
dir.create(here("output", "figures"), showWarnings = FALSE, recursive = TRUE)

grViz(prisma_dot) |>
  export_svg() |>
  charToRaw() |>
  rsvg_png(here("output", "figures", "supplementary_figure_s1.png"))
