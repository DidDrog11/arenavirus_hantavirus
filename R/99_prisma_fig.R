if (!require("pacman")) install.packages("pacman")

pkgs <- c("here", "tidyverse", "DiagrammeR", "DiagrammeRsvg", "rsvg")

pacman::p_load(pkgs, character.only = T)

# Load the consolidated citation data
combined_data <- read_rds(here("data", "data_cleaning", "01_01_output.rds"))

# 1. Calculate PRISMA flow values
# 2448 is derived from the manual de-duplication comment in 00_load_data.R
n_identified <- 2448 
n_sought <- nrow(combined_data$citations$all_citations)
n_excluded_ta <- n_identified - n_sought

# 'no_extractions' represents papers that were sought but neither excluded nor extracted
n_not_retrieved <- nrow(combined_data$citations$no_extractions)
n_assessed <- n_sought - n_not_retrieved
n_excluded_ft <- nrow(combined_data$citations$excluded)
n_included <- length(unique(combined_data$citations$extractions$full_text_id))

# 2. Dynamically format full-text exclusion reasons
# Clean and standardize exclusion reasons
combined_data$citations$excluded <- combined_data$citations$excluded |>
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
      
      # Catch-all (should be empty based on your current data)
      TRUE ~ "Other"
    )
  )

# 2. Dynamically format full-text exclusion reasons (Updated to use reason_clean)
exclusion_reasons <- combined_data$citations$excluded |> 
  count(reason_clean) |> 
  drop_na() |> 
  arrange(desc(n)) |> 
  mutate(label = paste0(reason_clean, " (n = ", n, ")")) |> 
  pull(label) |> 
  paste(collapse = "\n")

# 3. Construct the Graphviz DOT syntax
prisma_dot <- paste0('
  digraph PRISMA {
    graph [layout = dot, rankdir = TB]
    node [shape = box, style = rounded, fontname = Helvetica, margin = 0.2]
    edge [fontname = Helvetica]

    identified [label = "Records identified from databases\n(n = ', n_identified, ')"]
    screened [label = "Records screened\n(n = ', n_identified, ')"]
    excluded_ta [label = "Records excluded at title/abstract\n(n = ', n_excluded_ta, ')"]
    sought [label = "Reports sought for retrieval\n(n = ', n_sought, ')"]
    not_retrieved [label = "Reports not retrieved\n(n = ', n_not_retrieved, ')"]
    assessed [label = "Reports assessed for eligibility\n(n = ', n_assessed, ')"]
    excluded_ft [label = "Reports excluded:\n', exclusion_reasons, '"]
    included [label = "Studies included in review\n(n = ', n_included, ')"]

    # Define the vertical flow
    identified -> screened
    screened -> sought
    sought -> assessed
    assessed -> included

    # Define the horizontal exclusions
    screened -> excluded_ta
    sought -> not_retrieved
    assessed -> excluded_ft

    # Force horizontal alignment for exclusion boxes
    {rank = same; screened; excluded_ta}
    {rank = same; sought; not_retrieved}
    {rank = same; assessed; excluded_ft}
  }
')

# 4. Render and export
dir.create(here("output", "figures"), showWarnings = FALSE, recursive = TRUE)

grViz(prisma_dot) |> 
  export_svg() |> 
  charToRaw() |> 
  rsvg_png(here("output", "figures", "supplementary_figure_s1.png"))
