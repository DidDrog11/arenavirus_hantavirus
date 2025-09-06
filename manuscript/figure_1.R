library(here)
library(rnaturalearth)
library(rnaturalearthdata)
library(terra)
library(tidyverse)
library(tidyterra)

# read in data and unlist
data <- read_rds(here("data", "clean_data", "2025-03-05_data.rds"))

included_studies <- data$descriptive
citations <- data$citations |>
  filter(full_text_id %in% included_studies$full_text_id)
host <- data$host |>
  ungroup()
pathogen <- data$pathogen 


# Figure 1A - Research Effort ---------------------------------------------

citations_count <- citations |>
  group_by(`Publication Year`) |>
  summarise(n_citations = n())

detections_count <- included_studies |>
  ungroup() |>
  distinct(full_text_id, n_assays) |>
  left_join(citations |>
              select(`Publication Year`, full_text_id) |>
              distinct(full_text_id, `Publication Year`),
            by = "full_text_id") |>
  group_by(`Publication Year`) |>
  summarise(n_assays = sum(n_assays, na.rm = TRUE))

ylim_assay <- c(0, max(detections_count$n_assays, na.rm = TRUE))
ylim_citations <- c(0, max(citations_count$n_citations, na.rm = TRUE))
b <- diff(ylim_assay)/diff(ylim_citations)
a <- ylim_assay[1] - b*ylim_citations[1]

research_effort <- detections_count |>
  ggplot() +
  geom_col(aes(x = `Publication Year`, y = n_assays)) +
  geom_line(data = citations_count, aes(x = `Publication Year`, y = a + n_citations*b), colour = "blue") +
  scale_y_continuous("Number of assays", sec.axis = sec_axis(~ (. - a)/b, name = "Number of publications")) +
  theme_minimal()


# Figure 1B - Phylogeny of rodents and shrews with testing ----------------


# Figure 1C - Locations of sampling ---------------------------------------

world

pathogen 


