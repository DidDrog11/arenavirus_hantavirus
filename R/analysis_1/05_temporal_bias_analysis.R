# ---
# 05_temporal_bias_analysis.R
#
# Purpose: To quantify temporal biases in surveillance effort and sampled diversity.
# This script generates plots of sampling effort (hosts, studies) and
# species diversity (hosts, pathogens) over time.
# ---

# Setup -------------------------------------------------------------------
# Load necessary packages
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(scales)
library(ggrepel)
library(patchwork)
library(mgcv)

# Create output directory
output_dir <- here("output", "analysis_1")

# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-09-25.rds") # <-- Check this date
arha_db <- read_rds(db_path)

# Extract relevant tables
host_data <- arha_db$host
pathogen_data <- arha_db$pathogen
# We use citations to get a unique study count, as some citations have multiple host_data entries
citations_data <- arha_db$citations

# Data Preparation --------------------------------------------------------

# --- 1. Host Data Preparation ---
# We use sample_start_year as the primary temporal field.
# Filter for valid years and link to citations data to get publication year as a fallback
host_temporal_data <- host_data %>%
  select(study_id, host_record_id, start_date, number_of_hosts, host_species, iso3c) %>%
  left_join(citations_data %>% select(study_id, publication_year), by = "study_id") %>%
  # Coalesce sample_start_year with publication_year if sample year is missing
  mutate(year = coalesce(year(start_date), publication_year),
         decade = floor(year / 10) * 10  ) %>% # Create decade
  filter(year >= 1960, year <= 2025)

# --- 2. Pathogen Data Preparation ---
# We need to link pathogen data to host data to get the sample year
pathogen_temporal_data <- pathogen_data %>%
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) %>%
  select(host_record_id, pathogen_record_id, pathogen_family, pathogen_species_cleaned) %>%
  # Join with host_temporal_data to get the 'year'
  left_join(host_temporal_data %>% select(host_record_id, year),
            by = "host_record_id") %>%
  filter(!is.na(year))

# Analysis 1: Effort Over Time (Year & Decade) ----------------------------------

# --- 1a. Aggregate by Year ---
effort_by_year <- host_temporal_data %>%
  group_by(year) %>%
  summarise(n_hosts_sampled = sum(number_of_hosts, na.rm = TRUE),
            n_studies = n_distinct(study_id)) %>%
  ungroup()

# --- 1b. Aggregate by Decade ---
effort_by_decade <- host_temporal_data %>%
  group_by(decade) %>%
  summarise(n_hosts_sampled = sum(number_of_hosts, na.rm = TRUE),
            n_studies = n_distinct(study_id)) %>%
  ungroup()

# --- 1c. Plot: Effort by Year ---
# Create data for annotations
event_annotations <- tibble(year = c(1969, 1993),
                            label = c("Lassa virus first\ndescribed (1969)", "Sin Nombre virus\noutbreak (1993)"))

p_effort_year_hosts <- ggplot(effort_by_year, aes(x = year, y = n_hosts_sampled)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_vline(data = event_annotations, aes(xintercept = year), 
             linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_text_repel(data = event_annotations %>% mutate(y_pos = max(effort_by_year$n_hosts_sampled) * 0.9),
                  aes(x = year, y = y_pos, label = label),
                  direction = "y",
                  nudge_y = 5000,
                  nudge_x = c(-5, 5), # Nudge labels left/right
                  min.segment.length = 0.5) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(1960, 2020, 10)) +
  labs(title = "A)",
       x = "Sample Collection Year",
       y = "Total Individuals Sampled") +
  theme_minimal()

p_effort_year_studies <- ggplot(effort_by_year, aes(x = year, y = n_studies)) +
  geom_line(color = "firebrick", linewidth = 1) +
  geom_point(color = "firebrick", size = 2) +
  geom_vline(data = event_annotations, aes(xintercept = year), 
             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(1960, 2020, 10)) +
  labs(title = "B)",
       x = "Sample Collection Year",
       y = "Unique Studies") +
  theme_minimal()

# Combine and save the yearly effort plot
p_effort_combined_year <- p_effort_year_hosts / p_effort_year_studies

ggsave(filename = here(output_dir, "temporal_effort_by_year.png"),
       plot = p_effort_combined_year,
       width = 10, height = 8, dpi = 300)

# --- 1d. Plot: Effort by Decade ---
# Melt for faceted plotting
effort_by_decade_long <- effort_by_decade %>%
  pivot_longer(cols = c(n_hosts_sampled, n_studies),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric_label = factor(metric, levels = c("n_hosts_sampled", "n_studies"), labels = c("Individuals Sampled", "Unique Studies")))

p_effort_decade <- ggplot(effort_by_decade_long, aes(x = factor(decade), y = value, fill = metric_label)) +
  geom_col(position = "dodge", alpha = 0.9) +
  facet_wrap(~ metric_label, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c("Individuals Sampled" = "steelblue", "Unique Studies" = "firebrick")) +
  labs(title = "Total Sampling Effort by Decade",
       x = "Decade",
       y = "Total Count") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = here(output_dir, "temporal_effort_by_decade.png"),
       plot = p_effort_decade,
       width = 10, height = 5, dpi = 300)

# Analysis 2: Diversity Over Time -----------------------------------------
# --- 2a. Aggregate Diversity by Year ---
host_diversity_by_year <- host_temporal_data %>%
  group_by(year) %>%
  summarise(n_host_species = n_distinct(host_species)) %>%
  ungroup()

pathogen_diversity_by_year <- pathogen_temporal_data %>%
  group_by(year, pathogen_family) %>%
  summarise(n_pathogen_species = n_distinct(pathogen_species_cleaned),
            .groups = "drop")

# --- 2b. Plot: Diversity Over Time ---
p_diversity_time <- ggplot() +
  # Host species line
  geom_line(data = host_diversity_by_year,
            aes(x = year, y = n_host_species, color = "Host Species"),
            linewidth = 1.2) +
  # Pathogen species lines
  geom_line(data = pathogen_diversity_by_year,
            aes(x = year, y = n_pathogen_species, color = pathogen_family),
            linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1960, 2030, 10)) +
  scale_color_manual(name = "Taxa Sampled",
                     values = c("Host Species" = "black",
                                "Arenaviridae" = "#E69F00",
                                "Hantaviridae" = "#56B4E9")) +
  labs(x = "Sample Collection Year",
       y = "Number of Unique Species Sampled") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = here(output_dir, "temporal_diversity_over_time.png"),
       plot = p_diversity_time,
       width = 10, height = 6, dpi = 300)

# --- 3. Save Aggregated Data ---
temporal_aggregates <- list(effort_by_year = effort_by_year,
                            effort_by_decade = effort_by_decade,
                            host_diversity_by_year = host_diversity_by_year,
                            pathogen_diversity_by_year = pathogen_diversity_by_year)

write_rds(temporal_aggregates,  here("output", "analysis_1", "temporal_aggregates.rds"))


# Temporal Trends ---------------------------------------------------------
# --- 4. Continent Data ---
continents <- tibble(continent = countrycode::codelist$continent,
                     iso3c = countrycode::codelist$iso3c) %>%
  drop_na()

effort_by_year_continent <- host_temporal_data %>%
  drop_na(iso3c) %>%
  left_join(continents) %>%
  mutate(continent = as.factor(continent)) %>% 
  group_by(year, continent) %>%
  summarise(n_studies = n_distinct(study_id),
            n_hosts_sampled = sum(number_of_hosts, na.rm = TRUE)) %>%
  ungroup()

# --- 5. Fit GAM ---
gam_temporal_trends <- gam(n_hosts_sampled ~ s(year) + continent + s(year, by = continent),
                           data = effort_by_year_continent  %>% # Use pre-COVID data to build the "normal" trend
                             filter(year < 2020 & year >= 1960),
                           family = nb(),
                           method = "REML")

# --- 6. Interpretation ---
unique_continents <- effort_by_year_continent %>%
  filter(continent != "Oceania") %>%
  pull(continent) %>%
  unique()

prediction_grid <- expand.grid(year = seq(min(effort_by_year_continent$year), 
                                          max(effort_by_year_continent$year), by = 1),
                               continent = unique_continents)

predictions <- predict(gam_temporal_trends,
                       newdata = prediction_grid,
                       type = "link",
                       se.fit = TRUE)

plot_data <- data.frame(prediction_grid,
                        fit_link = predictions$fit,
                        se_link = predictions$se.fit)

invlink <- gam_temporal_trends$family$linkinv

plot_data <- plot_data %>%
  mutate(predicted_hosts = invlink(fit_link),
         conf.low_link = fit_link - (1.96 * se_link),
         conf.high_link = fit_link + (1.96 * se_link),
         conf.low = invlink(conf.low_link),
         conf.high = invlink(conf.high_link))

temporal_sampling_trends <- ggplot() +
  geom_point(data = effort_by_year_continent %>%
               filter(continent != "Oceania"),
             aes(x = year, y = n_hosts_sampled, color = continent),
             alpha = 0.6,
             size = 2) +
  geom_ribbon(data = plot_data,
              aes(x = year, ymin = conf.low, ymax = conf.high, fill = continent),
              alpha = 0.2) +
  geom_line(data = plot_data,
            aes(x = year, y = predicted_hosts, color = continent),
            linewidth = 1.2) +
  facet_wrap(~ continent, scales = "free_y") +
  labs(x = "Year",
       y = "Predicted Host Samples (Line) vs. Observed (Points)") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(filename = here(output_dir, "temporal_sampling_trends.png"),
       plot = temporal_sampling_trends,
       width = 10, height = 6, dpi = 300)
