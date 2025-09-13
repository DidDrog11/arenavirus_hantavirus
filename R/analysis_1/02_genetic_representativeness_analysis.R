# ---
# 02_genetic_representativeness_analysis.R
#
# Purpose: To assess the representativeness of the available genetic data
# (sequences) in relation to the overall geographic and host sampling coverage.
# This script identifies key gaps where sampling has occurred but from which
# sequence data is not available.
# ---

# --- 1. Setup ---
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(terra)
library(tidyterra)
library(biscale)
library(cowplot)
library(patchwork)
library(ggnewscale)
library(forcats)

# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-08-25.rds")
arha_db <- read_rds(db_path)
host_data <- arha_db$host
pathogen_data <- arha_db$pathogen
sequence_data <- arha_db$sequence

# --- 2. Geographic Representativeness Analysis ---

hosts_per_country <- host_data %>%
  filter(!is.na(iso3c)) %>%
  group_by(country, iso3c) %>%
  summarise(total_hosts_sampled = sum(number_of_hosts, na.rm = TRUE), .groups = "drop")

sequences_missing_country <- sequence_data %>%
  filter(is.na(iso3c)) %>%
  left_join(host_data, by = "host_record_id") %>%
  select(sequence_record_id, accession_primary, iso3c.x, iso3c.y, country.x, country.y, sequence_type) %>%
  mutate(iso3c = coalesce(iso3c.x, iso3c.y),
         country = coalesce(country.x, country.y)) %>%
  group_by(country, iso3c, sequence_type) %>%
  summarise(n_records_with_sequences = n_distinct(accession_primary), .groups = "drop") %>%
  tidyr::pivot_wider(
           names_from = sequence_type,
           values_from = n_records_with_sequences,
           values_fill = 0,
           names_prefix = "n_records_"
         )

sequences_per_country <- sequence_data %>%
  group_by(country, iso3c, sequence_type) %>%
  summarise(n_records_with_sequences = n_distinct(accession_primary), .groups = "drop") %>%
  # Pivot to have separate columns for Host and Pathogen sequences
  tidyr::pivot_wider(
    names_from = sequence_type,
    values_from = n_records_with_sequences,
    values_fill = 0,
    names_prefix = "n_records_"
  )

sequences_per_country_stratified <- bind_rows(sequences_missing_country, sequences_per_country) %>%
  group_by(country, iso3c) %>%
  summarise(n_records_Pathogen = sum(n_records_Pathogen, na.rm = TRUE),
            n_records_Host = sum(n_records_Host, na.rm = TRUE),
            n_records_NA = sum(n_records_NA, na.rm = TRUE), .groups = "drop")

geo_rep_data <- hosts_per_country %>%
  left_join(sequences_per_country_stratified, by = c("country", "iso3c")) %>%
  mutate(
    # Coalesce NA counts to 0
    n_records_Host = coalesce(n_records_Host, 0L),
    n_records_Pathogen = coalesce(n_records_Pathogen, 0L),
    # Calculate proportions
    prop_host_sequenced = n_records_Host / total_hosts_sampled,
    prop_host_sequenced = case_when(prop_host_sequenced >= 1 ~ 1,
                                    TRUE ~ prop_host_sequenced),
    prop_pathogen_sequenced = n_records_Pathogen / total_hosts_sampled,
    prop_pathogen_sequenced = case_when(prop_pathogen_sequenced >= 1 ~ 1,
                                        TRUE ~ prop_pathogen_sequenced)
  )

# Join with a world map shapefile
world_map <- ne_countries(scale = "medium", returnclass = "sf")
map_data <- vect(world_map) %>%
  left_join(geo_rep_data, by = c("iso_a3" = "iso3c"))

# Create a helper function to avoid repeating plot code
create_binned_map <- function(map_df, fill_var, title) {
  
  # Ensure the fill variable exists
  if (!fill_var %in% names(map_df)) {
    stop(paste("Fill variable", fill_var, "not found in data."))
  }
  
  # Define the refined breaks and labels
  breaks <- c(-Inf, 0, 0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.05, 0.1, 0.4, 0.999, Inf)
  labels <- c("0%", ">0 - 0.1%", "0.1% - 0.2%", "0.2% - 0.4%", "0.4% - 0.6%", 
              "0.6% - 0.8%", "0.8% - 1%", "1% - 5%", "5% - 10%", "10% - 40%", "40% - 100%", "100%")
  
  plot_data <- map_df %>%
    mutate(
      # Create a new binned factor variable from the continuous proportion
      completeness_bin_temp = cut(
        .data[[fill_var]],
        breaks = breaks,
        labels = labels,
        right = TRUE, 
        include.lowest = TRUE
      ),
      completeness_bin = factor(completeness_bin_temp, levels = labels)
    )

  num_gradient_colors <- length(labels) - 2 # We need 10 colours for the gradient
  color_ramp <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))
  palette <- color_ramp(num_gradient_colors)
  
  map_colors <- c(
    "0%" = "grey90",
    ">0 - 0.1%" = palette[1],
    "0.1% - 0.2%" = palette[2],
    "0.2% - 0.4%" = palette[3],
    "0.4% - 0.6%" = palette[4],
    "0.6% - 0.8%" = palette[5],
    "0.8% - 1%" = palette[6],
    "1% - 5%" = palette[7],
    "5% - 10%" = palette[8],
    "10% - 40%" = palette[9],
    "40% - 100%" = palette[10],
    "100%" = "gold"
  )
  
  p <- ggplot(data = plot_data) +
    geom_spatvector(aes(fill = completeness_bin), colour = "black", linewidth = 0.1) +
    scale_fill_manual(
      name = "Sequencing Completeness",
      values = map_colors,
      na.translate = FALSE,
      drop = FALSE
    ) +
    labs(title = title) +
    coord_sf(crs = "EPSG:8857") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

# Generate the two maps using the helper function
host_geo_map_final <- create_binned_map(map_data, "prop_host_sequenced", "Host Sequences")
pathogen_geo_map_final <- create_binned_map(map_data, "prop_pathogen_sequenced", "Pathogen Sequences")

shared_legend <- get_legend(
  host_geo_map_final + 
    guides(fill = guide_legend(nrow = 2)) # Ensure legend is a single row
)

host_geo_map_final <- host_geo_map_final + theme(legend.position = "none")
pathogen_geo_map_final <- pathogen_geo_map_final + theme(legend.position = "none")

plots_row <- host_geo_map_final + pathogen_geo_map_final

final_combined_plot <- plot_grid(
  plots_row, 
  shared_legend, 
  ncol = 1, 
  rel_heights = c(1, 0.1) # Give 90% of the height to the plots, 10% to the legend
)

final_combined_plot_with_title <- final_combined_plot +
  plot_annotation(
    title = "Geographic Representativeness of Genetic Data",
    caption = "Colour reflects the ratio of unique host records with sequences to the total number of hosts sampled per country."
  )

print(final_combined_plot_with_title)
ggsave(here("output", "analysis_1", "geographic_representativeness_sequencing.png"), plot = final_combined_plot_with_title, width = 16, height = 6, dpi = 300)

# --- 2.5. Visualize Bivariate Geographic Representativeness ---
breaks <- c(0, 0.01, 0.10, 0.50, 1.0)
labels <- c(">0% - 1%", "1% - 10%", "10% - 50%", "50% - 100%")

bivariate_manual_data <- map_data %>%
  mutate(
    # Create the binned factors for host and pathogen completeness
    host_bin = cut(prop_host_sequenced, breaks = breaks, labels = labels, include.lowest = TRUE, right = TRUE),
    pathogen_bin = cut(prop_pathogen_sequenced, breaks = breaks, labels = labels, include.lowest = TRUE, right = TRUE),
    
    # Create the final, comprehensive category for plotting
    map_category = case_when(
      is.na(total_hosts_sampled) ~ "No Sampling Data",
      prop_host_sequenced == 0 & prop_pathogen_sequenced == 0 ~ "Sampled, 0% Completeness",
      # For all other cases, create a combined factor level
      TRUE ~ as.character(interaction(host_bin, pathogen_bin, sep = " - "))
    )
  )

# Get the 16 color hex codes from the palette
bivariate_palette_colors <- bi_pal("BlueGold", dim = 4, preview = FALSE)

# Get the names of all 16 possible bivariate levels
bivariate_levels <- levels(interaction(factor(labels), factor(labels), sep = " - "))

# Name the palette colors with the bivariate levels
names(bivariate_palette_colors) <- bivariate_levels

# Create the final, complete named vector of colors for the scale
final_map_colors <- c(
  "No Sampling Data" = "white",
  "Sampled, 0% Completeness" = "black",
  bivariate_palette_colors
)

bivariate_map <- ggplot() +
  # Use a single geom_sf layer and map the fill to our new comprehensive category
  geom_spatvector(data = bivariate_manual_data, aes(fill = map_category), color = "black", linewidth = 0.1, show.legend = FALSE) +
  # Use the fully named color vector in scale_fill_manual
  scale_fill_manual(
    values = final_map_colors,
    na.value = "black" # Fallback for categories not in the data
  ) +
  labs(
    title = "Bivariate Geographic Representativeness of Genetic Data",
    subtitle = "Comparing sequencing completeness for hosts and pathogens",
    caption = "Countries in black contain host and pathogen sampling but no sequence data"
  ) +
  coord_sf(crs = "EPSG:8857") +
  theme_minimal()

legend_data <- expand.grid(x = labels, y = labels)

bivariate_legend <- ggplot() +
  geom_tile(data = legend_data, aes(x = x, y = y, fill = interaction(x, y, sep = " - ")), color = "white", linewidth = 0.5) +
  scale_fill_manual(values = bivariate_palette_colors, guide = "none") +
  labs(x = "Host Completeness ->", y = "Pathogen Completeness ->") +
  theme_void() +
  theme(
    axis.title = element_text(size = 10),
    axis.title.y = element_text(angle = 90),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  coord_fixed()

final_bivariate_plot <- ggdraw() +
  draw_plot(bivariate_map, 0, 0, 1, 1) +
  draw_plot(bivariate_legend, x = 0.1, y = 0.12, width = 0.25, height = 0.25)

ggsave(here("output", "analysis_1", "bivariate_plot_genetic_data.png"), plot = final_bivariate_plot, width = 12, height = 8, dpi = 300)

# --- 3. Taxonomic Representativeness Analysis ---
hosts_per_species <- host_data %>%
  drop_na(host_species) %>%
  group_by(host_species) %>%
  summarise(total_hosts_sampled = sum(number_of_hosts, na.rm = TRUE), .groups = "drop")

sequences_per_species <- sequence_data %>%
  left_join(host_data %>% select(host_record_id, host_species, number_of_hosts), by = "host_record_id") %>%
  filter(!is.na(host_species)) %>%
  group_by(host_species) %>%
  summarise(n_hosts_with_sequences = n(), .groups = "drop")

taxo_rep_data <- hosts_per_species %>%
  left_join(sequences_per_species, by = "host_species") %>%
  mutate(
    n_hosts_with_sequences = coalesce(n_hosts_with_sequences, 0L),
    proportion_sequenced = n_hosts_with_sequences / total_hosts_sampled
  ) %>%
  # Filter for the most-sampled hosts for a clear plot
  filter(total_hosts_sampled > 1000) %>%
  arrange(proportion_sequenced) %>%
  mutate(host_species = fct_inorder(host_species))

taxo_rep_plot <- ggplot(taxo_rep_data, aes(x = host_species, y = proportion_sequenced)) +
  geom_segment(aes(xend = host_species, yend = 0), linewidth = 0.8) +
  geom_point(aes(size = total_hosts_sampled), shape = 21, fill = "skyblue", color = "black") +
  scale_y_continuous(labels = scales::percent) +
  scale_size_continuous(name = "Total Hosts Sampled", labels = scales::comma) +
  coord_flip() +
  labs(
    title = "Taxonomic Representativeness of Genetic Data",
    subtitle = "For host species with >1000 individuals sampled globally",
    x = "Host Species",
    y = "Proportion of Sampled Hosts with Sequences"
  ) +
  theme_minimal()

print(taxo_rep_plot)
ggsave(here("output", "analysis_1", "taxonomic_representativeness_plot.png"), plot = taxo_rep_plot, width = 10, height = 8)

# --- 3b Joint Taxonomic Representativeness Analysis ---
hosts_per_species <- host_data %>%
  drop_na(host_species) %>%
  group_by(host_species) %>%
  summarise(total_hosts_sampled = sum(number_of_hosts, na.rm = TRUE), .groups = "drop")

sequences_per_species_stratified <- sequence_data %>%
  left_join(host_data %>% select(host_record_id, host_species, number_of_hosts), by = "host_record_id") %>%
  filter(!is.na(host_species), !is.na(sequence_type)) %>%
  group_by(host_species, sequence_type) %>%
  summarise(n_hosts_with_sequences = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = sequence_type,
    values_from = n_hosts_with_sequences,
    values_fill = 0
  )

taxo_rep_data_host <- hosts_per_species %>%
  left_join(sequences_per_species_stratified, by = "host_species") %>%
  mutate(
    n_sequence_host = coalesce(Host, 0L),
    n_sequence_pathogen = coalesce(Pathogen, 0L),
    prop_host_sequenced = n_sequence_host / total_hosts_sampled,
    prop_pathogen_sequenced = n_sequence_pathogen / total_hosts_sampled
  ) %>%
  filter(total_hosts_sampled > 1000) %>%
  mutate(host_species = paste0(host_species, " (N=", prettyNum(total_hosts_sampled, big.mark = ",", digits = 6), ")"),
                                                                     host_species = fct_reorder(host_species, total_hosts_sampled))

host_taxo_rep_plot <- ggplot(taxo_rep_data_host, aes(y = host_species)) +
  geom_segment(aes(x = prop_pathogen_sequenced, xend = prop_host_sequenced, yend = host_species), color = "grey", linewidth = 1) +
  geom_point(aes(x = prop_pathogen_sequenced, color = "Pathogen"), size = 4) +
  geom_point(aes(x = prop_host_sequenced, color = "Host"), size = 4) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(name = "Sequence Type", values = c("Pathogen" = "firebrick", "Host" = "darkgreen")) +
  labs(title = "Host Taxonomic Representativeness of Genetic Data", subtitle = "For host species with >1000 individuals sampled globally", x = "Proportion of Sampled Hosts with Sequences", y = "Host Species") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(here("output", "analysis_1", "host_taxonomic_representativeness_plot.png"), plot = host_taxo_rep_plot, width = 10, height = 8)

# --- 3.5. Visualize PATHOGEN Taxonomic Representativeness ---
pcr_positives_per_pathogen <- pathogen_data %>%
  filter(pathogen_family %in% c("Hantaviridae", "Arenaviridae")) %>%
  filter(taxonomic_level == "species") %>%
  filter(stringr::str_detect(assay, "PCR|Culture|Sequencing")) %>%
  filter(number_positive > 0) %>%
  group_by(pathogen_species_cleaned) %>%
  summarise(
    total_pcr_positives = sum(number_positive, na.rm = TRUE),
    .groups = "drop"
  )

sequences_per_pathogen <- sequence_data %>%
  filter(sequence_type == "Pathogen") %>%
  group_by(pathogen, pathogen_species_clean) %>%
  mutate(match_1 = case_when(!is.na(pathogen) & pathogen %in% detections_per_pathogen$pathogen_species_cleaned ~ TRUE,
                             TRUE ~ FALSE),
         match_2 = case_when(!is.na(pathogen_species_clean) & pathogen_species_clean %in% detections_per_pathogen$pathogen_species_cleaned ~ TRUE,
                             TRUE ~ FALSE),
         pathogen_species_cleaned = case_when(match_1 == TRUE ~ pathogen,
                                              match_2 == TRUE ~ pathogen_species_clean,
                                              TRUE ~ NA))  %>%
  group_by(pathogen) %>%
  arrange(pathogen_species_cleaned) %>%
  fill(pathogen_species_cleaned, .direction = "down") %>%
  drop_na(pathogen_species_cleaned) %>%
  group_by(pathogen_species_cleaned) %>%
  count(pathogen_species_cleaned, name = "n_sequences")

taxo_rep_data_pathogen <- pcr_positives_per_pathogen %>%
  filter(total_pcr_positives > 5) %>%
  full_join(sequences_per_pathogen, by = "pathogen_species_cleaned") %>%
  mutate(
    total_pcr_positives = coalesce(total_pcr_positives, 0L),
    n_sequences = coalesce(n_sequences, 0L),
    proportion_sequenced_raw = n_sequences / total_pcr_positives,
    proportion_sequenced = if_else(proportion_sequenced_raw > 1, 1, proportion_sequenced_raw)
  ) %>%
  mutate(pathogen_species_cleaned = fct_reorder(pathogen_species_cleaned, total_pcr_positives))

pathogen_taxo_rep_plot <- ggplot(taxo_rep_data_pathogen, aes(x = proportion_sequenced, y = pathogen_species_cleaned)) +
  geom_col(aes(fill = total_pcr_positives), width = 0.8) +
  geom_text(aes(label = scales::percent(proportion_sequenced, accuracy = 0.1)), 
            hjust = -0.1, size = 3) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.1), expand = c(0, 0)) +
  scale_fill_viridis_c(
    name = "Total PCR+ Detections", 
    trans = "log1p",
    labels = scales::comma
  ) +
  labs(
    title = "Pathogen Genetic Data Representativeness",
    subtitle = "Proportion of PCR-positive detections with at least one sequence in GenBank",
    x = "Proportion of PCR+ Detections Sequenced",
    y = "Pathogen Species"
  ) +
  theme_minimal()

print(pathogen_taxo_rep_plot)
ggsave(here("output", "analysis_1", "pathogen_taxonomic_representativeness_plot.png"), 
       plot = pathogen_taxo_rep_plot, width = 10, height = 8)
