# ---
# 03_geographic_temporal_bias_analysis.R
#
# Purpose: To quantify geographic and temporal biases in surveillance effort.
# This script generates maps of sampling locations, compares them to host ranges,
# correlates effort with socio-economic data, and plots temporal trends.
# ---

# --- 1. Setup ---
# Load necessary packages
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(stringr)
library(geodata)
library(rnaturalearth)
library(lubridate)
library(WDI) # For World Bank socio-economic data
library(taxize)

# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-08-25.rds") # <-- Update this date
host_traits <- read_rds(here("data", "external", "harmonised_species.rds"))
arha_db <- read_rds(db_path)
host_data <- arha_db$host
pathogen_data <- arha_db$pathogen

# --- 2. Load Additional Spatial Data ---
# a) Host Range Maps (from IUCN)
try({
  iucn_ranges <- vect(here("data", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
}, silent = TRUE)

iucn_subset <- iucn_ranges[iucn_ranges$`order_` %in% c("RODENTIA", "SORICOMORPHA")]

iucn_species_names <- tibble(iucn_species_name = iucn_subset$sci_name,
                             iucn_genus_name = iucn_subset$genus) %>%
  distinct() %>%
  left_join(host_traits %>%
              select(scientific_name, gbif_genus, gbif_species, species_gbifid, genus_gbifid),
            by = c("iucn_species_name" = "scientific_name")) %>%
  group_by(iucn_genus_name) %>%
  fill(gbif_genus, .direction= c("down")) %>%
  fill(genus_gbifid, .direction= c("down")) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(match_1 = case_when(iucn_genus_name == gbif_genus ~ TRUE,
                             TRUE ~ FALSE))

iucn_subset <- left_join(iucn_subset, iucn_species_names %>%
                           select(iucn_species_name, gbif_species, species_gbifid, gbif_genus, genus_gbifid, match_1),
                         by = c("sci_name" = "iucn_species_name"))

# b) Global ADM2 Shapefile
# This assumes you have a global ADM2 shapefile from GADM or another source
try({
  gadm_adm2 <- gadm(country = sort(unique(host_data %>% 
                                            filter(!is.na(gadm_adm2)) %>% 
                                                     pull(iso3c))),
                                          level = 2, path = here("data", "external"), version = "latest")
}, silent = TRUE)

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  select(iso_a3) %>%
  filter(!iso_a3 %in% gadm_adm2$GID_0) %>%
  vect() %>%
  project(project_crs)  %>%
  rename(GID_0 = iso_a3)

gadm_adm2 <- rbind(gadm_adm2, world_map)

# --- 3. Geographic Bias Analysis ---

# 3.1. Quantitative Assessment of Sampling Gaps in Host Ranges
if (exists("iucn_ranges") && exists("gadm_adm2")) {
  message("Performing quantitative analysis of host range sampling gaps using ADM2 units...")
  
  # Ensure all layers are in the same, valid projection for spatial operations
  project_crs <- crs("EPSG:8857")
  gadm_adm2_proj <- project(gadm_adm2, y = project_crs)
  
  # Get a unique list of all ADM2 units that have been sampled in our database
  sampled_adm2_units <- host_data %>%
    filter(!is.na(gadm_adm2) & !is.na(country)) %>%
    distinct(country, gadm_adm1, gadm_adm2)
  
  # --- 3.2.a. Analysis by Genus ---
  all_genera <- host_data %>%
    filter(host_order %in% c("Rodentia", "Soricomorpha")) %>%
    drop_na(host_genus) %>%
    group_by(host_genus) %>%
    summarise(n = sum(number_of_hosts, na.rm = TRUE)) %>%
    arrange(-n) %>%
    pull(host_genus)
  
  genus_range_coverage_summary <- purrr::map_dfr(all_genera[1], function(genus_name) {
    
    message(paste("... processing genus:", genus_name))
    
    genus_range <- iucn_subset %>%
      filter(str_detect(gbif_genus, genus_name)) %>%
      aggregate() %>%
      project(y = project_crs)
    
    cropped_adm2 <- gadm_adm2_proj %>%
      crop(ext(genus_range))
    
    inhabited_adm2 <- cropped_adm2 %>%
      relate(genus_range, relation = "intersects")
    
    inhabited_adm2_gadm <- cropped_adm2[inhabited_adm2]
    
    sampled_inhabited_adm2 <- host_data %>%
      filter(str_detect(host_genus, genus_name)) %>%
      group_by(gadm_adm2) %>%
      summarise(number_of_hosts = sum(number_of_hosts, na.rm = TRUE),
                trap_nights_clean = sum(trap_nights_clean, na.rm = TRUE)) %>%
      mutate(sampled = TRUE)
    
    inhabited_adm2_gadm <- left_join(inhabited_adm2_gadm, sampled_inhabited_adm2, by = c("GID_2" = "gadm_adm2"))
    
    total_inhabited_area_km2 <- sum(expanse(inhabited_adm2_gadm)) / 1e6
    total_sampled_area_km2 <- sum(expanse(inhabited_adm2_gadm %>%
                                            filter(sampled == TRUE))) / 1e6
    
    sampled_area <- tibble(genus = genus_name,
                           proportion_of_area_sampled = as.numeric(total_sampled_area_km2 / total_inhabited_area_km2))
    sampled_area_vect <- inhabited_adm2_gadm
    
    print("--- Proportion of Inhabited Geographic Area Sampled (by Genus) ---")
    print(sampled_area)
    
    return(list(sampled_area,
                sampled_area_vect))
    
  })
  
  write_csv(genus_range_coverage_summary$prop_area_sampled, here("R", "analysis_1", "outputs", "host_range_coverage_summary_genus.csv"))
  
}
  
  # --- 3.2.b. Analysis by Top 5 Sampled Species ---
  # First, identify the top 5 most sampled species from your database
  top_5_species <- host_data %>%
    group_by(host_species) %>%
    summarise(total_sampled = sum(number_of_hosts, na.rm = TRUE)) %>%
    slice_max(order_by = total_sampled, n = 5) %>%
    pull(host_species)
  
  message(paste("\nIdentified top 5 sampled species:", paste(top_5_species, collapse = ", ")))
  
  species_range_coverage_summary <- purrr::map_dfr(top_5_species, function(species_name) {
    
    message(paste("... processing species:", species_name))
    
    # Get the total range for the specific species
    species_range <- iucn_ranges %>%
      filter(sci_name == species_name, presence %in% c(1,2)) %>%
      st_union() %>%
      st_transform(crs = mollweide_crs)
    
    # Check if a range was found
    if (st_is_empty(species_range)) {
      warning(paste("No valid IUCN range found for", species_name))
      return(tibble(species = species_name, proportion_of_area_sampled = NA_real_))
    }
    
    inhabited_adm2 <- gadm_adm2_proj %>%
      st_filter(species_range, .predicate = st_intersects)
    
    sampled_inhabited_adm2 <- inhabited_adm2 %>%
      inner_join(sampled_adm2_units, by = c("NAME_0" = "country", "NAME_2" = "gadm_adm2"))
    
    total_inhabited_area_km2 <- sum(st_area(inhabited_adm2)) / 1e6
    total_sampled_area_km2 <- sum(st_area(sampled_inhabited_adm2)) / 1e6
    
    tibble(
      species = species_name,
      proportion_of_area_sampled = as.numeric(total_sampled_area_km2 / total_inhabited_area_km2)
    )
  })
  
  print("\n--- Proportion of Inhabited Geographic Area Sampled (by Top 5 Species) ---")
  print(species_range_coverage_summary)
  write_csv(species_range_coverage_summary, here("R", "analysis_1", "outputs", "host_range_coverage_summary_species.csv"))
}


# 3.3. Correlate sampling intensity with socio-economic data
# (This section is unchanged from the previous version)
try({
  socio_economic_data <- WDI(
    country = "all", indicator = "NY.GDP.PCAP.CD", 
    start = 2020, end = 2020, extra = TRUE
  ) %>% as_tibble() %>% filter(region != "Aggregates") %>%
    select(iso3c, gdp_per_capita = NY.GDP.PCAP.CD)
  
  country_sampling_intensity <- host_data %>%
    group_by(iso3c) %>%
    summarise(total_hosts_sampled = sum(number_of_hosts, na.rm = TRUE), .groups = "drop") %>%
    left_join(socio_economic_data, by = "iso3c") %>%
    filter(!is.na(gdp_per_capita), !is.na(total_hosts_sampled))
  
  gdp_model <- lm(log1p(total_hosts_sampled) ~ log10(gdp_per_capita), data = country_sampling_intensity)
  print(summary(gdp_model))
}, silent = TRUE)


# --- 4. Temporal Bias Analysis ---
# (This section is unchanged from the previous version)

# 4.1. Plot number of hosts sampled per year (with annotations)
temporal_summary <- host_data %>%
  filter(!is.na(start_date)) %>%
  mutate(year = year(start_date)) %>%
  group_by(year) %>%
  summarise(n_hosts_sampled = sum(number_of_hosts, na.rm = TRUE)) %>%
  filter(year >= 1960, year <= 2024) 

temporal_effort_plot <- ggplot(temporal_summary, aes(x = year, y = n_hosts_sampled)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1993, linetype = "dashed", color = "red") +
  annotate("text", x = 1993, y = max(temporal_summary$n_hosts_sampled) * 0.9,
           label = "Sin Nombre virus outbreak", angle = 90, vjust = -0.5, color = "red") +
  labs(
    title = "Total Number of Hosts Sampled Over Time",
    x = "Year of Sample Collection",
    y = "Number of Hosts Sampled"
  ) +
  theme_minimal()

print(temporal_effort_plot)

# 4.2. Analyze trends in the DIVERSITY of taxa sampled over time
first_year_sampled <- host_data %>%
  filter(!is.na(start_date), !is.na(host_species)) %>%
  mutate(year = year(start_date)) %>%
  group_by(host_species) %>%
  summarise(first_year = min(year), .groups = "drop")

first_year_detected <- pathogen_data %>%
  left_join(host_data %>% select(host_record_id, start_date), by = "host_record_id") %>%
  filter(!is.na(start_date), !is.na(pathogen_species_cleaned)) %>%
  mutate(year = year(start_date)) %>%
  group_by(pathogen_species_cleaned) %>%
  summarise(first_year = min(year), .groups = "drop")

discovery_curve_data <- tibble(year = 1960:2024) %>%
  rowwise() %>%
  mutate(
    cumulative_hosts = sum(first_year_sampled$first_year <= year),
    cumulative_pathogens = sum(first_year_detected$first_year <= year)
  ) %>%
  ungroup() %>%
  pivot_longer(
    cols = c(cumulative_hosts, cumulative_pathogens),
    names_to = "taxon_group",
    values_to = "cumulative_species"
  )

discovery_curve_plot <- ggplot(discovery_curve_data, aes(x = year, y = cumulative_species, color = taxon_group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1", name = "Taxon Group", labels = c("Hosts", "Pathogens")) +
  labs(
    title = "Species Discovery Curves Over Time",
    subtitle = "Cumulative number of unique host and pathogen species recorded in the dataset",
    x = "Year",
    y = "Cumulative Number of Unique Species"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(discovery_curve_plot)

