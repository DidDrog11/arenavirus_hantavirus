# ---
# 03_geographic_bias_analysis.R
#
# Purpose: To quantify geographic and temporal biases in surveillance effort.
# This script generates maps of sampling locations and compares them to host ranges.
# ---

# Setup -------------------------------------------------------------------
# Load necessary packages
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(tidyverse)
library(stringr)
library(geodata)
library(ggrepel)
library(rnaturalearth)
library(lubridate)
library(WDI)
library(taxize)
library(akima)
library(scales)
library(ggrepel)
library(metR)    

# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-09-25.rds") # <-- Update this date
host_traits <- read_rds(here("data", "external", "harmonised_species.rds"))
arha_db <- read_rds(db_path)
host_data <- arha_db$host
pathogen_data <- arha_db$pathogen

# Load harmonized traits data (needed for joining by gbif_id)
host_traits <- read_rds(here("data", "external", "harmonised_species.rds"))

mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

lookup_path <- here("data", "external", "adm2_species_lookup.rds")


# Create spatial lookup genus ---------------------------------------------
if (!file.exists(lookup_path)) {
  
  iucn_redlist <- vect(here("data", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
  
  iucn_subset <- iucn_redlist %>%
    filter(sci_name %in% host_traits$scientific_name) %>%
    filter(presence %in% c(1, 2, 3, 4)) %>%
    project(mollweide_crs) %>%
    select(scientific_name = sci_name)
  
  iucn_raw <- iucn_subset %>%
    arrange(scientific_name) %>%
    left_join(host_traits %>%
                select(scientific_name, gbif_genus), by = "scientific_name") %>%
    makeValid() %>%
    aggregate(by = "gbif_genus") %>%
    simplifyGeom(tolerance = 5000)
  
  gadm_adm2 <- vect(here("data", "gadm", "gadm_adm2_combined.shp"))
  
  gadm_adm2_proj <- project(gadm_adm2, y = mollweide_crs) %>%
    simplifyGeom(tolerance = 5000)
  
  gadm_adm2_proj$GID_2_AREA <- expanse(gadm_adm2_proj, unit = "km", transform = TRUE)
  
  countries_to_process <- unique(gadm_adm2_proj$GID_0)
  
  inhabited_adm2_lookup <- purrr::map_dfr(countries_to_process, function(country_iso) {
    
    message(paste("... processing", country_iso))
    
    # Filter for the current country's ADM2 polygons
    country_adm2 <- gadm_adm2_proj %>% filter(GID_0 == country_iso)
    
    message(paste("... dissolving", country_iso))
    country_border <- aggregate(country_adm2)
    message(paste("... dissolved", country_iso))
    
    iucn_cropped <- crop(iucn_raw, country_border)
    message(paste("... iucn cropped", country_iso))
    
    # If there are no relevant IUCN ranges, skip to the next country
    if (nrow(iucn_cropped) == 0) {
      return(NULL)
    }
    
    genus_per_iucn_poly <- as.data.frame(iucn_cropped) %>%
      pull(gbif_genus)
    
    intersection_matrix <- relate(country_adm2, iucn_cropped, relation = "intersects")
    
    genus_presence_tibble <- as_tibble(intersection_matrix, .name_repair = "minimal")
    colnames(genus_presence_tibble) <- genus_per_iucn_poly
    
    result <- as_tibble(country_adm2) %>%
      bind_cols(genus_presence_tibble)
    
    return(result)
    
  })
  
  
  lookup_table_wide <- inhabited_adm2_lookup
  
  lookup_table_long <- lookup_table_wide %>%
    pivot_longer(
      cols = any_of(na.omit(unique(host_traits$gbif_genus))), 
      names_to = "gbif_genus",
      values_to = "present"
    ) %>%
    filter(present == TRUE) %>%
    select(-present)
  
  lookup_combined <- list(
    lookup_table_wide = lookup_table_wide,
    lookup_table_long = lookup_table_long
  )
  
  # c) Save the lookup table to cache
  write_rds(lookup_combined, lookup_path)
  writeVector(gadm_adm2_proj, here("data", "gadm", "gadm_adm2_simplified.shp"), overwrite = TRUE)
  writeVector(iucn_raw, here("data", "external", "iucn_genus_subset.shp"), overwrite = TRUE)
  
} else {
  
  lookup_combined <- read_rds(lookup_path)
  gadm_adm2_proj <- vect(here("data", "gadm", "gadm_adm2_simplified.shp"))
  iucn_raw <- vect(here("data", "external", "iucn_genus_subset.shp"))
}

ne_world <- rnaturalearth::ne_countries(scale = 50) %>%
  st_transform(crs = mollweide_crs) %>%
  st_make_valid()

# --- 3. Geographic Bias Analysis ---

# 3.1. Quantitative Assessment of Sampling Gaps in Host Ranges

sampling_points <- host_data %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  group_by(longitude, latitude) %>%
  mutate(sampling_point_id = cur_group_id()) %>%
  ungroup() %>%
  vect(geom = c("longitude", "latitude"), crs = "EPSG:4326") %>%
  project(y = mollweide_crs)

# Get a unique list of all ADM2 units that have been sampled in our database
unique_sampling_points <- sampling_points %>%
  as_tibble() %>%
  filter(!is.na(gadm_adm2)) %>%
  group_by(gadm_adm2, host_genus) %>%
  summarise(number_detected = sum(number_of_hosts, na.rm = TRUE),
            .groups = "drop")

# --- 3.2.a. Analysis by Genus ---
all_genera <- host_data %>%
  filter(host_order %in% c("Rodentia", "Soricomorpha")) %>%
  drop_na(host_genus) %>%
  drop_na(gadm_adm2) %>%
  group_by(host_genus) %>%
  summarise(n = sum(number_of_hosts, na.rm = TRUE)) %>%
  arrange(-n) %>%
  pull(host_genus)

genera_to_process <- all_genera %>%
  purrr::set_names()

cache_dir <- here("output", "analysis_1", "coverage_cache")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
}

genus_range_coverage_summary <- purrr::map_dfr(genera_to_process, function(genus_name) {
  
  # Define the specific output file for this genus
  output_file_path <- file.path(cache_dir, paste0(genus_name, "_coverage.rds"))
  
  # Check if the result for this genus already exists
  if (file.exists(output_file_path)) {
    
    message(paste("... loading", genus_name, "from cache."))
    # If it exists, load it and return immediately
    cached_result <- read_rds(output_file_path)
    return(cached_result$summary)
    
  } else {
    
  message(paste("... processing genus:", genus_name, Sys.time()))
  
  genus_range <- iucn_raw %>%
    filter(gbif_genus == genus_name)
  
  message(paste("Associating sampling locations with adm2 ", Sys.time()))
  
  genus_sampling_locations <- as_tibble(unique_sampling_points) %>%
    filter(host_genus == genus_name) %>%
    pull(gadm_adm2)
  
  genus_sampling_locations_vect <- gadm_adm2_proj %>%
    filter(GID_2 %in% genus_sampling_locations) %>%
    left_join(as_tibble(unique_sampling_points) %>%
                filter(host_genus == genus_name), by = c("GID_2" = "gadm_adm2"))
  
  inhabited_adm2 <- lookup_table_long %>%
    filter(gbif_genus == genus_name)
  
  inhabited_adm2_vect <- gadm_adm2_proj %>%
    filter(GID_2 %in% inhabited_adm2$GID_2)
  
  message(paste("Sampled area defined ", Sys.time()))
  
  message(paste("Calculating areas ", Sys.time()))
  
  total_inhabited_area_km2 <- sum(inhabited_adm2_vect$GID_2_AREA, na.rm = TRUE)
  
  total_sampled_area_km2 <- sum(genus_sampling_locations_vect$GID_2_AREA, na.rm = TRUE)
  
  message(paste("Areas calculated ", Sys.time()))
  
  sampling_coverage <- tibble(
    genus = genus_name,
    range_area_km2 = as.numeric(total_inhabited_area_km2),
    sampled_area_km2 = as.numeric(total_sampled_area_km2),
    proportion_of_area_sampled = as.numeric(total_sampled_area_km2 / total_inhabited_area_km2)
  )
  
  message(paste("Producing map ", Sys.time()))
  
  bboxes_vec <- c(ext(genus_range), ext(genus_sampling_locations_vect))
  combined_extent_object <- c(
    xmin = min(bboxes_vec[[1]]$xmin, bboxes_vec[[2]]$xmin),
    ymin = min(bboxes_vec[[1]]$ymin, bboxes_vec[[2]]$ymin),
    xmax = max(bboxes_vec[[1]]$xmax, bboxes_vec[[2]]$xmax),
    ymax = max(bboxes_vec[[1]]$ymax, bboxes_vec[[2]]$ymax))
  
  sampling_map <- ggplot() +
    geom_sf(data = ne_world, lwd = 0.6) +
    geom_sf(data = genus_range, fill = "firebrick", alpha = 0.3) +
    geom_sf(data = genus_sampling_locations_vect, aes(fill = number_detected), alpha = 1, colour = "black", lwd = 0.4) +
    scale_fill_viridis_c(trans = scales::log10_trans()) +
    coord_sf(xlim = c(combined_extent_object[1], combined_extent_object[3]),
             ylim = c(combined_extent_object[2], combined_extent_object[4]), expand = TRUE) +
    labs(fill = "Number of individuals detected",
         title = paste0("Sampling of ", genus_name, " species"),
         caption = paste0("Red shaded area represents the IUCN range of species within the genus ", genus_name, "\n Administrative level 2 areas shaded grey are non-detections")) +
    theme_minimal()
  
  message(paste("Map produced ", Sys.time()))
  
  result_list <- list(
    summary = sampling_coverage,
    map = sampling_map
  )
  
  saveRDS(result_list, output_file_path)
  
  message(paste("Results for", genus_name, "saved to cache."))
  
  return(result_list$summary)
  
  }
  
})

plot_data <- genus_range_coverage_summary %>%
  filter(sampled_area_km2 > 0, range_area_km2 > 0)

# --- 1. Interpolate the data onto a regular grid ---
# We perform interpolation on the log10 scale to match the plot axes
plot_data <- genus_range_coverage_summary %>%
  filter(sampled_area_km2 > 0, range_area_km2 > 0) %>%
  # Identify the top 15 most widespread genera to label
  mutate(
    label = if_else(rank(desc(range_area_km2)) <= 15, genus, "")
  )

# --- Interpolate the data (no change here) ---
interpolated_surface <- interp(
  x = log10(plot_data$sampled_area_km2),
  y = log10(plot_data$range_area_km2),
  z = plot_data$proportion_of_area_sampled,
  linear = TRUE, extrap = FALSE, nx = 200, ny = 200
)

interpolated_df <- interp2xyz(interpolated_surface) %>%
  as_tibble() %>%
  mutate(x = 10^x, y = 10^y)

# --- Create the Plot ---
ggplot() +
  geom_raster(data = interpolated_df, aes(x = x, y = y, fill = z), alpha = 1) +
  scale_fill_viridis_c(
    name = "Interpolated\nSampling Proportion", 
    labels = label_percent(),
    na.value = "transparent"
  ) +
  geom_contour(
    data = interpolated_df, aes(x = x, y = y, z = z), 
    breaks = c(0.001, 0.01, 0.1, 1.0), 
    color = "white", linewidth = 0.6
  ) +
  geom_label_contour(
    data = interpolated_df, 
    aes(x = x, y = y, z = z), 
    breaks = c(0.001, 0.01, 0.1, 1.0),
    stroke = 0.2,
    label.placer = label_placer_fraction(frac = 0.5)
  ) +
  geom_point(data = plot_data, aes(x = sampled_area_km2, y = range_area_km2), 
             shape = 21, size = 2, fill = "white", color = "black", alpha = 0.7) +
  geom_label_repel(
    data = plot_data, 
    aes(x = sampled_area_km2, y = range_area_km2, label = label),
    size = 3, 
    max.overlaps = 10,
    min.segment.length = 0
  ) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  labs(
    title = "Sampling Effort vs. Total Inhabited Area by Genus",
    subtitle = "Background color shows interpolated sampling proportion",
    x = "Sampled Area within Range (km², log scale)",
    y = "Total Inhabited Area (km², log scale)"
  ) +
  theme_minimal()

# Load maps seperately
cache_dir <- here("output", "analysis_1", "coverage_cache")
all_files <- list.files(cache_dir, pattern = "\\.rds$", full.names = TRUE)

all_maps <- purrr::map(all_files, ~read_rds(.x)$map)

genus_names <- str_remove(basename(all_files), "_coverage\\.rds$")
all_maps <- set_names(all_maps, genus_names)

print(all_maps$Apodemus)
print(all_maps$Rattus)
print(all_maps$Mus)
print(all_maps$Peromyscus)
print(all_maps$Mastomys)
print(all_maps$Myodes)

# --- 3.2.b. Analysis by Top 5 Sampled Species ---
# --- Create and save a SPECIES-level lookup table ---
species_lookup_path <- here("data", "external", "adm2_species_lookup_species_level.rds")

# First, identify the top 10 most sampled species from the database
top_species <- host_data %>%
  drop_na(host_species) %>%
  group_by(host_species) %>%
  summarise(total_sampled = sum(number_of_hosts, na.rm = TRUE)) %>%
  slice_max(order_by = total_sampled, n = 10) %>%
  pull(host_species)

top_species[2] = "Clethrionomys glareolus"

if (!file.exists(species_lookup_path)) {
  
  iucn_redlist <- vect(here("data", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
  
  iucn_subset <- iucn_redlist %>%
    filter(sci_name %in% host_traits$scientific_name) %>%
    filter(presence %in% c(1, 2, 3, 4)) %>%
    project(mollweide_crs) %>%
    select(scientific_name = sci_name)
  
  iucn_species_subset <- iucn_subset %>%
    arrange(scientific_name) %>%
    left_join(host_traits %>%
                select(scientific_name, gbif_genus), by = "scientific_name") %>%
    filter(scientific_name %in% top_species) %>%
    aggregate(by = "scientific_name") %>%
    makeValid() %>%
    simplifyGeom(tolerance = 5000)
  
  simplified_adm2 <- gadm_adm2_proj %>%
    makeValid() %>%
    simplifyGeom(tolerance = 5000)
  
  countries_to_process <- unique(simplified_adm2$GID_0)
  
  cache_dir <- here("data", "external", "gadm_iucn_cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  inhabited_adm2_lookup <- purrr::map_dfr(countries_to_process, function(country_iso) {
    
    output_file_path <- file.path(cache_dir, paste0(country_iso, "_lookup.rds"))
    
    if (file.exists(output_file_path)) {
      
      message(paste("... loading", country_iso, "from cache."))
      return(read_rds(output_file_path))
      
    } else {
      
      if (is.na(country_iso)) {
        return(NULL)
      }
      
      message(paste(Sys.time(), ": processing", country_iso))
      
      country_adm2 <- simplified_adm2[simplified_adm2$GID_0 == country_iso, ]
      
      # Get the dissolved border of this country to quickly pre-filter species
      country_border <- aggregate(country_adm2)
      
      # Get IUCN ranges for ONLY species that overlap with this country's border
      intersects_country <- relate(iucn_species_subset, country_border, relation = "intersects")
      
      country_iucn_ranges <- iucn_species_subset[intersects_country, ]
      
      # If there are no relevant species ranges for this country, skip to the next
      if (nrow(country_iucn_ranges) == 0) {
        return(NULL)
      }
      
      # Perform the intersection on these two small, co-located datasets
      relation_matrix <- relate(country_adm2, country_iucn_ranges, relation = "intersects")
      
      intersection_indices <- as.data.frame(which(relation_matrix, arr.ind = TRUE))
      
      if (nrow(intersection_indices) == 0) {
        return(NULL)
      }
      
      # Build the final lookup table for this country by linking the indices back to the data
      result_tibble <- tibble(
        GID_2 = country_adm2$GID_2[intersection_indices$row],
        scientific_name = country_iucn_ranges$scientific_name[intersection_indices$col]
      )
      
      write_rds(result_tibble, output_file_path)
      
      message(paste(Sys.time(), ": completed", country_iso))
      
      return(result_tibble)
    }
  })
  
  # Process and save the long-format table
  lookup_table_long_species <- inhabited_adm2_lookup %>%
    filter(!is.na(scientific_name)) %>%
    select(GID_2, scientific_name) %>%
    distinct()
  
  write_rds(lookup_table_long_species, species_lookup_path)
  
  writeVector(iucn_species_subset, here("data", "external", "iucn_species_subset.shp"))
  
} else {
  
  lookup_table_long_species <- read_rds(species_lookup_path)
  iucn_species_subset <- vect(here("data", "external", "iucn_species_subset.shp"))
  
}

message(paste("\nIdentified top 9 sampled species:", paste(unique(lookup_table_long_species$scientific_name), collapse = ", ")))

# If this includes Clethrionomys glareolus it will need to be changed to Myodes glareolus to match with current taxonomy
top_species[2] = "Myodes glareolus"

species_sampling_points <- sampling_points %>%
  as_tibble() %>%
  filter(!is.na(gadm_adm2)) %>%
  filter(host_species %in% top_species) %>%
  group_by(gadm_adm2, host_species) %>%
  summarise(number_detected = sum(number_of_hosts, na.rm = TRUE),
            .groups = "drop")

species_range_coverage_summary <- purrr::map_dfr(top_species, function(species_name) {
  
  message(paste("... processing species:", species_name))
  
  # Get the total range for the specific species
  if(species_name == "Myodes glareolus") {
    
    iucn_species_range <- iucn_species_subset %>%
      filter(scientific_name == "Clethrionomys glareolus") %>%
      project(y = crs(simplified_adm2))
    
  } else {
   
    iucn_species_range <- iucn_species_subset %>%
      filter(scientific_name == species_name) %>%
      project(y = crs(simplified_adm2)) 
    
  }
  
  species_range <- lookup_table_long_species %>%
    filter(scientific_name == species_name)
  
  species_sampling <- species_sampling_points %>%
    filter(host_species == species_name)
  
  inhabited_adm2 <- simplified_adm2 %>%
    filter(GID_2 %in% species_range$GID_2)
  
  sampled_inhabited_adm2 <- simplified_adm2 %>%
    filter(GID_2 %in% species_sampling$gadm_adm2) %>%
    left_join(species_sampling, by = c("GID_2" = "gadm_adm2"))
  
  total_inhabited_area_km2 <- sum(expanse(iucn_species_range, unit = "km", transform = TRUE))
  total_sampled_area_km2 <- sum(expanse(sampled_inhabited_adm2, unit = "km", transform = TRUE))
  
  positive_sampling_points <- species_sampling %>%
    filter(number_detected > 0)
  
  extralimital_points <- simplified_adm2 %>%
    filter(GID_2 %in% positive_sampling_points$gadm_adm2) %>%
    left_join(species_sampling, by = c("GID_2" = "gadm_adm2")) %>%
    erase(iucn_species_range)
  
  species_sampling_summary <- tibble(
    species = species_name,
    number_detected = sum(species_sampling$number_detected, na.rm = TRUE),
    inhabited_area = total_inhabited_area_km2,
    sampled_area = total_sampled_area_km2,
    proportion_of_area_sampled = as.numeric(total_sampled_area_km2 / total_inhabited_area_km2),
    extralimital_n = sum(extralimital_points$number_detected, na.rm = TRUE),
    extralimital_area = sum(expanse(extralimital_points, unit = "km", transform = TRUE)),
    extralimital_proportion_of_total_area = extralimital_area / (sampled_area + extralimital_area),
    extralimital_proportion_of_detections = extralimital_n / (number_detected + extralimital_n))
  
  # a) Get the areas that were sampled but had zero detections of this species
  sampled_zeros <- sampled_inhabited_adm2 %>% 
    filter(is.na(number_detected) | number_detected == 0)
  
  # b) Get the areas that were sampled and had positive detections
  sampled_positives <- sampled_inhabited_adm2 %>% 
    filter(number_detected > 0)
  
  # c) Get the areas within the range that were never sampled at all
  unsampled_inhabited <- inhabited_adm2 %>%
    filter(!GID_2 %in% sampled_inhabited_adm2$GID_2)
  
  # d) Get the areas within the range that do not fall within the adm2 analysis
  iucn_species_range_outside_adm2 <- erase(iucn_species_range, inhabited_adm2)
  
  # Produce map of range and sampling
  bboxes_vec <- ext(inhabited_adm2) + ext(sampled_inhabited_adm2) + ext(iucn_species_range_outside_adm2)

  sampling_map <- ggplot() +
    geom_sf(data = iucn_species_range_outside_adm2, fill = "firebrick", alpha = 0.2) +
    geom_sf(data = crop(vect(ne_world), bboxes_vec), fill = "transparent", lwd = 0.3) +
    geom_sf(data = unsampled_inhabited, fill = "firebrick", alpha = 0.2, lwd = 0.2) +
    geom_sf(data = sampled_zeros, fill = "black", colour = "transparent", alpha = 1, lwd = 0.2) +
    geom_sf(data = sampled_positives, aes(fill = number_detected), colour = "transparent", alpha = 1, lwd = 0.2) +
    scale_fill_viridis_c(trans = scales::log10_trans(),
                         name = "Individuals Detected\n(log scale)") +
    guides(colour = "none") + 
    labs(title = paste0("Sampling of ", species_name),
         caption = paste0("Red shaded area represents the administrative level 2 regions within the IUCN range of ", species_name, "\n Administrative regions of the country outside of the IUCN range are not shown.\n Black administrative regions represent sampled but with no detections.")) +
    theme_minimal()
  
  cache_dir <- here("output", "analysis_1", "species_coverage_maps")
  ggsave(plot = sampling_map, filename = here(cache_dir, paste0(species_name, ".png")), width = 12, height = 8)
  
  return(species_sampling_summary)
})
