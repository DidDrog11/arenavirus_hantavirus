library(here)
library(terra)
library(tidyterra)
library(tidyverse)
library(rnaturalearth)
library(patchwork)

db_path <- here("data", "database", "Project_ArHa_database_2025-09-25.rds") # <-- Update this date
arha_db <- read_rds(db_path)

mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

combined <- arha_db$pathogen |> 
  left_join(arha_db$host, by = "host_record_id") |> 
  drop_na(host_species) |> 
  filter(taxonomic_level == "species")

combined_summary <- combined |> 
  group_by(pathogen_species_cleaned, host_species, country, iso3c) |> 
  summarise(number_detected = sum(number_of_hosts, na.rm = TRUE),
            number_tested = sum(number_tested, na.rm = TRUE),
            number_positive = sum(number_positive, na.rm = TRUE))

n_pathogens_in_host <- combined_summary |> 
  filter(number_positive != 0) |>
  group_by(host_species) |> 
  summarise(n_pathogens_in_species = n_distinct(pathogen_species_cleaned))

n_hosts_of_pathogen <- combined_summary |> 
  filter(number_positive != 0) |>
  group_by(pathogen_species_cleaned) |> 
  summarise(n_hosts_of_pathogen = n_distinct(host_species))

summary_information <- combined_summary |>
  left_join(n_pathogens_in_host, by = "host_species") |>
  left_join(n_hosts_of_pathogen, by = "pathogen_species_cleaned") |>
  group_by(country, iso3c) |>
  mutate(n_host_pathogen_records = n()) |>
  ungroup()

world <- ne_countries(scale = 10,  returnclass = "sv") |>
  project(y = mollweide_crs)

map_summary_1 <- merge(world |> 
                         select(adm0_a3),
                       summary_information |>
                         distinct(iso3c, n_host_pathogen_records),
                       by.x = "adm0_a3",
                       by.y = "iso3c",
                       all.x = TRUE)

map_summary_1 |>
  ggplot() +
  geom_spatvector(aes(fill = n_host_pathogen_records)) +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill = "Host-Pathogen records",
       caption = "Map shows the number of records (distinct pathogen species within different hosts).") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Specific examples
# China
country_iso = "BEL"

data <- combined |>
  filter(iso3c == country_iso) |>
  filter(coordinate_resolution_processed != "country") |>
  group_by(latitude, longitude, pathogen_species_cleaned, host_species) |>
  summarise(n_tested = sum(number_tested, na.rm = TRUE),
            n_positive = sum(number_positive, na.rm = TRUE),
            prop_positive = n_positive/n_tested)

coords <- data |>
  vect(geom = c("longitude", "latitude"), crs = "EPSG:4326") |>
  project(y = mollweide_crs) |>
  crds()

data_jittered <- data |>
  ungroup() |>
  select(-longitude, -latitude) |>
  cbind(coords) |>
  group_by(x, y) |>
  mutate(n_at_loc = n(),
         point_index = row_number())|>
  ungroup() |>
  mutate(jitter_radius = 10000, 
         angle = (point_index - 1) * (2 * pi / n_at_loc),
         final_x = if_else(n_at_loc > 1, x + cos(angle) * jitter_radius, x),
         final_y = if_else(n_at_loc > 1, y + sin(angle) * jitter_radius, y)) |>
  vect(geom = c("final_x", "final_y"), crs = mollweide_crs)

# Pathogen

pathogen_map <- ggplot() +
  geom_spatvector(data = world |>
                    filter(adm0_a3 == country_iso),
                  lwd = 1) +
  geom_spatvector(data = data_jittered,
                  aes(colour = pathogen_species_cleaned, size = prop_positive)) +
  scale_size_continuous(name = "Proportion Positive", range = c(0.5, 4)) +
  scale_color_viridis_d(name = "Pathogen Species") +
  labs(title = paste0("Sampling by pathogen in: ", country_iso)) +
  theme_minimal() +
  theme(legend.position = "bottom")

pathogen_barchart_tested <- data_jittered |>
  as_tibble() |>
  group_by(pathogen_species_cleaned, host_species) |>
  summarise(n_tested = sum(n_tested, na.rm = TRUE), .groups = "drop") |>
  filter(n_tested != 0) |>
  ggplot(aes(x = n_tested, y = pathogen_species_cleaned, fill = host_species)) +
  geom_col(colour = "black") +
  scale_fill_viridis_d(name = "Host Species") +
  labs(x = "Total Tested", y = NULL, title = "N Tested") +
  theme_minimal() +
  theme(legend.position = "none")

pathogen_barchart_prevalence <- data_jittered |>
  as_tibble() |>
  group_by(pathogen_species_cleaned, host_species) |>
  summarise(n_tested = sum(n_tested, na.rm = TRUE),
            n_positive = sum(n_positive, na.rm = TRUE),
            prop_positive = n_positive/n_tested,
            .groups = "drop") |>
  filter(n_tested != 0) |>
  ggplot(aes(x = prop_positive, y = pathogen_species_cleaned)) + 
  geom_point(aes(colour = host_species, fill = host_species), 
             shape = 21, size = 3,
             position = position_jitter(width = 0, height = 0.2)) +
  scale_color_viridis_d(name = "Host Species") +
  scale_fill_viridis_d(name = "Host Species") +
  labs(x = "Pathogen prevalence", y = NULL, title = "Prevalence") +
  theme_minimal()

pathogen_dashboard <- pathogen_map + pathogen_barchart_tested + pathogen_barchart_prevalence +
  plot_layout(widths = c(2, 1, 1),
              guides = "collect") +
  plot_annotation(title = paste0("Pathogen Surveillance Summary for: ", country_iso)) & 
  theme(legend.position = "bottom")

# Host

ggplot() +
  geom_spatvector(data = world |>
                    filter(adm0_a3 == country),
                  lwd = 1) +
  geom_spatvector(data = data_jittered,
                  aes(colour = host_species, size = prop_positive)) +
  scale_size_continuous(name = "Proportion Positive", range = c(0.5, 4)) +
  scale_color_viridis_d(name = "Host Species") +
  labs(title = paste0("Sampling by host in: ", country_iso)) +
  theme_minimal() +
  theme(legend.position = "bottom")
  
