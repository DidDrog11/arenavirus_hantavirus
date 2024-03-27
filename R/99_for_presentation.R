source(here::here("R", "00_load_data.R"))
combined_data <- read_rds(here("data", "clean_data", "2024-03-22_data.rds"))

# Additional packages
pkgs <- c(
  "cowplot",
  "terra",
  "geodata",
  "rnaturalearth",
  "sf"
)

pacman::p_load(pkgs, character.only = T)

# Read in clover data
clover_data <- read_rds(here("data", "clover_data.rds")) %>%
  mutate(pathogen_family = case_when(str_detect(Pathogen, "hantavir") ~ "Hantaviridae",
                                     str_detect(Pathogen, "arenavir") ~ "Arenaviridae",
                                     str_detect(Pathogen, "Juquitiba") ~ "Hantaviridae",
                                     str_detect(Pathogen, "Pampa") ~ "Arenaviridae"))

# Read in IUCN data
iucn_match <- read_rds(here("data", "iucn_match.rds"))
iucn_ranges <- read_rds(here("data", "iucn_ranges.rds"))

# Load geographic data
world_vect <- world(path = here("data"), resolution = 3)

# Country harmonisation
countrycodes <- country_codes()

# Map potential distribution of Hantaviridae
hanta_hosts <- clover_data %>%
  filter(pathogen_family == "Hantaviridae") %>%
  ungroup() %>%
  distinct(`Host species`)

hanta_host_ranges <- iucn_ranges %>%
  filter(SCI_NAME %in% hanta_hosts$`Host species`) %>%
  group_by(SCI_NAME) %>%
  group_split()

countries_containing_hanta_hosts <- lapply(hanta_host_ranges, function(x) {
  
  if(nrow(x) == 1) {
    
    range <- tibble(country = world_vect[st_overlaps(x, st_as_sf(world_vect))[[1]], ]$NAME_0) %>%
      mutate(host = x$SCI_NAME)
    
  }
  
  if(nrow(x) > 1) {
    
    range <- tibble(country = NA,
                    host = NA)
    
    for(i in 1:nrow(x)) {
      
      range <- bind_rows(range, 
                         tibble(country = world_vect[st_overlaps(x[i, ], st_as_sf(world_vect))[[1]], ]$NAME_0) %>%
                           mutate(host = x$SCI_NAME[i]))
      
    }
    
  }
  
  return(range)

}) %>%
  bind_rows() %>%
  drop_na(host)
  
n_hanta_hosts_country <- countries_containing_hanta_hosts %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  left_join(st_as_sf(world_vect), by = c("country" = "NAME_0")) %>%
  st_as_sf(crs = crs(st_as_sf(world_vect)))

# Map potential Arenaviridae

arena_hosts <- clover_data %>%
  filter(pathogen_family == "Arenaviridae") %>%
  ungroup() %>%
  distinct(`Host species`)

arena_host_ranges <- iucn_ranges %>%
  filter(SCI_NAME %in% arena_hosts$`Host species`) %>%
  group_by(SCI_NAME) %>%
  group_split()

countries_containing_arena_hosts <- lapply(arena_host_ranges, function(x) {
  
  if(nrow(x) == 1) {
    
    range <- tibble(country = world_vect[st_overlaps(x, st_as_sf(world_vect))[[1]], ]$NAME_0) %>%
      mutate(host = x$SCI_NAME)
    
  }
  
  if(nrow(x) > 1) {
    
    range <- tibble(country = NA,
                    host = NA)
    
    for(i in 1:nrow(x)) {
      
      range <- bind_rows(range, 
                         tibble(country = world_vect[st_overlaps(x[i, ], st_as_sf(world_vect))[[1]], ]$NAME_0) %>%
                           mutate(host = x$SCI_NAME[i]))
      
    }
    
  }
  
  return(range)
  
}) %>%
  bind_rows() %>%
  drop_na(host)


n_arena_hosts_country <- countries_containing_arena_hosts %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  left_join(st_as_sf(world_vect), by = c("country" = "NAME_0")) %>%
  st_as_sf(crs = crs(st_as_sf(world_vect)))


# Plots of hosts ----------------------------------------------------------

hanta_country_plot <- ggplot() +
  geom_sf(data = st_as_sf(world_vect)) +
  geom_sf(data = n_hanta_hosts_country, aes(fill = n)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Number of Hantavirus rodent host species within a country",
       fill = "N species")

save_plot(hanta_country_plot, filename = here("output", "hanta_country_hosts.png"), base_width = 12)

hanta_host_range <- ggplot() +
  geom_sf(data = st_as_sf(world_vect)) +
  geom_sf(data = iucn_ranges %>%
            filter(SCI_NAME %in% hanta_hosts$`Host species`),
          fill = "red",
          alpha = 0.2) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "IUCN ranges of Hantavirus rodent host species")

save_plot(hanta_host_range, filename = here("output", "hanta_hosts_range.png"), base_width = 12)

arena_country_plot <- ggplot() +
  geom_sf(data = st_as_sf(world_vect)) +
  geom_sf(data = n_arena_hosts_country, aes(fill = n)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Number of Arenavirus rodent host species within a country",
       fill = "N species")

save_plot(arena_country_plot, filename = here("output", "arena_country_hosts.png"), base_width = 12)

arena_host_range <- ggplot() +
  geom_sf(data = st_as_sf(world_vect)) +
  geom_sf(data = iucn_ranges %>%
  filter(SCI_NAME %in% arena_hosts$`Host species`),
  fill = "red",
  alpha = 0.2) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "IUCN ranges of Arenavirus rodent host species")

save_plot(arena_host_range, filename = here("output", "arena_hosts_range.png"), base_width = 12)


# Explore missing metadata GenBank ----------------------------------------

arena_genbank <- read_csv(here("data", "Arenas.csv"))

hanta_genbank <- read_csv(here("data", "Hantas.csv"))

incompleteness <- bind_rows(arena_genbank %>%
            mutate(family = "Arenaviridae (3,115)"), 
          hanta_genbank %>%
            mutate(family = "Hantaviridae (10,911)")) %>%
  mutate(country = if_else(is.na(Country), "Incomplete", "Complete"),
         host =  if_else(is.na(Host), "Incomplete", "Complete")) %>%
  group_by(family, country, host) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_col(aes(x = country, y = n, fill = host)) +
  facet_wrap(~ family) +
  theme_minimal() +
  labs(y = "N sequences",
       x = "Country data",
       fill = "Host data")

save_plot(plot = incompleteness,
            filename = here("output", "metadata_incompleteness.png"),
            base_width = 6,
            base_height = 10)

# Plot countries where studies were conducted -----------------------------

study_countries <- combined_data$citations %>%
  filter(!str_detect(decision, "Exclude")) %>%
  filter(!is.na(full_text_id)) %>%
  drop_na(Country) %>%
  select(Country) %>%
  separate_rows(Country, sep = ", ") %>%
  group_by(Country) %>%
  summarise(n = n())
  
study_countries$iso3 <- countrycode::countrycode(sourcevar = study_countries$Country, origin_regex = TRUE, origin = "country.name.en", destination = "iso3c")

study_map <- ggplot() + 
  geom_sf(data = st_as_sf(world_vect)) +
  geom_sf(data = st_as_sf(world_vect) %>%
  left_join(study_countries, by = c("GID_0" = "iso3")) %>%
  drop_na(n), aes(fill = n)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Countries with included studies (N = 810)",
       fill = "Number of publications")

save_plot(study_map, filename = here("output", "study_map.png"), base_width = 8, base_height = 8)


# H-P associations in data ------------------------------------------------

unique_host <- combined_data$pathogen %>%
  filter(!str_detect(associatedTaxa, "sp.|spp.")) %>%
  distinct(associatedTaxa) %>%
  drop_na() %>%
  arrange(associatedTaxa)

unique_path <- combined_data$pathogen %>%
  drop_na(scientificName) %>%
  distinct(scientificName)

library(taxize)
id <- get_gbifid(unique_host$associatedTaxa)
gbifid <- id2name(id)
gbif_df <- map_dfr(gbifid, as_tibble)
clean_host <- bind_cols(unique_host, gbif_df)

path_id <- get_uid(unique_path$scientificName)
taxize_options(ncbi_sleep = 0.6)
ncbiid <- id2name(path_id, db = "ncbi")
ncbi_df <- map_dfr(ncbiid, as_tibble)
clean_path <- bind_cols(unique_path, ncbi_df)

hp_df <- combined_data$pathogen %>%
  drop_na(scientificName) %>%
  filter(!str_detect(associatedTaxa, "sp.|spp.")) %>%
  filter(!str_detect(scientificName, "/")) %>%
  select(associatedTaxa, scientificName, family, occurrenceRemarks, organismQuantity, decimalLatitude, decimalLongitude) %>%
  left_join(clean_host %>%
              select(-value, -status) %>%
              rename("host_id" = id,
                     "host_clean" = name,
                     "host_rank" = rank), by = "associatedTaxa") %>%
  left_join(clean_path %>%
              select(-value, -status) %>%
              rename("path_id" = id,
                     "path_clean" = name,
                     "path_rank" = rank), by = "scientificName") %>%
  mutate(path_combined = coalesce(path_clean, scientificName)) %>%
  filter(host_rank == "species")

hp_summarised <- hp_df %>%
  group_by(family, host_clean, path_combined) %>%
  summarise(n_tested = sum(occurrenceRemarks),
            n_positive = sum(organismQuantity)) %>%
  filter(n_tested != 0) %>%
  group_by(host_clean, path_combined) %>%
  mutate(prevalence = round(n_positive/(n_tested), 3))

hp_list <- hp_summarised %>%
  group_by(family) %>%
  group_split()

library(tidyterra)

hanta_hp <- ggplot() +
  geom_point(data = hp_list[[1]], aes(y = path_combined, x = prevalence, colour = host_clean, size = n_tested), position = position_jitter(width = 0)) +
  guides(colour = "none") +
  theme_bw() +
  labs(title = paste(unique(hp_list[[1]]$family)),
       x = "Prevalence",
       y = element_blank(),
       size = "Number tested") +
  xlim(0, 1)

mammarena_hp <- ggplot() +
  geom_point(data = hp_list[[2]], aes(y = path_combined, x = prevalence, colour = host_clean, size = n_tested), position = position_jitter(width = 0)) +
  guides(colour = "none") +
  theme_bw() +
  labs(title = paste(unique(hp_list[[2]]$family)),
       x = "Prevalence",
       y = element_blank(),
       size = "Number tested") +
  xlim(0, 1)

ggsave(plot = hanta_hp, filename = here("output", "hantavirus_hp.png"), width = 8)
ggsave(plot = mammarena_hp, filename = here("output", "arenavirus_hp.png"), width = 8)

examples <- c("Orthohantavirus sinnombreense", "Mammarenavirus lassaense")

iucn_all <- vect(here("data", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))

map_hp <- hp_df %>%
  filter(path_combined %in% examples) %>%
  group_by(path_combined) %>%
  group_split()

hp_1_assays <- map_hp[[1]] %>%
  group_by("sci_name" = host_clean) %>%
  summarise(n_tested = sum(occurrenceRemarks),
            n_positive = sum(organismQuantity)) %>%
  arrange(-n_tested)

hp_1_ranges <- iucn_all[iucn_all$sci_name %in% hp_1_assays$sci_name] %>%
  as_tibble(geom = "WKT") %>%
  select(sci_name, geometry) %>%
  drop_na(geometry) %>%
  vect(geom = "geometry", crs = "EPSG:4326")

hp_1_points <- map_hp[[1]] %>%
  distinct("sci_name" = associatedTaxa, occurrenceRemarks, organismQuantity, decimalLatitude, decimalLongitude) %>%
  filter(organismQuantity >= 1) %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

plot_mammarena <- list()

for(i in 1:length(unique(hp_1_points$sci_name))){
  
  species <- sort(unique(hp_1_points$sci_name))[i]
  
  species_points <- hp_1_points %>%
    filter(sci_name == species)
  
  species_range <- hp_1_ranges %>%
    filter(sci_name == species)
  
  extent <- ext(species_range)
  
  plot_mammarena <- ggplot() +
    geom_spatvector(data = world_vect) +
    geom_spatvector(data = species_range, fill = "red", alpha = 0.2) +
    geom_spatvector(data = species_points, aes(colour = organismQuantity, size = occurrenceRemarks), alpha = 0.4) +
    scale_size_area() +
    scale_color_viridis_c() +
    theme_minimal() +
    labs(title = paste(species), 
         colour = "Positive",
         size = "Tested") +
    coord_sf(xlim = extent[1:2], ylim = extent[3:4])
  
  ggsave(plot = plot_mammarena, filename = here("output", "sample_maps", "arenaviridae", paste0(species, ".png")), width = 6)
  
}


hp_2_assays <- map_hp[[2]] %>%
  group_by("sci_name" = host_clean) %>%
  summarise(n_tested = sum(occurrenceRemarks),
            n_positive = sum(organismQuantity)) %>%
  arrange(-n_tested)

hp_2_ranges <- iucn_all[iucn_all$sci_name %in% hp_2_assays$sci_name] %>%
  as_tibble(geom = "WKT") %>%
  select(sci_name, geometry) %>%
  drop_na(geometry) %>%
  vect(geom = "geometry", crs = "EPSG:4326")

hp_2_points <- map_hp[[2]] %>%
  distinct("sci_name" = host_clean, occurrenceRemarks, organismQuantity, decimalLatitude, decimalLongitude) %>%
  filter(organismQuantity >= 1) %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

plot_hantavir <- list()

for(i in 1:length(unique(hp_2_points$sci_name))){
  
  species <- sort(unique(hp_2_points$sci_name))[i]
  
  species_points <- hp_2_points %>%
    filter(sci_name == species)
  
  species_range <- hp_2_ranges %>%
    filter(sci_name == species)
  
  extent <- ext(species_range)
  
  plot_hantavir <- ggplot() +
    geom_spatvector(data = world_vect) +
    geom_spatvector(data = species_range, fill = "red", alpha = 0.2) +
    geom_spatvector(data = species_points, aes(colour = organismQuantity, size = occurrenceRemarks), alpha = 0.4) +
    scale_size_area() +
    scale_color_viridis_c() +
    theme_minimal() +
    labs(title = paste(species), 
         colour = "Positive",
         size = "Tested") +
    coord_sf(xlim = extent[1:2], ylim = extent[3:4])
  
  ggsave(plot = plot_hantavir, filename = here("output", "sample_maps", "hantaviridae", paste0(species, ".png")), width = 6)
  
}
