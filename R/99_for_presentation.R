source(here::here("R", "00_load_data.R"))
combined_data <- read_rds(here("data", "clean_data", "2025-06-04_data.rds"))

# Additional packages
pkgs <- c(
  "cowplot",
  "terra",
  "geodata",
  "rnaturalearth",
  "sf"
)

pacman::p_load(pkgs, character.only = T)


# Dataset description -----------------------------------------------------
n_included <- combined_data$citations %>%
  distinct(study_id) %>%
  drop_na() %>%
  nrow()

n_rodent_records <- combined_data$host %>%
  distinct(rodent_record_id) %>%
  drop_na() %>%
  nrow()

n_rodents_sampled <- sum(combined_data$host$individualCount, na.rm = TRUE)

distinct_locations <- combined_data$host %>%
  distinct(decimalLatitude, decimalLongitude, eventDate) %>%
  drop_na() %>%
  nrow()

n_pathogen_records <- combined_data$pathogen %>%
  distinct(pathogen_record_id) %>%
  drop_na() %>%
  nrow()

n_assays <- sum(combined_data$pathogen$n_assayed, na.rm = TRUE)
n_positive <- sum(combined_data$pathogen$n_positive, na.rm = TRUE)
n_negative <- sum(combined_data$pathogen$n_negative, na.rm = TRUE)

n_sequences <- combined_data$sequences %>%
  filter(!is.na(associated_pathogen_record_id) | !is.na(associated_rodent_record_id))

n_path_host_sequences <- combined_data$sequences %>%
  filter(sequenceType == "Pathogen") %>%
  filter(!is.na(associated_pathogen_record_id) & !is.na(associated_rodent_record_id))

n_host_sequences <- combined_data$sequences %>%
  filter(sequenceType == "Host") %>%
  filter(!is.na(associated_pathogen_record_id) | !is.na(associated_rodent_record_id))

combined_data$descriptive %>%
  ungroup() %>%
  summarise(n_rodent_records = median(n_rodent_records, na.rm = TRUE),
            n_rodent_species = median(n_rodent_species, na.rm = TRUE),
            n_individuals = median(n_rodent_individuals, na.rm = TRUE),
            n_pathogen_records = median(n_pathogen_records, na.rm = TRUE),
            n_assays = median(n_assays, na.rm = TRUE),
            n_sequences = median(n_sequence_records[n_sequence_records > 0], na.rm = TRUE))

# Plots of samples --------------------------------------------------------

map_hanta_samples <- combined_data$pathogen %>%
  filter(n_assayed >= 1 &
           str_detect(family, "Hantaviridae")) %>%
  group_by(family, virus_clean, host_species, host_genus, host_family, host_order, decimalLatitude, decimalLongitude) %>%
  summarise(n_assayed = sum(n_assayed, na.rm = TRUE)) %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = project_crs) %>%
  ggplot() +
  geom_spatvector(data = world_shapefile, fill = "transparent") +
  geom_spatvector(aes(colour = log(n_assayed))) +
  scale_colour_viridis_c() +
  theme_minimal() +
  theme(title = element_text(size = 18)) +
  labs(title = "Hantaviridae",
       colour = "N assayed (log10)")

map_arena_samples <- combined_data$pathogen %>%
  filter(n_assayed >= 1 &
           str_detect(family, "Mammarenaviridae|Arenaviridae|Mammarenavirus")) %>%
  group_by(family, virus_clean, host_species, host_genus, host_family, host_order, decimalLatitude, decimalLongitude) %>%
  summarise(n_assayed = sum(n_assayed, na.rm = TRUE)) %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = project_crs) %>%
  ggplot() +
  geom_spatvector(data = world_shapefile, fill = "transparent") +
  geom_spatvector(aes(colour = log(n_assayed))) +
  scale_colour_viridis_c() +
  theme_minimal() +
  theme(title = element_text(size = 18)) +
  labs(title = "Arenaviridae",
       colour = "N assayed (log10)")

save_plot(map_hanta_samples, filename = here("output", "hanta_sample_locations.png"), base_width = 14, base_height = 10)
save_plot(map_arena_samples, filename = here("output", "arena_sample_locations.png"), base_width = 14, base_height = 10)

# Plots of hosts ----------------------------------------------------------

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

study_countries <- combined_data$host %>%
  ungroup() %>%
  mutate(iso3 = coalesce(iso3, GID_0)) %>%
  distinct(study_id, country, iso3) %>%
  mutate(iso3 = case_when(str_detect(country, "South Africa") ~ "ZAF",
                          str_detect(country, "Wales|England|Scotland") ~ "GBR",
                          str_detect(country, "Finland|Finald|Karelia") ~ "FIN",
                          str_detect(country, "US|U.S.|Texas|California|Nevada|Colorado") ~ "USA",
                          str_detect(country, "Ukraine") ~ "UKR",
                          str_detect(country, "Czech") ~ "CZE",
                          str_detect(country, "Phillipines") ~ "PHL",
                          str_detect(country, "Germany|Gemany") ~ "DEU",
                          str_detect(country, "Lithuania") ~ "LTU",
                          str_detect(country, "Nepal") ~ "NPL",
                          str_detect(country, "Venezuela|Venezuala") ~ "VEN",
                          str_detect(country, "Zambia") ~ "ZMB",
                          str_detect(country, "Solvenia|Slovenia") ~ "SVN",
                          str_detect(country, "Argentina") ~ "ARG",
                          str_detect(country, "France") ~ "FRA",
                          str_detect(country, "Korea") ~ "KOR",
                          str_detect(country, "Guinea") ~ "GIN",
                          str_detect(country, "Russia") ~ "RUS",
                          str_detect(country, "Uruaguay") ~ "URY",
                          str_detect(country, "Brazil") ~ "BRA",
                          str_detect(country, "Nambia") ~ "NAM",
                          TRUE ~ iso3
                          )) %>%
  group_by(iso3) %>%
  summarise(n = n())

study_map <- ggplot() + 
  geom_sf(data = st_as_sf(world_shapefile)) +
  geom_sf(data = st_as_sf(world_shapefile) %>%
            left_join(study_countries, by = c("GID_0" = "iso3")) %>%
            drop_na(n), aes(fill = n)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Countries with included studies",
       fill = "Number of publications")

save_plot(study_map, filename = here("output", "study_map.png"), base_width = 8, base_height = 8)


# Harmonising data --------------------------------------------------------
nrow(rodent_names)
n_species <- combined_data$host %>%
  ungroup() %>%
  drop_na(species) %>%
  distinct(species)
n_genera <- combined_data$host %>%
  ungroup() %>%
  drop_na(genus) %>%
  distinct(genus)

nrow(virus_mapping)
n_viruses <- combined_data$pathogen %>%
  ungroup() %>%
  filter(taxonomic_level == "species") %>%
  distinct(virus_clean)


# H-P associations in data ------------------------------------------------
library(ggrepel)
library(ggridges)

unique_host <- combined_data$pathogen %>%
  select(virus_clean, family, taxonomic_level, assay_clean, host_species, host_genus, host_family, host_order, n_assayed, n_positive) %>%
  mutate(family = case_when(str_detect(virus_clean, "Mammarenavirus allpahuayoense|Mammarenavirus lassaense|Mammarenavirus choriomeningitidis") ~ "Arenaviridae",
                            str_detect(virus_clean, "Orthohantavirus andesense|Orthohantavirus seoulense") ~ "Hantaviridae",
                            TRUE ~ family),
         family = case_when(str_detect(family, "Arenaviridae|Mammarenaviridae|Mammarenavirus") ~ "Arenaviridae",
                            TRUE ~ family))

# Limit to named pathogens and species for host
subset_hp <- unique_host %>%
  filter(str_detect(family, "Arena|Hanta")) %>%
  filter(taxonomic_level == "species") %>%
  filter(!str_detect(virus_clean, "Aza")) %>%
  filter(!is.na(host_species)) %>%
  filter(n_assayed >= 1) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(prop_pos = n_positive/n_assayed,
         assay_clean = case_when(str_detect(assay_clean, "PCR|Culture|Sequencing") ~ "Acute",
                                 str_detect(assay_clean, "Serology") ~ "Prior"))


# For Ricardo -------------------------------------------------------------
subset_hp_rr <- unique_host %>%
  filter(taxonomic_level == "species") %>%
  filter(!is.na(host_species)) %>%
  filter(n_assayed >= 1) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(virus_clean = factor(virus_clean),
         host_species = factor(host_species),
         prop_pos = n_positive/n_assayed,
         assay_clean = case_when(str_detect(assay_clean, "PCR|Culture|Sequencing") ~ "PCR, culture or sequencing",
                                 str_detect(assay_clean, "Serology") ~ "Serology"))

subset_hp_biodiv <- combined_data$pathogen %>%
  select(virus_clean, family, taxonomic_level, assay_clean, host_species, host_genus, host_family, host_order, n_assayed, n_positive, decimalLatitude, decimalLongitude) %>%
  filter(taxonomic_level == "species") %>%
  filter(!is.na(host_species)) %>%
  filter(n_assayed >= 1) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(virus_clean = factor(virus_clean),
         host_species = factor(host_species),
         prop_pos = n_positive/n_assayed,
         assay_clean = case_when(str_detect(assay_clean, "PCR|Culture|Sequencing") ~ "PCR, culture or sequencing",
                                 str_detect(assay_clean, "Serology") ~ "Serology"))

write_csv(subset_hp_rr, here("data", "data_outputs", "prevalence_2025-06-04.csv"))
write_csv(subset_hp_biodiv, here("data", "data_outputs", "prevalence_w_coords_2025-06-04.csv"))

# Continue HP -------------------------------------------------------------

subset_hp$virus_clean <- fct_rev(fct_infreq(subset_hp$virus_clean))
subset_hp$host_species <- fct_rev(fct_infreq(subset_hp$host_species))

labelled_hanta <- subset_hp %>%
  filter(family == "Hantaviridae") %>%
  group_by(virus_clean, host_species, assay_clean) %>%
  filter(n_assayed >= 25 &
           n_positive > 1) %>%
  arrange(-prop_pos) %>%
  slice(1)

hp_simple <- subset_hp %>%
  filter(str_detect(family, "Hanta|Arena")) %>%
  filter(n_positive >= 1)

acute_species <- hp_simple %>%
  filter(assay_clean == "Acute") %>%
  group_by(virus_clean, family) %>%
  summarise(acute_species = list(unique(host_species)),
            n_acute_species = n_distinct(host_species),
            n_acute_assayed = sum(n_assayed, na.rm = TRUE),
            .groups = "drop")

prior_species <- hp_simple %>%
  filter(assay_clean == "Prior") %>%
  left_join(acute_species %>% select(virus_clean, acute_species), by = "virus_clean") %>%
  rowwise() %>%
  filter(is.null(acute_species) || !(host_species %in% acute_species)) %>%
  ungroup() %>%
  group_by(virus_clean, family) %>%
  summarise(n_prior_species = n_distinct(host_species),
            n_prior_assayed = sum(n_assayed, na.rm = TRUE),
            .groups = "drop")

hp_summary <- acute_species %>%
  full_join(prior_species, by = c("virus_clean", "family")) %>%
  mutate(
    n_acute_species = replace_na(n_acute_species, 0),
    n_prior_species = replace_na(n_prior_species, 0),
    n_acute_assayed = replace_na(n_acute_assayed, 0),
    n_prior_assayed = replace_na(n_prior_assayed, 0),
    total_species = n_acute_species + n_prior_species,
    total_assayed = n_acute_assayed + n_prior_assayed
  )

hp_plots <- hp_summary %>%
  group_by(family) %>%
  group_split() %>%
  lapply(function(df) {
    ggplot(df) +
      # Acute segment
      geom_segment(aes(x = 0, xend = n_acute_species, y = virus_clean, yend = virus_clean, colour = log10(total_assayed)),
                   linewidth = 1.2, alpha = 1) +
      # Prior extension
      geom_segment(aes(x = n_acute_species, xend = total_species, y = virus_clean, yend = virus_clean, colour = log10(total_assayed)),
                   linewidth = 1.2, alpha = 0.4) +
  # Acute point (fully opaque)      
      geom_point(aes(x = n_acute_species, y = virus_clean, colour = log10(total_assayed)),
               size = 3, alpha = 1) +
      # Prior point (transparent)
      geom_point(aes(x = total_species, y = virus_clean, colour = log10(total_assayed)),
                 size = 3, alpha = 0.4) +
      scale_colour_viridis_c(name = "N samples assayed (log10)", direction = -1) +
      labs(x = "Number of host species with evidence of infection",
           y = NULL) +
      theme_minimal(base_size = 11) +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal")
    })

save_plot(plot_grid(plotlist = hp_plots, ncol = 2), filename = here("output", "simple_hp.png"), base_width = 12, base_height = 8, bg = "transparent")

hp_simple_acute <- subset_hp %>%
  filter(str_detect(family, "Hanta|Arena")) %>%
  filter(str_detect(assay_clean, "Acute")) %>%
  filter(n_positive >= 1) %>%
  group_by(virus_clean, family) %>%
  summarise(n_species = length(unique(host_species)),
            n_assayed = sum(n_assayed, na.rm = TRUE)) %>%
  ggplot() +
  geom_point(aes(x = n_species, y = virus_clean, colour = log(n_assayed))) +
  geom_segment(aes(x = 0, xend = n_species, y = virus_clean, yend = virus_clean, colour = log(n_assayed)),
               linewidth = 1) +
  scale_colour_viridis_c(direction = -1) +
  facet_grid(~family, scales = "free") +
  labs(y = element_blank(),
       x = "Number of species",
       colour = "N samples assayed (log10)") +
  theme_minimal()

save_plot(hp_simple_acute, filename = here("output", "simple_hp_acute.png"), base_width = 8, base_height = 8, bg = "transparent")

hanta_hp <- ggplot() +
  geom_point(data = subset_hp %>%
                   filter(family == "Hantaviridae") %>% drop_na(assay_clean),
                 aes(x = prop_pos, y = virus_clean, colour = host_genus, size = log10(n_assayed)),
             position = position_jitter(width = NULL)) +
  guides(colour = "none") +
  # geom_text_repel(data = labelled_hanta,
  #              aes(y = virus_clean, x = prop_pos, label = host_species),
  #              max.overlaps = 30,
  #              min.segment.length = 0,
  #              size = 3,
  #              force = 1.2) +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "black", linewidth = 0.2)) +
  facet_grid(~ assay_clean) +
  labs(title = "Hantaviruses",
       x = "Proportion positive",
       y = element_blank(),
       size = "Number tested (log10)",
       caption = "Labels only shown for >25 assayed and >1 positive") +
  expand_limits(x = c(0, max(subset_hp$prop_pos) * 1.1),
                y = c(0, nrow(subset_hp %>%
                                  filter(family == "Hantaviridae") %>%
                                  distinct(virus_clean)) * 1.1)) +
  coord_cartesian(xlim = c(0, 1))

ggsave(plot = hanta_hp, filename = here("output", "hantavirus_hp.png"), width = 12, height = 8)

labelled_arena <-  subset_hp %>%
  filter(family == "Arenaviridae") %>%
  group_by(virus_clean, host_species, assay_clean) %>%
  filter(n_assayed >= 25 &
           n_positive > 1) %>%
  arrange(-prop_pos) %>%
  slice(1)

mammarena_hp <- ggplot() +
  geom_point(data = subset_hp %>%
               filter(family == "Arenaviridae") %>%
               drop_na(assay_clean),
             aes(x = prop_pos, y = virus_clean, colour = host_genus, size = log10(n_assayed)),
             position = position_jitter(width = NULL)) +
  guides(colour = "none") +
  # geom_text_repel(data = labelled_arena,
  #                 aes(y = virus_clean, x = prop_pos, label = host_species),
  #                 max.overlaps = 30,
  #                 min.segment.length = 0,
  #                 size = 3,
  #                 force = 1.2) +
  theme_minimal() +
  theme(panel.grid.major.y = element_line(color = "black", linewidth = 0.2)) +
  facet_grid(~ assay_clean) +
  labs(title = "Arenaviridae",
       x = "Proportion positive",
       y = element_blank(),
       size = "Number tested (log10)",
       caption = "Labels only shown for >25 assayed and >1 positive") +
  expand_limits(x = c(0, max(subset_hp$prop_pos) * 1.1),
                y = c(0, nrow(subset_hp %>%
                                filter(family == "Arenaviridae") %>%
                                distinct(virus_clean)) * 1.1)) +
  coord_cartesian(xlim = c(0, 1))

ggsave(plot = mammarena_hp, filename = here("output", "arenavirus_hp.png"), width = 12, height = 8)


# Host ranges -------------------------------------------------------------

examples <- c("Orthohantavirus sinnombreense", "Mammarenavirus lassaense")

iucn_all <- vect(here("data", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))

map_hp <- combined_data$pathogen %>%
  filter(str_detect(assay_clean, "PCR|Sequencing|Culture")) %>%
  filter(virus_clean %in% examples) %>%
  filter(!study_id %in% c("grant_52")) %>%
  group_by(virus_clean) %>%
  group_split()

hp_1_assays <- map_hp[[1]] %>%
  group_by(assay_clean, host_species, host_genus, decimalLatitude, decimalLongitude) %>%
  summarise(n_assayed = sum(n_assayed, na.rm = TRUE),
            n_positive = sum(n_positive, na.rm = TRUE)) %>%
  arrange(-n_assayed)

hp_1_hosts <- hp_1_assays %>%
  filter(n_positive >= 1) %>%
  drop_na(host_species) %>%
  pull(host_species) %>%
  unique()

hp_1_ranges <- iucn_all[iucn_all$sci_name %in% hp_1_hosts] %>%
  as_tibble(geom = "WKT") %>%
  select(sci_name, geometry) %>%
  drop_na(geometry) %>%
  vect(geom = "geometry", crs = "EPSG:4326")

hp_1_points <- map_hp[[1]] %>%
  distinct(host_species, host_genus, n_assayed, n_positive, decimalLatitude, decimalLongitude) %>%
  filter(n_positive >= 1) %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

plot_mammarena <- list()

species <- sort(unique(hp_1_points$host_species))

species_points <- hp_1_points %>%
  filter(host_species %in% species)

species_range <- hp_1_ranges %>%
  filter(sci_name %in% species)

extent <- ext(species_range)

plot_mammarena <- ggplot() +
  geom_spatvector(data = world_shapefile) +
  geom_spatvector(data = species_range, fill = "red", alpha = 0.2) +
  geom_spatvector(data = species_points, aes(colour = ifelse(n_positive == 0, NA, n_positive), size = n_assayed), alpha = 0.4) +
  scale_size_area() +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Mammarenavirus lassaense", 
       colour = "Positive",
       size = "Tested")

ggsave(plot = plot_mammarena, filename = here("output", "sample_maps", "arenaviridae", "mammarenavirus_lassaense.png"), width = 6)


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
    geom_spatvector(data = world_shapefile) +
    geom_spatvector(data = species_range, fill = "red", alpha = 0.2) +
    geom_spatvector(data = species_points, aes(colour = n_positive, size = n_assayed), alpha = 0.4) +
    scale_size_area() +
    scale_color_viridis_c() +
    theme_minimal() +
    labs(title = paste(species), 
         colour = "Positive",
         size = "Tested") +
    coord_sf(xlim = extent[1:2], ylim = extent[3:4])
  
  ggsave(plot = plot_hantavir, filename = here("output", "sample_maps", "hantaviridae", paste0(species, ".png")), width = 6)
  
}


# Sequences ---------------------------------------------------------------

save_plot(
  combined_data$sequences %>%
  filter(sequenceType == "Pathogen") %>%
  filter(!str_detect(virus_clean, "Aza")) %>%
  distinct(accession_number, virus_clean) %>%
  group_by(virus_clean) %>%
  summarise(n = n()) %>%
  arrange(-n) %>%
  mutate(virus_clean = fct_inorder(virus_clean)) %>%
  mutate(virus_clean_rev = fct_rev(virus_clean)) %>%
  mutate(cumsum_n = cumsum(n)) %>%
  ggplot(aes(x = 1, y = n, fill = virus_clean_rev, alpha = n)) +
  geom_col(position = "stack") +
  geom_label(aes(label = n, y = cumsum_n - n / 2), fill = "white") +
  scale_fill_viridis_d(direction = -1) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  guides(
    fill = guide_legend(),    # keep this
    alpha = "none"            # remove alpha legend
  ) +
  labs(x = NULL, y = NULL, fill = "Virus name") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")     
  ),
  filename = here("output", "sequences_available.png"), base_height = 8, base_width = 12)

