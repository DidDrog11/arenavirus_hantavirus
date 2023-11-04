source(here::here("R", "00_load_data.R"))

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
  labs(title = "Countries with included studies (N = 790)",
       fill = "Number of publications")

save_plot(study_map, filename = here("output", "study_map.png"), base_width = 8, base_height = 8)
