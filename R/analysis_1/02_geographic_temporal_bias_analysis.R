# ---
# 02_geographic_temporal_bias_analysis.R
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
library(rnaturalearth)
library(lubridate)
library(WDI) # For World Bank socio-economic data

# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-08-21.rds") # <-- Update this date
arha_db <- read_rds(db_path)
host_data <- arha_db$host

# --- 2. Load Additional Data ---
# a) Host Range Maps (e.g., from IUCN)
# This should be a shapefile or geopackage.
try({
  # Example for a single genus; this would need to be looped or combined for multiple taxa
  apodemus_range <- st_read(here("data", "external", "iucn_ranges", "apodemus.shp"))
}, silent = TRUE)

# b) Socio-economic Data (from World Bank)
# This downloads data directly.
try({
  socio_economic_data <- WDI(
    country = "all",
    indicator = c("NY.GDP.PCAP.CD", "GB.XPD.RSDV.GD.ZS"), # GDP per capita, R&D expenditure
    start = 2020, end = 2020,
    extra = TRUE
  ) %>%
    as_tibble() %>%
    filter(region != "Aggregates") %>%
    select(iso3c, gdp_per_capita = NY.GDP.PCAP.CD, rd_expenditure = GB.XPD.RSDV.GD.ZS)
}, silent = TRUE)


# --- 3. Geographic Bias Analysis ---

# 3.1. Generate map of all sampling locations
world_map <- ne_countries(scale = "medium", returnclass = "sf")
sampling_points <- host_data %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

ggplot() +
  geom_sf(data = world_map, fill = "grey90") +
  geom_sf(data = sampling_points, size = 1, alpha = 0.5, color = "darkred") +
  labs(title = "Global Distribution of Sampling Locations") +
  theme_void()

ggsave(here("R", "analysis_1", "outputs", "sampling_locations_map.png"), width = 12, height = 7)

# 3.2. Overlay sampling points on a host range map (example)
if (exists("apodemus_range")) {
  ggplot() +
    geom_sf(data = apodemus_range, fill = "lightblue", alpha = 0.6) +
    geom_sf(data = world_map, fill = NA, color = "grey50") +
    geom_sf(data = sampling_points %>% filter(host_genus == "Apodemus"), color = "darkred") +
    labs(title = "Sampling Locations for Apodemus vs. Known Geographic Range") +
    theme_void()
  
  ggsave(here("R", "analysis_1", "outputs", "apodemus_range_vs_sampling.png"), width = 12, height = 7)
}

# 3.3. Correlate sampling intensity with socio-economic data
if (exists("socio_economic_data")) {
  country_sampling_intensity <- host_data %>%
    group_by(iso3c) %>%
    summarise(total_hosts_sampled = sum(number_of_hosts, na.rm = TRUE), .groups = "drop") %>%
    left_join(socio_economic_data, by = "iso3c") %>%
    filter(!is.na(gdp_per_capita), !is.na(total_hosts_sampled))
  
  # Fit a simple linear model
  gdp_model <- lm(log1p(total_hosts_sampled) ~ log10(gdp_per_capita), data = country_sampling_intensity)
  
  summary(gdp_model)
  
  # Plot the correlation
  ggplot(country_sampling_intensity, aes(x = gdp_per_capita, y = total_hosts_sampled)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_x_log10(labels = scales::dollar) +
    scale_y_log10() +
    labs(
      title = "Correlation between Sampling Intensity and GDP",
      x = "GDP per Capita (log scale)",
      y = "Total Hosts Sampled (log scale)"
    ) +
    theme_minimal()
  
  ggsave(here("R", "analysis_1", "outputs", "gdp_vs_sampling_plot.png"), width = 8, height = 6)
}


# --- 4. Temporal Bias Analysis ---

# 4.1. Plot number of hosts sampled per year
temporal_summary <- host_data %>%
  filter(!is.na(start_date)) %>%
  mutate(year = year(start_date)) %>%
  group_by(year) %>%
  summarise(n_hosts_sampled = sum(number_of_hosts, na.rm = TRUE)) %>%
  filter(year > 1950) # Filter out potential date parsing errors

ggplot(temporal_summary, aes(x = year, y = n_hosts_sampled)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Total Number of Hosts Sampled Over Time",
    x = "Year of Sample Collection",
    y = "Number of Hosts Sampled"
  ) +
  theme_minimal()

ggsave(here("R", "analysis_1", "outputs", "temporal_sampling_trend.png"), width = 10, height = 6)
