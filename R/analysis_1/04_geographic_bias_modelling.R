# ---
# 04_geographic_bias_modeling.R
#
# Purpose: To fit hierarchical models testing the drivers of geographic sampling bias
# at both the country and subnational (ADM2) scales.
# ---

# --- 1. Setup ---
library(here)
library(readr)
library(dplyr)
library(sf)
library(terra)
library(tidyterra)
library(WDI)
library(lme4)
library(countrycode)
library(MASS)
library(brms)
library(wopr)
library(performance)
library(scales)
library(tidyverse)
library(tidybayes)
library(marginaleffects)
library(patchwork)
library(biscale)
library(cowplot)
library(rnaturalearth)

conflicted::conflict_prefer_all("MASS", quiet = TRUE)
conflicted::conflict_prefer_all("brms", quiet = TRUE)
conflicted::conflict_prefer_all("purrr", quiet = TRUE)
conflicted::conflict_prefer_all("terra", quiet = TRUE)
conflicted::conflict_prefer_all("readr", quiet = TRUE)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
conflicted::conflict_prefer_all("tidyr", quiet = TRUE)

# --- 2. Load Pre-processed Data ---
# Main database object
db_path <- here("data", "database", "Project_ArHa_database_2025-09-25.rds")
arha_db <- read_rds(db_path)
host_data <- arha_db$host

# Harmonized host traits table (contains gbif_id, scientific_name, etc.)
host_traits <- read_rds(here("data", "external", "harmonised_species.rds"))

# Simplified GADM ADM2 shapefile (created in the previous script)
gadm_adm2_proj <- vect(here("data", "gadm", "gadm_adm2_simplified.shp"))

# IUCN range data (needed to create the host richness raster)
iucn_ranges_proj <- vect(here("data", "external", "iucn_subset.shp"))

# --- 3. Country-Level Analysis ---
# 3.1. Prepare country-level response variables
country_summary <- host_data %>%
  drop_na(gbif_id) %>%
  drop_na(iso3c) %>%
  group_by(iso3c) %>%
  summarise(
    # Model 1 Response: Sampling Intensity
    total_hosts_sampled = sum(number_of_hosts, na.rm = TRUE),
    # Model 2 Response: Sampling Breadth
    n_species_sampled = n_distinct(gbif_id),
    .groups = "drop"
  ) %>%
  filter(total_hosts_sampled > 0)

# 3.2. Get country-level predictor (GDP)
median_year_per_country <- host_data %>%
  drop_na(iso3c) %>%
  mutate(start_year = lubridate::year(start_date),
         end_year = lubridate::year(end_date)) %>%
  distinct(iso3c, start_year, end_year) %>%
  rowwise() %>%
  mutate(mid_year = round(mean(c(start_year, end_year), na.rm = TRUE), 0)) %>%
  ungroup() %>%
  group_by(iso3c) %>%
  summarise(median_year = as.integer(median(mid_year, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(iso2c = countrycode(iso3c, origin = "iso3c", destination = "iso2c"))

gdp_data <- purrr::map_dfr(1:nrow(median_year_per_country),
                           function(i) {
                             iso2_code <- median_year_per_country$iso2c[i] 
                             iso3_code <- median_year_per_country$iso3c[i] 
                             year <- median_year_per_country$median_year[i]
                             result <- try(WDI(country = iso2_code,
                                               indicator = c("gdp_per_capita" = "NY.GDP.PCAP.KD",
                                                             "population" = "SP.POP.TOTL"),
                                               start = year,
                                               end = year), silent = TRUE)
                             if (inherits(result, "try-error")) {
                               return(tibble(iso2c = iso2_code, iso3c = iso3_code, status = "failed"))
                             } else {
                               return(result)
                             }
                           })


failed_countries <- gdp_data %>%
  filter(is.na(gdp_per_capita) |
           status == "failed") %>%
  distinct(iso2c, iso3c)

if(nrow(failed_countries) > 0) {
  
  gdp_data_retry <- WDI(country = failed_countries$iso2c,
                        indicator = c("gdp_per_capita" = "NY.GDP.PCAP.KD",
                                      "population" = "SP.POP.TOTL"),
                        start = 1980, 
                        end = 2025,
                        extra = TRUE) %>%
    as_tibble() %>%
    filter(!is.na(gdp_per_capita)) %>%
    left_join(median_year_per_country %>% 
                select(iso3c, median_year), by = "iso3c") %>%
    group_by(iso3c) %>%
    mutate(year_diff = abs(year - median_year)) %>%
    slice_min(order_by = year_diff, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  gdp_data_final <- gdp_data %>%
      left_join(gdp_data_retry %>% 
                  select(country, year, iso3c, gdp_per_capita_retry = gdp_per_capita, population_retry = population), 
      by = "iso3c") %>%
    mutate(gdp_per_capita = coalesce(gdp_per_capita, gdp_per_capita_retry),
           population = coalesce(population, population_retry),
           country = coalesce(country.x, country.y),
           year = coalesce(year.y, year.x)) %>%
    select(country, year, iso2c, iso3c, gdp_per_capita, population, status)
} else {
  gdp_data_final <- gdp_data
}

france_gdp <- gdp_data_final %>%
  filter(iso3c == "FRA") %>%
  pull(gdp_per_capita)

# Manually added values from https://unctadstat.unctad.org/datacentre/dataviewer/US.GDPTotal
gdp_data_imputed <- gdp_data_final %>%
  mutate(gdp_per_capita = case_when(iso3c %in% c("GUF", "MYT", "REU") & is.na(gdp_per_capita) ~ france_gdp, #2005
                                    iso3c == "PRK" ~ 632, #2010
                                    iso3c == "TWN" ~ 10131, #1994
                                    iso3c == "VEN" ~ 10719, #1996
                                    TRUE ~ gdp_per_capita),
         population = case_when(iso3c == "GUF" ~ 205954, #2005
                                iso3c == "MYT" ~ 179601, #2005
                                iso3c == "REU" ~ 781962, #2005
                                iso3c == "PRK" ~ 24987000, #2010
                                iso3c == "TWN" ~ 21357000, #1994
                                iso3c == "VEN" ~ 22662000, #1996
                                TRUE ~ gdp_per_capita)) %>%
  mutate(country = case_when(iso3c == "GUF" ~ "French Guiana",
                             iso3c == "MYT" ~ "Mayotte",
                             iso3c == "REU" ~ "Reunion",
                             iso3c == "TWN" ~ "Taiwan",
                             TRUE ~ country),
         year = case_when(iso3c %in% c("GUF", "MYT", "REU") ~ 2005,
                          iso3c == "PRK" ~ 2010,
                          iso3c == "TWN" ~ 1994,
                          iso3c == "VEN" ~ 1996,
                          TRUE ~ year),
         iso2c = if_else(iso3c == "NAM", "NA", iso2c)) %>%
  select(-status) %>%
  filter(!is.na(gdp_per_capita))

# 3.3. Create analytical dataset and fit models
country_analysis_data <- country_summary %>%
  left_join(gdp_data_imputed %>%
              select(iso3c, gdp_per_capita, population, year), by = "iso3c")

country_analysis_data_scaled <- country_analysis_data %>%
  mutate(
    gdp_per_capita_scaled = as.numeric(scale(log10(gdp_per_capita))),
    population_scaled = as.numeric(scale(log10(population)))
  )

# Model 1: GDP vs. Sampling Intensity
model_intensity_bayesian <- brm(formula = log10(total_hosts_sampled) ~ gdp_per_capita_scaled + offset(log10(population)),
                                data = country_analysis_data_scaled,
                                family = gaussian(),
                                prior = c(prior(normal(0, 1), class = "b"),
                                          prior(normal(0, 5), class = "Intercept"),
                                          prior(student_t(3, 0, 2.5), class = "sigma")),
                                control = list(adapt_delta = 0.95),
                                chains = 4, iter = 4000, warmup = 1000,  cores = 4, seed = 12345)
summary(model_intensity_bayesian)
fixef(model_intensity_bayesian)

# Model 2: GDP vs. Sampling Breadth
model_breadth_bayesian <- brm(formula = n_species_sampled ~ gdp_per_capita_scaled + offset(log10(population)),
                              data = country_analysis_data_scaled,
                              family = negbinomial(),
                              prior = c(
                                prior(normal(0, 1), class = "b"),
                                prior(normal(0, 5), class = "Intercept"),
                                prior(gamma(2, 0.1), class = "shape")),
                              control = list(adapt_delta = 0.95),
                              chains = 4, iter = 4000, warmup = 1000, cores = 4, seed = 12345)
summary(model_breadth_bayesian)
fixef(model_breadth_bayesian)

# --- 4. Subnational (ADM2) Analysis ---
# 4.1. Prepare ADM2-level predictors (Zonal Statistics)
adm2_predictors_path <- here("data", "external", "adm2_predictors.rds")

if (!file.exists(adm2_predictors_path)) {
  
  template_raster <- rast(ext(gadm_adm2_proj), crs = crs(gadm_adm2_proj), resolution = 25000)
  
  # a) Load global raster datasets (assumes you have downloaded them)
  #    - WorldPop: https://www.worldpop.org/
  #    - VIIRS Nighttime Lights: https://eogdata.mines.edu/nighttime_light/annual/v22/2024/
  #    - Travel Time to Cities: https://malariaatlas.org/research-project/accessibility-to-cities/
  pop_density_raw <- rast(here("data", "external", "rasters", "worldpop_2025_1km.tif"))
  pop_density_25km <- project(pop_density_raw, template_raster, method = "average")
  
  night_lights_raw <- rast(here("data", "external", "rasters", "viirs_2024.tif"))
  night_lights_25km <- project(night_lights_raw, template_raster, method = "average")
  
  travel_time_raw <- rast(here("data", "external", "rasters", "travel_time_to_cities.tif"))
  travel_time_25km <- project(travel_time_raw, template_raster, method = "average")
  
  # write these smaller rasters and remove previous
  writeRaster(pop_density_25km, filename = here("data", "external", "rasters", "worldpop_2025_25km.tif"))
  writeRaster(night_lights_25km, filename = here("data", "external", "rasters", "viirs_2024_25km.tif"))
  writeRaster(travel_time_25km, filename = here("data", "external", "rasters", "travel_time_25km.tif"))
  
  # b) Create Host Species Richness Raster
  iucn_for_raster <- iucn_ranges_proj %>% mutate(species_count = 1)
  host_richness_25km <- rasterize(iucn_for_raster, template_raster, field = "species_count", fun = "sum")
  
  # c) Extract values for each predictor (zonal statistics)
  adm2_predictors <- gadm_adm2_proj %>%
    tidyterra::select(GID_0, GID_2) %>%
    mutate(area_km2 = as.numeric(expanse(.) / 1e6),
           mean_pop_density = terra::extract(pop_density_25km, ., fun = "mean", ID = FALSE, na.rm = TRUE)[[1]],
           mean_night_lights = terra::extract(night_lights_25km, ., fun = "mean", ID = FALSE, na.rm = TRUE)[[1]],
           mean_travel_time = terra::extract(travel_time_25km, ., fun = "mean", ID = FALSE, na.rm = TRUE)[[1]],
           mean_host_richness = terra::extract(host_richness_25km, ., fun = "mean", ID = FALSE, na.rm = TRUE)[[1]])
  
  adm2_predictors_df <- as_tibble(adm2_predictors) %>%
    mutate(mean_night_lights = if_else(mean_night_lights < 0, 0, mean_night_lights),
           mean_pop_density = coalesce(mean_pop_density, 0),
           mean_night_lights = coalesce(mean_night_lights, 0),
           mean_host_richness = coalesce(mean_host_richness, 0))
  
  write_rds(adm2_predictors_df, adm2_predictors_path)

  } else {
  message("Loading cached ADM2-level predictors...")
  adm2_predictors <- read_rds(adm2_predictors_path)
}

# 4.2. Prepare ADM2-level response variables
adm2_response_vars <- host_data %>%
  drop_na(gadm_adm2) %>%
  group_by(gadm_adm2) %>%
  summarise(
    n_hosts_sampled = sum(number_of_hosts, na.rm = TRUE),
    n_species_sampled = n_distinct(gbif_id)
  )

# 4.3. Create final analytical dataset
adm2_analysis_data <- adm2_predictors %>%
  left_join(adm2_response_vars, by = c("GID_2" = "gadm_adm2")) %>%
  mutate(n_hosts_sampled = coalesce(n_hosts_sampled, 0L),
         n_species_sampled = coalesce(n_species_sampled, 0L)) %>%
  drop_na(mean_pop_density) %>%
  drop_na(mean_night_lights) %>%
  drop_na(mean_travel_time) %>%
  drop_na(mean_host_richness)

# Check collinearity among predictors
vif_check_model <- lm(
  n_hosts_sampled ~ scale(log1p(mean_pop_density)) + 
    scale(log1p(mean_night_lights)) + 
    scale(log1p(mean_travel_time)) + 
    scale(mean_host_richness),
  data = adm2_analysis_data
)

# Use the performance::check_collinearity() function to get VIF scores
collinearity_check <- check_collinearity(vif_check_model)
print(collinearity_check)

# 4.4. Fit Subnational Bayesian GLMMs
priors <- c(
  prior(normal(0, 2), class = "b"),      # Prior for scaled fixed effects (slopes)
  prior(normal(0, 5), class = "Intercept"),
  prior(student_t(3, 0, 2.5), class = "sd") # Robust prior for the random effect standard deviation
)

# Model A: Sampling Intensity
# glmm_intensity_nb <- brm(
#   formula = n_hosts_sampled ~ scale(log1p(mean_pop_density)) + scale(log1p(mean_night_lights)) +
#     scale(log1p(mean_travel_time)) + scale(mean_host_richness) +
#     offset(log(area_km2)) + (1 | GID_0),
#   data = adm2_analysis_data,
#   family = negbinomial(),
#   prior = priors,
#   chains = 4, iter = 4000, warmup = 1000, cores = 4, seed = 12345,
#   control = list(adapt_delta = 0.95),
#   save_pars = save_pars(all = TRUE)
# )
# 
# print(summary(glmm_intensity_brm))
# pp_check(glmm_intensity_brm, ndraws = 100) + 
#   scale_x_log10(
#     name = "Number of Hosts Sampled (log scale)",
#     breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
#     labels = scales::comma
#   )
# pp_check(glmm_intensity_brm, type = "rootogram", style = "hanging", ndraws = 100) + 
#   scale_x_log10(
#     name = "Number of Hosts Sampled (log scale)",
#     breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
#     labels = scales::comma
#   )
# 
# write_rds(glmm_intensity_nb, here("output", "analysis_1", "glmm_intensity_model.rds"))
# glmm_intensity_nb <- read_rds(here("output", "analysis_1", "glmm_intensity_model.rds"))

# Model A2: Sampling Intensity (with zero-inflation)
glmm_intensity_zinb <- brm(
  formula = bf(
    n_hosts_sampled ~ scale(log1p(mean_pop_density)) + scale(log1p(mean_night_lights)) + 
      scale(log1p(mean_travel_time)) + scale(mean_host_richness) + 
      offset(log(area_km2)) + (1 | GID_0),
    zi ~ scale(log1p(mean_pop_density)) + scale(log1p(mean_night_lights)) + 
      scale(log1p(mean_travel_time)) + scale(mean_host_richness)
  ),
  data = adm2_analysis_data,
  family = zero_inflated_negbinomial(), 
  prior = priors,
  chains = 4, iter = 4000, warmup = 1000, 
  cores = 4, seed = 12345,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE)
)
# 
# pp_check(glmm_intensity_zinb, ndraws = 100) +
#   scale_x_log10(
#     name = "Number of Hosts Sampled (log scale)",
#     breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
#     labels = scales::comma
#   )
# 
# pp_check(glmm_intensity_zinb, type = "rootogram", style = "hanging", ndraws = 100) +
#   scale_x_log10(
#     name = "Number of Hosts Sampled (log scale)",
#     breaks = c(0, 1, 10, 100, 1000, 10000, 100000),
#     labels = scales::comma
#   )
# 
write_rds(glmm_intensity_zinb, here("output", "analysis_1", "glmm_intensity_model_zinb.rds"))
glmm_intensity_zinb <- read_rds(here("output", "analysis_1", "glmm_intensity_model_zinb.rds"))

# Model A comparison
# loo_nb <- loo(glmm_intensity_nb, cores = 4, moment_match = TRUE)
# loo_zinb <- loo(glmm_intensity_zinb, cores = 4, moment_match = TRUE)
# 
# model_comparison <- loo_compare(loo_nb, loo_zinb)

# Model B1: Sampling Breadth
# glmm_breadth_brm <- brm(
#   formula = n_species_sampled ~ scale(log1p(mean_pop_density)) + scale(log1p(mean_night_lights)) + 
#     scale(log1p(mean_travel_time)) + scale(mean_host_richness) + 
#     offset(log(area_km2)) + (1 | GID_0),
#   data = adm2_analysis_data,
#   family = negbinomial(), # NB
#   prior = priors,
#   chains = 4, iter = 4000, warmup = 1000, cores = 4, seed = 12345,
#   control = list(adapt_delta = 0.95),
#   save_pars = save_pars(all = TRUE)
# )
# 
# write_rds(glmm_breadth_brm, here("output", "analysis_1", "glmm_breadth_model.rds"))

# Model B2
glmm_breadth_zinb <- brm(
  formula = bf(
    n_species_sampled ~ scale(log1p(mean_pop_density)) + scale(log1p(mean_night_lights)) + 
      scale(log1p(mean_travel_time)) + scale(mean_host_richness) + 
      offset(log(area_km2)) + (1 | GID_0),
    zi ~ scale(log1p(mean_pop_density)) + scale(log1p(mean_night_lights)) + 
      scale(log1p(mean_travel_time)) + scale(mean_host_richness)
  ),
  data = adm2_analysis_data,
  family = zero_inflated_negbinomial(),
  prior = priors,
  chains = 4, iter = 4000, warmup = 1000, 
  cores = 4, seed = 12345,
  control = list(adapt_delta = 0.95),
  save_pars = save_pars(all = TRUE)
)

write_rds(glmm_breadth_zinb, here("output", "analysis_1", "glmm_breadth_model_zinb.rds"))
glmm_breadth_zinb <- read_rds(here("output", "analysis_1", "glmm_breadth_model_zinb.rds"))

# loo_nb_breadth <- loo(glmm_breadth_brm, cores = 4, moment_match = TRUE)
# loo_zinb_breadth <- loo(glmm_breadth_zinb, cores = 4, moment_match = TRUE)

# Compare the two models
# model_comparison_breadth <- loo_compare(loo_nb_breadth, loo_zinb_breadth)

# --- 5. Subnational visualisation ---
# Produces a map for intensity of sampling
# Create the map
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = mollweide_crs) %>%
  st_make_valid()

# 5.1 Filter for only the ADM2 units that were actually sampled
data_sampled <- adm2_analysis_data %>%
  filter(n_hosts_sampled > 0)

# 5.2 Get the predicted count for these sampled units
pred_sampled <- posterior_epred(glmm_intensity_zinb,
                                newdata = data_sampled,
                                re_formula = NA,
                                dpar = "mu") %>%
  apply(2, mean)

# 5.3 Calculate residuals on the log scale (Observed - Predicted)
residual_data <- data_sampled %>%
  mutate(pred_intensity = pred_sampled,
         log_residual = log10(n_hosts_sampled) - log10(pred_intensity),
         log_n_hosts_sampled = log10(n_hosts_sampled)) %>%
  dplyr::select(GID_0, GID_2, n_hosts_sampled, log_n_hosts_sampled, log_residual)

# 5.4 Join residuals to the spatial data
map_data_resid <- gadm_adm2_proj %>%
  filter(GID_2 %in% residual_data$GID_2) %>%
  left_join(residual_data, by = "GID_2") %>%
  filter(!is.na(log_residual))

max_res <- max(abs(map_data_resid$log_residual), na.rm = TRUE)

# 5.5 Plot the global map
p_map_residuals <- ggplot() +
  geom_sf(data = world_map, fill = "grey50", colour = "black") +
  geom_sf(data = map_data_resid, aes(fill = log_residual), colour = NA, show.legend = TRUE) +
  scale_fill_gradient2(name = "Surveillance Gap\n(log10[Obs/Pred])",
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       limits = c(-max_res, max_res)) +
  coord_sf(crs = mollweide_crs) +
  labs(title = "A)") +
  theme_minimal()

# 5.6 Generate insets
# Identify centroids of sampled regions
sampled_centroids <- map_data_resid %>% 
  centroids() %>%
  crds()

# K-means clustering to find 20 main clusters
set.seed(123)
kmeans_result <- kmeans(sampled_centroids, centers = 20, nstart = 25)
sampled_centroids <- as_tibble(sampled_centroids) %>%
  mutate(cluster = kmeans_result$cluster)

# Get the bounding box of a cluster
roi_list <- sampled_centroids %>%
  group_by(cluster) %>%
  vect(geom = c("x", "y"), crs = crs(map_data_resid)) %>%
  aggregate(by = "cluster")

# Function to produce maps
regional_map <- function(roi_list, tag) {
  ggplot() +
    geom_spatvector(data = vect(world_map), fill = "grey50", colour = "black") +
    geom_spatvector(data = map_data_resid, aes(fill = log_residual), colour = NA, show.legend = TRUE) +
    scale_fill_gradient2(name = "Surveillance Gap\n(log10[Obs/Pred])",
                         low = "blue", mid = "white", high = "red",
                         midpoint = 0,
                         limits = c(-max_res, max_res)) +
    coord_sf(xlim = c(ext(roi_list)[1], ext(roi_list)[2]), 
             ylim = c(ext(roi_list)[3], ext(roi_list)[4]), 
             expand = FALSE) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = tag)
}

inset_B <- regional_map(roi_list[3], "B)")
inset_C <- regional_map(roi_list[12], "C)")
inset_D <- regional_map(roi_list[16], "D)")
inset_E <- regional_map(roi_list[17], "E)")

labels <- tibble(label = c("B)", "C)", "D)", "E)"),
                 x = c(ext(roi_list[3])[1], ext(roi_list[12])[2], ext(roi_list[16])[2], ext(roi_list[17])[2]),
                 y = c(ext(roi_list[3])[3], ext(roi_list[12])[3], ext(roi_list[16])[4], ext(roi_list[17])[3])) |>
  vect(geom = c("x", "y"), crs = mollweide_crs)

# Add insets to global map
p_map_residuals_updated <- p_map_residuals +
  geom_spatvector(data = vect(ext(roi_list[3]), crs = mollweide_crs), fill = NA, lwd = 0.8, colour = "darkgreen") +
  geom_spatvector(data = vect(ext(roi_list[12]), crs = mollweide_crs), fill = NA, lwd = 0.8, colour = "darkgreen") +
  geom_spatvector(data = vect(ext(roi_list[16]), crs = mollweide_crs), fill = NA, lwd = 0.8, colour = "darkgreen") +
  geom_spatvector(data = vect(ext(roi_list[17]), crs = mollweide_crs), fill = NA, lwd = 0.8, colour = "darkgreen") +
  geom_spatvector_label(data = labels, aes(label = label), size = 3) +
  theme_minimal()

inset_grid <- (inset_B | inset_C) / 
  (inset_D | inset_E)

shared_legend <- get_legend(p_map_residuals_updated + theme(legend.position = "bottom"))

final_regional_dashboard <- ((p_map_residuals_updated + theme(legend.position = "none", axis.title = element_blank())) / inset_grid) / shared_legend +
  plot_layout(heights = c(1, 2, 0.5))

ggsave(here("output", "analysis_1", "regional_surveillance_gaps_dashboard.png"), plot = final_regional_dashboard, width = 8, height = 10, dpi = 300)

# Reproduces the map for the breadth of sampling
# Get the predicted count for these sampled units
pred_sampled_br <- posterior_epred(glmm_breadth_zinb,
                                   newdata = data_sampled,
                                   re_formula = NA,
                                   dpar = "mu") %>%
  apply(2, mean)

# Calculate residuals on the log scale (Observed - Predicted)
residual_data_br <- data_sampled %>%
  mutate(pred_intensity_br = pred_sampled_br,
         log_residual = log10(n_species_sampled) - log10(pred_sampled_br),
         log_species_sampled = log10(n_species_sampled)) %>%
  dplyr::select(GID_0, GID_2, n_species_sampled, log_species_sampled, log_residual)

# Join residuals to the spatial data
map_data_resid_br <- gadm_adm2_proj %>%
  filter(GID_2 %in% residual_data_br$GID_2) %>%
  left_join(residual_data_br, by = "GID_2") %>%
  filter(!is.na(log_residual))

max_res_br <- max(abs(map_data_resid_br$log_residual), na.rm = TRUE)

# Plot the global map
p_map_residuals_br <- ggplot() +
  geom_sf(data = world_map, fill = "grey50", colour = "black") +
  geom_sf(data = map_data_resid_br, aes(fill = log_residual), colour = NA, show.legend = TRUE) +
  scale_fill_gradient2(name = "Surveillance Gap\n(log10[Obs/Pred])",
                       low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       limits = c(-max_res_br, max_res_br)) +
  coord_sf(crs = mollweide_crs) +
  labs(title = "A)") +
  theme_minimal()

# Function to produce maps
regional_map_br <- function(roi_list, tag) {
  ggplot() +
    geom_spatvector(data = vect(world_map), fill = "grey50", colour = "black") +
    geom_spatvector(data = map_data_resid_br, aes(fill = log_residual), colour = NA, show.legend = TRUE) +
    scale_fill_gradient2(name = "Surveillance Gap\n(log10[Obs/Pred])",
                         low = "blue", mid = "white", high = "red",
                         midpoint = 0,
                         limits = c(-max_res_br, max_res_br)) +
    coord_sf(xlim = c(ext(roi_list)[1], ext(roi_list)[2]), 
             ylim = c(ext(roi_list)[3], ext(roi_list)[4]), 
             expand = FALSE) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = tag)
}

inset_B_br <- regional_map_br(roi_list[3], "B)")
inset_C_br <- regional_map_br(roi_list[12], "C)")
inset_D_br <- regional_map_br(roi_list[16], "D)")
inset_E_br <- regional_map_br(roi_list[17], "E)")

# Add insets to global map
p_map_residuals_updated_br <- p_map_residuals_br +
  geom_spatvector(data = vect(ext(roi_list[3]), crs = mollweide_crs), fill = NA, lwd = 0.8, colour = "darkgreen") +
  geom_spatvector(data = vect(ext(roi_list[12]), crs = mollweide_crs), fill = NA, lwd = 0.8, colour = "darkgreen") +
  geom_spatvector(data = vect(ext(roi_list[16]), crs = mollweide_crs), fill = NA, lwd = 0.8, colour = "darkgreen") +
  geom_spatvector(data = vect(ext(roi_list[17]), crs = mollweide_crs), fill = NA, lwd = 0.8, colour = "darkgreen") +
  geom_spatvector_label(data = labels, aes(label = label), size = 3) +
  theme_minimal()

inset_grid_br <- (inset_B_br | inset_C_br) / 
  (inset_D_br | inset_E_br)

shared_legend_br <- get_legend(p_map_residuals_updated_br + theme(legend.position = "bottom"))

final_regional_dashboard_br <- ((p_map_residuals_updated_br + theme(legend.position = "none", axis.title = element_blank())) / inset_grid_br) / shared_legend_br +
  plot_layout(heights = c(1, 2, 0.5))

ggsave(here("output", "analysis_1", "regional_surveillance_gaps_dashboard_breadth.png"), plot = final_regional_dashboard_br, width = 8, height = 10, dpi = 300)

# Tabular format of hotspots and coldspots
country_summary_intensity <- residual_data %>%
  group_by(GID_0) %>%
  summarise(country_mean_residual = mean(log_residual, na.rm = TRUE),
            country_residual_range = paste0("[", round(min(log_residual, na.rm = TRUE), 2), " - ", round(max(log_residual, na.rm = TRUE), 2), "]"),
            adm2_sampled = n())

country_summary <- as_tibble(gadm_adm2_proj) %>% 
  group_by(GID_0) %>% 
  summarise(adm2_total = n())

hs_cs_intensity <- residual_data %>%
  mutate(data = "Sampling intensity (number of individuals)") %>%
  arrange(-log_residual) %>%
  mutate(row_number = row_number()) %>%
  filter(row_number %in% c(1:10, 1329:1338)) %>%
  mutate(classification = case_when(row_number %in% 1:10 ~ "Hotspot",
                                    TRUE ~ "Coldspot")) %>%
  left_join(as_tibble(gadm_adm2_proj) %>%
              dplyr::select(GID_2, COUNTRY, NAME_2)) %>%
  mutate(NAME_2 = case_when(GID_2 == "GBR.1.25_1" ~ "Derbyshire",
                            TRUE ~ NAME_2),
         observed = n_hosts_sampled,
         predicted = exp(log_n_hosts_sampled)) %>%
  left_join(country_summary) %>%
  left_join(country_summary_intensity) %>%
  rowwise() %>%
  mutate(adm2_sampled_total = paste0(adm2_sampled, "/", adm2_total, " (", round((adm2_sampled/adm2_total) * 100, 1), "%)"),
         country_residual = paste0(round(country_mean_residual, 2), " ", country_residual_range)) %>%
  ungroup() %>%
  dplyr::select(data, classification, country = COUNTRY, adm2_sampled_total, adm2_name = NAME_2, observed, predicted, log_residual, country_residual)
  
hs_cs_breadth <- residual_data_br %>%
  mutate(data = "Sampling breadth (number of species)") %>%
  arrange(-log_residual) %>%
  mutate(row_number = row_number()) %>%
  filter(row_number %in% c(1:10, 1329:1338)) %>%
  mutate(classification = case_when(row_number %in% 1:10 ~ "Hotspot",
                                    TRUE ~ "Coldspot")) %>%
  left_join(as_tibble(gadm_adm2_proj) %>%
              dplyr::select(GID_2, COUNTRY, NAME_2)) %>%
  mutate(observed = n_species_sampled,
         predicted = exp(log_species_sampled)) %>%
  left_join(country_summary) %>%
  left_join(country_summary_intensity) %>%
  rowwise() %>%
  mutate(adm2_sampled_total = paste0(adm2_sampled, "/", adm2_total, " (", round((adm2_sampled/adm2_total) * 100, 1), "%)"),
         country_residual = paste0(round(country_mean_residual, 2), " ", country_residual_range)) %>%
  ungroup() %>%
  dplyr::select(data, classification, country = COUNTRY, adm2_sampled_total, adm2_name = NAME_2, observed, predicted, log_residual, country_residual)

write_rds(list(intensity = hs_cs_intensity,
               breadth = hs_cs_breadth),
          here("output", "analysis_1", "hot_spot.rds"))

# --- 6. Marginal Effects --- #
# Interpret covariate marginal effects

create_mfx_plot <- function(model, predictor, dpar_component,use_log_scale = FALSE) {
  
  y_label <- if (dpar_component == "mu") "Predicted Count (mu)" else "Zero-Inflation Prob. (zi)"
  
  # Apply log scale if requested
  if (use_log_scale) {
    
    observed_data <- tibble(id = 1:length(model$data[[predictor]] |>
                                            unique()),
                            predictor = model$data[[predictor]] |>
                              unique() + 0.001)
    
    x_scale_range <- summary(log10(observed_data$predictor + 0.001))
    
    p <- plot_predictions(model,
                          condition = predictor,
                          type = "response",
                          dpar = dpar_component) +
      geom_rug(data = observed_data, aes(y = NULL, x = predictor), sides = "b", alpha = 0.1, inherit.aes = FALSE) +
      theme_minimal() +
      scale_x_log10() +
      coord_cartesian(xlim = c(x_scale_range[1], x_scale_range[6])) +
      labs(x = paste0(str_to_title(str_replace_all(predictor, "_", " ")), " log()"), y = y_label)
    
  } else {
    
    observed_data <- tibble(id = 1:length(model$data[[predictor]] |>
                                            unique()),
                            predictor = model$data[[predictor]] |>
                              unique())
    
    x_scale_range <- summary(observed_data$predictor)
    
    p <- plot_predictions(model,
                          condition = predictor,
                          type = "response",
                          dpar = dpar_component) +
      geom_rug(data = observed_data, aes(y = NULL, x = predictor), sides = "b", alpha = 0.1, inherit.aes = FALSE) +
      theme_minimal() +
      coord_cartesian(xlim = c(x_scale_range[1], x_scale_range[6])) +
      labs(x = str_to_title(str_replace_all(predictor, "_", " ")), y = y_label)
    
  }
  return(p)
}

# --- Intensity Model (glmm_intensity_zinb) ---

# Pop Density (log scale)
p_int_pop_mu <- create_mfx_plot(glmm_intensity_zinb, "mean_pop_density", "mu", use_log_scale = TRUE)
p_int_pop_zi <- create_mfx_plot(glmm_intensity_zinb, "mean_pop_density", "zi", use_log_scale = TRUE)

# Night Lights (linear scale)
p_int_nlight_mu <- create_mfx_plot(glmm_intensity_zinb, "mean_night_lights", "mu")
p_int_nlight_zi <- create_mfx_plot(glmm_intensity_zinb, "mean_night_lights", "zi")

# Travel Time (log scale)
p_int_travel_mu <- create_mfx_plot(glmm_intensity_zinb, "mean_travel_time", "mu", use_log_scale = TRUE)
p_int_travel_zi <- create_mfx_plot(glmm_intensity_zinb, "mean_travel_time", "zi", use_log_scale = TRUE)

# Host Richness (linear scale)
p_int_rich_mu <- create_mfx_plot(glmm_intensity_zinb, "mean_host_richness", "mu")
p_int_rich_zi <- create_mfx_plot(glmm_intensity_zinb, "mean_host_richness", "zi")

# --- Breadth Model (glmm_breadth_zinb) ---

# Pop Density (log scale)
p_brd_pop_mu <- create_mfx_plot(glmm_breadth_zinb, "mean_pop_density", "mu", use_log_scale = TRUE)
p_brd_pop_zi <- create_mfx_plot(glmm_breadth_zinb, "mean_pop_density", "zi", use_log_scale = TRUE)

# Night Lights (linear scale)
p_brd_nlight_mu <- create_mfx_plot(glmm_breadth_zinb, "mean_night_lights", "mu")
p_brd_nlight_zi <- create_mfx_plot(glmm_breadth_zinb, "mean_night_lights", "zi")

# Travel Time (log scale)
p_brd_travel_mu <- create_mfx_plot(glmm_breadth_zinb, "mean_travel_time", "mu", use_log_scale = TRUE)
p_brd_travel_zi <- create_mfx_plot(glmm_breadth_zinb, "mean_travel_time", "zi", use_log_scale = TRUE)

# Host Richness (linear scale)
p_brd_rich_mu <- create_mfx_plot(glmm_breadth_zinb, "mean_host_richness", "mu")
p_brd_rich_zi <- create_mfx_plot(glmm_breadth_zinb, "mean_host_richness", "zi")


# --- Figure 1: Sampling Intensity Model ---
# A 4-row, 2-column grid
plot_intensity_grid <- 
  (p_int_pop_mu + p_int_pop_zi) /
  (p_int_nlight_mu + p_int_nlight_zi) /
  (p_int_travel_mu + p_int_travel_zi) /
  (p_int_rich_mu + p_int_rich_zi) +
  plot_annotation(title = "Conditional Predictions: Sampling Intensity (n_hosts_sampled)",
                  subtitle = "Left: Predicted Count (mu). Right: Zero-Inflation Probability (zi). \nAll plots hold other covariates at their mean/mode.")

# Save the plot
ggsave(here("output", "analysis_1", "conditional_predictions_intensity.png"), 
       plot = plot_intensity_grid, width = 10, height = 14, dpi = 300)

# --- Figure 2: Sampling Breadth Model ---
# A 4-row, 2-column grid
plot_breadth_grid <-
  (p_brd_pop_mu + p_brd_pop_zi) /
  (p_brd_nlight_mu + p_brd_nlight_zi) /
  (p_brd_travel_mu + p_brd_travel_zi) /
  (p_brd_rich_mu + p_brd_rich_zi) +
  plot_annotation(title = "Conditional Predictions: Sampling Breadth (n_species_sampled)",
                  subtitle = "Left: Predicted Count (mu). Right: Zero-Inflation Probability (zi). \nAll plots hold other covariates at their mean/mode.")

# Save the plot
ggsave(here("output", "analysis_1", "conditional_predictions_breadth.png"), 
       plot = plot_breadth_grid, width = 10, height = 14, dpi = 300)

# --- 7. Additional Model Checks ---
# --- 7.1 Posterior Predictive Checks (PPCs) ---
# --- Intensity Model ---
# Check 1: Does the model's simulated data distribution (light blue)
# match the observed data (dark blue)?
pp_check(glmm_intensity_zinb, type = "hist", ndraws = 50) +
  scale_x_continuous(trans = 'log1p') + 
  coord_cartesian(xlim = c(0, 2e5), ylim = c(0, 800)) +
  ggtitle("Intensity: Histogram (Log1p Scale)")

# Check 2: Does the model predict the correct *proportion* of zeros?
# This is the single most important check for a ZINB model.
prop_zero <- function(y) {mean(y == 0)}
pp_check(glmm_intensity_zinb, type = "stat", stat = "prop_zero", ndraws = 100) +
  ggtitle("Intensity: Proportion of Zeros")

# --- Breadth Model ---
pp_check(glmm_breadth_zinb, type = "hist", ndraws = 50) +
  scale_x_continuous(trans = 'log1p') + 
  coord_cartesian(xlim = c(0, 2e5), ylim = c(0, 800)) +
  ggtitle("Breadth: Histogram of Observed vs. Predicted")

pp_check(glmm_breadth_zinb, type = "stat", stat = "prop_zero", ndraws = 100) +
  ggtitle("Breadth: Proportion of Zeros")

# --- 3. DHARMa Randomized Quantile Residuals ---
library(DHARMa)
# n = 250 (i.e., 250 posterior draws) for a fast-yet-reliable check.
n_sims <- 250
y_obs_int <- glmm_intensity_zinb$data$n_hosts_sampled
y_sim_int <- posterior_predict(glmm_intensity_zinb,
                               ndraws = n_sims,
                               re_formula = NA)
y_pred_int_draws <- posterior_epred(glmm_intensity_zinb,
                                    ndraws = n_sims,
                                    re_formula = NA)
y_pred_int_mean <- apply(y_pred_int_draws, 2, mean)

dharma_res_int <- createDHARMa(simulatedResponse = t(y_sim_int),
                               observedResponse = y_obs_int,
                               fittedPredictedResponse = y_pred_int_mean,
                               integerResponse = TRUE)

plot(dharma_res_int)

y_obs_brd <- glmm_breadth_zinb$data$n_species_sampled

y_sim_brd_raw <- posterior_predict(glmm_breadth_zinb,
                                   ndraws = n_sims,
                                   re_formula = NA)

y_pred_brd_draws <- posterior_epred(glmm_breadth_zinb,
                                    ndraws = n_sims,
                                    re_formula = NA)

y_pred_brd_mean <- apply(y_pred_brd_draws, 2, mean)

# THE FIX: Transpose the simulated response matrix
dharma_res_brd <- createDHARMa(simulatedResponse = t(y_sim_brd_raw),
                               observedResponse = y_obs_brd,
                               fittedPredictedResponse = y_pred_brd_mean,
                               integerResponse = TRUE)

plot(dharma_res_brd)

# --- 4. DHARMa Specific Pattern Checks ---
# We can use the DHARMa residuals to test for your two specific concerns:
# 1. Residuals vs. Predictors (Systematic error?)
# 2. Spatial Autocorrelation (Missing spatial component?)

# --- 4a. Residuals vs. Predictors ---
# Looking for *no patterns* (i.e., flat red lines).
model_data <- glmm_intensity_zinb$data

# 2. Create the predictor vectors as they are in the formula
pred_pop_density <- as.numeric(scale(log1p(model_data$mean_pop_density)))
pred_night_lights <- as.numeric(scale(log1p(model_data$mean_night_lights)))
pred_travel_time <- as.numeric(scale(log1p(model_data$mean_travel_time)))
pred_host_richness <- as.numeric(scale(model_data$mean_host_richness))

# 3. Plot residuals against each predictor vector one by one
#plotResiduals(dharma_res_int, predictor = pred_pop_density, quantreg = TRUE)
#plotResiduals(dharma_res_int, predictor = pred_night_lights, quantreg = TRUE)
#plotResiduals(dharma_res_int, predictor = pred_travel_time, quantreg = TRUE)
#plotResiduals(dharma_res_int, predictor = pred_host_richness, quantreg = TRUE)

# --- 4b. Spatial Autocorrelation ---
# We need coordinates for the *full* dataset (all 38,897 rows)
# Assuming 'gadm_adm2_proj' has the geometries for all GID_2s in your model data.

# Create a centroid list for *all* ADM2 units in the model
full_model_data <- adm2_analysis_data %>%
  dplyr::select(GID_2) %>%
  left_join(gadm_adm2_proj %>% 
              dplyr::select(GID_2),
            .,
            by = "GID_2") %>%
  aggregate(by = "GID_2") %>%
  centroids()

# Extract coordinates in the *exact same order* as the model data
coords <- tibble(x = crds(full_model_data)[, 1],
                 y = crds(full_model_data)[, 2],
                 GID_2 = c(full_model_data$GID_2)) %>%
  filter(GID_2 %in% adm2_analysis_data$GID_2)

# testSpatialAutocorrelation(simulationOutput = dharma_res_int,
#                            x = coords$x,
#                            y = coords$y)
# 
# testSpatialAutocorrelation(simulationOutput = dharma_res_brd,
#                            x = coords$x,
#                            y = coords$y)