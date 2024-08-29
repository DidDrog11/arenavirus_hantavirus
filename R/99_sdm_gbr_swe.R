
# Setup -------------------------------------------------------------------
source(here::here("R", "00_load_data.R"))

pkgs <- c(
  "cowplot",
  "dbarts",
  "devtools",
  "geodata",
  "here",
  "reshape",
  "rnaturalearth",
  "sf",
  "terra",
  "tidyterra",
  "tidyverse"
)

pacman::p_load(pkgs, character.only = T)

devtools::install_github('cjcarlson/embarcadero')
# then load that specific package
library(embarcadero)

project_crs <- "EPSG:4326"
europe_crs <- "EPSG:3035"


# Load data ---------------------------------------------------------------
# Ricardo's priorities
rr_priority <- read_csv(here("data", "misc", "Hanta_Hosts_Extraction.csv")) %>%
  filter(str_detect(Species, "puumal|seoul")) %>%
  group_by(Species) %>%
  group_split()
names(rr_priority) <- lapply(rr_priority, function(x) unique(x$Species))

# Read in IUCN data
iucn_match <- read_rds(here("data", "iucn_match.rds"))
iucn_ranges <- read_rds(here("data", "iucn_ranges.rds"))

# Load geographic data
world_vect <- world(path = here("data"), resolution = 3)

puumala_hosts <- rr_priority$`Orthohantavirus puumalaense` %>%
  distinct(Host) %>%
  filter(!str_detect(Host, "Homo|befordiae"))

puumala_host_ranges <- iucn_match %>%
  filter(`Host species` %in% puumala_hosts$Host) %>%
  left_join(iucn_ranges, by = c("IUCN_name" = "SCI_NAME")) %>%
  st_as_sf()

# Puumala host occurrences manually downloaded from GBIF
# GBIF.org (14 August 2024) GBIF Occurrence Download  https://doi.org/10.15468/dl.9cdyb2
puumala_host_occurrences <- read_tsv(here("data", "misc", "hanta_host_occurrences.csv")) %>%
  group_by(species)

puumala_hosts <- puumala_host_occurrences %>%
  group_split()
names(puumala_hosts) <- group_keys(puumala_host_occurrences) %>%
  pull(species)

# Can also bring in the data from ArHa
arha <- read_rds(here("data", "clean_data", "2024-08-14_data.rds"))
puumala_host_occurrence_arha <- arha$host %>%
  drop_na(species) %>%
  filter(species %in% group_keys(puumala_host_occurrences)$species) %>%
  group_by(species)
puumala_hosts_arha <- puumala_host_occurrence_arha %>%
  group_split()
names(puumala_hosts_arha) <- group_keys(puumala_host_occurrence_arha) %>%
  pull(species)

# Myodes glareolus --------------------------------------------------------

myodes_gl_iucn <- iucn_ranges %>%
  filter(SCI_NAME == "Clethrionomys glareolus") %>%
  vect()

myodes_gl_point <- bind_rows(
  puumala_hosts$`Myodes glareolus` %>%
    dplyr::select(species, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, eventDate, day, month, year) %>%
    mutate(source = "GBIF",
           individualCount = 1),
  puumala_hosts_arha$`Myodes glareolus` %>%
    dplyr::select(species, decimalLatitude, decimalLongitude, individualCount, eventDate) %>%
    mutate(source = "ArHa"))

myodes_gl_pres <- myodes_gl_point %>%
  filter(individualCount >= 1) %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = project_crs)

myodes_gl_abs <- myodes_gl_point %>%
  filter(individualCount == 0) %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = project_crs)

myodes_gl_iucn_ext <- vect(ext(myodes_gl_iucn), crs = project_crs) 
myodes_gl_area <- buffer(myodes_gl_iucn_ext, 500000)

# Produce covariate raster ------------------------------------------------
# 
# trees <- landcover("trees", path = here("data", "misc"))
# water <- landcover("water", path = here("data", "misc"))
# cropland <- landcover("cropland", path = here("data", "misc"))
# grassland <- landcover("grassland", path = here("data", "misc"))
# bioclim_data <- worldclim_global(var = c("bio"), path = here("data", "misc"), res = 0.5) %>%
#   tidyterra::select("wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_5", "wc2.1_30s_bio_6", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_bio_15")
# elevation <- elevation_global(res = 0.5, path = here("data", "misc"))
# anth_footprint <- footprint(year = 2009, path = here("data", "misc"))
# anth_travel_time <- travel_time(to = "city", size = 9, up = TRUE, path = here("data", "misc"))
# 
# raster_list <- list(trees,
#                     water,
#                     cropland,
#                     grassland,
#                     bioclim_data,
#                     elevation,
#                     anth_footprint,
#                     anth_travel_time)
# 
# crop_to_ext <- function(raster_list, host_ranges = puumala_host_ranges) {
#   c_rast <- lapply(raster_list, function(x) x %>%
#                      crop(ext(vect(host_ranges))))
#   combined_raster <- rast(c_rast)
#   return(combined_raster)
# }
# 
# puumala_covs <- crop_to_ext(raster_list, host_ranges = puumala_host_ranges)
# names(puumala_covs) <- c("trees", "water", "cropland", "grassland", "bio_1", "bio_2", "bio_5", "bio_6", "bio_12", "bio_13", "bio_14", "bio_15", "elevation", "anth_footprint", "travel_cities_1")
# Downsample
# writeRaster(puumala_covs, filename = here("data", "misc", "covs_rast2.tiff"))
# puumala_covs <- rast(here("data", "misc", "covs_rast2.tiff"))
# ds_res <- 0.05
# res_factor <- ds_res/res(puumala_covs)
# ds_puumala_covs <- terra::aggregate(puumala_covs, fact = res_factor, fun = "median", cores = 6)
# Downsampled to 17.5 million raster cells ~5km at equator
# writeRaster(ds_puumala_covs, filename = here("data", "misc", "covs_rast3.tiff"))
# puumala_covs <- rast(here("data", "misc", "covs_rast3.tiff"))


# Map occurrence data -----------------------------------------------------

myodes_gl_loc <- ggplot() +
  geom_spatvector(data = world_vect, fill = "transparent") +
  geom_spatvector(data = myodes_gl_iucn, fill = "red", alpha = 0.2) +
  geom_spatvector(data = myodes_gl_pres, aes(colour = source), alpha = 0.4) +
  geom_spatvector(data = myodes_gl_abs, aes(colour = source), alpha = 1) +
  geom_spatvector(data = myodes_gl_area, fill = "transparent")

ggsave(plot = myodes_gl_loc, filename = here("data", "misc", "p1.png"), width = 12, height = 10)

# I'm not too sure about the accuracy of the occurrences in the US. I assume the ones around the equator are also miscoded.
# Produce a buffer of ~500km around the IUCN range and limit observations to this area
# Limit the covs to this area too
puumala_covs_myodes_gl <- crop(puumala_covs, myodes_gl_area)
# Limit the occurrences to this area too
myodes_gl_pres <- crop(myodes_gl_pres, myodes_gl_area)


# Thin occurrence data ----------------------------------------------------

# Thin occurrence data to the cells of the covs raster
myodes_gl_rast <- rasterize(myodes_gl_pres, puumala_covs_myodes_gl, field = "individualCount", fun = "count", background = 0, tolerance = 0.0001)
names(myodes_gl_rast) <- "occurrence"
# Check rasterization count is similar to point count
sum(values(myodes_gl_rast))
# Convert to binary presence and filter to occurrence only
presence_myodes_gl_rast <- ifel(myodes_gl_rast >= 1, 1, myodes_gl_rast) %>%
  filter(occurrence == 1)
# Add back to raster for later
puumala_covs_myodes_gl$myodes_gl <- presence_myodes_gl_rast
# Convert back to points
presence_myodes_gl_vect <- as.points(presence_myodes_gl_rast, values = TRUE)
# Extract the values from the raster
covs_myodes_gl_df <- extract(puumala_covs_myodes_gl, presence_myodes_gl_vect, cells = TRUE, xy = TRUE, method = "simple")
# These are cells that contain occurrences
covs_myodes_gl_df$myodes_gl <- 1

# Check which contain NA for any of the covs
na_covs <- covs_myodes_gl_df %>%
  filter(if_any(everything(), is.na))
na_covs_vect <- vect(na_covs, geom = c("x", "y"), crs = project_crs)


# Figure out how to deal with NAs later -----------------------------------
# 
# # Focal window (adjust size as needed)
# w <- matrix(1, nrow = 3, ncol = 3)
# # Use focal to impute missing values from nearby cells to reduce missing
# for(var in names(puumala_covs_myodes_gl)) {
#   
#   missing_in_var <- na_covs_vect[is.na(values(na_covs_vect[var]))]
#   
#   layer <- puumala_covs_myodes_gl[var]
#   
#   focal_layer <- focal(layer, w = w, fun = median, na.rm = TRUE)
#   
#   a <- extract(focal_layer, missing_in_var, ID = TRUE, xy = TRUE, method = "simple")
#   b <- extract(focal_layer, missing_in_var, ID = TRUE, xy = TRUE, method = "bilinear")
#   
#   imputed_values <- bind_cols(a,
#                              b %>%
#                                dplyr::select(bilinear_median = focal_median)) %>%
#     mutate(value = coalesce(focal_median, bilinear_median)) %>%
#     dplyr::select(value)
#   
#   imputed_values$cell <- missing_in_var$cell
#   
#   puumala_covs_myodes_gl[var][imputed_values$cell] <- imputed_values$value
#   
# }
# 
# 
# # NAs for bioclim variables are typically in areas of 100% water. It is also unlikely that voles exist within these locations so they will be removed.
# # This will remove 239 NAs
# ID_water <- na_covs %>%
#   filter(water == 1) %>%
#   pull(ID)
# # Elevation is also missing from a lot of these
# remaining_NA <- na_covs %>%
#   filter(!ID %in% ID_water)
# 
# remaining_NA_vect <- vect(remaining_NA, geom = c("x", "y"), crs = project_crs)
# Locations of NA
ggplot() +
  geom_spatvector(data = na_covs_vect) +
  geom_spatvector(data = crop(world_vect, na_covs_vect), fill = "transparent")


# Remove NA, ideally we don't want to do this -----------------------------

covs_myodes_gl_df <- covs_myodes_gl_df %>%
  drop_na()

# Exploring variables -----------------------------------------------------

# Check the variability of these variables
variable_plots <- list()
for(var in names(covs_myodes_gl_df)[2:16])  {
  variable_plots[[var]] <- ggplot(data = covs_myodes_gl_df,
                                  aes(x = !!sym(var))) +
    geom_histogram(fill = "lightblue", color = "black", binwidth = 0.1) +
    labs(title = paste("Histogram of", var), x = var, y = "Frequency") +
    theme_minimal()
}

# Create pseudoabsences ---------------------------------------------------
# We may want to consider using non-detections from the ArHa data but for simplicity here I'm just using randomPoints from the raster
n_presence <- nrow(covs_myodes_gl_df)
pseudoabsence <- spatSample(puumala_covs_myodes_gl %>%
                              tidyterra::select(-myodes_gl), size = n_presence, method = "random", replace = FALSE, cells = TRUE, xy = TRUE, as.df = TRUE, values = TRUE, na.rm = TRUE)
pseudoabsence$myodes_gl <- 0
# Make sure none of the pseudoabsences occur in the same cells as presences
all_cov_myodes_gl <- bind_rows(covs_myodes_gl_df,
                               pseudoabsence %>%
                                 filter(!cell %in% covs_myodes_gl_df$cell)) %>%
  dplyr::select(-ID) %>%
  drop_na()
table(all_cov_myodes_gl$myodes_gl)

# Run models --------------------------------------------------------------
xvars <- names(puumala_covs_myodes_gl)[1:15]
yvar <- "myodes_gl"




# Subset to UK and Sweden -------------------------------------------------
### GBR Model -------------------------------------------------------------
GBR <- world_vect %>% filter(GID_0 == "GBR")

GBR_rast <- crop(puumala_covs_myodes_gl, GBR)
n_presence_GBR <- sum(values(GBR_rast$myodes_gl), na.rm = TRUE)
pseudoabsence_GBR_pool <- as.data.frame(GBR_rast, cells = TRUE) %>%
  drop_na(bio_1, trees) %>%
  filter(is.na(myodes_gl)) %>%
  pull(cell)
pseudoabsence_GBR_cell <- sample(pseudoabsence_GBR_pool, size = n_presence_GBR, replace = FALSE)
GBR_rast$myodes_gl[pseudoabsence_GBR_cell] <- 0

GBR_myodes_df <- as.data.frame(GBR_rast, cells = TRUE, xy = TRUE) %>%
  mutate(myodes_gl = as.integer(myodes_gl)) %>%
  drop_na()
GBR_myodes_vect <- vect(GBR_myodes_df, geom = c("x", "y"), crs = project_crs)

# Running as model.0 works, including all variables
x.data = GBR_myodes_df[, xvars]
y.data = GBR_myodes_df[yvar]
GBR_model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                         ri.data = NULL, n.trees = 200)
summary(GBR_model.0)
# Check predict works
GBR_pred.0 <- predict(GBR_model.0, raster::stack(GBR_rast %>% tidyterra::select(any_of(xvars))))
ggplot() + geom_spatraster(data = rast(GBR_pred.0, crs = project_crs))

# Now try variable selection to reduce variables based on RMSE
GBR_model_var <- variable.step2(x.data, y.data, ri.data = NULL, n.trees = 10, iter = 50, quiet = FALSE)

GBR_model_final <- bart.flex(x.data = x.data[, GBR_model_var], y.data = y.data, ri.data = NULL, n.trees = 200)
GBR_model_summary <- summary(GBR_model_final)
GBR_model_varimp <- varimp(GBR_model_final)

# Do the spatial prediction
GBR_pred.final <- predict(GBR_model_final,
                          raster::stack(GBR_rast %>% tidyterra::select(any_of(xvars))),
                          quantiles = c(0.025, 0.975))

GBR_pred_rast <- rast(GBR_pred.final, crs = project_crs)
GBR_bin_pred_rast <- GBR_pred_rast >= 0.467
GBR_bin_pred_rast <- as.numeric(GBR_bin_pred_rast)

### SWE Model -------------------------------------------------------------
SWE <- world_vect %>% filter(GID_0 == "SWE")

SWE_rast <- crop(puumala_covs_myodes_gl, SWE)
n_presence_SWE <- sum(values(SWE_rast$myodes_gl), na.rm = TRUE)
pseudoabsence_SWE_pool <- as.data.frame(SWE_rast, cells = TRUE) %>%
  drop_na(bio_1, trees) %>%
  filter(is.na(myodes_gl)) %>%
  pull(cell)
pseudoabsence_SWE_cell <- sample(pseudoabsence_SWE_pool, size = n_presence_SWE, replace = FALSE)
SWE_rast$myodes_gl[pseudoabsence_SWE_cell] <- 0

SWE_myodes_df <- as.data.frame(SWE_rast, cells = TRUE, xy = TRUE) %>%
  mutate(myodes_gl = as.integer(myodes_gl)) %>%
  drop_na()
SWE_myodes_vect <- vect(SWE_myodes_df, geom = c("x", "y"), crs = project_crs)

# Running as model.0 works, including all variables
x.data = SWE_myodes_df[, xvars]
y.data = SWE_myodes_df[yvar]
SWE_model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                         ri.data = NULL, n.trees = 200)
summary(SWE_model.0)
# Check predict works
SWE_pred.0 <- predict(SWE_model.0, raster::stack(SWE_rast %>% tidyterra::select(any_of(xvars))))
ggplot() + geom_spatraster(data = rast(SWE_pred.0, crs = project_crs))

# Now try variable selection to reduce variables based on RMSE
SWE_model_var <- variable.step2(x.data, y.data, ri.data = NULL, n.trees = 10, iter = 50, quiet = FALSE)

SWE_model_final <- bart.flex(x.data = x.data[, SWE_model_var], y.data = y.data, ri.data = NULL, n.trees = 200)
SWE_model_summary <- summary(SWE_model_final)
SWE_model_varimp <- varimp(SWE_model_final)

# Do the spatial prediction
SWE_pred.final <- predict(SWE_model_final,
                          raster::stack(SWE_rast %>% tidyterra::select(any_of(xvars))),
                          quantiles = c(0.025, 0.975))

SWE_pred_rast <- rast(SWE_pred.final, crs = project_crs)
SWE_bin_pred_rast <- SWE_pred_rast >= 0.424
SWE_bin_pred_rast <- as.numeric(SWE_bin_pred_rast)

## Explore the uncertainty in these predictions ----------------------------
GBR_central_estimate <- ggplot() +
  geom_spatraster(data = GBR_pred_rast[[1]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

GBR_lower_conf <- ggplot() +
  geom_spatraster(data = GBR_pred_rast[[2]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

GBR_upper_conf <- ggplot() +
  geom_spatraster(data = GBR_pred_rast[[3]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

plot_grid(GBR_lower_conf, GBR_central_estimate, GBR_upper_conf, nrow = 1)

SWE_central_estimate <- ggplot() +
  geom_spatraster(data = SWE_pred_rast[[1]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

SWE_lower_conf <- ggplot() +
  geom_spatraster(data = SWE_pred_rast[[2]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

SWE_upper_conf <- ggplot() +
  geom_spatraster(data = SWE_pred_rast[[3]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

plot_grid(SWE_lower_conf, SWE_central_estimate, SWE_upper_conf, nrow = 1)

## Explore uncertainty in binarised predictions ----------------------------
GBR_bin_central <- ggplot() +
  geom_spatraster(data = GBR_bin_pred_rast[[1]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  geom_spatvector(data = GBR_myodes_vect %>%
                    tidyterra::select(myodes_gl) %>%
                    filter(myodes_gl == 1),
                  colour = "orange",
                  size = 0.2) +
  theme_minimal()

GBR_bin_lower <- ggplot() +
  geom_spatraster(data = GBR_bin_pred_rast[[2]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

GBR_bin_upper <- ggplot() +
  geom_spatraster(data = GBR_bin_pred_rast[[3]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

plot_grid(GBR_bin_lower, GBR_bin_central, GBR_bin_upper, nrow = 1)

SWE_bin_central <- ggplot() +
  geom_spatraster(data = SWE_bin_pred_rast[[1]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  geom_spatvector(data = SWE_myodes_vect %>%
                    tidyterra::select(myodes_gl) %>%
                    filter(myodes_gl == 1),
                  colour = "orange",
                  size = 0.2) +
  theme_minimal()

SWE_bin_lower <- ggplot() +
  geom_spatraster(data = SWE_bin_pred_rast[[2]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

SWE_bin_upper <- ggplot() +
  geom_spatraster(data = SWE_bin_pred_rast[[3]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

plot_grid(SWE_bin_lower, SWE_bin_central, SWE_bin_upper, nrow = 1)