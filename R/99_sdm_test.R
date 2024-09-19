
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
# Convert to equal area projection for myodes gl analysis in Europe
# covs_myodes_gl_crop <- crop(puumala_covs, myodes_gl_area)
# # Calculate the projected extent of your area of interest in the new CRS
# projected_extent <- ext(project(myodes_gl_area, europe_crs))
# # Create an empty template raster with 1km resolution
# template_raster <- rast(
#   ext = projected_extent,     # Use the projected extent
#   resolution = 1000,          # 1 km resolution
#   crs = europe_crs            # Target CRS
# )
# covs_myodes_gl <- project(covs_myodes_gl_crop, template_raster, method = "bilinear", res = 1000)
# writeRaster(covs_myodes_gl, filename = here("data", "misc", "covs_rast4.tiff"))
covs_myodes_gl <- rast(here("data", "misc", "covs_rast4.tiff"))
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
# Limit the occurrences to this area too
myodes_gl_area <- project(myodes_gl_area, europe_crs)
myodes_gl_pres <- project(myodes_gl_pres, europe_crs)
myodes_gl_pres <- crop(myodes_gl_pres, myodes_gl_area)


# Thin occurrence data ----------------------------------------------------

# Thin occurrence data to the cells of the covs raster
myodes_gl_rast <- rasterize(myodes_gl_pres, covs_myodes_gl, field = "individualCount", fun = "count", background = 0, tolerance = 0.0001)
names(myodes_gl_rast) <- "occurrence"
# Check rasterization count is similar to point count
sum(values(myodes_gl_rast))
# Convert to binary presence and filter to occurrence only
presence_myodes_gl_rast <- ifel(myodes_gl_rast >= 1, 1, myodes_gl_rast) %>%
  filter(occurrence == 1)
# Add back to raster for later
covs_myodes_gl$myodes_gl <- presence_myodes_gl_rast
# Convert back to points
presence_myodes_gl_vect <- as.points(presence_myodes_gl_rast, values = TRUE)
# Extract the values from the raster
covs_myodes_gl_df <- extract(covs_myodes_gl, presence_myodes_gl_vect, cells = TRUE, xy = TRUE, method = "simple")
# These are cells that contain occurrences
covs_myodes_gl_df$myodes_gl <- 1

# Check which contain NA for any of the covs
na_covs <- covs_myodes_gl_df %>%
  filter(if_any(everything(), is.na))
na_covs_vect <- vect(na_covs, geom = c("x", "y"), crs = europe_crs)

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
  geom_spatvector(data = crop(project(world_vect, europe_crs), na_covs_vect), fill = "transparent") +
  geom_spatvector(data = na_covs_vect, colour = "black", size = 1)
  
  
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
pseudoabsence <- spatSample(covs_myodes_gl %>%
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
xvars <- names(covs_myodes_gl)[1:15]
yvar <- "myodes_gl"


# Full data model ---------------------------------------------------------
# Running as model.0 including all variables
x.data = all_cov_myodes_gl[, xvars]
y.data = all_cov_myodes_gl[yvar]
model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                     ri.data = NULL, n.trees = 200)
summary(model.0)
# Check predict works
# Predict to a reasonable resolution raster
temp_proj <- project(x = covs_myodes_gl, y = "EPSG:3857")

# 10km cells
temp_proj_lr <- aggregate(temp_proj, fact = 10, fun = "median", na.rm = TRUE, cores = 4)
new_data <- project(x = temp_proj_lr, y = "EPSG:3035")
new_data <- raster::stack(new_data %>% tidyterra::select(any_of(xvars)))

pred.0 <- predict(model.0, new_data, splitby = 50, quiet = FALSE)
ggplot() + geom_spatraster(data = rast(pred.0, crs = europe_crs))

# Now try variable selection to reduce variables based on RMSE
model_var <- variable.step2(x.data, y.data, ri.data = NULL, n.trees = 10, iter = 50, quiet = FALSE)

model_final <- bart.flex(x.data = x.data[, model_var], y.data = y.data, ri.data = NULL, n.trees = 200)
model_summary <- summary(model_final)
model_varimp <- varimp(model_final)

# Do the spatial prediction
pred_final <- predict(model_final,
                      new_data,
                      quantiles = c(0.025, 0.975), 
                      splitby = 50, 
                      quiet = FALSE)
crs(pred_final) <- europe_crs

central_estimate <- ggplot() +
  geom_spatraster(data = pred_final[[1]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

overlay_occurrence <- ggplot() +
  geom_spatraster(data = pred_final[[1]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  geom_spatvector(data = myodes_gl_pres, colour = "orange", size = 0.2) +
  theme_minimal()

lower_conf <- ggplot() +
  geom_spatraster(data = pred_final[[2]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

upper_conf <- ggplot() +
  geom_spatraster(data = pred_final[[3]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

plot_grid(lower_conf, central_estimate, upper_conf, nrow = 1)
plot_grid(central_estimate, overlay_occurrence, nrow = 1)

bin_pred_rast <- pred_final >= 0.4601
bin_pred_rast <- as.numeric(bin_pred_rast)

bin_central <- ggplot() +
  geom_spatraster(data = bin_pred_rast[[1]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

overlay_occurrence_bin <- ggplot() +
  geom_spatraster(data = bin_pred_rast[[1]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  geom_spatvector(data = myodes_gl_pres, colour = "orange", size = 0.2) +
  theme_minimal()

bin_lower <- ggplot() +
  geom_spatraster(data = bin_pred_rast[[2]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

bin_upper <- ggplot() +
  geom_spatraster(data = bin_pred_rast[[3]], na.rm = TRUE) + 
  scale_fill_viridis_c(na.value = NA, limits = c(0, 1)) +
  theme_minimal()

plot_grid(bin_lower, bin_central, bin_upper, nrow = 1)
plot_grid(bin_central, overlay_occurrence_bin, nrow = 1)

write_rds(model_final, here("data", "misc", "myodes_gl_final.rds"))

# Updated bart.step functions ---------------------------------------------
variable.step2 <- function (x.data, y.data, ri.data = NULL, n.trees = 10, iter = 50, 
                            quiet = FALSE) 
{
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }
  comp <- complete.cases(x.data)
  if (length(comp) < (nrow(x.data))) {
    message("Some rows with NA's have been automatically dropped. \n")
  }
  x.data <- x.data[comp, ]
  y.data <- y.data[comp, ]
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               ri.data = ri.data, n.trees = 200))
  if (class(model.0) == "rbart") {
    fitobj <- model.0$fit[[1]]
  }
  if (class(model.0) == "bart") {
    fitobj <- model.0$fit
  }
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x, 
                                                                                "drop")) == FALSE)))]
  if (length(dropnames) > 0) {
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames, collapse = " "), " \n")
  }
  x.data <- x.data %>% dplyr::select(-any_of(dropnames))
  nvars <- ncol(x.data)
  varnums <- c(1:nvars)
  varlist.orig <- varlist <- colnames(x.data)
  rmses <- data.frame(Variable.number = c(), RMSE = c())
  dropped.varlist <- c()
  for (var.j in c(nvars:3)) {
    print(noquote(paste("Number of variables included:", 
                        var.j)))
    print(noquote("Dropped:"))
    print(if (length(dropped.varlist) == 0) {
      noquote("")
    }
    else {
      noquote(dropped.varlist)
    })
    rmse.list <- c()
    if (!quiet) {
      pb <- txtProgressBar(min = 0, max = iter, style = 3)
    }
    for (index in 1:iter) {
      quietly(model.j <- bart.flex(x.data = x.data[, varnums], 
                                   y.data = y.data, ri.data = ri.data, n.trees = n.trees))
      quietly(vi.j <- varimp(model.j))
      if (index == 1) {
        vi.j.df <- vi.j
      }
      else {
        vi.j.df[, index + 1] <- vi.j[, 2]
      }
      pred.p <- colMeans(pnorm(model.j$yhat.train))[y.data == 
                                                      1]
      pred.a <- colMeans(pnorm(model.j$yhat.train))[y.data == 
                                                      0]
      pred.c <- c(pred.p, pred.a)
      true.c <- c(rep(1, length(pred.p)), rep(0, length(pred.a)))
      rmsej.i <- Metrics::rmse(true.c, pred.c)
      rmse.list <- c(rmse.list, rmsej.i)
      if (!quiet) {
        setTxtProgressBar(pb, index)
      }
    }
    vi.j <- data.frame(vi.j.df[, 1], rowMeans(vi.j.df[, 
                                                      -1]))
    vi.j <- vi.j[order(vi.j[, 2]), ]
    drop.var <- vi.j[1, 1]
    dropped.varlist <- c(dropped.varlist, as.character(drop.var))
    rmsej <- mean(rmse.list)
    rmses <- rbind(rmses, c(nvars - var.j, rmsej))
    colnames(rmses) <- c("VarsDropped", "RMSE")
    varnums <- varnums[!(varnums == which(varlist.orig == 
                                            drop.var))]
    varlist <- varlist.orig[varnums]
    print(noquote("---------------------------------------"))
  }
  g1 <- ggplot2::ggplot(rmses, aes(y = RMSE, x = VarsDropped)) + 
    geom_line(color = "black") + geom_point(size = 3) + 
    theme_bw() + ylab("RMSE of model\n") + xlab("\nVariables dropped") + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, 
                                                                         face = "bold")) + scale_x_discrete(limits = c(0:(nrow(rmses))))
  print(g1)
  print(noquote("---------------------------------------"))
  print(noquote("Final recommended variable list"))
  varlist.final <- varlist.orig[!(varlist.orig %in% dropped.varlist[0:(which(rmses$RMSE == 
                                                                               min(rmses$RMSE)) - 1)])]
  print(noquote(varlist.final))
  invisible(varlist.final)
}

bart.step2 <- function (x.data, y.data, ri.data = NULL, iter.step = 100, tree.step = 10, 
                        iter.plot = 100, full = FALSE, quiet = FALSE) 
{
  quietly <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }
  quietly(model.0 <- bart.flex(x.data = x.data, y.data = y.data, 
                               ri.data = ri.data, n.trees = 200))
  if (class(model.0) == "rbart") {
    fitobj <- model.0$fit[[1]]
  }
  if (class(model.0) == "bart") {
    fitobj <- model.0$fit
  }
  dropnames <- colnames(x.data)[!(colnames(x.data) %in% names(which(unlist(attr(fitobj$data@x, 
                                                                                "drop")) == FALSE)))]
  if (length(dropnames) == 0) {
  }
  else {
    message("Some of your variables have been automatically dropped by dbarts.")
    message("(This could be because they're characters, homogenous, etc.)")
    message("It is strongly recommended that you remove these from the raw data:")
    message(paste(dropnames, collapse = " "), " \n")
  }
  x.data <- x.data %>% dplyr::select(-any_of(dropnames))
  quiet2 <- quiet
  if (full == TRUE) {
    varimp.diag(x.data, y.data, ri.data, iter = iter.plot, 
                quiet = quiet2)
  }
  vs <- variable.step2(x.data, y.data, ri.data, n.trees = tree.step, 
                       iter = iter.step, quiet = quiet2)
  invisible(best.model <- bart.flex(x.data = x.data[, vs], 
                                    y.data = y.data, ri.data = ri.data, n.trees = 200))
  if (full == TRUE) {
    varimp(best.model, plots = TRUE)
  }
  if (full == TRUE) {
    p <- summary(best.model, plots = TRUE)
    print(p)
  }
  else {
    p <- summary(best.model, plots = FALSE)
    print(p)
  }
  invisible(best.model)
}


