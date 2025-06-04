pkgs <- c(
  "countrycode",
  "geodata",
  "ggpattern",
  "googledrive",
  "googlesheets4",
  "here",
  "readxl",
  "magick",
  "stringdist",
  "terra",
  "taxize",
  "tidyverse",
  "tidyterra",
  "colorspace",
  "htmltools",
  "leaflet",
  "plotly",
  "sf"
)

pacman::p_load(pkgs, character.only = T)

project_crs <- "EPSG:4326"

world_shapefile <- world(path = here("data"))

# Change to most recent
analysis_date <- "2025-05-20"
arha_data <- read_rds(here("data", "clean_data", paste0(analysis_date, "_data.rds")))

# Combine with rodent all pathogen WA dataset
wa_rodent_data_all <- readRDS("~/Downloads/rodent_spatial.rds")

wa_rodent_data <- as_tibble(wa_rodent_data_all) |>
  mutate(decimalLongitude = st_coordinates(wa_rodent_data_all)[ , 1],
         decimalLatitude = st_coordinates(wa_rodent_data_all)[ , 2],
         country = as.character(country)) |>
  select(rodent_record_id = record_id, country, decimalLongitude, decimalLatitude) |>
  distinct(rodent_record_id, country, decimalLongitude, decimalLatitude, .keep_all = TRUE) |>
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = project_crs) |>
  mutate(project = "WA rodent pathogen")

africa <- codelist |>
  filter(continent == "Africa")

africa_shapefile <- world_shapefile |>
  filter(GID_0 %in% africa$iso3c)

sampling_sites <- arha_data$host |>
  select(rodent_record_id, decimalLatitude, decimalLongitude) |>
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE) |>
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = project_crs) |>
  mutate(project = "ArHa")

combined_sampling_sites <- rbind(wa_rodent_data, sampling_sites) |>
  as_tibble(geom = c("XY")) |>
  distinct(x, y, .keep_all = TRUE) |>
  select(rodent_record_id, country, project, x, y) |>
  vect(geom = c("x", "y"), crs = project_crs)

sampling_sites_africa <- terra::intersect(combined_sampling_sites, africa_shapefile)

mpox <- c("Angola", "Burundi", "Cameroon", "Central African Republic", "Democratic Republic of the Congo", "Gabon", "Ghana", "Guinea", "Ivory Coast", "Kenya", "Liberia", "Malawi", "Morocco","Mozambique", "Nigeria", "Republic of the Congo", "Rwanda", "Sierra Lene","South Africa", "Tanzania", "Zambia", "Zimbabwe")
ebola <- c("Democratic Republic of the Congo", "Gabon", "Guinea", "Liberia", "Mali", "Nigeria", "Republic of the Congo", "Senegal", "Sierra Leone", "South Sudan", "Uganda")
lassa <- c("Benin", "Burkina Faso", "Ghana", "Guinea", "Liberia", "Mali", "Ivory Coast", "Nigeria", "Sierra Leone", "Togo")
marburg <- c("Angola", "Equatorial Guinea", "Guinea", "Ghana", "Kenya", "Rwanda", "South Africa", "Uganda", "Tanzania", "Zimbabwe")

focal_viruses <- tibble(disease = c(rep("Ebola", length(ebola)),
                                    rep("Lassa", length(lassa)),
                                    rep("Marburg", length(marburg)),
                                    rep("Mpox", length(mpox))),
                        country = c(ebola, lassa, marburg, mpox)) |>
  pivot_wider(names_from = disease,
              values_from = disease) |>
  mutate(outbreak_observed = fct_rev(fct(case_when(Ebola == "Ebola" & Lassa == "Lassa" & Marburg == "Marburg" & Mpox == "Mpox" ~ "Ebola, Lassa, Marburg and Mpox",
                                                   Ebola == "Ebola" & Lassa == "Lassa" & Marburg == "Marburg" ~ "Ebola, Lassa and Marburg",
                                                   Ebola == "Ebola" & Mpox == "Mpox" ~ "Ebola and Mpox",
                                                   Ebola == "Ebola" & Lassa == "Lassa" ~ "Ebola and Lassa",
                                                   Lassa == "Lassa" & Marburg == "Marburg" & Mpox == "Mpox" ~ "Lassa, Marburg and Mpox",
                                                   Lassa == "Lassa" & Mpox == "Mpox" ~ "Lassa and Mpox",
                                                   Marburg == "Marburg" & Mpox == "Mpox" ~ "Marburg",
                                                   Ebola == "Ebola" ~ "Ebola",
                                                   Lassa == "Lassa" ~ "Lassa",
                                                   Marburg == "Marburg" ~ "Marburg",
                                                   Mpox == "Mpox" ~ "Mpox"),
                                         levels = c("Ebola, Lassa, Marburg and Mpox",
                                                    "Lassa, Marburg and Mpox",
                                                    "Ebola and Lassa",
                                                    "Ebola and Mpox",
                                                    "Lassa and Mpox",
                                                    "Ebola",
                                                    "Lassa",
                                                    "Marburg",
                                                    "Mpox")))) |>
  rowwise() |>
  mutate(virus_count = fct(as.character(sum(!is.na(Ebola), !is.na(Lassa), !is.na(Marburg), !is.na(Mpox))), levels = c("1", "2", "3", "4")),
         virus_count_name = case_when(virus_count == "1" ~ "A",
                                      virus_count == "2" ~ "B",
                                      virus_count == "3" ~ "C",
                                      virus_count == "4" ~ "D")) |>
  select(country, outbreak_observed, virus_count, virus_count_name) %>%
  left_join(x = africa_shapefile, y = ., by = c("NAME_0" = "country"), copy = TRUE) |>
  drop_na(outbreak_observed)

target_countries <- africa_shapefile |>
  filter(NAME_0 %in% c("Democratic Republic of the Congo", "Sierra Leone", "Tanzania"))

sites <- ggplot() +
  geom_spatvector(data = africa_shapefile) +
  geom_spatvector(data = focal_viruses, aes(fill = outbreak_observed)) +
  geom_spatvector(data = target_countries, fill = NA, colour = "white", lwd = 1) +
  geom_spatvector(data = sampling_sites_africa, size = 0.6, aes(colour = project, shape = project)) +
  theme_minimal() +
  scale_colour_manual(values = c("red", "cyan")) +
  labs(colour = "Dataset",
       shape = "Dataset",
       fill = "Observed outbreak",
       caption = "Countries with associatied outbreaks may include imported cases.")


# ArHa data ---------------------------------------------------------------

arha_pathogen <- arha_data$pathogen

jittered_points <- arha_pathogen |>
  drop_na(decimalLongitude, decimalLatitude) |>
  mutate(
    decimalLongitude = jitter(decimalLongitude, amount = 0.001),
    decimalLatitude = jitter(decimalLatitude, amount = 0.001),
    positive_binary = case_when(n_assayed >= 1 & n_positive >= 1 ~ TRUE,
                                n_assayed >= 1 & n_positive == 0 ~ FALSE,
                                TRUE ~ NA),
    log_tested = 1 + log(n_assayed)) |>
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = project_crs) |>
  select(virus_clean, taxonomic_level, assay_clean, n_assayed, n_positive, host_species, host_genus, eventDate, locality, coordinate_resolution, positive_binary, log_tested)

# Host
myodes_glareolus <- mapview::mapview(jittered_points |>
                                       filter(str_detect(host_species, "Myodes glareolus")),
                                     zcol = "virus_clean",
                                     cex = "log_tested",
                                     layer.name = "Myodes glareolus")

mastomys_natalensis <- mapview::mapview(jittered_points |>
                                          filter(str_detect(host_species, "Mastomys natalensis")),
                                        zcol = "virus_clean",
                                        cex = "log_tested",
                                        layer.name = "Mastomys natalensis")

# Pathogen
# Puumala
puumala_map <- mapview::mapview(jittered_points |>
                                  filter(str_detect(virus_clean, "puum")),
                                zcol = "positive_binary",
                                cex = "log_tested",
                                layer.name = "Puumala")

# Dobrova
dobrova_map <- mapview::mapview(jittered_points |>
                                  filter(str_detect(virus_clean, "dobr")),
                                zcol = "positive_binary",
                                cex = "log_tested",
                                layer.name = "Dobrova")

# Lassa
lassa_map <- mapview::mapview(jittered_points |>
                                  filter(str_detect(virus_clean, "lassa")),
                                zcol = "positive_binary",
                                cex = "log_tested",
                                layer.name = "Lassa")
