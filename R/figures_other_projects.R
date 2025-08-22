pkgs <- c(
  "countrycode",
  "cowplot",
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
analysis_date <- "2025-06-04"
arha_data <- read_rds(here("data", "clean_data", paste0(analysis_date, "_data.rds")))

# Combine with rodent all pathogen WA dataset
wa_rodent_data_all <- readRDS("C:/Users/ucbtds4/R_Repositories/scoping_review/data_clean/rodent_spatial.rds")

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

mpox <- c("Angola", "Burundi", "Cameroon", "Central African Republic", "Democratic Republic of the Congo", "Gabon", "Ghana", "Guinea", "Côte d'Ivoire", "Kenya", "Liberia", "Malawi", "Morocco","Mozambique", "Nigeria", "Congo", "Rwanda", "Sierra Leone","South Africa", "Tanzania", "Zambia", "Zimbabwe")
ebola <- c("Democratic Republic of the Congo", "Gabon", "Guinea", "Liberia", "Mali", "Nigeria", "Congo", "Senegal", "Sierra Leone", "South Sudan", "Uganda")
lassa <- c("Benin", "Burkina Faso", "Ghana", "Guinea", "Liberia", "Mali", "Côte d'Ivoire", "Nigeria", "Sierra Leone", "Togo")
marburg <- c("Angola", "Democratic Republic of the Congo", "Equatorial Guinea", "Guinea", "Ghana", "Kenya", "Rwanda", "South Africa", "Uganda", "Tanzania", "Zimbabwe")

focal_viruses <- tibble(disease = c(rep("Ebola", length(ebola)),
                                    rep("Lassa", length(lassa)),
                                    rep("Marburg", length(marburg)),
                                    rep("Mpox", length(mpox))),
                        country = c(ebola, lassa, marburg, mpox)) |>
  pivot_wider(names_from = disease,
              values_from = disease) |>
  mutate(outbreak_observed = fct_rev(fct(case_when(Ebola == "Ebola" & Lassa == "Lassa" & Marburg == "Marburg" & Mpox == "Mpox" ~ "Ebola, Lassa, Marburg and Mpox",
                                                   Ebola == "Ebola" & Lassa == "Lassa" & Marburg == "Marburg" ~ "Ebola, Lassa and Marburg",
                                                   Ebola == "Ebola" & Lassa == "Lassa" & Mpox == "Mpox" ~ "Ebola, Lassa and Mpox",
                                                   Ebola == "Ebola" & Marburg == "Marburg" & Mpox == "Mpox" ~ "Ebola, Marburg and Mpox",
                                                   Ebola == "Ebola" & Mpox == "Mpox" ~ "Ebola and Mpox",
                                                   Ebola == "Ebola" & Lassa == "Lassa" ~ "Ebola and Lassa",
                                                   Lassa == "Lassa" & Marburg == "Marburg" & Mpox == "Mpox" ~ "Lassa, Marburg and Mpox",
                                                   Ebola == "Ebola" & Marburg == "Marburg" ~ "Ebola and Marburg",
                                                   Lassa == "Lassa" & Mpox == "Mpox" ~ "Lassa and Mpox",
                                                   Marburg == "Marburg" & Mpox == "Mpox" ~ "Marburg and Mpox",
                                                   Ebola == "Ebola" ~ "Ebola",
                                                   Lassa == "Lassa" ~ "Lassa",
                                                   Marburg == "Marburg" ~ "Marburg",
                                                   Mpox == "Mpox" ~ "Mpox"),
                                         levels = c("Ebola, Lassa, Marburg and Mpox",
                                                    "Ebola, Lassa and Mpox",
                                                    "Ebola, Marburg and Mpox",
                                                    "Lassa, Marburg and Mpox",
                                                    "Ebola and Lassa",
                                                    "Ebola and Marburg",
                                                    "Ebola and Mpox",
                                                    "Lassa and Mpox",
                                                    "Marburg and Mpox",
                                                    "Ebola",
                                                    "Lassa",
                                                    "Marburg",
                                                    "Mpox")))) |>
  rowwise() |>
  mutate(virus_count = factor(rowSums(across(c(Ebola, Lassa, Marburg, Mpox), ~ !is.na(.x))), levels = c("1", "2", "3", "4"))) |>
  select(country, outbreak_observed, virus_count)

focal_viruses <- africa_shapefile |>
  left_join(y = focal_viruses, by = c("NAME_0" = "country"), copy = TRUE) |>
  drop_na(outbreak_observed)

target_countries <- africa_shapefile |>
  filter(NAME_0 %in% c("Democratic Republic of the Congo", "Sierra Leone", "Tanzania"))

library(MetBrewer)
fill_palette <- setNames(met.brewer(n = 13, name = "Signac"), rev(levels(focal_viruses$outbreak_observed)))

library(magick)
sites_main <- ggplot() +
  geom_spatvector(data = africa_shapefile, fill = NA) +
  geom_spatvector(data = focal_viruses, aes(fill = outbreak_observed)) +
  geom_sf_pattern(data = focal_viruses, aes(pattern_fill = outbreak_observed, pattern_type = virus_count),
                  pattern = "polygon_tiling",
                  pattern_spacing = 0.03) +
  geom_spatvector(data = sampling_sites_africa, size = 1, aes(shape = project), colour = "black") +
  geom_spatvector(data = target_countries, fill = NA, colour = "white", lwd = 1) +
  theme_minimal() +
  scale_fill_manual(values = fill_palette) +
  scale_pattern_fill_manual(values = fill_palette) +
  scale_pattern_type_manual(values = c("square", "triangular", "hexagonal", "rhombille")) +
  guides(
    fill = guide_legend(
      override.aes = list(fill = rev(fill_palette), 
                          colour = NA)
      ),
    pattern_fill = "none",
    pattern_type = guide_legend(),
    pattern_alpha = "none",       
    colour = "none"               
  ) +
  labs(shape = "Dataset",
       fill = "Observed outbreak",
       pattern_type = "Number of distinct viruses")

bbox_sierra_leone <- ext(focal_viruses |>
                           filter(GID_0 == "SLE"))
bbox_tanzania <- ext(focal_viruses |>
                       filter(GID_0 == "TZA"))

sl_zoom_virus <- crop(focal_viruses, ext(bbox_sierra_leone))
sl_zoom_sampling <- crop(sampling_sites_africa, ext(bbox_sierra_leone))
tz_zoom_virus <- crop(focal_viruses, ext(bbox_tanzania))
tz_zoom_sampling <- crop(sampling_sites_africa, ext(bbox_tanzania))

sl_map <- ggplot() +
  geom_spatvector(data = crop(africa_shapefile, ext(bbox_sierra_leone)), fill = NA) +
  geom_sf_pattern(data = sl_zoom_virus, aes(pattern_fill = outbreak_observed, pattern_type = virus_count),
                  pattern = "polygon_tiling",
                  pattern_spacing = 0.2) +
  geom_spatvector(data = sl_zoom_sampling, size = 2, aes(shape = project), colour = "black") +
  geom_spatvector(data = target_countries |>
                    filter(GID_0 == "SLE"), fill = NA, colour = "white", lwd = 1) +
  theme_void() +
  scale_pattern_fill_manual(values = fill_palette) +
  scale_pattern_type_manual(values = c("square", "triangular", "hexagonal", "rhombille")) +
  scale_pattern_alpha_manual(values = c("1" = 0.4, "2" = 0.6, "3" = 0.8, "4" = 1)) +
  guides(pattern_fill = "none",
         pattern_type = "none",
         shape = "none",
         colour = "none") +
  labs(title = "A)")

tz_map <- ggplot() +
  geom_spatvector(data = crop(africa_shapefile, ext(bbox_tanzania)), fill = NA) +
  geom_sf_pattern(data = tz_zoom_virus, aes(pattern_fill = outbreak_observed, pattern_type = virus_count),
                  pattern = "polygon_tiling",
                  pattern_spacing = 0.2) +
  geom_spatvector(data = tz_zoom_sampling, size = 2, aes(shape = project), colour = "black") +
  geom_spatvector(data = target_countries |>
                    filter(GID_0 == "TZA"), fill = NA, colour = "white", lwd = 1) +
  theme_void() +
  scale_pattern_fill_manual(values = fill_palette) +
  scale_pattern_type_manual(values = c("square", "triangular", "hexagonal", "rhombille")) +
  scale_pattern_alpha_manual(values = c("1" = 0.4, "2" = 0.6, "3" = 0.8, "4" = 1)) +
  guides(pattern_fill = "none",
         pattern_type = "none",
         shape = "none",
         colour = "none") +
  labs(title = "B)")
  
library(cowplot)
final_map <- ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(sites_main, 0, 0, 1, 1) +  # full map
  draw_plot(sl_map, x = 0.04, y = 0.4, width = 0.18, height = 0.18) +  # Sierra Leone inset (moved inward)
  draw_plot(tz_map, x = 0.05, y = 0.25, width = 0.16, height = 0.16)    # Tanzania inset (moved inward)

save_plot(
  plot = final_map,
  filename = here("misc", "ASSET_fig.png"),
  base_height = 10,
  base_width = 9,
  dpi = 300
)

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



# ESPIDAM poster ----------------------------------------------------------
c <- country_codes()

a <- arha_data$pathogen |>
  ungroup() |>
  left_join(arha_data$host |>
              ungroup() |>
              select(rodent_record_id, iso3, iso3_composite),
            by  = c("associated_rodent_record_id" = "rodent_record_id")) |>
  mutate(iso3 = case_when(str_detect(country, "Argentina") ~ "ARG", 
                          str_detect(country, "Brazil") ~ "BRA", 
                          str_detect(country, "California U.S.|Colorado|Nevada|California|Texas|USA") ~ "USA", 
                          str_detect(country, "Czech") ~ "CZE", 
                          str_detect(country, "England|Scotland|Wales") ~ "GBR", 
                          str_detect(country, "Finland|Finald") ~ "FIN",
                          str_detect(country, "France") ~ "FRA",
                          str_detect(country, "Gemany|Germany") ~ "DEU",
                          str_detect(country, "Guinea") ~ "GIN",
                          str_detect(country, "Karelia|Russa|Russia|Siberia") ~ "RUS",
                          str_detect(country, "Korea") ~ "KOR",
                          str_detect(country, "Lithuania") ~ "LTU",
                          str_detect(country, "Nambia") ~ "NAM",
                          str_detect(country, "Nepal") ~ "NPL",
                          str_detect(country, "Phillipines") ~ "PHL",
                          str_detect(country, "Slovenia|Solvenia") ~ "SVN",
                          str_detect(country, "South Africa") ~ "ZAF",
                          str_detect(country, "Ukraine") ~ "UKR",
                          str_detect(country, "Uruaguay") ~ "URY",
                          str_detect(country, "Venezeula|Venezuala|Venezuela") ~ "VEN",
                          str_detect(country, "Zambia") ~ "ZMB",
                          str_detect(locality, "Lusaka|Livingstone") ~ "ZMB",
                          TRUE ~ iso3
                          )) |>
  mutate(
    latest_year = map_int(eventDate, ~ {
      # extract all 4‑digit runs
      yrs_chr <- str_extract_all(.x, "\\b\\d{4}\\b")[[1]]
      # turn to integer
      yrs_int <- as.integer(yrs_chr)
      # if none found, return NA, else the max year
      if (length(yrs_int) == 0L) NA_integer_ else max(yrs_int)
    })
  ) |>
  filter(latest_year <= 2025) |>
  left_join(c |>
              select(continent, ISO3),
            by = c("iso3" = "ISO3")) |>
  drop_na(continent) |>
  group_by(continent, latest_year) |>
  summarise(n_assayed = sum(n_assayed, na.rm = TRUE)) |>
  ungroup() |>
  complete(
    continent,
    latest_year = full_seq(1966:2025, 1),
    fill = list(n_assayed = 0)
  ) |>
  arrange(continent, latest_year) |>
  group_by(continent) |>
  mutate(cum_assayed = cumsum(n_assayed)) |>
  ungroup()

library(scales) 
library(RColorBrewer)

cumulative_effort <- ggplot(a, aes(x = latest_year, y = cum_assayed, colour = continent)) +
  geom_line(linewidth = 1.2) +
  
  # Better axis formatting
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    breaks = seq(1970, 2025, by = 5),
    limits = c(1965, 2025)
  ) +
  
  scale_colour_brewer(palette = "Set2") +
  
  # Add axis and legend labels
  labs(
    title = "Cumulative Assay Effort by Continent (1966–2025)",
    x = "Year",
    y = "Cumulative assays",
    colour = "Continent"
  ) +
  
  # Improve theme for posters
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

save_plot(filename = here("misc", "sampling_effort.png"), plot = cumulative_effort, base_width = 7, base_height = 5)
