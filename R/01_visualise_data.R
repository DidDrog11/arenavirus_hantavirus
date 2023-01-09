source(here("R", "00_load_data.R"))

pkgs <- c(
  "colorspace",
  "htmltools",
  "leaflet",
  "sf"
)

pacman::p_load(pkgs, character.only = T)

project_crs <- "EPSG:4326"

rodent_locations <- combined_data$rodent %>%
  select(species, latitude, longitude, number) %>%
  mutate(number = as.numeric(number)) %>%
  drop_na(number) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = project_crs) %>%
  group_by(geometry, species) %>%
  mutate(number = sum(number)) %>%
  ungroup() %>%
  distinct() %>%
  rowwise() %>%
  mutate(labels = HTML(paste0("Species: ", species, "<br>", "Number detected: ", number)))

pal <- colorFactor(diverging_hcl(n = length(sort(unique(rodent_locations$species))), palette = "Berlin"), domain = sort(unique(rodent_locations$species)))

leaflet(data = rodent_locations) %>%
  addTiles() %>%
  addCircleMarkers(radius = ~ifelse(sqrt(number) == 0, 1, sqrt(number)),
                   color = ~pal(species),
                   stroke = FALSE,
                   fillOpacity = 0.8,
                   clusterOptions = markerClusterOptions(spiderfyOnMaxZoom = TRUE,
                                                         spiderLegPolylineOptions = list(length = 6)),
                   label = ~labels,
                   popup = ~labels)

                   