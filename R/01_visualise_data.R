source(here::here("R", "00_load_data.R"))

pkgs <- c(
  "colorspace",
  "htmltools",
  "leaflet",
  "plotly",
  "sf"
)

pacman::p_load(pkgs, character.only = T)

project_crs <- "EPSG:4326"


# Rodents -----------------------------------------------------------------

rodent_locations <- combined_data$host %>%
  select(genus, scientificName, decimalLatitude, decimalLongitude, individualCount) %>%
  mutate(number = as.numeric(individualCount)) %>%
  drop_na(number, decimalLongitude, decimalLatitude) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = project_crs) %>%
  group_by(geometry, scientificName) %>%
  mutate(number = sum(number, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct() %>%
  rowwise() %>%
  mutate(labels = HTML(paste0("Species: ", scientificName, "<br>", "Number detected: ", number)))

pal <- colorFactor(diverging_hcl(n = length(sort(unique(rodent_locations$genus))), palette = "Berlin"),
                   domain = sort(unique(rodent_locations$genus)))

leaflet(data = rodent_locations) %>%
  addTiles() %>%
  addCircleMarkers(radius = ~ifelse(sqrt(number) == 0, 1, sqrt(number)),
                   color = ~pal(genus),
                   stroke = FALSE,
                   fillOpacity = 0.8,
                   clusterOptions = markerClusterOptions(spiderfyOnMaxZoom = TRUE,
                                                         spiderLegPolylineOptions = list(length = 6),
                                                         spiderfyDistanceMultiplier = 2),
                   label = ~labels,
                   popup = ~labels) %>%
  addLegend("bottomright", pal = pal, values = ~genus, title = "Host genus", opacity = 0.8)


# Pathogens ---------------------------------------------------------------

pathogen_locations <- combined_data$pathogen %>%
  select(identificationRemarks, family, scientificName, associatedTaxa, decimalLatitude, decimalLongitude, occurrenceRemarks, organismQuantity) %>%
  mutate(number_tested = as.numeric(occurrenceRemarks),
         number_positive = as.numeric(organismQuantity)) %>%
  drop_na(number_tested, number_positive) %>%
  drop_na(decimalLongitude, decimalLatitude) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = project_crs) %>%
  group_by(geometry, identificationRemarks, family, scientificName, associatedTaxa) %>%
  mutate(number_tested = sum(number_tested),
         number_positive = sum(number_positive)) %>%
  ungroup() %>%
  distinct() %>%
  mutate(percent_positive = case_when(is.nan(number_positive/number_tested * 100) ~ "NA",
                                      number_positive == 0 ~ "0%",
                                      number_positive/number_tested * 100 < 1 ~ "< 1%",
                                      TRUE ~ paste0(as.character(round(number_positive/number_tested * 100, 2)), "%")),
         detection_method = case_when(identificationRemarks == "pcr" ~ "Viral",
                                      identificationRemarks == "antibody" ~ "Antibody"),
         detection = case_when(number_positive >= 1 ~ "Detected",
                               TRUE ~ "Not detected")) %>%
  rowwise() %>%
  mutate(labels = HTML(paste0("Pathogen family: ", family, "<br>",
                              "Pathogen species: ", scientificName, "<br>",
                              "Host species: ", associatedTaxa, "<br>",
                              "Number tested: ", number_tested, "<br>",
                              "Number positive (%): ", number_positive, " (", percent_positive, ")", "<br>")))

path_pal <- colorFactor(diverging_hcl(n = length(unique(pathogen_locations$scientificName)), palette = "Berlin"),
                        domain = sort(unique(pathogen_locations$scientificName)))
stroke_pal <- colorFactor(palette = c("white", "gray"), levels = c("Detected", "Not detected"))

leaflet(data = pathogen_locations) %>%
  addTiles() %>%
  addCircleMarkers(data = pathogen_locations %>%
                     filter(detection_method == "Viral") %>%
                     filter(number_tested != 0),
                   radius = ~ifelse(sqrt(number_tested) == 0, 1, sqrt(number_tested)),
                   fill = TRUE,
                   fillColor = ~path_pal(scientificName),
                   stroke = TRUE,
                   color = ~stroke_pal(detection),
                   weight = 2,
                   opacity = 1,
                   fillOpacity = 0.8,
                   clusterOptions = markerClusterOptions(spiderfyOnMaxZoom = TRUE,
                                                         spiderLegPolylineOptions = list(length = 6),
                                                         spiderfyDistanceMultiplier = 2),
                   label = ~labels,
                   popup = ~labels,
                   group = "Viral detection") %>%
  addCircleMarkers(data = pathogen_locations %>%
                     filter(detection_method == "Antibody") %>%
                     filter(number_tested != 0),
                   radius = ~ifelse(sqrt(number_tested) == 0, 1, sqrt(number_tested)),
                   fill = TRUE,
                   fillColor = ~path_pal(scientificName),
                   stroke = TRUE,
                   color = ~stroke_pal(detection),
                   weight = 2,
                   opacity = 1,
                   fillOpacity = 0.8,
                   clusterOptions = markerClusterOptions(spiderfyOnMaxZoom = TRUE,
                                                         spiderLegPolylineOptions = list(length = 6),
                                                         spiderfyDistanceMultiplier = 2),
                   label = ~labels,
                   popup = ~labels,
                   group = "Antibody detection") %>%
  addLayersControl(overlayGroups = c("Viral detection", "Antibody detection"),
                   options = layersControlOptions(collapsed = FALSE)) %>%
  addLegend("bottomright", pal = path_pal, values = ~scientificName, title = "Pathogen species", opacity = 0.8)

# Host-pathogen associations ----------------------------------------------

pathogen_summary <- combined_data$pathogen %>%
  select(genus = host_genus,
         species = host_species,
         pathogen_family, pathogen_species,
         number_tested, number_positive,
         detection_method) %>%
  mutate(number_tested = as.numeric(number_tested)) %>%
  filter(number_tested >= 1) %>%
  left_join(combined_data$zoonoses %>%
              select(pathogen_species, known_zoonosis),
            by = c("pathogen_species"))

hp_associations <- combined_data$rodent %>%
  select(genus, species) %>%
  left_join(pathogen_summary %>%
              group_by(detection_method, species, pathogen_species, known_zoonosis) %>%
              mutate(number_tested = sum(number_tested),
                     number_positive = sum(number_positive)), by = c("genus", "species")) %>%
  drop_na(number_tested) %>%
  distinct() %>%
  filter(!pathogen_species %in% c("Multiple", "Puumala virus/Saaremaa virus")) %>%
  arrange(pathogen_family, pathogen_species, genus, species) %>%
  mutate(association = case_when(number_positive >= 1 ~ TRUE,
                                 TRUE ~ FALSE),
         detection_method = case_when(detection_method == "pcr" ~ "Viral",
                                      detection_method == "antibody" ~ "Antibody"))

scaled_n <- preProcess(as.data.frame(log10(hp_associations$number_tested)), method = "range")
hp_associations$alpha <- predict(scaled_n, as.data.frame(log10(hp_associations$number_tested)))[, 1]

host_pathogen_plot <- hp_associations %>%
  ggplot(aes(x = pathogen_species, y = species)) +
  geom_tile(aes(fill = association, colour = known_zoonosis, alpha = alpha,
                text = paste0("Pathogen species: ", pathogen_species, "<br>",
                              "Host species: ", species, "<br>",
                              "Number tested: ", number_tested, "<br>",
                              "Number positive: ", number_positive, "<br>",
                              "Known zoonosis: ", known_zoonosis))) +
  scale_fill_viridis_d() +
  scale_colour_manual(values = c("white", "black")) +
  facet_wrap(~ detection_method) +
  guides(alpha = "none",
         colour = "none") +
  labs(fill = "H-P association",
       x = element_blank(),
       y = element_blank()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplotly(host_pathogen_plot, tooltip = "text")
