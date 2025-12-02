# ---
# 01_taxonomic_bias_analysis.R
#
# Purpose: To quantify taxonomic biases in the sampling of hosts and pathogens.
# This script generates summary tables of sampling effort, compares sampled
# species to a comprehensive checklist, and fits a GAM to test for predictors
# of host sampling intensity.
# ---

# --- 1. Setup ---
# Load necessary packages
library(here)
library(countrycode)
library(data.table)
library(ggplot2)
library(patchwork)
library(readr)
library(readxl)
library(dplyr)
library(janitor)
library(stringr)
library(taxize)
library(terra)
library(tidyterra)
library(tidyverse)
library(traitdata)
library(mgcv) # For GAMs

# Update external data
update_external = FALSE

# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-09-25.rds") # <-- Update this date
arha_db <- read_rds(db_path)
host_data <- arha_db$host
pathogen_data <- arha_db$pathogen

# --- 2. Load Additional Data ---
# a) Comprehensive Host Checklist
try({
  if(update_external == TRUE) {
    mdd_github_url <- "https://raw.githubusercontent.com/mammaldiversity/mammaldiversity.github.io/master/_data/mdd.csv"
    mdd_checklist <- readr::read_csv(mdd_github_url)
    comprehensive_checklist <- mdd_checklist %>%
      clean_names() %>%
      dplyr::filter(order %in% c("Rodentia", "Eulipotyphla")) %>%
      mutate(
        # Clean scientific name by replacing underscores
        sci_name = str_replace_all(sci_name, "_", " "),
        # Parse the country distribution string, convert to ISO3c codes, and collapse back to a string
        iso3c_distribution = purrr::map_chr(str_split(country_distribution, "\\|"), ~{
          # The `countrycode` function handles the conversion
          codes <- countrycode(.x, "country.name", "iso3c", nomatch = NA)
          # Collapse the non-missing codes into a single, comma-separated string
          paste(stats::na.omit(codes), collapse = ", ")
        }),
        
        # Convert the biogeographic realm to a factor
        biogeographic_realm = forcats::as_factor(biogeographic_realm)
      ) %>%
      # Select and rename the final, clean columns
      select(
        scientific_name = sci_name,
        order,
        family,
        genus,
        iucn_status,
        country_distribution,
        iso3c_distribution,
        biogeographic_realm
      )
    
    write_rds(comprehensive_checklist, here("data", "external", "mdd_checklist.rds"))
    
  } else {
    
    comprehensive_checklist <- read_rds(here("data", "external", "mdd_checklist.rds"))
    
  }
}, silent = TRUE)

# b) Host Trait Data (EltonTraits)
try({
  if(update_external == TRUE) {
    data("elton_mammals")
    
    elton_traits <- elton_mammals %>%
      clean_names() %>%
      filter(scientific_name_std %in% comprehensive_checklist$scientific_name) %>%
      select(
        scientific_name = scientific_name_std,
        adult_body_mass_g = body_mass_value,
        diet_invert = diet_inv,
        diet_vertebrate = diet_vend,
        diet_fruit = diet_fruit,
        diet_seed = diet_seed,
        foraging_stratum = for_strat_value,
        nocturnal = activity_nocturnal,
        crepuscular = activity_crepuscular,
        diurnal = activity_diurnal
      )
    
    write_rds(elton_traits, here("data", "external", "elton_v1.rds")) 
  } else {
    elton_traits <- read_rds(here("data", "external", "elton_v1.rds"))
  }
}, silent = TRUE)

# c) Host Ranges (IUCN RedList)
try({
  if(update_external == TRUE) {
    iucn_redlist <- vect(here("data", "iucn", "MAMMALS_TERRESTRIAL_ONLY.shp"))
    
    mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
    
    iucn_subset <- iucn_redlist %>%
      filter(sci_name %in% comprehensive_checklist$scientific_name) %>%
      filter(presence %in% c(1, 2, 3, 4)) %>%
      aggregate(by = "sci_name", fun = "sum") %>%
      project(mollweide_crs) %>%
      mutate(range_size_km2 = expanse(.) / 1e6) %>%
      select(scientific_name = sci_name,
             origin = origin,
             range_size_km2)
    
    writeVector(iucn_subset, here("data", "external", "iucn_subset.shp"), overwrite = TRUE) 
  } else {
    iucn_subset <- vect(here("data", "external", "iucn_subset.shp"))
  }
}, silent = TRUE)

iucn_range_sizes <- as.data.frame(iucn_subset) %>%
  as_tibble()

# d) Rodent synanthropy
if(update_external == TRUE) {
  
  synanthropy <- read_xlsx(here("data", "external", "rodent_synanthropy.xlsx")) %>%
    select(Rodents, Synanthropic, `S-index`) %>%
    mutate(scientific_name = str_to_sentence(str_replace(Rodents, "_", " ")),
           synanthropy = Synanthropic) %>%
    group_by(scientific_name) %>%
    summarise(s_index =  median(`S-index`, na.rm = TRUE),
              synanthropy = paste(unique(synanthropy), collapse = ", "),
              .groups = "drop")
  
  write_rds(synanthropy, here("data", "external", "synanthropy.rds"))
  
} else {
  
  synanthropy <- read_rds(here("data", "external", "synanthropy.rds"))
  
}

# Join the EltonTraits and IUCN range sizes to create a single, comprehensive traits table
if(update_external == TRUE) {
  
  host_traits <- comprehensive_checklist %>%
    select(-country_distribution) %>%
    full_join(elton_traits, by = "scientific_name") %>%
    full_join(iucn_range_sizes, by = "scientific_name") %>%
    full_join(synanthropy, by = "scientific_name") %>%
    distinct()
  
  resolved_ids <- taxize::get_gbifid(sort(unique(host_traits$scientific_name)), ask = FALSE)
  
  resolved_ids_classification <- classification(resolved_ids, db = "gbif") 
  names(resolved_ids_classification) <- sort(unique(host_traits$scientific_name))
  
  valid_classifications <- resolved_ids_classification[sapply(resolved_ids_classification, is.data.frame)]
  
  resolved_ids_df <- data.table::rbindlist(valid_classifications, idcol = "query") %>%
    as_tibble() %>%
    filter(rank %in% c("species", "genus", "family")) %>%
    pivot_wider(
      id_cols = query,
      names_from = rank,
      values_from = c(name, id),
      names_sep = "_"
    )
  
  unresolved <- host_traits$scientific_name[!host_traits$scientific_name %in% resolved_ids_df$query]
  
  unresolved <- tibble(unresolved_names = unresolved,
                       unresolved_gbif = NA) %>%
    distinct() %>%
    arrange(unresolved_names)
  
  # Loop through each unresolved name
  for (i in 131:length(unresolved$unresolved_names)) {
    
    message(paste("\n--- Resolving:", i, "---"))
    
    name = unresolved$unresolved_names[i]
    
    # Use tryCatch to prevent the loop from stopping if one name fails
    result <- taxize::get_gbifid(name, ask = TRUE)
    
    # Store the result in the list, named by the query
    unresolved$unresolved_gbif[i] <- result

    Sys.sleep(0.5)
  }
  
  unresolved_classifications <- classification(as.gbifid(unresolved$unresolved_gbif), db = "gbif")
  names(unresolved_classifications) = unresolved$unresolved_names
  
  valid_classifications <- unresolved_classifications[sapply(unresolved_classifications, is.data.frame)]
  
  unresolved_ids_df <- data.table::rbindlist(valid_classifications, idcol = "query") %>%
    as_tibble() %>%
    filter(rank %in% c("species", "genus", "family")) %>%
    pivot_wider(
      id_cols = query,
      names_from = rank,
      values_from = c(name, id),
      names_sep = "_"
    )
  
  gbif_harmonised <- bind_rows(resolved_ids_df, unresolved_ids_df) %>%
    arrange(query) %>%
    distinct() %>%
    rename(gbif_species = name_species,
           gbif_genus = name_genus,
           gbif_family = name_family,
           species_gbifid = id_species,
           genus_gbifid = id_genus,
           family_gbifid = id_family)
  
  host_traits <- host_traits %>%
    left_join(gbif_harmonised, by = c("scientific_name" = "query")) %>%
    distinct()
  
  write_rds(host_traits, here("data", "external", "harmonised_species.rds"))
  
} else {
  
  host_traits <- read_rds(here("data", "external", "harmonised_species.rds"))
  
}

# --- 3. Host Bias Analysis ---
# Helper function for calculating the mean, returning NA for all-NA input
safe_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  } else {
    return(mean(x, na.rm = TRUE))
  }
}

# Helper function for calculating the max, returning NA for all-NA input
safe_max <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  } else {
    return(max(x, na.rm = TRUE))
  }
}

# Helper function for calculating the mode (most common value), returning NA for all-NA input
safe_mode <- function(x) {
  x_clean <- stats::na.omit(x)
  if (length(x_clean) == 0) {
    return(NA_character_)
  } else {
    ux <- unique(x_clean)
    return(ux[which.max(tabulate(match(x_clean, ux)))])
  }
}

host_traits_cleaned <- host_traits %>%
  mutate(iucn_status = case_when(str_detect(iucn_status, "EX") ~ "Extinct",
                                 str_detect(iucn_status, "EW") ~ "Extinct in the wild",
                                 str_detect(iucn_status, "CR") ~ "Critically endangered",
                                 str_detect(iucn_status, "EN") ~ "Endangered",
                                 str_detect(iucn_status, "VU") ~ "Vulnerable",
                                 str_detect(iucn_status, "NT") ~ "Near threatened",
                                 str_detect(iucn_status, "LC") ~ "Least concern",
                                 str_detect(iucn_status, "DD") ~ "Data deficient",
                                 str_detect(iucn_status, "NE") ~ "Not evaluated",
                                 TRUE ~ NA),
         iucn_status = factor(iucn_status, levels = c("Extinct", "Extinct in the wild", "Critically endangered", "Endangered", "Vulnerable", "Near threatened", "Least concern", "Data deficient", "Not evaluated")),
         synanthropy = factor(synanthropy, levels = c("TS", "TS, NA", "O", "N, O", "N")),
         species = case_when(scientific_name == gbif_species ~ gbif_species,
                             scientific_name != gbif_species ~ gbif_species,
                             is.na(scientific_name) ~ NA),
         species_match = case_when(scientific_name == gbif_species ~ TRUE,
                                   scientific_name != gbif_species ~ FALSE,
                                   TRUE ~ NA),
         family = case_when(family == gbif_family ~ gbif_family,
                            family != gbif_family ~ gbif_family,
                            is.na(family) ~ NA),
         family_match = case_when(family == gbif_family ~ TRUE,
                                     family != gbif_family ~ FALSE,
                                     is.na(family) | is.na(gbif_family) ~ NA),
         genus = case_when(genus == gbif_genus ~ gbif_genus,
                           genus != gbif_genus ~ gbif_genus,
                            is.na(genus) ~ NA),
         genus_match = case_when(genus == gbif_genus ~ TRUE,
                                 genus != gbif_genus ~ FALSE,
                                  is.na(genus) | is.na(gbif_genus) ~ NA)
         ) %>%
  filter(!is.na(species_gbifid)) %>%
  group_by(species_gbifid, genus_gbifid, family_gbifid) %>%
  mutate(n = n()) %>%
  summarise(
    n = unique(n),
    # Keep all original scientific names that were aggregated
    synonym_names = paste(unique(scientific_name), collapse = ", "),
    scientific_name = first(species),
    genus = first(genus),
    family = first(family),
    order = first(order),
    # Use helpers for aggregation
    adult_body_mass_g = safe_mean(adult_body_mass_g),
    range_size_km2 = safe_mean(range_size_km2),
    s_index = safe_mean(s_index),
    iucn_status = safe_mode(iucn_status),
    synanthropy = safe_mode(synanthropy),
    nocturnal = safe_max(nocturnal),
    crepuscular = safe_max(crepuscular),
    diurnal = safe_max(diurnal),
    # Combine unique biogeographic realms
    biogeographic_realm = paste(unique(na.omit(biogeographic_realm)), collapse = "|"),
    .groups = "drop"
  ) %>%
  mutate(synanthropy = case_when(str_detect(synanthropy, "TS") ~ "Totally Synanthropic",
                                 synanthropy == "O" ~ "Occasionally Synanthropic",
                                 synanthropy == "N" ~ "Not Synanthropic",
                                 TRUE ~ NA)) %>%
  arrange(scientific_name)

write_rds(host_traits_cleaned, here("data", "external", "host_traits_cleaned.rds"))

# 3.1. Create a summary table of sampling effort per host species
host_sampling_summary <- host_data %>%
  drop_na(host_species) %>%
  left_join(
    pathogen_data %>% 
      group_by(host_record_id) %>%
      summarise(pathogen_samples = sum(number_tested, na.rm = TRUE),
                positive_samples = sum(number_positive, na.rm = TRUE)),
    by = "host_record_id"
  ) %>%
  group_by(host_species, host_genus, host_family, gbif_id) %>%
  summarise(
    n_studies = n_distinct(study_id),
    n_pathogen_samples = sum(pathogen_samples, na.rm = TRUE),
    n_positive_samples = sum(positive_samples > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pathogen_samples))

# 3.2. Compare sampled species to the comprehensive checklist
sampled_species <- unique(host_sampling_summary$gbif_id)
sampled_in_checklist <- host_traits_cleaned %>%
  filter(species_gbifid %in% sampled_species)
missing_from_sampling <- host_traits %>%
  filter(!species_gbifid %in% sampled_species)

comparison_data <- host_traits_cleaned %>%
  drop_na(species_gbifid) %>%
  mutate(sampling_status = if_else(species_gbifid %in% sampled_in_checklist$species_gbifid, "Sampled", "Not Sampled"))

# a) Taxonomic Comparison (by Family)
family_comparison <- comparison_data %>%
  drop_na(family) %>%
  count(order, family, sampling_status) %>%
  pivot_wider(names_from = sampling_status, values_from = n, values_fill = 0) %>%
  mutate(proportion_sampled = Sampled / (Sampled + `Not Sampled`)) %>%
  arrange(desc(proportion_sampled))

plot_data_long <- family_comparison %>%
  mutate(total_species = Sampled + `Not Sampled`) %>%
  mutate(family_label = forcats::fct_reorder2(
    paste0(family, " (n=", total_species, ")"), 
    order, 
    -total_species
  )) %>%
  select(family_label, order, Sampled, `Not Sampled`, proportion_sampled) %>%
  tidyr::pivot_longer(
    cols = c(Sampled, `Not Sampled`),
    names_to = "sampling_status",
    values_to = "count"
  ) %>%
  mutate(
    sampling_status = factor(sampling_status, levels = c("Not Sampled", "Sampled")),
    proportion_sampled = if_else(sampling_status == "Not Sampled", NA_real_, proportion_sampled)
  ) %>%
  arrange(family_label, sampling_status)

p_combined <- ggplot(plot_data_long, aes(x = family_label, y = count, fill = proportion_sampled)) +
  geom_col(color = "white") + 
  coord_polar() +
  scale_y_log10() +
  facet_wrap(~ order, scales = "free_x") +
  scale_fill_viridis_c(
    name = "Proportion Sampled",
    labels = scales::percent,
    limits = c(0, 1),
    na.value = "grey80" # This colors the "Not Sampled" bars grey
  ) +
  labs(
    title = "Taxonomic Sampling Bias by Family and Order",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

ggsave(here("output", "analysis_1","family_bias_circular_stacked_plot.png"), 
       plot = p_combined, width = 20, height = 10, dpi = 300)

# b) Biogeographic Realm Comparison
realm_comparison <- comparison_data %>%
  separate_rows(biogeographic_realm, sep = "\\|") %>%
  filter(biogeographic_realm != "") %>%
  mutate(biogeographic_realm = str_remove_all(biogeographic_realm, "\\?")) %>%
  count(biogeographic_realm, sampling_status) %>%
  pivot_wider(names_from = sampling_status, values_from = n, values_fill = 0) %>%
  mutate(proportion_sampled = Sampled / (Sampled + `Not Sampled`)) %>%
  arrange(desc(proportion_sampled))

realm_plot_data_long <- realm_comparison %>%
  mutate(
    total_species = Sampled + `Not Sampled`,
    realm_label = forcats::fct_reorder(paste0(biogeographic_realm, " (n=", total_species, ")"), -total_species)
  ) %>%
  select(realm_label, Sampled, `Not Sampled`, proportion_sampled) %>%
  tidyr::pivot_longer(
    cols = c(Sampled, `Not Sampled`),
    names_to = "sampling_status",
    values_to = "count"
  ) %>%
  mutate(
    sampling_status = factor(sampling_status, levels = c("Not Sampled", "Sampled")),
    proportion_sampled = if_else(sampling_status == "Not Sampled", NA_real_, proportion_sampled)
  ) %>%
  arrange(realm_label, sampling_status)

# Create the final plot
realm_bias_plot <- ggplot(realm_plot_data_long, aes(x = realm_label, y = count, fill = proportion_sampled)) +
  geom_col(color = "white") + 
  coord_polar() +
  scale_y_log10() +
  scale_fill_viridis_c(
    name = "Proportion Sampled",
    labels = scales::percent,
    limits = c(0, 1),
    na.value = "grey80"
  ) +
  labs(
    title = "Geographic Sampling Bias by Biogeographic Realm",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(vjust = 3.5, size = 10),
    axis.text.y = element_blank()
  )

ggsave(here("output", "analysis_1","realm_bias_circular_stacked_plot.png"), 
       plot = realm_bias_plot, width = 15, height = 15, dpi = 300)

# c) IUCN Status Comparison (Visualization)
iucn_plot <- comparison_data %>%
  filter(!is.na(iucn_status)) %>%
  count(iucn_status, sampling_status) %>%
  ggplot(aes(x = iucn_status, y = n, fill = sampling_status)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportion of Species Sampled by IUCN Conservation Status",
    x = "IUCN Status",
    y = "Proportion of Species",
    fill = "Sampling Status"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("output", "analysis_1","iucn_sampling.png"), 
       plot = iucn_plot, width = 15, height = 15, dpi = 300)

# 3.3. Fit a GAM to test predictors of sampling intensity
# Create the base analytical dataset
analysis_data_base <- comparison_data %>%
  # Create a single activity pattern factor for modeling
  mutate(
    activity_pattern = case_when(
      nocturnal == 1 & diurnal == 0 & crepuscular == 0 ~ "Nocturnal",
      nocturnal == 0 & diurnal == 1 & crepuscular == 0 ~ "Diurnal",
      nocturnal == 1 | diurnal == 1 | crepuscular == 1 ~ "Cathemeral", # Any overlap/crepuscular is cathemeral
      TRUE ~ NA_character_
    ),
    activity_pattern = as.factor(activity_pattern),
    synanthropy_status = as.factor(synanthropy)
  )  %>%
  left_join(host_sampling_summary %>%
              select(gbif_id, n_studies, n_pathogen_samples, n_positive_samples), by = c("species_gbifid" = "gbif_id")) %>%
  filter(iucn_status != "Extinct") %>%
  mutate(biogeographic_realm = as.character(biogeographic_realm)) %>%
  separate_rows(biogeographic_realm, sep = "\\|") %>%
  filter(biogeographic_realm != "") %>%
  distinct(scientific_name, species_gbifid, n_pathogen_samples, adult_body_mass_g, range_size_km2, 
           synanthropy_status, activity_pattern, biogeographic_realm) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols = c(scientific_name, species_gbifid, n_pathogen_samples, adult_body_mass_g, range_size_km2, synanthropy_status, activity_pattern),
    names_from = biogeographic_realm,
    values_from = value,
    values_fill = 0
  ) %>%
  mutate(n_pathogen_samples = replace_na(n_pathogen_samples, 0)) %>%
  clean_names()

realm_predictors <- names(analysis_data_base)[names(analysis_data_base) %in% tolower(levels(comprehensive_checklist$biogeographic_realm))]

# GAMS on full dataset

# 3.3.1 Realms
model_1_formula <- as.formula(paste("n_pathogen_samples ~", paste(realm_predictors, collapse = " + ")))
full_data <- analysis_data_base %>% filter(!is.na(range_size_km2) &
                                             !is.na(adult_body_mass_g) &
                                             range_size_km2 > 0 & 
                                             adult_body_mass_g > 0 &
                                             !is.na(activity_pattern))
gam_model_1 <- gam(model_1_formula, data = full_data, family = nb())

summary(gam_model_1)

# 3.3.2 Realms and range size
model_2_formula <- as.formula(paste("n_pathogen_samples ~ s(log10(range_size_km2), k = 64) +", paste(realm_predictors, collapse = " + ")))

gam_model_2 <- gam(model_2_formula, data = full_data, family = nb())

summary(gam_model_2)

# 3.3.3 Realms, range size and body mass
model_3_formula <- as.formula(paste("n_pathogen_samples ~ s(log10(range_size_km2), k = 64) + s(log10(adult_body_mass_g), k = 32) +", paste(realm_predictors, collapse = " + ")))

gam_model_3 <- gam(model_3_formula, data = full_data, family = nb())

summary(gam_model_3)

# 3.3.4 Realms, range size, body mass and activity
model_4_formula <- as.formula(paste("n_pathogen_samples ~ s(log10(range_size_km2), k = 64) + s(log10(adult_body_mass_g), k = 32) + activity_pattern +", paste(realm_predictors, collapse = " + ")))

gam_model_4 <- gam(model_4_formula, data = full_data, family = nb())

summary(gam_model_4)

write_rds(gam_model_4, here("output", "analysis_1", "models", "sampling_intensity_full_dataset.rds"))

# 3.3.5 Model comparison
tibble(model = c(deparse1(model_1_formula),
                 deparse1(model_2_formula),
                 deparse1(model_3_formula),
                 deparse1(model_4_formula)),
       AIC = c(AIC(gam_model_1),
               AIC(gam_model_2),
               AIC(gam_model_3),
               AIC(gam_model_4)))

# GAM on dataset with synanthropy
# 3.3.6 Realms
data_subset <- analysis_data_base %>% filter(!is.na(range_size_km2) &
                                               !is.na(adult_body_mass_g) &
                                               range_size_km2 > 0 & 
                                               adult_body_mass_g > 0 &
                                               !is.na(activity_pattern) &
                                               !is.na(synanthropy_status))
gam_model_5 <- gam(model_1_formula, data = data_subset, family = nb())

summary(gam_model_5)

# 3.3.7 Realms and range size
gam_model_6 <- gam(model_2_formula, data = data_subset, family = nb())

summary(gam_model_6)

# 3.3.8 Realms, range size and body mass
gam_model_7 <- gam(model_3_formula, data = data_subset, family = nb())

summary(gam_model_7)

# 3.3.9 Realms, range size, body mass and activity
gam_model_8 <- gam(model_4_formula, data = data_subset, family = nb())

summary(gam_model_8)

# 3.3.10 Realms, range size, body mass, activity and synanthropy status
model_9_formula <- as.formula(paste("n_pathogen_samples ~ s(log10(range_size_km2), k = 64) + s(log10(adult_body_mass_g), k = 32) + activity_pattern + synanthropy_status +",
                                    paste(realm_predictors, collapse = " + ")))
gam_model_9 <- gam(model_9_formula, data = data_subset, family = nb())

summary(gam_model_9)

write_rds(gam_model_9, here("output", "analysis_1", "models", "sampling_intensity_dataset_subset.rds"))

# 3.3.11 Model comparison
tibble(model = c(deparse1(model_1_formula),
                 deparse1(model_2_formula),
                 deparse1(model_3_formula),
                 deparse1(model_4_formula),
                 deparse1(model_9_formula)),
       AIC = c(AIC(gam_model_5),
               AIC(gam_model_6),
               AIC(gam_model_7),
               AIC(gam_model_8),
               AIC(gam_model_9)))

# --- 4. Pathogen Bias Analysis ---

# 4.1. Create a summary table of detection frequency per virus species
pathogen_detection_summary <- pathogen_data %>%
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) %>%
  filter(!str_detect(pathogen_species_cleaned, ",")) %>%
  group_by(pathogen_species_cleaned, pathogen_family) %>%
  summarise(
    n_positive = sum(number_positive, na.rm = TRUE),
    n_tested = sum(number_tested, na.rm = TRUE),
    n_studies = n_distinct(study_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_tested), desc(n_positive)) %>%
  mutate(pathogen_species_cleaned = fct_rev(fct_inorder(pathogen_species_cleaned)))

arha_pathogen_summary <- ggplot(pathogen_detection_summary) +
  geom_col(aes(y = pathogen_species_cleaned, x = n_tested, fill = n_studies)) +
  facet_wrap(~ pathogen_family, scales = "free_y") +
  scale_x_log10(labels = scales::label_log()) +
  scale_fill_viridis_c() +
  theme(legend.position = "bottom") +
  labs(y = "Pathogen",
       x = "Number of samples tested",
       fill = "Number of included studies") +
  theme_bw()

ggsave(here("output", "analysis_1","arha_pathogen_sampling.png"), 
       plot = arha_pathogen_summary, width = 15, height = 15, dpi = 300)

# 4.2. Analyze sampling effort of different pathogens within a species
pathogens_tested_per_host <- pathogen_data %>%
  left_join(host_data %>% select(host_record_id, host_species), by = "host_record_id") %>%
  filter(!is.na(host_species), !is.na(pathogen_species_cleaned)) %>%
  separate_rows(pathogen_species_cleaned, sep = ",\\s*")

surveillance_breadth_summary <- pathogens_tested_per_host %>%
  group_by(host_species) %>%
  summarise(n_viruses_tested = n_distinct(pathogen_species_cleaned), .groups = "drop") %>%
  arrange(desc(n_viruses_tested))

surveillance_breadth_plot <- surveillance_breadth_summary %>%
  count(n_viruses_tested) %>%
  ggplot(aes(x = n_viruses_tested, y = n)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(
    title = "Distribution of Surveillance Breadth",
    x = "Number of Unique Virus Species Tested For (per host species)",
    y = "Count of Host Species"
  ) +
  theme_minimal()

ggsave(here("output", "analysis_1", "surveillance_breadth_plot.png"), 
       plot = surveillance_breadth_plot, width = 8, height = 6)

# 4.3 Co-detection

co_detection_summary <- pathogen_data %>%
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) %>%
  left_join(host_data %>% select(host_record_id, host_species), by = "host_record_id") %>%
  drop_na(host_species) %>%
  drop_na(pathogen_species_cleaned) %>%
  drop_na(ncbi_id) %>%
  # Separate rows with multiple comma-separated pathogens
  separate_rows(pathogen_species_cleaned, sep = ", ") %>%
  # Now, group and summarize the testing effort for each unique pair
  group_by(host_species, pathogen_species_cleaned, pathogen_family) %>%
  summarise(
    number_tested = sum(number_tested, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(number_tested > 0)

# Pivot the data into a wide host x pathogen matrix
surveillance_matrix <- co_detection_summary %>%
  pivot_wider(
    id_cols = host_species,
    names_from = pathogen_species_cleaned,
    values_from = number_tested,
    values_fill = 0
  )

# Separate the host names from the numeric data for clustering
host_names <- surveillance_matrix$host_species
numeric_matrix <- as.matrix(surveillance_matrix[, -1])
rownames(numeric_matrix) <- host_names

# Perform hierarchical clustering on both hosts (rows) and pathogens (columns)
host_clustering <- hclust(dist(numeric_matrix))
pathogen_clustering <- hclust(dist(t(numeric_matrix)))

# Get the ordered labels from the clustering results
ordered_hosts <- rownames(numeric_matrix)[host_clustering$order]
ordered_pathogens <- colnames(numeric_matrix)[pathogen_clustering$order]

# Convert the matrix back to a long format for ggplot
heatmap_data <- co_detection_summary %>%
  mutate(
    host_species = factor(host_species, levels = ordered_hosts),
    pathogen_species_cleaned = factor(pathogen_species_cleaned, levels = ordered_pathogens)
  )

co_surveillance_heatmap <- ggplot(heatmap_data, aes(x = pathogen_species_cleaned, y = host_species, fill = log1p(number_tested))) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "magma", name = "log(N Tested + 1)") +
  facet_wrap(~ pathogen_family, scales = "free") +
  labs(
    title = "Co-surveillance Landscape of Arenaviruses and Hantaviruses",
    subtitle = "Hosts (y-axis) and pathogens (x-axis) are clustered by sampling similarity",
    x = "Pathogen Species",
    y = "Host Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

co_detection_summary <- pathogen_data %>%
  filter(pathogen_family %in% c("Arenaviridae", "Hantaviridae")) %>%
  left_join(host_data, by = c("host_record_id")) %>%
  select(pathogen_record_id, host_record_id, pathogen_species_cleaned, pathogen_species_ncbi, ncbi_id, assay,
         number_tested, number_positive, host_species, host_genus, gbif_id) %>%
  drop_na(host_species, host_genus, gbif_id, pathogen_species_cleaned, pathogen_species_ncbi) %>%
  group_by(host_species, host_genus, gbif_id, pathogen_species_cleaned, pathogen_species_ncbi, ncbi_id) %>%
  summarise(number_tested = sum(number_tested, na.rm = TRUE),
            number_positive = sum(number_positive, na.rm = TRUE)) %>%
  summarise(n_viruses_tested = n_distinct(pathogen_species_cleaned), .groups = "drop")

# Plot the distribution
ggplot(co_detection_summary, aes(x = n_viruses_tested)) +
  geom_bar() +
  labs(
    title = "Frequency of Co-detection",
    x = "Number of Unique Virus Species Tested per Host Record",
    y = "Count of Host Records"
  ) +
  theme_minimal()

ggsave(here("R", "analysis_1", "outputs", "co_detection_plot.png"), width = 8, height = 6)
