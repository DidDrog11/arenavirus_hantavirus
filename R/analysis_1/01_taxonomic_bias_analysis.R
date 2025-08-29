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
library(readr)
library(readxl)
library(dplyr)
library(janitor)
library(stringr)
library(taxize)
library(terra)
library(tidyterra)
library(tidyr)
library(traitdata)
library(mgcv) # For GAMs

# Update external data
update_external = TRUE

# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-08-25.rds") # <-- Update this date
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
# This should have columns like 'scientificName', 'adult_body_mass_g', 'range_size_km2', etc.
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
             range_size_km2)
    
    writeVector(iucn_subset, here("data", "external", "iucn_subset.shp")) 
  } else {
    iucn_subset <- vect(here("data", "external", "iucn_subset.shp"))
  }
}, silent = TRUE)

iucn_range_sizes <- as.data.frame(iucn_subset) %>%
  as_tibble() %>%
  rename(scientific_name = scientific)

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
  
  write_rds(synanthropy_combined, here("data", "external", "synanthropy.rds"))
  
} else {
  
  synanthropy <- read_rds(here("data", "external", "synanthropy.rds"))
  
}

# Join the EltonTraits and IUCN range sizes to create a single, comprehensive traits table
if(update_external == TRUE) {
  
  host_traits <- comprehensive_checklist %>%
    select(-country_distribution) %>%
    full_join(elton_traits, by = "scientific_name") %>%
    full_join(iucn_range_sizes, by = "scientific_name") %>%
    full_join(synanthropy, by = "scientific_name")
  
  resolved_ids <- taxize::get_gbifid(sort(unique(host_traits$scientific_name)), ask = FALSE)
  
  resolved_ids_classification <- classification(resolved_ids, db = "gbif") 
  names(resolved_ids_classification) <- sort(unique(host_traits$scientific_name))
  
  valid_classifications <- resolved_ids_classification[sapply(resolved_ids_classification, is.data.frame)]
  
  resolved_ids_df <- data.table::rbindlist(valid_classifications, idcol = "query") %>%
    as_tibble() %>%
    filter(rank == "species")
  
  unresolved <- host_traits$scientific_name[!host_traits$scientific_name %in% resolved_ids_df$query]
  
  ---
  
  unresolved <- tibble(unresolved_names = unresolved,
                       unresolved_gbif = NA)
  
  # Loop through each unresolved name
  for (i in 1:length(unresolved$unresolved_names)) {
    
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
    filter(rank == "species")
  
  gbif_harmonised <- bind_rows(resolved_ids_df, unresolved_ids_df) %>%
    arrange(query) %>%
    distinct()
  
  host_traits <- host_traits %>%
    left_join(gbif_harmonised %>%
                rename(scientific_name = query, resolved_name = name, gbif_id = id), by = "scientific_name")
  
  write_rds(host_traits, here("data", "external", "harmonised_species.rds"))
  
} else {
  
  host_traits <- read_rds(here("data", "external", "harmonised_species.rds"))
  
}

# --- 3. Host Bias Analysis ---

host_traits <- host_traits %>%
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
         synanthropy = factor(synanthropy, levels = c("TS", "TS, NA", "O", "N, O", "N"))) %>%
  filter(!is.na(gbif_id)) %>%
  group_by(gbif_id) %>%
  summarise(
    # Keep the first scientific name as the representative for the group
    scientific_name = paste(scientific_name, collapse = ", "),
    # Take the mean for continuous traits, ignoring NAs
    adult_body_mass_g = mean(adult_body_mass_g, na.rm = TRUE),
    range_size_km2 = mean(range_size_km2, na.rm = TRUE),
    s_index = mean(s_index, na.rm = TRUE),
    # For categorical traits, take the most common value (mode)
    iucn_status = names(which.max(table(iucn_status))),
    synanthropy = names(which.max(table(synanthropy))),
    nocturnal = max(nocturnal, na.rm = TRUE),
    crepuscular = max(crepuscular, na.rm = TRUE),
    diurnal = max(diurnal, na.rm = TRUE),
    biogeographic_realm = paste(unique(na.omit(biogeographic_realm)), collapse = "|"),
    .groups = "drop"
  )


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
sampled_in_checklist <- host_traits %>%
  filter(gbif_id %in% sampled_species)
missing_from_sampling <- host_traits %>%
  filter(!gbif_id %in% sampled_species)

comparison_data <- host_traits %>%
  drop_na(gbif_id) %>%
  mutate(sampling_status = if_else(gbif_id %in% sampled_in_checklist$gbif_id, "Sampled", "Not Sampled"))

# a) Taxonomic Comparison (by Family)
family_comparison <- comparison_data %>%
  count(family, sampling_status) %>%
  pivot_wider(names_from = sampling_status, values_from = n, values_fill = 0) %>%
  mutate(proportion_sampled = Sampled / (Sampled + `Not Sampled`)) %>%
  arrange(desc(proportion_sampled))

# b) Biogeographic Realm Comparison
realm_comparison <- comparison_data %>%
  separate_rows(biogeographic_realm, sep = "\\|") %>%
  filter(biogeographic_realm != "") %>%
  mutate(biogeographic_realm = str_remove_all(biogeographic_realm, "\\?")) %>%
  count(biogeographic_realm, sampling_status) %>%
  pivot_wider(names_from = sampling_status, values_from = n, values_fill = 0) %>%
  mutate(proportion_sampled = Sampled / (Sampled + `Not Sampled`)) %>%
  arrange(desc(proportion_sampled))

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
              select(gbif_id, n_studies, n_pathogen_samples, n_positive_samples), by = "gbif_id") %>%
  filter(iucn_status != "Extinct") %>%
  mutate(biogeographic_realm = as.character(biogeographic_realm)) %>%
  separate_rows(biogeographic_realm, sep = "\\|") %>%
  filter(biogeographic_realm != "") %>%
  distinct(scientific_name, gbif_id, n_pathogen_samples, adult_body_mass_g, range_size_km2, 
           synanthropy_status, activity_pattern, foraging_stratum, biogeographic_realm) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols = c(scientific_name, gbif_id, n_pathogen_samples, adult_body_mass_g, range_size_km2, synanthropy_status, activity_pattern, foraging_stratum),
    names_from = biogeographic_realm,
    values_from = value,
    values_fill = 0
  ) %>%
  clean_names()

realm_predictors <- names(analysis_data_base)[names(analysis_data_base) %in% tolower(levels(comprehensive_checklist$biogeographic_realm))]

model_1_formula <- as.formula(paste("n_pathogen_samples ~", paste(realm_predictors, collapse = " + ")))
model_1_data <- analysis_data_base %>% select(n_pathogen_samples, all_of(realm_predictors))

# --- 4. Pathogen Bias Analysis ---

# 4.1. Create a summary table of detection frequency per virus species
pathogen_detection_summary <- pathogen_data %>%
  group_by(pathogen_species_cleaned, pathogen_family) %>%
  summarise(
    n_positive_records = n(),
    n_hosts_tested_on = n_distinct(host_record_id),
    n_studies_detected_in = n_distinct(study_id),
    .groups = "drop"
  ) %>%
  arrange(desc(n_positive_records))

write_csv(pathogen_detection_summary, here("R", "analysis_1", "outputs", "pathogen_detection_summary.csv"))

# 4.2. Analyze co-detection frequency (number of viruses tested per host record)
co_detection_summary <- pathogen_data %>%
  group_by(host_record_id) %>%
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
