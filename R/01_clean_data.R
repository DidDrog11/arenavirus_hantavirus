source(here::here("R", "00_load_data.R"))

# Cleaned data will be consolidated as combined_data
combined_data <- list()

# Clean citation sheet ----------------------------------------------------
clean_citations <- combined_data_v2$citations %>%
  left_join(bind_rows(
    combined_data_v2$studies %>%
      select(study_id, full_text_id),
    combined_data_v3$studies %>%
      select(study_id, full_text_id)), by = "full_text_id") %>%
  relocate(study_id, .before = 1)


# Clean Rodent sheet ------------------------------------------------------

# Check rodent IDs are unique
rodent_uid <- tibble(source = c(rep("v2", times = nrow(combined_data_v2$host)), rep("v3", times = nrow(combined_data_v3$host))),
                      row_number = c(seq(from = 1, to = nrow(combined_data_v2$host)), seq(from = 1, to = nrow(combined_data_v3$host))),
                      rodent_id = c(combined_data_v2$host$rodent_record_id, combined_data_v3$host$rodent_record_id)) %>%
  group_by(rodent_id) %>%
  mutate(n = n())

table(id_count = rodent_uid$n)

# Clean rodent names

rodent_names <- tibble(rodent_names = sort(c(unique(combined_data_v2$host$scientificName), unique(combined_data_v3$host$scientificName))),
                       species_level = case_when(str_detect(rodent_names, "sp\\.?$|spp\\.?$") ~ FALSE,
                                                 TRUE ~ TRUE)) %>%
  distinct()
rodent_genus <- tibble(rodent_genus = sort(c(unique(combined_data_v2$host$genus), str_split(unique(combined_data_v3$host$scientificName), " ", simplify = TRUE)[, 1])),
                       species_level = FALSE) %>%
  distinct()

# Species level resolution
species_names <- read_rds(here("data", "raw_data", "gbif_species_names.rds"))
species_hierarchy <- read_rds(here("data", "raw_data", "gbif_species_hierarchy.rds"))
# Compare the species names previously resolved with names in the data
# If these do not match it will access GBIF and ask for those with multiple matches to match to the closest species
if(!nrow(species_names) == length(rodent_names$rodent_names[rodent_names$species_level == TRUE])) {
  
  gbif_species <- get_gbifid(rodent_names$rodent_names[rodent_names$species_level == TRUE])
  resolve_species <- classification(gbif_species, db = "gbif")
  species_hierarchy <- rbind(resolve_species) %>% 
    group_by(query) %>% 
    distinct() %>% 
    ungroup() %>% 
    pivot_wider(id_cols = "query", names_from = rank, values_from = name) %>% 
    left_join(rbind(resolve_species) %>% 
                group_by(query) %>% 
                distinct() %>%
                summarise(gbif_id = id[which.max(rank == "species")]), by = "query") %>%
    arrange(kingdom, phylum, class, order, family, genus, species, subspecies) %>%
    drop_na(species)
  species_names <- tibble(current_name = rodent_names$rodent_names[rodent_names$species_level == TRUE],
                          query = as.character(gbif_species)) %>%
    left_join(species_hierarchy %>%
                distinct(species, query, gbif_id), by = c("query")) %>%
    rename("resolved_name" = species)
  write_rds(species_names, here("data", "raw_data", "gbif_species_names.rds"))
  write_rds(species_hierarchy, here("data", "raw_data", "gbif_species_hierarchy.rds"))
  
}
# Some species names remain unmatched
# Need to think about the processes to incorporate or rename these
unmatched_species_names <- species_names %>%
  filter(is.na(gbif_id)) %>%
  pull(current_name)
#gnr_species <- gnr_resolve(unmatched_species_names)

# Genus level resolution
genus_names <- read_rds(here("data", "raw_data", "gbif_genus_names.rds"))
genus_hierarchy <- read_rds(here("data", "raw_data", "gbif_genus_hierarchy.rds"))
# Compare the genus names previously resolved with names in the data
# If these do not match it will access GBIF and ask for those with multiple matches to match to the closest genus
if(!nrow(genus_names) == length(unique(sort(c(rodent_names$rodent_names[rodent_names$species_level == FALSE], rodent_genus$rodent_genus))))) {
  
  gbif_genus <- get_gbifid(unique(sort(c(rodent_names$rodent_names[rodent_names$species_level == FALSE], rodent_genus$rodent_genus))))
  resolve_genus <- classification(gbif_genus, db = "gbif")
  genus_hierarchy <- rbind(resolve_genus) %>%
    group_by(query) %>%
    distinct() %>% 
    ungroup() %>% 
    pivot_wider(id_cols = "query", names_from = rank, values_from = name) %>% 
    left_join(rbind(resolve_genus) %>% 
                group_by(query) %>% 
                distinct() %>%
                summarise(gbif_id = id[which.max(rank == "genus")]), by = "query") %>% 
    arrange(kingdom, phylum, class, order, family, genus) %>%
    drop_na(genus)
  genus_names <- tibble(current_name = unique(sort(c(rodent_names$rodent_names[rodent_names$species_level == FALSE], rodent_genus$rodent_genus))),
                        query = as.character(gbif_genus)) %>%
    left_join(genus_hierarchy %>%
                distinct(genus, query, gbif_id), by = c("query"))
  write_rds(genus_names, here("data", "raw_data", "gbif_genus_names.rds"))
  write_rds(genus_hierarchy, here("data", "raw_data", "gbif_genus_hierarchy.rds"))
  
}

combined_hierarchy <- bind_rows(species_hierarchy,
                                genus_hierarchy) %>%
  distinct() %>%
  left_join(
    bind_rows(species_names,
              genus_names) %>%
      select(current_name, query) %>%
      distinct(),
    by = "query"
  )

clean_v2_host <- combined_data_v2$host %>%
  mutate(scientificName = case_when(is.na(scientificName) ~ genus,
                                    TRUE ~ scientificName)) %>%
  left_join(combined_hierarchy, by = c("scientificName" = "current_name")) %>%
  rename("reported_name" = scientificName)

table(species_id = !is.na(clean_v2_host$species), genus_id = !is.na(clean_v2_host$genus.y))

clean_v3_host <- combined_data_v3$host %>%
  left_join(combined_hierarchy, by = c("scientificName" = "current_name")) %>%
  rename("reported_name" = scientificName)

table(species_id = !is.na(clean_v3_host$species), genus_id = !is.na(clean_v3_host$genus))

# impute non-detected species for studies with summarised detections
studies_to_impute <- combined_data_v3$studies %>%
  filter(!str_detect(data_access, "individual"))

list_v3_host <- clean_v3_host %>%
  group_by(study_id) %>%
  group_split()

expanded_v3_host <- lapply(list_v3_host, function(x) {
  
  if(unique(x$study_id) %in% studies_to_impute$study_id) {
    
    all_species <- unique(x$gbif_id)
    
    site_session_list <- x %>%
      group_by(eventDate, locality, country, verbatimLocality, decimalLatitude, decimalLongitude) %>%
      group_split()
    
    site_session_list <- map(site_session_list, ~ {
      missing_species <- setdiff(all_species, unique(.x$gbif_id))
      if(length(missing_species) > 0) {
        missing_rows <- tibble(
          gbif_id = missing_species,
          # You can add other columns here if needed
          individualCount = 0,
          imputed = TRUE
        )
        bind_rows(.x, missing_rows) %>%
          fill(study_id, eventDate, locality, country, verbatimLocality, coordinate_resolution,
               decimalLatitude, decimalLongitude, .direction = "downup")
      } else {
        .x
      }
    }) %>%
      bind_rows() %>%
      group_by(gbif_id) %>%
      fill(kingdom, phylum, class, order, family, genus, species, .direction = "downup") %>%
      mutate(reported_name = coalesce(reported_name, species, genus)) %>%
      ungroup()
    
  } else {
    
    x
    
  }
}) %>%
  bind_rows()

clean_host <- clean_v2_host %>%
  select(-genus.x) %>%
  rename("genus" = genus.y) %>%
  mutate(version = "v2") %>%
  bind_rows(expanded_v3_host %>%
              mutate(version = "v3"))

# Clean Pathogen sheet ----------------------------------------------------
virus_names <- tibble(virus = sort(unique(c(unique(combined_data_v2$pathogen$scientificName), unique(combined_data_v3$pathogen$scientificName),
                                            unique(combined_data_v2$sequence_data$scientificName), unique(combined_data_v3$sequence_data$scientificName))))) %>%
  # Using NCBI Taxonomy Browser
  mutate(virus_clean = case_when(str_detect(tolower(virus), "amga") ~ "Amga orthohantavirus",
                                 str_detect(tolower(virus), "andes") ~ "Orthohantavirus andesense",
                                 str_detect(tolower(virus), "araraqua") ~ "Araraquara virus",
                                 str_detect(tolower(virus), "asama") ~ "Orthohantavirus asamaense",
                                 str_detect(tolower(virus), "asikkala") ~ "Orthohantavirus asikkalaense",
                                 str_detect(tolower(virus), "bayou|bayoui") ~ "Orthohantavirus bayoui",
                                 str_detect(tolower(virus), "bitu") ~ "Bitu mammarenavirus",
                                 str_detect(tolower(virus), "chapare") ~ "Mammarenavirus chapareense",
                                 str_detect(tolower(virus), "dobrova|dobrava|belgrade|dobravaense") ~ "Orthohantavirus dobravaense",
                                 str_detect(tolower(virus), "hantaan|hantann|hantanense") ~ "Orthohantavirus hantanense",
                                 str_detect(tolower(virus), "jabora") ~ "Jabora hantavirus",
                                 str_detect(tolower(virus), "junin") ~ "Mammarenavirus juninense",
                                 str_detect(tolower(virus), "juquitiba") & !str_detect(tolower(virus), "jabora") ~ "Juquitiba virus",
                                 str_detect(tolower(virus), "kielder")~ "Kielder hantavirus",
                                 str_detect(tolower(virus), "kitale") ~ "Mammarenavirus kitaleense",
                                 str_detect(tolower(virus), "kurkino") ~ "Kurkino virus",
                                 str_detect(tolower(virus), "kwanza") ~ "Kwanza mammarenavirus",
                                 str_detect(tolower(virus), "lassa") ~ "Mammarenavirus lassaense",
                                 str_detect(tolower(virus), "latino") ~ "Mammarenavirus latinum",
                                 str_detect(tolower(virus), "lcmv|lymphocytic") ~ "Mammarenavirus choriomeningitidis",
                                 str_detect(tolower(virus), "limestone") ~ "Limestone Canyon orthohantavirus",
                                 str_detect(tolower(virus), "mayotte") ~ "Mayotte virus",
                                 str_detect(tolower(virus), "mobola") ~ "Mobola mammarenavirus",
                                 str_detect(tolower(virus), "morogoro") ~ "Morogoro mammarenavirus",
                                 str_detect(tolower(virus), "muleshoe") ~ "Muleshoe hantavirus",
                                 str_detect(tolower(virus), "puumala") & !str_detect(tolower(virus), "saaremaa") ~ "Orthohantavirus puumalaense",
                                 str_detect(tolower(virus), "rusne") ~ "Rusne orthohantavirus",
                                 str_detect(tolower(virus), "saaremaa") & !str_detect(tolower(virus), "puumala") ~ "Kurkino virus",
                                 str_detect(tolower(virus), "seewis") ~ "Orthohantavirus seewisense",
                                 str_detect(tolower(virus), "seoul") ~ "Orthohantavirus seoulense",
                                 str_detect(tolower(virus), "sin nombre") ~ "Orthohantavirus sinnombreense",
                                 str_detect(tolower(virus), "thailand") ~ "Orthohantavirus thailandense",
                                 str_detect(tolower(virus), "tigray") ~ "Orthohantavirus tigrayense",
                                 str_detect(tolower(virus), "tula") ~ "Orthohantavirus tulaense",
                                 str_detect(tolower(virus), "wenzhou") ~ "Mammarenavirus wenzhouense",
                                 TRUE ~ virus))

# To-do generate taxonomy from looking up on NCBI or equivalent for pathogens

assay_clean <- tibble(assay = unique(c(unique(combined_data_v2$pathogen$identificationRemarks), unique(combined_data_v3$pathogen$assay)))) %>%
  mutate(assay_clean = case_when(str_detect(tolower(assay), "elisa|antibody|antigen|ifa|immuno|serology|frnt") ~ "Serology",
                                 str_detect(tolower(assay), "pcr|rt-pcr") ~ "PCR"))
# Check pathogen IDs are unique
pathogen_uid <- tibble(source = c(rep("v2", times = nrow(combined_data_v2$pathogen)), rep("v3", times = nrow(combined_data_v3$pathogen))),
                     row_number = c(seq(from = 1, to = nrow(combined_data_v2$pathogen)), seq(from = 1, to = nrow(combined_data_v3$pathogen))),
                     pathogen_id = c(combined_data_v2$pathogen$pathogen_record_id, combined_data_v3$pathogen$pathogen_record_id)) %>%
  group_by(pathogen_id) %>%
  mutate(n = n())

table(id_count = pathogen_uid$n)

# v2 pathogen processing --------------------------------------------------

clean_v2_path <- combined_data_v2$pathogen %>%
  mutate(family = case_when(str_detect(tolower(family), "hanta") ~ "Hantaviridae",
                            str_detect(tolower(family), "arena") ~ "Mammarenaviridae")) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  rename("original_name" = scientificName) %>%
  mutate(virus_clean = coalesce(virus_clean, family),
         original_name = coalesce(original_name, family)) %>%
  left_join(assay_clean, by = c("identificationRemarks" = "assay")) %>%
  rename("original_assay" = identificationRemarks) %>%
  mutate(n_assayed = as.numeric(occurrenceRemarks),
         n_negative = as.numeric(number_negative),
         n_positive = as.numeric(organismQuantity),
         n_inconclusive = as.numeric(number_inconclusive)) %>%
  mutate(n_negative = replace_na(n_negative, 0),
         n_positive = replace_na(n_positive, 0),
         n_inconclusive = replace_na(n_inconclusive, 0)) %>%
  select(pathogen_record_id, associated_rodent_record_id, study_id, family, virus_clean, assay_clean,
         n_assayed, n_positive, n_negative, n_inconclusive, original_name, original_assay, note) %>%
  left_join(clean_host %>%
              select(rodent_record_id, study_id, eventDate, host = reported_name,
                     locality, country, verbatimLocality, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id"))

# Some records are missing an associated rodent record id
orphaned_pathogen_v2 <- clean_v2_path %>%
  filter(is.na(associated_rodent_record_id))

# v3 pathogen processing --------------------------------------------------

clean_v3_path <- combined_data_v3$pathogen %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  rename("original_name" = scientificName) %>%
  mutate(virus_clean = coalesce(virus_clean, family),
         original_name = coalesce(original_name, family)) %>%
  left_join(assay_clean, by = "assay") %>%
  mutate(n_assayed = as.numeric(tested),
         n_negative = as.numeric(negative),
         n_positive = as.numeric(positive),
         n_inconclusive = as.numeric(number_inconclusive)) %>%
  mutate(n_negative = replace_na(n_negative, 0),
         n_positive = replace_na(n_positive, 0),
         n_inconclusive = replace_na(n_inconclusive, 0)) %>%
  select(pathogen_record_id, associated_rodent_record_id, study_id, family, virus_clean, assay_clean,
         n_assayed, n_positive, n_negative, n_inconclusive, original_name, note) %>%
  left_join(clean_host %>%
              select(rodent_record_id, study_id, eventDate, host = reported_name,
                     locality, country, verbatimLocality, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id"))

orphaned_pathogen_v3 <- clean_v3_path %>%
  filter(is.na(associated_rodent_record_id))

orphaned_pathogen <- bind_rows(orphaned_pathogen_v2, orphaned_pathogen_v3)

clean_pathogen <- clean_v2_path %>%
  mutate(version = "v2") %>%
  bind_rows(clean_v3_path %>%
              mutate(version = "v3"))

# Clean sequence_data sheet -----------------------------------------------

pathogen_sequences_v2 <- combined_data_v2$sequence_data %>%
  filter(str_detect(sequenceType, "Pathogen")) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, virus_clean, original_name = scientificName, accession_number) %>%
  left_join(clean_host %>%
              select(rodent_record_id, study_id, eventDate, species, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id")) %>%
  left_join(clean_pathogen %>%
              select(pathogen_record_id, study_id, n_assayed, n_positive, n_negative, host),
            by = c("associated_pathogen_record_id" = "pathogen_record_id", "study_id"))

host_sequences_v2 <-  combined_data_v2$sequence_data %>%
  filter(str_detect(sequenceType, "Host")) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, virus_clean, original_name = scientificName, accession_number) %>%
  left_join(clean_host %>%
              select(rodent_record_id, study_id, eventDate, species, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id")) %>%
  left_join(clean_pathogen %>%
              select(pathogen_record_id, study_id, n_assayed, n_positive, n_negative, host),
            by = c("associated_pathogen_record_id" = "pathogen_record_id", "study_id"))

other_sequences_v2 <-  combined_data_v2$sequence_data %>%
  filter(!str_detect(sequenceType, "Host|Pathogen")) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, virus_clean, original_name = scientificName, accession_number) %>%
  left_join(clean_host %>%
              select(rodent_record_id, study_id, eventDate, species, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id")) %>%
  left_join(clean_pathogen %>%
              select(pathogen_record_id, study_id, n_assayed, n_positive, n_negative, host),
            by = c("associated_pathogen_record_id" = "pathogen_record_id", "study_id"))

clean_sequences_v2 <- bind_rows(pathogen_sequences_v2,
                                host_sequences_v2,
                                other_sequences_v2) %>%
  mutate(version = "v2")

pathogen_sequences_v3 <- combined_data_v3$sequence_data %>%
  filter(str_detect(sequenceType, "Pathogen")) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, virus_clean, original_name = scientificName, accession_number) %>%
  left_join(clean_host %>%
              select(rodent_record_id, study_id, eventDate, species, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id")) %>%
  left_join(clean_pathogen %>%
              select(pathogen_record_id, study_id, n_assayed, n_positive, n_negative, host),
            by = c("associated_pathogen_record_id" = "pathogen_record_id", "study_id"))

host_sequences_v3 <-  combined_data_v3$sequence_data %>%
  filter(str_detect(sequenceType, "Host")) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, virus_clean, original_name = scientificName, accession_number) %>%
  left_join(clean_host %>%
              select(rodent_record_id, study_id, eventDate, species, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id")) %>%
  left_join(clean_pathogen %>%
              select(pathogen_record_id, study_id, n_assayed, n_positive, n_negative, host),
            by = c("associated_pathogen_record_id" = "pathogen_record_id", "study_id"))

other_sequences_v3 <-  combined_data_v3$sequence_data %>%
  filter(!str_detect(sequenceType, "Host|Pathogen")) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, virus_clean, original_name = scientificName, accession_number) %>%
  left_join(clean_host %>%
              select(rodent_record_id, study_id, eventDate, species, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id")) %>%
  left_join(clean_pathogen %>%
              select(pathogen_record_id, study_id, n_assayed, n_positive, n_negative, host),
            by = c("associated_pathogen_record_id" = "pathogen_record_id", "study_id"))

clean_sequences_v3 <- bind_rows(pathogen_sequences_v3,
                                host_sequences_v3,
                                other_sequences_v3) %>%
  mutate(version = "v3")

clean_sequences <- bind_rows(clean_sequences_v2,
                             clean_sequences_v3)


# Clean data --------------------------------------------------------------
combined_data <- list(citations = clean_citations,
                      host = clean_host,
                      pathogen = clean_pathogen,
                      sequences = clean_sequences)


write_rds(combined_data, here("data", "clean_data", paste0(Sys.Date(), "_data.rds")))
