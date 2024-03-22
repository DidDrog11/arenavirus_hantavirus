source(here::here("R", "00_load_data.R"))


# Clean citation sheet ----------------------------------------------------
combined_data$citations <- combined_data$citations %>%
  left_join(combined_data$studies %>%
              select(study_id, full_text_id), by = "full_text_id") %>%
  relocate(study_id, .before = 1)


# Clean Pathogen sheet ----------------------------------------------------
virus_names <- tibble(virus = sort(unique(combined_data$pathogen$scientificName))) %>%
  mutate(virus_clean = case_when(str_detect(tolower(virus), "andes") ~ "Andes orthohantavirus",
                                 str_detect(tolower(virus), "araraqua") ~ "Araraquara orthohantavirus",
                                 str_detect(tolower(virus), "asama") ~ "Asama orthohantavirus",
                                 str_detect(tolower(virus), "bayou|bayoui") ~ "Bayou orthohantavirus",
                                 str_detect(tolower(virus), "chapare") ~ "Chapare mannarenavirus",
                                 str_detect(tolower(virus), "dobrova|dobrava|belgrade|dobravaense") ~ "Dobrava-Belgrade orthohantavirus",
                                 str_detect(tolower(virus), "hantaan|hantann") ~ "Hantaan orthohantavirus",
                                 str_detect(tolower(virus), "junin") ~ "Junin mammarenavirus",
                                 str_detect(tolower(virus), "juquitiba") & !str_detect(tolower(virus), "jabora") ~ "Juquitiba orthohantavirus",
                                 str_detect(tolower(virus), "lassa") ~ "Lassa mammarenavirus",
                                 str_detect(tolower(virus), "lcmv|lymphocytic") ~ "Lymphocytic choriomeningitis mammarenavirus",
                                 str_detect(tolower(virus), "limestone") ~ "Limestone Canyon orthohantavirus",
                                 str_detect(tolower(virus), "muleshoe") ~ "Muleshoe orthohantavirus",
                                 str_detect(tolower(virus), "puumala") & !str_detect(tolower(virus), "saaremaa") ~ "Puumala orthohantavirus",
                                 str_detect(tolower(virus), "rusne") ~ "Rusne orthohantavirus",
                                 str_detect(tolower(virus), "seewis") ~ "Seewis orthohantavirus",
                                 str_detect(tolower(virus), "seoul") ~ "Seoul orthohantavirus",
                                 str_detect(tolower(virus), "sin nombre") ~ "Sin Nombre orthohantavirus",
                                 str_detect(tolower(virus), "thailand") ~ "Thailand orthohantavirus",
                                 str_detect(tolower(virus), "tula") ~ "Tula orthohantavirus",
                                 str_detect(tolower(virus), "wenzhou") ~ "Wenzhou mammarenavirus",
                                 TRUE ~ virus))

assay_clean <- tibble(assay = sort(unique(combined_data$pathogen$identificationRemarks))) %>%
  mutate(assay_clean = case_when(str_detect(tolower(assay), "elisa|antibody|antigen|ifa|immuno|serology") ~ "Serology",
                                 str_detect(tolower(assay), "pcr") ~ "PCR"))

combined_data$pathogen <- combined_data$pathogen %>%
  mutate(host_genus = str_to_sentence(host_genus),
         associatedTaxa = str_to_sentence(associatedTaxa),
         decimalLatitude = round(decimalLatitude, 2),
         decimalLongitude = round(decimalLongitude, 2),
         family = case_when(str_detect(tolower(family), "hanta") ~ "Hantaviridae",
                            str_detect(tolower(family), "arena") ~ "Mammarenaviridae")) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  select(-scientificName, scientificName = virus_clean) %>%
  left_join(assay_clean, by = c("identificationRemarks" = "assay")) %>%
  select(-identificationRemarks, identificationRemarks = assay_clean) %>%
  mutate(occurrenceRemarks = as.numeric(occurrenceRemarks),
         number_negative = as.numeric(number_negative),
         organismQuantity = as.numeric(organismQuantity),
         number_inconclusive = as.numeric(number_inconclusive)) %>%
  mutate(number_negative = replace_na(number_negative, 0),
         organismQuantity = replace_na(organismQuantity, 0),
         number_inconclusive = replace_na(number_inconclusive, 0)) %>%
  select(any_of(names(combined_data$pathogen)))


# Clean sequence_data sheet -----------------------------------------------
virus_names <- tibble(virus = sort(unique(combined_data$sequence_data$scientificName))) %>%
  mutate(virus_clean = case_when(str_detect(tolower(virus), "amga") ~ "Amga orthohantavirus",
                                 str_detect(tolower(virus), "andes") ~ "Andes orthohantavirus",
                                 str_detect(tolower(virus), "araraqua") ~ "Araraquara orthohantavirus",
                                 str_detect(tolower(virus), "asama") ~ "Asama orthohantavirus",
                                 str_detect(tolower(virus), "bayou|bayoui") ~ "Bayou orthohantavirus",
                                 str_detect(tolower(virus), "bitu") ~ "Bitu mammarenavirus",
                                 str_detect(tolower(virus), "chapare") ~ "Chapare mannarenavirus",
                                 str_detect(tolower(virus), "dobrova|dobrava|belgrade|dobravaense") ~ "Dobrava-Belgrade orthohantavirus",
                                 str_detect(tolower(virus), "hantaan|hantann") ~ "Hantaan orthohantavirus",
                                 str_detect(tolower(virus), "jabora") ~ "Jabora orthohantavirus",
                                 str_detect(tolower(virus), "junin") ~ "Junin mammarenavirus",
                                 str_detect(tolower(virus), "juquitiba") & !str_detect(tolower(virus), "jabora") ~ "Juquitiba orthohantavirus",
                                 str_detect(tolower(virus), "kitale") ~ "Kitale mammarenavirus",
                                 str_detect(tolower(virus), "kwanza") ~ "Kwanza mammarenavirus",
                                 str_detect(tolower(virus), "lassa") ~ "Lassa mammarenavirus",
                                 str_detect(tolower(virus), "lcmv|lymphocytic") ~ "Lymphocytic choriomeningitis mammarenavirus",
                                 str_detect(tolower(virus), "limestone") ~ "Limestone Canyon orthohantavirus",
                                 str_detect(tolower(virus), "muleshoe") ~ "Muleshoe orthohantavirus",
                                 str_detect(tolower(virus), "puumala") & !str_detect(tolower(virus), "saaremaa") ~ "Puumala orthohantavirus",
                                 str_detect(tolower(virus), "rusne") ~ "Rusne orthohantavirus",
                                 str_detect(tolower(virus), "seewis") ~ "Seewis orthohantavirus",
                                 str_detect(tolower(virus), "seoul") ~ "Seoul orthohantavirus",
                                 str_detect(tolower(virus), "sin nombre") ~ "Sin Nombre orthohantavirus",
                                 str_detect(tolower(virus), "thailand") ~ "Thailand orthohantavirus",
                                 str_detect(tolower(virus), "tula") ~ "Tula orthohantavirus",
                                 str_detect(tolower(virus), "wenzhou") ~ "Wenzhou mammarenavirus",
                                 TRUE ~ virus))


combined_data$sequence_data <- combined_data$sequence_data %>%
  mutate(study_id = as.numeric(study_id)) %>%
  left_join(virus_names, by = c("scientificName" = "virus")) %>%
  left_join(combined_data$host %>%
              select(rodent_record_id, eventDate, coordinate_resolution, decimalLatitude, decimalLongitude, individualCount),
            by = c("associated_rodent_record_id" = "rodent_record_id")) %>%
  left_join(combined_data$pathogen %>%
              select(pathogen_record_id, coordinate_resolution, decimalLatitude, decimalLongitude, family),
            by = c("associated_pathogen_record_id" = "pathogen_record_id")) %>%
  mutate(coordinate_resolution = coalesce(coordinate_resolution.y, coordinate_resolution.x),
         decimalLatitude = coalesce(decimalLatitude.y, decimalLatitude.x),
         decimalLongitude = coalesce(decimalLongitude.y, decimalLongitude.x),
         family = case_when(sequenceType == "Host" ~ "",
                            str_detect(virus_clean, "hanta|HV|Hanta") ~ "Hantaviridae",
                            str_detect(virus_clean, "arena|AreV|Arena") ~ "Arenaviridae")) %>%
  select(sequence_record_id, associated_rodent_record_id, associated_pathogen_record_id, eventDate, study_id, host_genus, associatedTaxa,
         sequenceType, family, virus_clean, coordinate_resolution, decimalLatitude, decimalLongitude, accession_number, method, note, date_sampled, sample_location)

write_rds(combined_data, here("data", "clean_data", paste0(Sys.Date(), "_data.rds")))
