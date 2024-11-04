library(tidyverse)
library(readxl)
library(here)
library(stringr)
library(tidygeocoder)

rodent_data %>%
  filter(study_id == 15) %>%
  group_by(scientificName, country) %>%
  summarise(n = sum(as.numeric(individualCount), na.rm = TRUE)) %>%
  arrange(country, scientificName) %>%
  mutate(genus = str_split(scientificName, " ", simplify = TRUE)[1]) %>%
  select(n) %>%
  write.table(., "clipboard", sep="\t", row.names=FALSE)


# Extracting data from ft_12 ----------------------------------------------
ft_12_t2 <- read_xlsx(here("data", "data_to_extract", "ft_12_mod_t2.xlsx"))
names(ft_12_t2) <- c("site", "year", "species", "trapped", "leptospira", "leptospira_type", "hantavirus")

ft_12_locations <- tibble(site = unique(ft_12_t2$site[!is.na(ft_12_t2$site)]),
                          latitude = c(56.0782, 55.5019, 54.1886, 56.002, 55.7296, 55.2416, 55.543, 56.0219, 55.5556, 54.3979, 55.8452, 56.01776, 56.1891, 54.5435, 55.3013, 55.3913, 54.7929, 54.6373, 55.6065, 55.6386, 54.4184, 55.9774, 55.7475),
                          longitude = c(24.3937, 25.0324, 24.6701, 24.6288, 24.408, 22.4707, 21.1213, 25.3896, 21.9929, 23.9312, 26.19353, 24.7191, 24.6967, 23.9616, 21.3650, 24.7556, 25.5985, 24.9334, 23.1757, 25.5809, 23.9707, 24.7577, 25.7674))

ft_12_t2 <- ft_12_t2 %>%
  filter(species != "all") %>%
  fill(site, .direction = "down") %>%
  left_join(ft_12_locations, by = "site") %>%
  mutate(site = str_remove_all(site, pattern = " \\(\\d+\\)"),
         species = str_replace_all(species, pattern = "Apodemus flavicolis", replacement = "Apodemus flavicollis")) %>%
  fill(year, .direction = "down") %>%
  select(year, species, site, latitude, longitude, trapped)

ft_12_t2_rodent <- ft_12_t2 %>%
  distinct(year, site, latitude, longitude) %>%
  mutate(site = fct_inorder(site)) %>%
  cross_join(ft_12_t2 %>%
               distinct(species) %>%
               select(species)) %>%
  arrange(site, year, species) %>%
  left_join(ft_12_t2 %>%
              select(-latitude, -longitude), by = c("year", "site", "species")) %>%
  mutate(genus = str_split(species, pattern = " ", simplify = TRUE)[,1],
         trapped = replace_na(trapped, 0)) %>%
  fill(year, .direction = "down") %>%
  select(year, site, genus, species, latitude, longitude, trapped) %>%
  write.table(., "clipboard", sep="\t", row.names=FALSE)

# Extracting data from ft_19 ----------------------------------------------

ft_19_supp_2 <- read_xlsx(here("data", "data_to_extract", "ft_19_supp_2.xlsx"))

# Function to convert degrees, minutes, seconds in string format to decimal degrees
dms_str_to_dd <- function(coord_str) {
  # Extract the first two numbers as degrees
  degrees <- as.numeric(substr(coord_str, 1, 2))
  # Remove the first two numbers and trim the remaining string
  remaining_str <- str_trim(substr(coord_str, 3, nchar(coord_str)), side = "both")
  # Replace non-numeric characters with "." if surrounded by numbers
  remaining_str <- gsub("(?<=[0-9])[^0-9]+(?=[0-9])", ".", remaining_str, perl = TRUE)
  # Remove the last non-numeric character if it is not a number
  if (!grepl("[0-9]$", substr(remaining_str, nchar(remaining_str), nchar(remaining_str)))) {
    remaining_str <- gsub("[^0-9]+$", "", remaining_str)
  }
  # Convert to numeric
  numeric_value <- as.numeric(remaining_str)
  # Calculate decimal degrees
  dd <- degrees + numeric_value/60
  
  return(dd)
}

# Coordinates in string format
converted_coords <- ft_19_supp_2 %>%
  mutate(Latitude = case_when(Latitude == "14 07.641'" ~ "24 07.641'", # correct a coordinate
                              TRUE ~ Latitude)) %>%
  rowwise() %>%
  mutate(dd_lat = round(dms_str_to_dd(Latitude) * -1, 4),
         dd_lon = round(dms_str_to_dd(Longitude) * -1, 4)) %>%
  select(Genus, Species, Date, RNA, dd_lat, dd_lon) # keep these to match to NCBI manually

converted_coords %>%
  select(dd_lat, dd_lon) %>%
  write.table(., "clipboard", sep="\t", row.names=FALSE)

sf::st_as_sf(converted_coords, coords = c("dd_lon", "dd_lat"), crs = "EPSG:4326") %>%
  mapview::mapview()

ft_19_supp_2 %>%
  mutate(PCR = case_when(RNA == 0 ~ 0,
                         RNA == 1 ~ 1,
                         TRUE ~ NA),
         Antibody = case_when(`Sero status` == 0 ~ 0,
                                   `Sero status` == 1 ~ 1,
                                   TRUE ~ NA),
         IFA_note = case_when(`IFA Titer` == 0 ~ "Negative titer < 32",
                              `IFA Titer` == 1 ~ "Postive titer (32-64)",
                              `IFA Titer` == 2 ~ "Postive titer (128-256)",
                              `IFA Titer` == 3 ~ "Postive titer (512-1024)",
                              `IFA Titer` == 4 ~ "Postive titer (2048-5096)",
                              `IFA Titer` == 5 ~ "Postive titer (>8192)",
                              TRUE ~ NA),
         IFA = case_when(`IFA Titer` == 0 ~ 0,
                         `IFA Titer` %in% c(1, 2, 3, 4) ~ 1,
                         TRUE ~ NA)) %>%
  select(PCR, Antibody, IFA, IFA_note) %>%
  write.table(., "clipboard", sep="\t", row.names=FALSE)


# ft_20 centre of study region --------------------------------------------
bw <- tibble(site = c("BWF1", "BWF2", "BWF3", "BWG1", "BWG2", "BWG3"),
            lat = c(48.829307, 48.829981, 48.844419, 48.841850, 48.844470, 48.810801),
            lon = c(8.966542, 8.962033, 8.957213, 8.951865, 8.945196, 8.968887)) %>%
  sf::st_as_sf(., coords = c("lon", "lat"), crs = "EPSG:4326")

bw %>%
  summarize(geometry = sf::st_union(geometry)) %>% 
  sf::st_centroid()

nw <- tibble(site = c("NWF1", "NWF2", "NWF3", "NWG1", "NWG2", "NWG3"),
             lat = c(51.998752, 51.978863, 51.993606, 51.993511, 51.987140, 51.980575),
             lon = c(7.335144, 7.329545, 7.316983, 7.314114, 7.325719, 7.321710)) %>%
  sf::st_as_sf(., coords = c("lon", "lat"), crs = "EPSG:4326")

nw %>%
  summarize(geometry = sf::st_union(geometry)) %>% 
  sf::st_centroid()

# ft_21 formatting individual level data ----------------------------------

ft_21_supp_1 <-  read_xlsx(here("data", "data_to_extract", "ft_21_supp_1.xlsx"))

ft_21_supp_1 %>%
  mutate(date = as.Date(`Sample Collection Date`),
         genus = str_split(`Genus & Species`, " ", simplify = TRUE)[, 1],
         species = `Genus & Species`,
         serology = as.numeric(`Sero Status`),
         ifa_titre = case_when(`IFA Titer` == 0 ~ "Negative titer < 32",
                               `IFA Titer` == 1 ~ "Postive titer (32-64)",
                               `IFA Titer` == 2 ~ "Postive titer (128-256)",
                               `IFA Titer` == 3 ~ "Postive titer (512-1024)",
                               `IFA Titer` == 4 ~ "Postive titer (2048-5096)",
                               `IFA Titer` == 5 ~ "Postive titer (>8192)",
                               TRUE ~ NA),
         verbatim_locality = Grid) %>%
  select(date, genus, species, verbatim_locality, serology, ifa_titre) %>%
  write.table(., here("data", "data_to_extract", "ft_21_processed.txt"), sep="\t", row.names=FALSE)


# ft_35 formatting sequence data ------------------------------------------

ft_35_supp_1 <- read_xlsx(here("data", "data_to_extract", "ft_35_supp_1.xlsx"))

ft_35_supp_1 %>%
  arrange(`Accession number`) %>%
  mutate(`Date captured` = as.Date(`Date captured`, "%d.%m.%Y")) %>%
  select(`Accession number`, `Date captured`, Area) %>%
  write.table("clipboard", sep="\t", row.names=FALSE)

# ft_40 formatting --------------------------------------------------------

ft_40_supp_4 <- read_xlsx(here("data", "data_to_extract", "ft_40_supp_4.xlsx"))

ft_40_locations <- tibble(Districts = unique(ft_40_supp_4$Districts),
                          coordinate_resolution = c("village", "village", "village", "village", "city", "village", "village"),
                          lat = c(-1.6309, -1.6443, -1.6505, -1.6229, -1.6324, -1.6267, -1.5995),
                          lon = c(13.5475,  13.5711, 13.6081, 13.6029, 13.5881, 13.6137, 13.5996))

species_locations <- tibble(Districts = unique(ft_40_supp_4$Districts),
                            Species = unique(ft_40_supp_4$Species))

# Creating all possible combinations
species_locations <- expand(species_locations, Districts, Species)

ft_40_supp_4 %>%
  mutate(n = 1) %>%
  full_join(species_locations, .) %>%
  mutate(genus = str_split(Species, " ", simplify = TRUE)[ , 1],
         taxon = case_when(str_detect(Species, "sp$") ~ "genus",
                           TRUE ~ "species"),
         country = "Gabon") %>%
  left_join(ft_40_locations) %>%
  select(taxon, genus, Species, Districts, country, coordinate_resolution, lat, lon, n, Hantavirus, Mammarenavirus, LCMV) %>%
  mutate(n = replace_na(n, 0),
         Hantavirus = replace_na(Hantavirus, 0),
         Mammarenavirus = replace_na(Mammarenavirus, 0),
         LCMV = replace_na(LCMV, 0)) %>%
  write.table("clipboard", sep="\t", row.names=FALSE)

# ft_44 formatting --------------------------------------------------------

ft_44_t1 <- read_xlsx(here("data", "data_to_extract", "ft_44_table_1.xlsx"))

ft_44_species <- sort(unique(ft_44_t1$`Animal host`))


# ft_45 formatting --------------------------------------------------------
ft_45 <- read_xlsx(here("data", "data_to_extract", "ft_45_supp.xlsx"))

ft_45_locations <- tibble(site = c("Dingzhou", "Luan", "Xinji", "Xian", "Guantao", "Kuancheng"),
                          lat = c(38.516, 39.5, 37.943, 38.19, 36.533, 40.611),
                          lon = c(114.99, 118.7, 115.218, 116.123, 115.3, 118.485))
ft_45 %>%
  select(-`No. rats`, -`Other species`, -Density) %>%
  fill(year, .direction = "down") %>%
  pivot_longer(cols = c(`Rattus norvegicus`, `Mus musculus`, `Tscherskia triton`, `Apodemus agrarius`),
               names_to = "species",
               values_to = "n") %>%
  left_join(ft_45_locations, by = "site") %>%
  mutate(genus = str_split(species, " ", simplify = TRUE)[, 1]) %>%
  select(year, genus, species, site, setting, lat, lon, n, trap_effort = `No. rattraps`, n_positive = `Virus-carrying Rate`) %>%
  arrange(year, site, setting, genus, species) %>%
  mutate(n_positive = n_positive*10) %>%
  select(n_positive) %>%
  write.table("clipboard", sep="\t", row.names=FALSE)


# ft_48 formatting --------------------------------------------------------

ft_48_1 <- read_xlsx(here("data", "data_to_extract", "ft_48_supp_1.xlsx"), sheet = 1)
ft_48_2 <- read_xlsx(here("data", "data_to_extract", "ft_48_supp_1.xlsx"), sheet = 2)

species_names_48 <- tibble(names = names(ft_48_1 %>%
                                           select(-Total, -c("Locality", "Sampling time", "Coordinates", "Sample type"))),
                           genus = str_split(names, pattern = " ", simplify = TRUE)[, 1]) %>%
  mutate(genus = case_when(str_detect(genus, "\\.") ~ as.character(NA),
                           TRUE ~ genus)) %>%
  fill(genus, .direction = "downup") %>%
  mutate(taxonomic_level = case_when(str_detect(names, " sp.") ~ "genus",
                                     TRUE ~ "species")) %>%
  mutate(species = case_when(taxonomic_level == "species" ~ paste(genus, str_extract(names, "\\S+$")),
                             TRUE ~ names)) %>%
  select(taxonomic_level, genus, species, names)

ft_48_1_processed <- ft_48_1 %>%
  fill(Locality, .direction = "down") %>%
  fill(`Sampling time`, .direction = "down") %>%
  fill(Coordinates, .direction = "down") %>%
  mutate(lat = str_trim(str_split(Coordinates, ",", simplify = TRUE)[, 1]),
         lon = str_trim(str_split(Coordinates, ",", simplify = TRUE)[, 2])) %>%
  select(-Coordinates, -`Sample type`, -Total) %>%
  rename("date" = `Sampling time`) %>%
  mutate(start_date = str_split(date, "-", simplify = TRUE)[, 1],
         end_date = str_split(date, "-", simplify = TRUE)[, 2]) %>%
  select(-date) %>%
  mutate(start_date = parse_date(start_date, "%B %Y"),
         end_date = parse_date(end_date, "%B %Y"),
         date = case_when(!is.na(end_date) ~ paste0(start_date, "/", end_date),
                          TRUE ~ as.character(start_date))) %>%
  select(-start_date, -end_date) %>%
  relocate(Locality, date, lat, lon, .before = 1) %>%
  pivot_longer(cols = 5:59, names_to = "names", values_to = "n") %>%
  mutate(n_individuals = str_split(n, "/", simplify = TRUE)[, 2],
         n_individuals = case_when(is.na(n) ~ "0",
                                   TRUE ~ n_individuals),
         n_positive =  str_split(n, "/", simplify = TRUE)[, 1],
         n_positive = replace_na(n_positive, "0")) %>%
  left_join(species_names_48) %>%
  mutate(rodent_record_id = seq(8302, length.out = nrow(.))) %>%
  select(rodent_record_id, date, taxonomic_level, genus, species, Locality, lat, lon, n_individuals, n_positive)

write.table(ft_48_1_processed, here("data", "data_to_extract", "ft_48_processed.txt"), sep="\t", row.names=FALSE)

clean_sp_name <- function(names) {
  gen_sp <- str_split(names, " ", simplify = TRUE)[, 1:2]
  apply(gen_sp, 1, function(row) paste(row, collapse = " "))
}

ft_48_2_processed <- ft_48_2 %>%
  mutate(species = clean_sp_name(ft_48_2$`cytb identification host`)) %>%
  filter(`GenBank AN cytb` != "-") %>%
  filter(`Mamm-arenavirus screening` == "+") %>%
  select(species, Locality, Latitude, Longitude, accession = `GenBank AN cytb`, `Mamm-arenavirus screening`) %>%
  mutate(lat = as.numeric(Latitude),
         lon = as.numeric(Longitude)) %>%
  left_join(ft_48_1_processed %>%
              filter(n_individuals >= 1) %>%
              mutate(lat = as.numeric(lat),
                     lon = as.numeric(lon)) %>%
              select(rodent_record_id, species, Locality, lat, lon, n_individuals),
            by = c("species", "Locality")) %>%
  distinct(rodent_record_id, species, Locality, accession) %>%
  mutate(species = case_when(species == "Acomys " ~ "Acomys percivali",
                             TRUE ~ species)) %>%
  mutate(genus = str_split(species, " ", simplify = TRUE)[, 1]) %>%
  write.table("clipboard", sep="\t", row.names=FALSE)
# needs further work to associated each id with a single genbank record. For now just matched on species, not location
# currently only matching for virus positive


# ft_64 -------------------------------------------------------------------
library(stringdist)
ft_64_locations <- read_xlsx(here("data", "data_to_extract", "ft_64_locations.xlsx")) %>%
  mutate(locality = str_split(`Locality; district (region)`, ";", simplify = TRUE)[, 1],
         region = str_trim(str_split(`Locality; district (region)`, ";", simplify = TRUE)[, 2]))
ft_64_samples <- read_xlsx(here("data", "data_to_extract", "ft_64_samples.xlsx"))  %>%
  mutate(event_date = `Year of collection`,
         genus = case_when(str_detect(`Host speci es`, "AF|AS|AA") ~ "Apodemus",
                           str_detect(`Host speci es`, "MA") ~ "Microtus",
                           str_detect(`Host speci es`, "CG") ~ "Clethrionomys",
                           str_detect(`Host speci es`, "SM|SA") ~ "Sorex",
                           str_detect(`Host speci es`, "CL|CS") ~ "Crocidura",
                           str_detect(`Host speci es`, "NF") ~ "Neomys",
                           TRUE ~ "error"),
         species = case_when(str_detect(`Host speci es`, "AF") ~ "Apodemus flavicollis",
                             str_detect(`Host speci es`, "AS") ~ "Apodemus sylvaticus",
                             str_detect(`Host speci es`, "AA") ~ "Apodemus agrarius",
                             str_detect(`Host speci es`, "MA") ~ "Microtus arvalis",
                             str_detect(`Host speci es`, "CG") ~ "Clethrionomys glareolus",
                             str_detect(`Host speci es`, "SM") ~ "Sorex minutus",
                             str_detect(`Host speci es`, "SA") ~ "Sorex araneus",
                             str_detect(`Host speci es`, "CS") ~ "Crocidura leucodon",
                             str_detect(`Host speci es`, "CL") ~ "Crocidura suaveolens",
                             str_detect(`Host speci es`, "NF") ~ "Neomys fodiens",
                             TRUE ~ "error"),
         `Kurkino virus` = case_when(str_detect(`Virus detected, positive tissue`, "KURV") ~ 1,
                          TRUE ~ 0),
         `Tula virus` = case_when(str_detect(`Virus detected, positive tissue`, "TULV") ~ 1,
                          TRUE ~ 0),
         `Asikkala virus` = case_when(str_detect(`Virus detected, positive tissue`, "ASIV") ~ 1,
                          TRUE ~ 0),
         `Seewis virus` = case_when(str_detect(`Virus detected, positive tissue`, "SWSV") ~ 1,
                          TRUE ~ 0)) %>%
  pivot_longer(cols = c("Kurkino virus", "Tula virus", "Asikkala virus", "Seewis virus"), names_to = "Pathogen", values_to = "Positive") 

# Unique localities from the first dataset
localities_ft_64 <- unique(ft_64_locations$locality)

# Unique localities from the second dataset
localities_ft_64_samples <- unique(ft_64_samples$Locality)

# Initialize a vector to store matched localities
matched_localities <- character(length(localities_ft_64))

# Set a threshold for the maximum distance
threshold <- 2  # You can adjust this threshold as needed

# Loop through each locality in the first dataset
for (i in seq_along(localities_ft_64)) {
  # Compute string distances between the current locality and all localities in the second dataset
  distances <- stringdist(localities_ft_64[i], localities_ft_64_samples)
  
  # Check if any of the distances fall below the threshold
  if (any(distances <= threshold)) {
    # If a match is found, store the matched locality
    matched_localities[i] <- localities_ft_64_samples[which.min(distances)]
  } else {
    # If no match is found, store NA
    matched_localities[i] <- NA
  }
}

# Create a data frame with original localities and matched localities
matched_df <- data.frame(original_locality = localities_ft_64, matched_locality = matched_localities)

# Merge the matched localities dataframe with ft_64_locations to get the coordinates
merged_data <- merge(matched_df, ft_64_locations, by.x = "original_locality", by.y = "locality", all.x = TRUE)

# Merge the merged_data with ft_64_samples to add the coordinates to ft_64_samples
ft_64_final <- merge(ft_64_samples, merged_data, by.x = "Locality", by.y = "matched_locality", all.x = TRUE) %>%
  mutate(locality = paste0(Locality, ", ", region)) %>%
  select(-Locality, -region) %>%
  select(event_date, genus, species, locality, verbatimLocality = `Locality type`, Latitude, Longitude, Pathogen, Positive, `Sample code`) %>%
  arrange(Pathogen, event_date, locality, species)

writeClipboard(ft_64_final$`Sample code`)
  
write.table(ft_64_final %>%
              select(Pathogen, Positive), here("data", "data_to_extract", "ft_64_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")


# ft_66 -------------------------------------------------------------------
ft_66 <- read_csv(here("data", "data_to_extract", "ft_66_supp.csv")) %>%
  select(-`...15`)

ft_66_coords <- read_csv(here("data", "data_to_extract", "ft_66_locations.csv"))

ft_66 %>%
  rename("site" = 1) %>%
  pivot_longer(cols = !contains("site"), names_to = "species", values_to = "n") %>%
  mutate(n = replace_na(n, 0),
         genus = str_split(species, " ", simplify = TRUE)[, 1]) %>%
  left_join(ft_66_coords) %>%
  select(genus, species, site, lat, lon, n) %>%
  write.table(here("data", "data_to_extract", "ft_66_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# ft_67 -------------------------------------------------------------------

ft_67_rodent <- read_csv(here("data", "data_to_extract", "ft_67_supp_1.csv"))

ft_67_path <- read_csv(here("data", "data_to_extract", "ft_67_supp_2.csv"))

# ft_67_rodent %>%
#   mutate(season = factor(season, levels = c("spring", "summer", "autumn")),
#          site = factor(site, levels = c("Weissach", "Jeeser", "Billerbeck", "Gotha"))) %>%
#   arrange(site, year, season, locality) %>%
#   select(year, season, species, site, locality, lat, lon) %>%
#   write_csv(here("data", "data_to_extract", "ft_67_supp_1.csv"))

summarised_supp_1 <- ft_67_rodent %>%
  group_by(year, season, species, site, lat, lon, assay, pathogen) %>%
  summarise(n = sum(lt, st),
            tested = sum(tested),
            positive = sum(positive),
            te = sum(tn))

rodent_data_67 <- summarised_supp_1 %>%
  ungroup() %>%
  mutate(season = factor(season, levels = c("spring", "summer", "autumn"))) %>%
  select(year, season, species, site, lat, lon, n, te) %>%
  bind_rows(ft_67_path %>%
              ungroup() %>%
              mutate(season = factor(season, levels = c("spring", "summer", "autumn"))) %>%
              select(year, season, species, site, lat, lon, n, te) %>%
              filter(species == "Microtus agrestis") %>%
              distinct()) %>%
  arrange(species, year, season, site) %>%
  mutate(rodent_record_id = 19292 + (row_number(.) - 1)) %>%
  mutate(year = case_when(season == "spring" ~ paste0(year, "-04"),
                          season == "summer" ~ paste0(year, "-07"),
                          season == "autumn" ~ paste0(year, "-10")))

write.table(rodent_data_67, here("data", "data_to_extract", "ft_67_rodent_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

pathogen_67 <- summarised_supp_1 %>%
  ungroup() %>%
  select(year, season, species, site, lat, lon, assay, pathogen, tested, positive) %>%
  bind_rows(ft_67_path %>%
              select(year, season, species, site, lat, lon, assay, pathogen, tested, positive)) %>%
  left_join(rodent_data_67 %>%
              select(rodent_record_id, year, season, species, site) %>%
              mutate(year = as.numeric(str_split(year, "-", simplify = TRUE)[, 1])))  %>%
  mutate(year = case_when(season == "spring" ~ paste0(year, "-04"),
                          season == "summer" ~ paste0(year, "-07"),
                          season == "autumn" ~ paste0(year, "-10"))) %>%
  select(-season) %>%
  group_by(assay, pathogen) %>%
  arrange(assay, pathogen, rodent_record_id) %>%
  ungroup() %>%
  mutate(pathogen_record_id = 26909 + (row_number(.) - 1)) %>%
  relocate(pathogen_record_id, rodent_record_id, .before = 1)
  
write.table(pathogen_67, here("data", "data_to_extract", "ft_67_pathogen_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

ft_67_sequences <- read_csv(here("data", "data_to_extract", "ft_67_supp_3.csv")) %>%
  separate(`Collection\rdate`, into = c("year", "season"), sep = "\\r", remove = FALSE) %>%
  mutate(S = str_sub(S, end = 8),
         M = str_sub(M, end = 8),
         L = str_sub(L, end = 8)) %>%
  pivot_longer(cols = c(S, M, L), names_to = "segment", values_to = "accession") %>%
  filter(str_length(accession) >= 8) %>%
  select(year, 
         season,
         "species" = Host,
         "site" = `Trapping\rlocation`,
         accession,
         Isolate)  %>%
  mutate(year = case_when(season == "spring" ~ paste0(year, "-04"),
                          season == "summer" ~ paste0(year, "-07"),
                          season == "fall" ~ paste0(year, "-10"))) %>%
  select(-season) %>%
  left_join(pathogen_67 %>%
              filter(assay == "PCR" & pathogen == "TULV") %>%
              select(pathogen_record_id, rodent_record_id, year, species, site) %>%
              distinct()) %>%
  select(pathogen_record_id, rodent_record_id, species, accession)

write.table(ft_67_sequences, here("data", "data_to_extract", "ft_67_sequences_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")


# ft_71 -------------------------------------------------------------------

ft_71 <- read_csv(here("data", "data_to_extract", "ft_71_supp_1.csv")) %>%
  rename("lab_id" = 1,
         "field_id" = 2,
         "museum_id" = 3,
         "date" = 4,
         "species" = 5,
         "province" = 6,
         "sex" = 7,
         "samples" = 8,
         "arenavirus_result" = 9,
         "hantavirus_result" = 10) %>%
  mutate(lat = -24.5914,
         long = 27.6258,
         date = as_date(date, format = "%d-%b-%y"),
         clean_species = case_when(
           str_detect(species, "Mus\\s*\\(Nannomys\\)\\ssp\\.?$") ~ "Mus spp.",
           str_detect(species, "Mus\\s*\\(Nannomys\\)\\s*sp\\.?$") ~ "Mus spp.",
           str_detect(species, "Mus (Nannomys) sp..") ~ "Mus spp.",
           str_detect(species, "Mus\\s*\\(Nannomys\\)\\s*minutoides\\ss\\.l\\.?$") ~ "Mus minutoides",
           str_detect(species, "s\\.l\\.?$") ~ str_remove(species, "\\s*s\\.l\\.?"),
           str_detect(species, "Elephantulus\\s*brachyrhynchus") ~ "Elephantulus brachyrhynchus",
           str_detect(species, "Aethomys\\s*/\\s*Micaelamys\\s*sp\\.?|Aethomys\\ssp\\.?$") ~ "Aeothomys spp.",
           str_detect(species, "Aethomys ineptus s.l..") ~ "Aethomys ineptus",
           str_detect(species, "Gerbilliscus\\s*leucogaster\\s*\\.$") ~ "Gerbilliscus leucogaster",
           str_detect(species, "Steatomys sp.") ~ "Steatomys spp.",
           TRUE ~ species
           ),
         taxa = case_when(str_detect(clean_species, "spp.") ~ "genus",
                          TRUE ~ "species"),
         genus = str_split(clean_species, " ", simplify = TRUE)[, 1], 
         arenavirus_result = case_when(str_detect(arenavirus_result, "Neg") ~ 0,
                                       str_detect(arenavirus_result, "Pos") ~ 1),
         hantavirus_result = case_when(str_detect(hantavirus_result, "Neg") ~ 0,
                                       str_detect(hantavirus_result, "Pos") ~ 1),
         rodent_record_id = 19393 + row_number() - 1) %>%
  pivot_longer(cols = c("arenavirus_result", "hantavirus_result"), names_to = "family", values_to = "result") %>%
  mutate(family = case_when(family == "arenavirus_result" ~ "Arenaviridae",
                            family == "hantavirus_result" ~ "Hantaviridae")) %>%
  group_by(family) %>%
  arrange(family, rodent_record_id) %>%
  select(rodent_record_id, date, taxa, genus, clean_species, province, lat, long, family, result, lab_id)

write.table(ft_71, here("data", "data_to_extract", "ft_71_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")


# ft_m_8 ------------------------------------------------------------------

ft_m_8_1 <- read_csv(here("data", "data_to_extract", "ft_m_8_supp_1.csv"))
ft_m_8_2 <- read_csv(here("data", "data_to_extract", "ft_m_8_supp_2.csv"))
ft_m_8_locations <- read_csv(here("data", "data_to_extract", "ft_m_8_locations.csv"))

ft_m_8_c <- bind_rows(ft_m_8_1, ft_m_8_2) %>%
  mutate(Session = ym(Session)) %>%
  pivot_longer(cols = 2:37, 
               names_to = c("site", ".value"), 
               names_pattern = "(.*)\\.(.*)", 
               values_to = "value") %>%
  drop_na() %>%
  left_join(ft_m_8_locations, by = "site") %>%
  mutate(genus = "Peromyscus",
         species = "Peromyscus maniculatus",
         name = paste0(name, ", ", state)) %>%
  select("event_date" = Session, genus, species, name, lat, lon, MNA, MNI, traps) %>%
  group_by(name) %>%
  arrange(name, event_date)

write.table(ft_m_8_c, here("data", "data_to_extract", "ft_m_8_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")


# ft_105 ------------------------------------------------------------------
ft_105 <- read_xlsx(here("data", "data_to_extract", "ft_105_supp.xlsx")) %>%
  rowwise() %>%
  mutate(event_date = paste0(Year, if(Season == 1) "-06" else "-09"),
         scientificName = "Myodes glareolus",
         locality = "",
         country = "Sweden",
         verbatim_locality = paste(Provyta.1, if(Lokaltyp == "Brand") "fire area" else if(Lokaltyp == "Hygge") "clear cut" else "unburned forest"),
         Provyta.1 = as.numeric(Provyta.1),
         individual_count = 1,
         serostatus = PUUV) %>%
  arrange(Provyta.1, Year, Season) %>%
  select(event_date, scientificName, locality, country, verbatim_locality, individual_count, serostatus)

write.table(ft_105, here("data", "data_to_extract", "ft_105_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")


# ft_114 ------------------------------------------------------------------

ft_114 <- read_csv(here("data", "data_to_extract", "ft_114.csv")) %>%
  fill(Year, .direction = "down") %>%
  fill(Provinces, .direction = "down") %>%
  mutate(location = paste0(City, ", ", Provinces))

write.table(ft_114, here("data", "data_to_extract", "ft_114_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")


# ft_115 ------------------------------------------------------------------

ft_115 <- read_xlsx(here("data", "data_to_extract", "ft_115_supp.xlsx")) %>%
  rename("ID" = 1,
         "scientificName" = 2,
         "cyt_b_accession" = 4,
         "co1_accession" = 5,
         "rag1_accession" = 6,
         "locality" = 9)  %>%
  filter(!apply(., 1, function(row) any(str_detect(row, "\\*\\*")))) %>%
  filter(!str_detect(locality, "Experimenta")) %>%
  mutate(eventDate = as.character(as.Date(as.numeric(Collection_date), origin = "1899-12-30")),
         eventDate = coalesce(eventDate, Collection_date),
         decimalLatitude = as.numeric(str_extract(`Lat,Lon`, "[0-9.]+(?=[NS])")) * ifelse(str_detect(`Lat,Lon`, "S"), -1, 1),
         decimalLongitude = as.numeric(str_extract(`Lat,Lon`, "[0-9.]+(?=[EW])")) * ifelse(str_detect(`Lat,Lon`, "W"), -1, 1),
         across(everything(), ~ str_replace_all(., "\\*", "")),
         rodent_record_id = row_number() + 35068,
         individualCount = 1)

ft_115_rodent <- ft_115 %>%
  select(rodent_record_id, scientificName, eventDate, locality, country = Country, decimalLatitude, decimalLongitude, individualCount)

write.table(ft_115_rodent, here("data", "data_to_extract", "ft_115_rodent_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

ft_115_sequence <- ft_115 %>%
  select(rodent_record_id, cyt_b_accession, co1_accession, rag1_accession) %>%
  pivot_longer(cols = c("cyt_b_accession", "co1_accession", "rag1_accession"),
               names_to = "note",
               values_to = "accession_number") %>%
  arrange(note, rodent_record_id, accession_number) %>%
  filter(!str_detect(accession_number, "ND|NR")) %>%
  select(rodent_record_id, accession_number, note)

write.table(ft_115_sequence, here("data", "data_to_extract", "ft_115_sequence_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")


# ft_124 ------------------------------------------------------------------

ft_124_traps <- read_csv(here("data", "data_to_extract", "ft_124_supp_1.csv"))

ft_124_rodents <- read_csv(here("data", "data_to_extract", "ft_124_supp_2.csv")) %>%
  mutate(verbatimLocality = paste0(`4_Trap_Line`, str_pad(`5_Trap_Number`, pad = 0, width = 4)),
         scientificName = paste0(`8_Genus`, " ", `9_Species`),
         eventDate = as_date(`3_Date`, format = "%d/%m/%Y"),
         IFA_value = as.character(`24_IFA_SCORE`),
         PCR = `25_RTPCR_Status`,
         rodent_id = `1_Specimen_Number`) %>%
  select(scientificName, eventDate, verbatimLocality, IFA_value, PCR, rodent_id) %>%
  left_join(ft_124_traps %>%
              select(`Trap number`, decimalLatitude = Latitude, decimalLongitude = Longitude),
            by = c("verbatimLocality" = "Trap number")) %>%
  select(scientificName, eventDate, verbatimLocality, decimalLatitude, decimalLongitude, IFA_value, PCR, rodent_id) %>%
  mutate(rodent_record_id = 35330 + row_number()) %>%
  relocate(rodent_record_id, .before = 1)

ft_124_pathogen <- ft_124_rodents %>%
  pivot_longer(cols = c("IFA_value", "PCR"),
               names_to = "assay",
               values_to = "result") %>%
  mutate(scientificName = "",
         assay = case_when(str_detect(assay, "PCR") ~ "PCR",
                           str_detect(assay, "IFA") ~ "IFA"),
         note = case_when(as.numeric(result) > 0 ~ as.numeric(result),
                          TRUE ~ NA),
         tested = case_when(assay == "IFA" ~ 1,
                            assay == "PCR" & result == "Not Tested" ~ 0,
                            TRUE ~ 1),
         positive = case_when(as.numeric(result) > 0 ~ 1,
                              str_detect(result, "POSITIVE") ~ 1,
                              TRUE ~ 0)) %>%
  arrange(assay, rodent_record_id) %>%
  mutate(pathogen_record_id = 43694 + row_number()) %>%
  select(pathogen_record_id, associated_rodent_record_id = rodent_record_id, scientificName, assay, tested, positive, rodent_id)

ft_124_sequences <- read_xlsx(here("data", "data_to_extract", "ft_124_sequences.xlsx")) %>%
  left_join(ft_124_pathogen %>%
              filter(assay == "PCR") %>%
              select(pathogen_record_id, associated_rodent_record_id, tested, positive, rodent_id),
            by = c("isolate_id" = "rodent_id"))

ft_124_pathogen$tested[ft_124_pathogen$pathogen_record_id %in% ft_124_sequences$pathogen_record_id] <- 1
ft_124_pathogen$positive[ft_124_pathogen$pathogen_record_id %in% ft_124_sequences$pathogen_record_id] <- 1

write.table(ft_124_rodents %>%
              select(rodent_record_id, scientificName, eventDate, verbatimLocality, decimalLatitude, decimalLongitude),
            here("data", "data_to_extract", "ft_124_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ft_124_pathogen %>%
              select(pathogen_record_id, associated_rodent_record_id, scientificName, assay, tested, positive),
            here("data", "data_to_extract", "ft_124_pathogen_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ft_124_sequences %>%
              select(pathogen_record_id, associated_rodent_record_id, accession),
            here("data", "data_to_extract", "ft_124_sequence_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# ft_169 ------------------------------------------------------------------

ft_169_rodents <- read_csv(here("data", "data_to_extract", "ft_169_supp_1.csv")) %>%
  select(site, location, country, species, n, pathogen, assay, tested, positive) %>%
  mutate(combined_location = paste0(location, ", ", country)) %>%
  geocode(address = combined_location, method = "osm", lat = latitude, long = longitude)

unmatched_locations <- ft_169_rodents %>%
  filter(is.na(latitude)) %>%
  distinct(combined_location) %>%
  bind_cols(latitude = c(54.17, 52.4, 51.63, 52.96, 52.04, 50.92, 51.30, 48.69),
            longitude = c(13.26, 11.29, 14.01, 11.95, 12.57, 14.8, 13.02, 12.21))

ft_169_rodents_processed <- ft_169_rodents %>%
  left_join(unmatched_locations, by = "combined_location", suffix = c("", "_unmatched")) %>%
  mutate(latitude = coalesce(latitude, latitude_unmatched),
         longitude = coalesce(longitude, longitude_unmatched),
         rodent_record_id = rep(36093:36329, 2)) %>%
  select(-latitude_unmatched, -longitude_unmatched)

ft_169_pathogen_processed <- ft_169_rodents_processed %>%
  mutate(pathogen_record_id = 44940 + row_number())

ft_169_sequences <- read_csv(here("data", "data_to_extract", "ft_169_supp_2.csv")) %>%
  mutate(species = case_when(str_detect(sequence, "Marv") ~ "Microtus arvalis",
                             str_detect(sequence, "Magr") ~ "Microtus agrestis",
                             str_detect(sequence, "Arv") ~ "Arvicola",
                             TRUE ~ NA)) %>%
  left_join(ft_169_pathogen_processed %>%
              filter(assay == "PCR") %>%
              select(pathogen_record_id, site, species), by = c("site", "species"))

ft_169_sequences_2 <- read_csv(here("data", "data_to_extract", "ft_169_supp_3.csv"))

write.table(ft_169_rodents_processed %>%
              select(rodent_record_id, scientificName = species, locality = location, country, verbatimLocality = site, decimalLatitude = latitude, decimalLongitude = longitude, individualCount = n),
            here("data", "data_to_extract", "ft_169_rodent_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

write.table(ft_169_pathogen_processed %>%
              select(pathogen_record_id, rodent_record_id, scientificName = pathogen, assay, tested, positive),
            here("data", "data_to_extract", "ft_169_pathogen_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

write.table(ft_169_sequences %>%
              select(pathogen_record_id, accession) %>%
              bind_rows(ft_169_sequences_2 %>%
                          select(accession)),
            here("data", "data_to_extract", "ft_169_sequences_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# ft_177 ------------------------------------------------------------------

ft_177_raw <- read_csv(here("data", "data_to_extract", "ft_177_t1.csv")) %>%
  select(1:9) %>%
  drop_na(species)

ft_177_separated <- ft_177_raw %>%
  # Using `separate` to split each `fields` and `residential` column
  separate(fields_2008, into = c("fields_2008_tested", "fields_2008_positive"), sep = "/", convert = TRUE) %>%
  separate(residential_2008, into = c("residential_2008_tested", "residential_2008_positive"), sep = "/", convert = TRUE) %>%
  separate(fields_2009, into = c("fields_2009_tested", "fields_2009_positive"), sep = "/", convert = TRUE) %>%
  separate(residential_2009, into = c("residential_2009_tested", "residential_2009_positive"), sep = "/", convert = TRUE) %>%
  separate(fields_2010, into = c("fields_2010_tested", "fields_2010_positive"), sep = "/", convert = TRUE) %>%
  separate(residential_2010, into = c("residential_2010_tested", "residential_2010_positive"), sep = "/", convert = TRUE) %>%
  separate(fields_2011, into = c("fields_2011_tested", "fields_2011_positive"), sep = "/", convert = TRUE) %>%
  separate(residential_2011, into = c("residential_2011_tested", "residential_2011_positive"), sep = "/", convert = TRUE)

# Step 2: Pivot the data to a long format
ft_177_long <- ft_177_separated %>%
  # Pivoting longer to create the `location`, `eventDate`, `tested`, and `positive` columns
  pivot_longer(
    cols = -species,  # All columns except `species`
    names_to = c("location", "eventDate", ".value"),
    names_pattern = "(fields|residential)_(\\d{4})_(tested|positive)"
  ) %>%
  # Renaming the `eventDate` column to align with the example structure
  rename(location = location, eventDate = eventDate, tested = tested, positive = positive) %>%
  mutate(species = case_when(str_detect(species, "R. ") ~ str_replace(species, "R. ", "Rattus "),
                             str_detect(species, "N. ") ~ str_replace(species, "N. ", "Niviventer "),
                             TRUE ~ species),
         eventDate = as.integer(eventDate)) %>%
  # Arrange for better readability (optional)
  arrange(eventDate, location, species) %>%
  mutate(rodent_record_id = 36445 + row_number(),
         pathogen_record_id = 45572 + row_number())

ft_177_sequences <- read_csv(here("data", "data_to_extract", "ft_177_s2.csv")) %>%
  mutate(s = str_trim(str_remove_all(s, "\\(partial\\)")),
         m = str_trim(str_remove_all(m, "\\(partial\\)")),
         location = "Longquan",
         virus = str_remove_all(virus, "Longquan"),
         species = str_extract(virus, "^[^-]+"),
         species = case_when(str_detect(species, "Aa") ~ "Apodemus agrarius",
                             str_detect(species, "Hu") ~ "Homo sapiens",
                             str_detect(species, "Mf") ~ "Microtus fortis",
                             str_detect(species, "Rf") ~ "Rattus flavipectus",
                             str_detect(species, "Rn") ~ "Rattus norvegicus"),
         year = if_else(
           str_detect(virus, "-\\d{2}-"),                 # Check if formatted like Aa-09-98
           str_extract(virus, "(?<=-)(\\d{2})(?=-)"),     # Extract 2 digits between `-`
           str_extract(virus, "\\d{4}")                   # Otherwise, extract the first 4 digits
         ),
         year = case_when(str_detect(year, "08") ~ as.integer(2008),
                          str_detect(year, "09") ~ as.integer(2009),
                          str_detect(year, "10") ~ as.integer(2010),
                          str_detect(year, "11") ~ as.integer(2011),
                          TRUE ~ as.integer(year))) %>%
  left_join(ft_177_long %>%
              filter(positive >= 1),
            by = c("species", "year" = "eventDate"))

write.table(ft_177_long %>%
              select(rodent_record_id, pathogen_record_id, species, eventDate, location, tested, positive),
            here("data", "data_to_extract", "ft_177_rodent_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(ft_177_sequences %>%
              select(rodent_record_id, pathogen_record_id, species, year, s, m),
            here("data", "data_to_extract", "ft_177_sequences_processed.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
