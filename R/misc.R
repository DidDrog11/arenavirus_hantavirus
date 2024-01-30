rodent_data %>%
  filter(study_id == 15) %>%
  group_by(scientificName, country) %>%
  summarise(n = sum(as.numeric(individualCount), na.rm = TRUE)) %>%
  arrange(country, scientificName) %>%
  mutate(genus = str_split(scientificName, " ", simplify = TRUE)[1]) %>%
  select(n) %>%
  write.table(., "clipboard", sep="\t", row.names=FALSE)


# Extracting data from ft_12 ----------------------------------------------
library(readxl)
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
  