source(here::here("R", "00_load_data.R"))
combined_data <- read_rds(here("data", "clean_data", "2025-03-05_data.rds"))

# Additional packages
pkgs <- c(
  "cowplot",
  "terra",
  "geodata",
  "rnaturalearth",
  "sf"
)

pacman::p_load(pkgs, character.only = T)

a <- combined_data$descriptive %>%
  mutate(data_extractor = str_split(study_id, "_", simplify = TRUE)[, 1])

table(a$data_extractor)

studies_for_grant_to_check <- a %>%
  filter(!str_detect(data_extractor, "grant")) %>%
  mutate(weight = n_rodent_individuals + n_pathogen_records) %>%
  filter(!is.na(weight), weight > 0) %>%
  filter(!study_id %in% c("anna_133", "anna_110", "anna_200", "anna_205", "anna_207", "anna_211", "anna_212", "anna_221", "anna_223",
                          "harry_103", "harry_102", "harry_100", "harry_104", "harry_109", "harry_111", "harry_116", "harry_119", "harry_122", "harry_126")) %>%
  group_by(data_extractor) %>%
  slice_sample(prop = 0.10, weight_by = weight) %>%
  select(study_id, full_text_id)

write_csv(studies_for_grant_to_check, here("misc", "studies_for_double_extraction_grant.csv"))
