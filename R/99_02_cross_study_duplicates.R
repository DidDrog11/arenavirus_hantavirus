# Project ArHa: cross-study duplicate check
# 99_02_cross_study_duplicates.R
# Re-runs the count-fingerprint check from data/duplicate_check/ (sent 2026-09-10,
# probably run on an older release) against the current database, tightened to
# require the same country and overlapping sampling dates, so that each flagged
# pair can be judged against the source papers. Read-only: writes a pair table to
# data/temp/ and a review list to data/duplicate_check/, and changes nothing in
# the pipeline.

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(here)

db_file <- sort(list.files(here("data", "database"), pattern = "\\.rds$", full.names = TRUE)) |> tail(1)
arha_db <- read_rds(db_file)
message("Database: ", basename(db_file))

flagged <- read_csv(here("data", "duplicate_check", "ArHa_duplicate_records.csv"), show_col_types = FALSE)

# 1. Fingerprint each assay record ---------------------------------------------
# The collaborator's definition (host species, pathogen, number positive, number
# tested) plus country (ISO3, as country names are not harmonised). Host species,
# country and dates come from the host table, joined on host_record_id.
# Zero-denominator records are excluded, as they carry no counts.
fingerprints <- arha_db$pathogen |>
  filter(number_tested > 0) |>
  left_join(arha_db$host |> select(host_record_id, host_species, iso3c, start_date, end_date), by = "host_record_id") |>
  # Unresolved host taxa would otherwise match each other across unrelated genera
  filter(!is.na(host_species)) |>
  # study_id is a factor in the database; character so pairs can be ordered below
  mutate(study_id = as.character(study_id),
         fp_key = paste(iso3c, host_species, pathogen_species_cleaned, number_positive, number_tested, sep = "|")) |>
  distinct(study_id, fp_key, host_species, number_positive, number_tested, start_date, end_date)

n_fp_per_study <- fingerprints |> group_by(study_id) |> summarise(n_fp = n_distinct(fp_key))

# 2. Pairs of studies sharing at least one fingerprint -------------------------
# A match also needs the two records' sampling windows to overlap. Where either
# record has no dates the match is kept, as it cannot be ruled out.
shared <- fingerprints |>
  inner_join(fingerprints |> select(-host_species, -number_positive, -number_tested),
             by = "fp_key", suffix = c("_1", "_2"), relationship = "many-to-many") |>
  filter(study_id_1 < study_id_2) |>
  filter(is.na(start_date_1) | is.na(end_date_1) | is.na(start_date_2) | is.na(end_date_2) |
           (start_date_1 <= end_date_2 & start_date_2 <= end_date_1)) |>
  distinct(study_id_1, study_id_2, fp_key, host_species, number_positive, number_tested)

pairs <- shared |>
  group_by(study_id_1, study_id_2) |>
  summarise(shared_records = n(),
            shared_positive = sum(number_positive > 0),
            shared_animals = sum(number_tested),
            largest_shared = max(number_tested),
            example = paste0(host_species, " ", number_positive, "/", number_tested)[which.max(number_tested)],
            .groups = "drop") |>
  left_join(n_fp_per_study |> rename(study_id_1 = study_id, n_fp_1 = n_fp), by = "study_id_1") |>
  left_join(n_fp_per_study |> rename(study_id_2 = study_id, n_fp_2 = n_fp), by = "study_id_2") |>
  mutate(match_type = if_else(shared_records == n_fp_1 & shared_records == n_fp_2, "EXACT", "PARTIAL"),
         # The collaborator's TRIVIAL rule: one shared record of one or two animals
         trivial = shared_records == 1 & largest_shared <= 2) |>
  arrange(match_type, desc(shared_animals))

# Add citation details for both sides
cites <- arha_db$citations |>
  distinct(study_id, .keep_all = TRUE) |>
  transmute(study_id = as.character(study_id), full_text_id = as.character(full_text_id), doi, first_author = str_extract(author, "^[^,;]+"),
            publication_year, title = str_trunc(title, 70))

pairs <- pairs |>
  left_join(cites |> rename_with(~ paste0(.x, "_1")), by = "study_id_1") |>
  left_join(cites |> rename_with(~ paste0(.x, "_2")), by = "study_id_2") |>
  mutate(same_full_text = full_text_id_1 == full_text_id_2,
         same_doi = !is.na(doi_1) & doi_1 == doi_2)

message("Pairs sharing >= 1 fingerprint: ", nrow(pairs),
        " | EXACT: ", sum(pairs$match_type == "EXACT"),
        " | PARTIAL non-trivial: ", sum(pairs$match_type == "PARTIAL" & !pairs$trivial))

# Pairs to check against the papers: at least one shared record with positives
# and more than two animals, or three or more shared records of any kind. Matches
# on single all-negative records are common by chance and are left out.
review_pairs <- pairs |>
  filter((shared_positive >= 1 & largest_shared > 2) | shared_records >= 3) |>
  arrange(desc(shared_positive), desc(shared_animals))

message("Pairs for review: ", nrow(review_pairs))

review_pairs |>
  select(study_id_1, study_id_2, match_type, shared_records, shared_positive, n_fp_1, n_fp_2, shared_animals, example,
         first_author_1, publication_year_1, first_author_2, publication_year_2, same_doi) |>
  print(n = Inf, width = Inf)

# 3. Status of the collaborator's 19 pairs in the current database -------------
# A pair can be missing now because one copy was dropped as a non-canonical
# extraction of the same paper (02_01_clean_descriptive.R, added 2026-08-17).
current_ids <- unique(as.character(arha_db$host$study_id))

flagged_status <- flagged |>
  select(study_1, study_2, match_type_sent = match_type, same_doi) |>
  mutate(in_db_1 = study_1 %in% current_ids,
         in_db_2 = study_2 %in% current_ids,
         key_1 = pmin(study_1, study_2),
         key_2 = pmax(study_1, study_2)) |>
  left_join(pairs |> select(study_id_1, study_id_2, match_type_now = match_type, shared_records),
            by = c("key_1" = "study_id_1", "key_2" = "study_id_2")) |>
  select(-key_1, -key_2)

cat("\nCollaborator's pairs, then and now:\n")
print(flagged_status, n = Inf, width = Inf)

# 4. Same publication still entered under two study_ids ------------------------
# Should be empty if the canonical-extraction step is working.
cat("\nfull_text_ids with data under more than one study_id:\n")
cites |>
  filter(study_id %in% current_ids) |>
  add_count(full_text_id) |>
  filter(n > 1) |>
  print(n = Inf)

# 5. Side-by-side summary of each non-trivial pair -----------------------------
# Dates, places and assays are what separate "same animals entered twice" from
# "the same team reporting different sessions with coincident counts".
linked <- read_rds(here("data", "data_cleaning", "02_01_output.rds"))$descriptives_cleaned |>
  select(study_id, linked_manuscripts_raw)

study_summary <- function(id) {
  h <- arha_db$host |> filter(study_id == id)
  p <- arha_db$pathogen |> filter(study_id == id)
  cite <- cites |> filter(study_id == id)
  tibble(study_id = id,
         citation = paste(first(cite$first_author), first(cite$publication_year), first(cite$doi)),
         countries = paste(sort(unique(h$country)), collapse = "; "),
         dates = paste(min(h$start_date, na.rm = TRUE), "to", max(h$end_date, na.rm = TRUE)),
         n_localities = n_distinct(h$locality_raw),
         lat_range = paste(round(range(h$latitude, na.rm = TRUE), 2), collapse = " to "),
         lon_range = paste(round(range(h$longitude, na.rm = TRUE), 2), collapse = " to "),
         host_records = nrow(h),
         animals = sum(h$number_of_hosts, na.rm = TRUE),
         pathogen_families = paste(sort(unique(p$pathogen_family)), collapse = "; "),
         assays = paste(sort(unique(p$assay)), collapse = "; "),
         linked = paste(linked$linked_manuscripts_raw[linked$study_id == id], collapse = "; "))
}

for (i in seq_len(nrow(review_pairs))) {
  cat("\n====", review_pairs$study_id_1[i], "vs", review_pairs$study_id_2[i],
      "|", review_pairs$match_type[i], "|", review_pairs$shared_records[i], "shared records\n")
  bind_rows(study_summary(review_pairs$study_id_1[i]), study_summary(review_pairs$study_id_2[i])) |>
    pivot_longer(-study_id, values_transform = as.character) |>
    pivot_wider(names_from = study_id, values_from = value) |>
    print(width = Inf)
}

# 6. Specific follow-ups from the collaborator's notes -------------------------
# Londono et al. 2011 (Colombia): its citation once carried the Czech Tula data.
# Check whether its own data were ever extracted.
cat("\nLondono 2011 (10.1089/vbz.2010.0129):\n")
londono <- arha_db$citations |> filter(str_detect(tolower(doi), "vbz\\.2010\\.0129"))
print(londono |> select(study_id, full_text_id, doi, decision, extractor))
cat("Host rows under these study_ids:", sum(arha_db$host$study_id %in% londono$study_id), "\n")

# grant_189 is cited to an arenavirus review chapter but holds a Sin Nombre record.
cat("\ngrant_189 pathogen records:\n")
arha_db$pathogen |>
  filter(study_id == "grant_189") |>
  count(pathogen_family, pathogen_species_cleaned, number_positive, number_tested) |>
  print()

# 7. Matching study totals -----------------------------------------------------
# The record-level check misses the same animals reported at a different
# breakdown: Klein 2012 reports Firing Points 10 and 60 by year, O'Guinn 2008 by
# site, and no single record matches although every species total does. Here
# records are summed per study, country, host species and pathogen, and studies
# are paired where those totals match and the study windows overlap.
totals <- arha_db$pathogen |>
  filter(number_tested > 0) |>
  left_join(arha_db$host |> select(host_record_id, host_species, iso3c, start_date, end_date), by = "host_record_id") |>
  filter(!is.na(host_species)) |>
  mutate(study_id = as.character(study_id)) |>
  group_by(study_id, iso3c, host_species, pathogen_species_cleaned) |>
  summarise(total_positive = sum(number_positive), total_tested = sum(number_tested),
            start = min(start_date, na.rm = TRUE), end = max(end_date, na.rm = TRUE), .groups = "drop") |>
  mutate(fp_key = paste(iso3c, host_species, pathogen_species_cleaned, total_positive, total_tested, sep = "|"))

n_totals_per_study <- totals |> count(study_id, name = "n_totals")

shared_totals <- totals |>
  inner_join(totals |> select(study_id, fp_key, start, end), by = "fp_key", suffix = c("_1", "_2"),
             relationship = "many-to-many") |>
  filter(study_id_1 < study_id_2) |>
  filter(is.infinite(start_1) | is.infinite(start_2) | (start_1 <= end_2 & start_2 <= end_1))

totals_pairs <- shared_totals |>
  group_by(study_id_1, study_id_2) |>
  summarise(shared_totals = n(),
            shared_positive = sum(total_positive > 0),
            shared_animals = sum(total_tested),
            largest_shared = max(total_tested),
            example = paste0(host_species, " ", total_positive, "/", total_tested)[which.max(total_tested)],
            .groups = "drop") |>
  left_join(n_totals_per_study |> rename(study_id_1 = study_id, n_totals_1 = n_totals), by = "study_id_1") |>
  left_join(n_totals_per_study |> rename(study_id_2 = study_id, n_totals_2 = n_totals), by = "study_id_2") |>
  filter((shared_positive >= 1 & largest_shared > 2) | shared_totals >= 3) |>
  # Flag pairs already on the record-level review list
  left_join(review_pairs |> transmute(study_id_1, study_id_2, on_record_list = TRUE), by = c("study_id_1", "study_id_2")) |>
  mutate(on_record_list = coalesce(on_record_list, FALSE)) |>
  left_join(cites |> select(study_id, full_text_id, first_author, publication_year) |> rename_with(~ paste0(.x, "_1")),
            by = "study_id_1") |>
  left_join(cites |> select(study_id, full_text_id, first_author, publication_year) |> rename_with(~ paste0(.x, "_2")),
            by = "study_id_2") |>
  arrange(on_record_list, desc(shared_positive), desc(shared_animals))

message("Pairs with matching study totals: ", nrow(totals_pairs),
        " | not already on the record-level list: ", sum(!totals_pairs$on_record_list))

totals_pairs |>
  filter(!on_record_list) |>
  select(study_id_1, study_id_2, shared_totals, shared_positive, n_totals_1, n_totals_2, shared_animals, example,
         first_author_1, publication_year_1, first_author_2, publication_year_2) |>
  print(n = Inf, width = Inf)

totals_pairs |>
  select(study_id_1, full_text_id_1, first_author_1, publication_year_1,
         study_id_2, full_text_id_2, first_author_2, publication_year_2,
         on_record_list, shared_totals, shared_positive, n_totals_1, n_totals_2, shared_animals, example) |>
  mutate(verdict = NA_character_, keep = NA_character_, notes = NA_character_) |>
  write_csv(here("data", "duplicate_check", "review_pairs_totals.csv"), na = "")

write_csv(pairs, here("data", "temp", "cross_study_duplicate_pairs.csv"))

# Review list, with full-text IDs so the papers can be pulled, and empty columns
# for recording the decision on each pair
review_pairs |>
  select(study_id_1, full_text_id_1, first_author_1, publication_year_1,
         study_id_2, full_text_id_2, first_author_2, publication_year_2,
         match_type, shared_records, shared_positive, n_fp_1, n_fp_2, shared_animals, example) |>
  mutate(verdict = NA_character_, keep = NA_character_, notes = NA_character_) |>
  write_csv(here("data", "duplicate_check", "review_pairs.csv"), na = "")
message("Written: data/temp/cross_study_duplicate_pairs.csv, data/duplicate_check/review_pairs.csv and review_pairs_totals.csv")
