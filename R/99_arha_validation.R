# ArHa release validation
# Run against the recompiled parquet files, outside the Shiny app.
# Reports the figures that go into the data descriptor and checks the
# referential-integrity claims made in Technical Validation.

library(arrow)
library(dplyr)
library(stringr)
library(tidyr)

parquet_dir <- here::here("arha_app", "data", "parquet")

sampling_events        <- read_parquet(file.path(parquet_dir, "sampling_events.parquet"))
host_occurrences       <- read_parquet(file.path(parquet_dir, "host_occurrences.parquet"))
pathogen_mof           <- read_parquet(file.path(parquet_dir, "pathogen_mof.parquet"))
resource_relationships <- read_parquet(file.path(parquet_dir, "resource_relationships.parquet"))

rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")


# 1. TABLE COUNTS -------------------------------------------------------------
rule("1. Table row counts")

tibble(
  table = c("sampling_events", "host_occurrences", "pathogen_mof", "resource_relationships"),
  rows  = c(nrow(sampling_events), nrow(host_occurrences),
            nrow(pathogen_mof), nrow(resource_relationships))
) |> print()

cat("\nIndividual animals (sum individualCount):",
    format(sum(host_occurrences$individualCount, na.rm = TRUE), big.mark = ","), "\n")
cat("Individuals screened (sum number_tested):",
    format(sum(pathogen_mof$number_tested, na.rm = TRUE), big.mark = ","), "\n")


# 2. REFERENTIAL INTEGRITY ----------------------------------------------------
rule("2. Referential integrity")

null_key <- pathogen_mof |> filter(is.na(occurrenceID))
broken_key <- pathogen_mof |>
  filter(!is.na(occurrenceID)) |>
  anti_join(host_occurrences, by = "occurrenceID")

cat("Assay records with NULL host key:", nrow(null_key), "\n")
if (nrow(null_key) > 0) print(null_key |> count(eventID, sort = TRUE))

cat("\nAssay records with BROKEN host key (expected 0):", nrow(broken_key), "\n")
if (nrow(broken_key) > 0) print(broken_key |> count(eventID, sort = TRUE))

# Guard against a literal "_NA" string surviving the recode
cat("\nLiteral '*_NA' strings still present (expected 0):",
    sum(str_detect(pathogen_mof$occurrenceID, "_NA$"), na.rm = TRUE), "\n")

# Compound key -- catches a valid occurrenceID paired with a mismatched eventID
cat("Broken on compound key (expected 0):",
    pathogen_mof |> filter(!is.na(occurrenceID)) |>
      anti_join(host_occurrences, by = c("eventID", "occurrenceID")) |> nrow(), "\n")

# Sequences: resourceID is polymorphic (host OR assay)
valid_ids <- c(host_occurrences$occurrenceID, pathogen_mof$measurementID)
cat("Sequences with unresolvable resourceID (expected 0):",
    sum(!resource_relationships$resourceID %in% valid_ids), "\n")

cat("Resource relationships: total rows vs distinct accessions:",
    nrow(resource_relationships), "vs",
    n_distinct(resource_relationships$accession_primary), "\n")

# Host -> event
cat("Host records with unresolvable eventID (expected 0):",
    host_occurrences |> anti_join(sampling_events, by = "eventID") |> nrow(), "\n")

cat("Non-integer test/positive counts:", sum(pathogen_mof$non_integer, na.rm = TRUE), "\n")
cat("Tested exceeds host individualCount:", sum(pathogen_mof$tested_detected, na.rm = TRUE), "\n")

bad_strings <- c("Chaeotdipus", "Peromycus eremicus", "Baimoys taylori")
cat("Uncorrected taxonomy strings still present:",
    sum(host_occurrences$genus %in% bad_strings | host_occurrences$scientificName %in% bad_strings), "\n")

# 3. THE CUYPERS FIX ----------------------------------------------------------
# EXPECTED: 1,226 tested = 1,225 study-level screen + 1 site-level L. striatus
#           6 positives (G. surdaster 2, M. minutoides 3, L. rosalia 1)
rule("3. event_0043 (Cuypers et al. 2023)")

pathogen_mof |>
  filter(eventID == "event_0043") |>
  summarise(
    assay_records = n(),
    tested        = sum(number_tested, na.rm = TRUE),
    positives     = sum(measurementValue, na.rm = TRUE),
    null_key      = sum(is.na(occurrenceID))
  ) |> print()

cat("\nPositives by species (expect surdaster 2, minutoides 3, rosalia 1):\n")
pathogen_mof |>
  filter(eventID == "event_0043", measurementValue > 0) |>
  left_join(host_occurrences |> select(occurrenceID, scientificName), by = "occurrenceID") |>
  select(scientificName, number_tested, measurementValue) |> print()

# 4. ZERO-DENOMINATOR RECORDS -------------------------------------------------
# Recompute for the Data Records paragraph -- the 21,640 / 79.7% figures move.
rule("4. Zero-denominator assay records")

zero_den <- pathogen_mof |> filter(number_tested == 0)
cat("Records with number_tested == 0:", nrow(zero_den),
    sprintf("(%.1f%% of assay records)\n", 100 * nrow(zero_den) / nrow(pathogen_mof)))

zero_den |>
  left_join(host_occurrences |> select(occurrenceID, individualCount), by = "occurrenceID") |>
  mutate(parent = case_when(
    is.na(occurrenceID)     ~ "no host key",
    is.na(individualCount)  ~ "host individualCount NA",
    individualCount == 0    ~ "host individualCount = 0 (nothing caught)",
    TRUE                    ~ "host individualCount > 0 (caught, not screened)"
  )) |>
  count(parent, sort = TRUE) |>
  mutate(pct = round(100 * n / sum(n), 1)) |> print()

# Logical impossibility -- positives without a denominator
cat("\nRecords with 0 tested but >0 positive (expected 0):",
    pathogen_mof |> filter(number_tested == 0, measurementValue > 0) |> nrow(), "\n")
cat("Records with positives > tested:",
    pathogen_mof |> filter(measurementValue > number_tested) |> nrow(), "\n")


# 5. QUALITY METRICS ----------------------------------------------------------
# Denominator is HOST RECORDS, not the joined table -- these are the figures
# that belong in the paper and that the app's value boxes now report.
rule("5. Quality metrics (denominator = host records)")

n_hosts <- nrow(host_occurrences)
host_occurrences |>
  summarise(
    host_records     = n_hosts,
    precise_spatial  = sum(coordinate_resolution_processed == "site", na.rm = TRUE),
    precise_temporal = sum(temporal_resolution %in% c("full_date", "day_range_resolution"), na.rm = TRUE),
    species_level    = sum(taxonRank == "species", na.rm = TRUE)
  ) |>
  mutate(across(-host_records, ~ sprintf("%s (%.1f%%)", format(.x, big.mark = ","), 100 * .x / n_hosts))) |>
  glimpse()

cat("\nCoordinate resolution breakdown:\n")
host_occurrences |> count(coordinate_resolution_processed, sort = TRUE) |> print()

# The 55 re-keyed Cuypers records must NOT be counted as site-level precision
cat("\nCoordinate resolution for event_0043 host records:\n")
host_occurrences |> filter(eventID == "event_0043") |>
  count(coordinate_resolution_processed) |> print()

# 6. DUPLICATE-BLOCK SWEEP ----------------------------------------------------
# Looks for the Uzia signature: a long species list at one site-date where most
# taxa appear both as a zero row and a non-zero row. Repeat trapping sessions
# vary between sessions; a near-total zero/non-zero pairing does not.
rule("6. Possible duplicated site blocks")

host_occurrences |>
  group_by(eventID, decimalLatitude, decimalLongitude, eventDate, verbatimLocality, scientificName) |>
  summarise(paired = n() > 1 &&
              any(individualCount == 0, na.rm = TRUE) &&
              any(individualCount > 0, na.rm = TRUE),
            .groups = "drop") |>
  group_by(eventID, decimalLatitude, decimalLongitude, eventDate) |>
  summarise(n_species = n(), prop_paired = mean(paired), .groups = "drop") |>
  filter(n_species >= 10, prop_paired > 0.8) |>
  arrange(desc(n_species)) |>
  print(n = 20)

cat("(Empty = no candidates. Any hit needs checking against the source paper.)\n")


# 7. number_negative DERIVATION ----------------------------------------------
rule("7. Inconclusive results")

cat("Records with number_inconclusive > 0:",
    pathogen_mof |> filter(number_inconclusive > 0) |> nrow(), "\n")
cat("...of which negatives are overcounted:",
    pathogen_mof |> filter(number_inconclusive > 0,
                           number_negative != number_tested - measurementValue - number_inconclusive) |> nrow(), "\n")


# 8. TOP CONTRIBUTORS ---------------------------------------------------------
# For manual verification against source manuscripts.
rule("8. Largest contributing studies")

host_occurrences |>
  count(eventID, name = "host_records", sort = TRUE) |>
  head(20) |>
  left_join(
    pathogen_mof |> group_by(eventID) |>
      summarise(tested = sum(number_tested, na.rm = TRUE),
                positives = sum(measurementValue, na.rm = TRUE), .groups = "drop"),
    by = "eventID") |>
  left_join(host_occurrences |> group_by(eventID) |>
              summarise(individuals = sum(individualCount, na.rm = TRUE), .groups = "drop"),
            by = "eventID") |>
  left_join(sampling_events |> select(eventID, associatedReferences), by = "eventID") |>
  mutate(reference = str_trunc(associatedReferences, 60)) |>
  select(-associatedReferences) |>
  head(20) |>
  slice_sample(n = 4)

rule("Done")
