library(DBI)
library(duckdb)
library(dplyr)
library(countrycode)

con <- dbConnect(duckdb(), dbdir = ":memory:")
DBI::dbExecute(con, "SET enable_progress_bar = false;")

dbExecute(con, "CREATE VIEW events AS SELECT * FROM read_parquet('data/parquet/sampling_events.parquet')")
dbExecute(con, "CREATE VIEW hosts AS SELECT * FROM read_parquet('data/parquet/host_occurrences.parquet')")
dbExecute(con, "CREATE VIEW pathogens AS SELECT * FROM read_parquet('data/parquet/pathogen_mof.parquet')")
dbExecute(con, "CREATE VIEW sequences AS SELECT * FROM read_parquet('data/parquet/resource_relationships.parquet')")

tbl_events    <- tbl(con, "events")
tbl_hosts_raw <- tbl(con, "hosts")
tbl_pathogens <- tbl(con, "pathogens")
tbl_sequences <- tbl(con, "sequences")

# --- Country name harmonisation ---
# host_occurrences$country is derived from country_processed, which is not
# guaranteed consistent. countryCode (ISO3C) is consistent
country_lookup <- tbl_hosts_raw |>
  distinct(countryCode) |>
  collect() |>
  filter(!is.na(countryCode)) |>
  mutate(country_clean = countrycode::countrycode(countryCode, origin = "iso3c", destination = "country.name"))

dbWriteTable(con, "country_lookup", country_lookup, overwrite = TRUE)

tbl_hosts <- tbl_hosts_raw |>
  left_join(tbl(con, "country_lookup"), by = "countryCode") |>
  mutate(country = coalesce(country_clean, country)) |>
  select(-country_clean)