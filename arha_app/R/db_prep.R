# ArHa_App/R/db_prep.R
library(DBI)
library(duckdb)
library(dplyr)

# Create a transient DuckDB connection
con <- dbConnect(duckdb(), dbdir = ":memory:")

DBI::dbExecute(con, "SET enable_progress_bar = false;")

# Register Parquet files as DuckDB views
dbExecute(con, "CREATE VIEW events AS SELECT * FROM read_parquet('data/parquet/sampling_events.parquet')")
dbExecute(con, "CREATE VIEW hosts AS SELECT * FROM read_parquet('data/parquet/host_occurrences.parquet')")
dbExecute(con, "CREATE VIEW pathogens AS SELECT * FROM read_parquet('data/parquet/pathogen_mof.parquet')")
dbExecute(con, "CREATE VIEW sequences AS SELECT * FROM read_parquet('data/parquet/resource_relationships.parquet')")

# dbplyr lazy tables
tbl_events <- tbl(con, "events")
tbl_hosts <- tbl(con, "hosts")
tbl_pathogens <- tbl(con, "pathogens")
tbl_sequences <- tbl(con, "sequences")
