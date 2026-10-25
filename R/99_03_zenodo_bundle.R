# Assemble the files for a Zenodo release in data/data_release/zenodo/<version>/.
#
# The files are uploaded by hand as a new version of the existing Zenodo record
# (concept DOI 10.5281/zenodo.21978982). Run after 07 and 08, once the release commit
# has been tagged with `version`. Tag with `git tag` and push the tag; do not create a
# GitHub release, which would also trigger the Zenodo-GitHub integration and add a
# second, automatic version.
#
# Files written:
#   ArHa_<version>_tables_parquet.zip   the four relational tables as Parquet
#   ArHa_<version>_tables_csv.zip       the same tables as CSV
#   ArHa_<version>_dwca.zip             Darwin Core Archive
#   ArHa_<version>_data_dictionary.csv  field definitions
#   ArHa_<version>_pipeline.zip         scripts, app, corrections and matching tables at the
#                                       tagged commit, plus the extraction data snapshot

library(here)

version <- "v1.2"               # git tag of the release commit
ref <- version                  # set to "HEAD" only for a draft bundle before tagging
analysis_date <- "2026-10-20"   # date of the extraction data snapshot, as in 00_load_data.R

release_dir <- here("data", "data_release")
out_dir <- here("data", "data_release", "zenodo", version)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- function(name) file.path(out_dir, paste0("ArHa_", version, "_", name))

# The deposited data must match the tagged code: stop if the tag is missing or the
# tracked release files differ from it
git <- function(...) system2("git", c("-C", shQuote(here()), ...), stdout = FALSE, stderr = FALSE)
if (ref != "HEAD" && git("rev-parse", "--verify", "--quiet", paste0("refs/tags/", ref)) != 0) {
  stop("Tag ", ref, " does not exist. Commit, tag the release commit, then rerun.")
}
if (git("diff", "--quiet", ref, "--", "data/data_release", "arha_app/data") != 0) {
  stop("Release files in data/data_release or arha_app/data differ from ", ref, ".")
}

# Tables, Darwin Core Archive and data dictionary
tables <- c("sampling_events", "host_occurrences", "pathogen_mof", "resource_relationships")
zip::zip(out_file("tables_parquet.zip"), paste0(tables, ".parquet"), root = file.path(release_dir, "parquet"))
zip::zip(out_file("tables_csv.zip"), paste0(tables, ".csv"), root = file.path(release_dir, "csv"))
file.copy(file.path(release_dir, "ArHa_dwca.zip"), out_file("dwca.zip"), overwrite = TRUE)
file.copy(file.path(release_dir, "data_dictionary.csv"), out_file("data_dictionary.csv"), overwrite = TRUE)

# Pipeline: the files needed to rebuild the database, taken from the tagged commit.
# Old database versions and intermediate files in the repository are left out.
pipeline_paths <- c("R", "arha_app", "data/corrections", "data/matching", "data/duplicate_resolution.csv",
                    "data/data_release/data_dictionary.csv", "README.md", "LICENSE", "CITATION.cff", ".here")
prefix <- paste0("arenavirus_hantavirus-", version)
pipeline_zip <- out_file("pipeline.zip")
if (file.exists(pipeline_zip)) file.remove(pipeline_zip)
status <- system2("git", c("-C", shQuote(here()), "archive", "--format=zip", paste0("--prefix=", prefix, "/"),
                           "-o", shQuote(pipeline_zip), ref, "--", pipeline_paths))
if (status != 0) stop("git archive failed")

# The extraction data snapshot is not tracked in git, so it is added to the pipeline zip
# from data/raw_data/ (one CSV per sheet, plus the column types)
snapshot <- list.files(here("data", "raw_data"), pattern = paste0("^", analysis_date, "_.*[.]csv$"))
staging <- file.path(tempdir(), "zenodo_snapshot")
dir.create(file.path(staging, prefix, "data", "raw_data"), recursive = TRUE, showWarnings = FALSE)
file.copy(here("data", "raw_data", snapshot), file.path(staging, prefix, "data", "raw_data"), overwrite = TRUE)
zip::zip_append(pipeline_zip, file.path(prefix, "data", "raw_data", snapshot), root = staging)

# Summary of what was written, with checksums to compare against Zenodo after upload
written <- list.files(out_dir, full.names = TRUE)
data.frame(file = basename(written), MB = round(file.size(written) / 1e6, 1),
           md5 = unname(tools::md5sum(written)))
