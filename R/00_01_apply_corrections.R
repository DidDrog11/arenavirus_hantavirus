# Project ArHa: Apply corrections to the extraction sheets
# 00_01_apply_corrections.R
# Sourced at the end of 00_load_data.R. Corrections to extracted records are
# never made in the Google Sheets or the downloaded snapshot; they are listed
# in tracked tables in data/corrections/ and applied here, so every change is
# recorded with its source in the paper. The corrected copies are saved as
# <analysis_date>_v2_data_corrected.rds and _v3_data_corrected.rds, which the
# cleaning scripts read.
#
# data/corrections/field_corrections.csv
#   One row per value: version (v2/v3), sheet, id_column, record_id, field,
#   current_value (as in the sheet; the run stops if it no longer matches),
#   corrected_value, source (paper and table/page), note.
# data/corrections/row_deletions.csv
#   Rows removed: version, sheet, id_column, record_id, source, note. Used for
#   rows replaced by a re-extraction and for rows that are not records (e.g.
#   0/0 cells where nothing was caught).
# data/corrections/rows_added/<version>_<sheet>.csv
#   Re-extracted rows in the sheet's own columns, with IDs already prefixed
#   (e.g. "david_49001").

corrections_dir <- here("data", "corrections")

field_corrections <- read_csv(file.path(corrections_dir, "field_corrections.csv"),
                              col_types = cols(.default = col_character()))
row_deletions     <- read_csv(file.path(corrections_dir, "row_deletions.csv"),
                              col_types = cols(.default = col_character()))

# Sheet values compared as text; NA and "" are treated as the same
as_text <- function(x) if_else(is.na(x), "", as.character(x))

# --- 1. Field corrections ------------------------------------------------------
apply_field_corrections <- function(data, version) {
  fc <- field_corrections |> filter(.data$version == !!version)
  for (i in seq_len(nrow(fc))) {
    r  <- fc[i, ]
    df <- data[[r$sheet]]
    rows <- which(as.character(df[[r$id_column]]) == r$record_id)
    if (length(rows) != 1) stop("Correction ", i, ": ", r$record_id, " found ", length(rows), " times in ", version, " ", r$sheet)
    if (as_text(df[[r$field]][rows]) != as_text(r$current_value)) {
      stop("Correction ", i, ": ", r$record_id, " ", r$field, " is '", df[[r$field]][rows],
           "', expected '", r$current_value, "'. The sheet has changed; review the correction.")
    }
    # Keep the column's type: numeric fields take numeric values; ID columns are
    # factors, and a value that is not yet a level would become NA, so the
    # factor is rebuilt in order of appearance
    if (is.numeric(df[[r$field]])) {
      df[[r$field]][rows] <- as.numeric(r$corrected_value)
    } else if (is.factor(df[[r$field]])) {
      col <- as.character(df[[r$field]])
      col[rows] <- r$corrected_value
      df[[r$field]] <- fct_inorder(col)
    } else {
      df[[r$field]][rows] <- r$corrected_value
    }
    data[[r$sheet]] <- df
  }
  message(version, ": ", nrow(fc), " field corrections applied")
  data
}

# --- 2. Row deletions ----------------------------------------------------------
apply_row_deletions <- function(data, version) {
  rd <- row_deletions |> filter(.data$version == !!version)
  for (s in unique(rd$sheet)) {
    ids  <- rd$record_id[rd$sheet == s]
    idc  <- unique(rd$id_column[rd$sheet == s])
    df   <- data[[s]]
    missing <- setdiff(ids, as.character(df[[idc]]))
    if (length(missing) > 0) stop(version, " ", s, ": rows to delete not found: ", paste(missing, collapse = ", "))
    data[[s]] <- df[!as.character(df[[idc]]) %in% ids, ]
  }
  message(version, ": ", nrow(rd), " rows deleted")
  data
}

# --- 3. Added rows -------------------------------------------------------------
# IDs in the sheets are factors in order of appearance; added rows are appended
# and the factors rebuilt so that order is kept
apply_added_rows <- function(data, version) {
  files <- list.files(file.path(corrections_dir, "rows_added"), pattern = paste0("^", version, "_.*\\.csv$"), full.names = TRUE)
  for (f in files) {
    s   <- sub(paste0("^", version, "_(.*)\\.csv$"), "\\1", basename(f))
    df  <- data[[s]]
    add <- read_csv(f, col_types = cols(.default = col_character()))
    if (!all(names(add) %in% names(df))) stop(basename(f), ": columns not in the ", s, " sheet")
    # Match the sheet's column types before binding
    for (col in names(add)) if (is.numeric(df[[col]])) add[[col]] <- as.numeric(add[[col]])
    id_cols <- names(df)[vapply(df, is.factor, logical(1))]
    out <- bind_rows(df |> mutate(across(all_of(id_cols), as.character)), add) |>
      mutate(across(all_of(id_cols), fct_inorder))
    data[[s]] <- out
    message(version, " ", s, ": ", nrow(add), " rows added from ", basename(f))
  }
  data
}

combined_data_v2 <- combined_data_v2 |> apply_row_deletions("v2") |> apply_added_rows("v2") |> apply_field_corrections("v2")
combined_data_v3 <- combined_data_v3 |> apply_row_deletions("v3") |> apply_added_rows("v3") |> apply_field_corrections("v3")

write_rds(combined_data_v2, here("data", "raw_data", paste0(analysis_date, "_v2_data_corrected.rds")))
write_rds(combined_data_v3, here("data", "raw_data", paste0(analysis_date, "_v3_data_corrected.rds")))
