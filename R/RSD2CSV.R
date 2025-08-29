# RDS2CSV.r
# Minimal utility to export a named table from an .rds bundle to CSV (UTF-8)
# Usage:
#   Rscript extract_db.R /path/input.rds table_name /path/out.csv
#   Rscript extract_db.R /path/input.rds --list      # lists available names

extract_db <- function(pathToInput, databaseToExtract, outputPath) {
  # Dependencies (readr only)
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required. Install with install.packages('readr').")
  }

  # 1) Load the .rds object
  obj <- readRDS(pathToInput)

  # 2) Allow simple bundles: either a list of tables or a single table
  if (identical(databaseToExtract, "--list")) {
    if (is.list(obj) && !is.null(names(obj))) {
      cat("Available tables:\n")
      cat(paste0(" - ", names(obj)), sep = "\n")
    } else {
      cat("The .rds does not contain a named list; it contains a single object of class: ",
          paste(class(obj), collapse = ","), "\n", sep = "")
    }
    return(invisible(NULL))
  }

  if (!is.list(obj)) {
    # Single object: allow extracting it by using any token (or "*")
    if (!(databaseToExtract %in% c("*", "1", ""))) {
      stop("Input .rds contains a single object (class: ",
           paste(class(obj), collapse=","), "). ",
           "Use databaseToExtract='*' or pass '--list' to inspect.")
    }
    df <- obj
  } else {
    # Named list: pick requested element
    if (is.null(names(obj)) || !databaseToExtract %in% names(obj)) {
      stop("Database '", databaseToExtract, "' not found. ",
           if (!is.null(names(obj))) paste0("Available: ",
           paste(names(obj), collapse = ", ")))
    }
    df <- obj[[databaseToExtract]]
  }

  # 3) Coerce to a rectanglular table if possible
  if (inherits(df, "data.table")) df <- as.data.frame(df)
  if (!inherits(df, "data.frame")) {
    stop("Selected object is not a data.frame-like table (class: ",
         paste(class(df), collapse=","), ").")
  }

  # 4) Ensure output directory exists
  outdir <- dirname(outputPath)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # 5) Write UTF-8 CSV (no BOM)
  readr::write_csv(df, file = outputPath)

  # If you need Excel-friendly CSV (UTF-8 with BOM), use:
  # readr::write_excel_csv(df, file = outputPath)

  invisible(TRUE)
}

# ---- CLI wrapper ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2 && args[2] == "--list") {
  # List available tables in the bundle
  extract_db(args[1], "--list", outputPath = tempfile())
} else if (length(args) == 3) {
  ok <- extract_db(args[1], args[2], args[3])
  if (isTRUE(ok)) {
    message("Wrote CSV to: ", args[3])
  }
} else {
  cat(
"Usage:
  Rscript extract_db.R /path/input.rds table_name /path/out.csv
  Rscript extract_db.R /path/input.rds --list

Notes:
- Output encoding is UTF-8 (no BOM). For Excel-compatible UTF-8 with BOM, use write_excel_csv().
- 'table_name' must match a named element in the .rds if it is a list; use '--list' to inspect.\n"
  )
}