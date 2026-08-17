# ArHa_App/R/utils.R
library(dplyr)
library(tibble)

# Single source of truth mapping each filter input (as named in mod_filters.R)
# to the database column it filters on. Used by BOTH get_faceted_counts() and
# mod_filters.R's filtered_data(), so the two can't drift out of sync (they
# previously duplicated this mapping separately, which is how the facet counts
# ended up silently ignoring the year filter).
filter_spec <- list(
  country  = "country",
  genus    = "genus",
  species  = "scientificName",
  fam      = "family_pathogen",
  path     = "scientificName_pathogen",
  modality = "measurementMethod"
)

# Apply the full set of active filters to a lazy dbplyr query.
# - exclude_col: DB column name to skip (used when computing facet counts for
#   that same column, so a filter doesn't restrict its own dropdown's counts)
# - apply_year: whether to also apply the year range filter (inputs$year)
apply_filters <- function(data, inputs, spec = filter_spec, exclude_col = NULL, apply_year = TRUE) {
  q <- data
  for (input_name in names(spec)) {
    col <- spec[[input_name]]
    if (!identical(col, exclude_col) && length(inputs[[input_name]]) > 0) {
      q <- q |> filter(.data[[col]] %in% !!inputs[[input_name]])
    }
  }
  if (apply_year && !is.null(inputs$year)) {
    q <- q |> filter(year >= !!inputs$year[1] & year <= !!inputs$year[2])
  }
  q
}

# Helper to calculate available options for a column, applying every OTHER
# active filter (including year) so the displayed counts always match what
# filtered_data() would actually return.
get_faceted_counts <- function(data, column_name, current_inputs, all_possible_values) {
  
  q <- apply_filters(data, current_inputs, exclude_col = column_name)
  
  # Calculate counts in SQL
  active_counts <- q |>
    filter(!is.na(!!sym(column_name))) |>
    group_by(!!sym(column_name)) |>
    summarise(n = n(), .groups = 'drop') |>
    collect()
  
  # Join and order by n (descending)
  result <- tibble(!!column_name := all_possible_values) |>
    left_join(active_counts, by = column_name) |>
    mutate(n = coalesce(n, 0L)) |>
    arrange(desc(n), !!sym(column_name))
  
  return(result)
}

build_publications_summary <- function(filtered_data) {
  q <- filtered_data()
  
  host_level <- q |>
    filter(!is.na(occurrenceID)) |>
    distinct(occurrenceID, eventID, country, decimalLatitude, decimalLongitude,
             eventDate, scientificName, individualCount) |>
    collect()
  
  path_level  <- q |> filter(!is.na(measurementID)) |> collect()
  
  n_studies     <- n_distinct(c(host_level$eventID, path_level$eventID))
  n_countries   <- n_distinct(host_level$country, na.rm = TRUE)
  n_sites       <- host_level |> distinct(decimalLatitude, decimalLongitude, eventDate) |> nrow()
  n_individuals <- sum(host_level$individualCount, na.rm = TRUE)
  n_species     <- n_distinct(host_level$scientificName[host_level$individualCount > 0], na.rm = TRUE)
  
  sampling_years <- as.integer(unlist(stringr::str_extract_all(host_level$eventDate, "\\d{4}")))
  yr_range <- if (length(sampling_years) > 0) range(sampling_years, na.rm = TRUE) else c(NA, NA)
  
  n_pathogens <- n_distinct(path_level$scientificName_pathogen, na.rm = TRUE)
  n_positive  <- sum(path_level$measurementValue, na.rm = TRUE)
  
  pl <- function(n, singular, plural = paste0(singular, "s")) if (n == 1) singular else plural
  year_txt <- if (is.na(yr_range[1])) "an unspecified period"
  else if (yr_range[1] == yr_range[2]) as.character(yr_range[1])
  else paste0(yr_range[1], "\u2013", yr_range[2])
  
  glue::glue(
    "{format(n_studies, big.mark=',')} {pl(n_studies,'study')}, in {n_countries} {pl(n_countries,'country')} ",
    "conducted between {year_txt}. Trapping occurred at {format(n_sites, big.mark=',')} site-session ",
    "{pl(n_sites,'location')}, detecting {format(n_individuals, big.mark=',')} individuals across {n_species} species. ",
    "Sampling of these individuals for {n_pathogens} pathogens returned {format(n_positive, big.mark=',')} positive results."
  )
}