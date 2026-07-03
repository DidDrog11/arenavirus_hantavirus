# ArHa_App/R/utils.R
library(dplyr)
library(tibble)

# Helper to calculate available options for a column, 
get_faceted_counts <- function(data, column_name, current_inputs, all_possible_values) {
  
  inputs_to_apply <- current_inputs
  inputs_to_apply[[column_name]] <- NULL 
  
  q <- data
  
  # Apply filters that aren't the current column
  if (length(inputs_to_apply$country) > 0)   q <- q |> filter(country %in% !!inputs_to_apply$country)
  if (length(inputs_to_apply$genus) > 0)     q <- q |> filter(genus %in% !!inputs_to_apply$genus)
  if (length(inputs_to_apply$species) > 0 && column_name != "scientificName") q <- q |> filter(scientificName %in% !!inputs_to_apply$species)
  if (length(inputs_to_apply$fam) > 0 && column_name != "family_pathogen")     q <- q |> filter(family_pathogen %in% !!inputs_to_apply$fam)
  if (length(inputs_to_apply$path) > 0 && column_name != "scientificName_pathogen") q <- q |> filter(scientificName_pathogen %in% !!inputs_to_apply$path)
  
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