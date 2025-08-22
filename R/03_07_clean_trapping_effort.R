# Project ArHa: Clean Trapping Effort Data
# 03_07_clean_trapping_effort.R
# Purpose: This script parses and standardizes the `trapEffort` and
# `trapEffortResolution` columns into a clean, structured format.

combined_data <- read_rds(here("data", "data_cleaning", "03_06_output.rds"))

numbers_regex <- "\\d+(?:,\\d+)*"

trap_effort_cleaning <- combined_data$host %>%
  select(rodent_record_id, study_id, locality_raw, verbatimLocality_raw, coordinate_resolution_raw, coordinate_resolution_processed, trapEffort, trapEffortResolution, event_date_raw, temporal_resolution, start_date, end_date) %>%
  mutate(
    # Use a single, comprehensive regex to find all numbers and get the max
    trap_nights_clean = map_dbl(str_extract_all(trapEffort, "\\d+(?:,\\d+)*"), ~ {
      # Handle cases where no number is found
      if (length(.x) == 0) {
        return(NA_real_)
      }
      # Otherwise, get the max number from the list
      suppressWarnings(max(as.numeric(str_remove_all(.x, ",")), na.rm = TRUE))
    }),
    trap_nights_clean = case_when(trap_nights_clean == -Inf ~ NA_real_,
                                  TRUE ~ trap_nights_clean),
    
    # Create the status flag based on the presence of a number and keywords
    trap_nights_status = case_when(
      is.na(trap_nights_clean) ~ "not_reported",
      trap_nights_clean < 10 ~ "nights",
      str_detect(tolower(trapEffort), "trap-?nights?|trap night") | trap_nights_clean >= 10 ~ "trapnight_reported",
      str_detect(trapEffort, "\\d+-\\d+") ~ "trapnight_range",
      TRUE ~ "unspecified"
    ),
    trap_nights_status = fct(trap_nights_status, levels = c("trapnight_reported", "trap_night_range", "nights", "unspecified", "not_reported")),
    
    # Clean trapEffortResolution
    trap_effort_resolution_clean = case_when(
      str_detect(str_to_lower(trapEffortResolution), "site|village|verbatim") & str_detect(str_to_lower(trapEffortResolution), "session|season|visit") ~ "site-session",
      str_detect(str_to_lower(trapEffortResolution), "site|verbatim|region") ~ "site",
      str_detect(str_to_lower(trapEffortResolution), "session|season|visit") ~ "session",
      str_detect(str_to_lower(trapEffortResolution), "study") ~ "study",
      TRUE ~ NA
    ),
    trap_effort_resolution_clean = fct(trap_effort_resolution_clean, levels = c("site-session", "site", "session", "study"))
  )

combined_data$host <- combined_data$host %>%
  select(-trapEffort, -trapEffortResolution) %>%
  left_join(trap_effort_cleaning %>%
              select(rodent_record_id, trap_nights_status, trap_nights_clean, trap_effort_raw = trapEffort, trap_effort_resolution_clean, trap_effort_resolution_raw = trapEffortResolution), by = "rodent_record_id") %>%
  relocate(individual_count = individualCount, .before = country_raw)

write_rds(combined_data, here("data", "data_cleaning", "03_07_output.rds"))
