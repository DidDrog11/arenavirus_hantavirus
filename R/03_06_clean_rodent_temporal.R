# Project ArHa: Clean Rodent Temporal Data
# 03_05_clean_rodent_temporal.R
# Purpose: This script parses and standardizes the `eventDate` column into
# a clean, structured format for temporal analysis.

combined_data <- read_rds(here("data", "data_cleaning", "03_05_output.rds"))

event_dates <- combined_data$host %>%
  select(rodent_record_id, study_id, eventDate) %>%
  group_by(eventDate) %>%
  mutate(date_id = cur_group_id()) %>%
  ungroup()

clean_event_dates <- event_dates %>%
  distinct(date_id, .keep_all = TRUE) %>%
    mutate(
      temporal_resolution = case_when(
        # High-priority: Day ranges and full dates
        str_detect(eventDate, "\\d{4}-\\d{2}-\\d{2}/\\d{4}-\\d{2}-\\d{2}") ~ "day_range_resolution",
        str_detect(eventDate, "\\d{4}-\\d{2}-\\d{2}/\\d{1,2}") ~ "day_range_resolution",
        str_detect(eventDate, "\\d{1,2}/\\d{1,2}/\\d{4}") ~ "full_date",
        str_detect(eventDate, "\\d{4}-\\d{2}-\\d{2}") ~ "full_date",
        
        # High-priority: Month and year ranges
        str_detect(eventDate, "\\d{4}-\\d{2}/\\d{4}-\\d{2}") ~ "month_range_resolution",
        str_detect(eventDate, "\\d{4}/\\d{2}-\\d{4}/\\d{2}") ~ "month_range_resolution",
        str_detect(eventDate, "\\d{1,2}/\\d{4}-\\d{1,2}/\\d{4}") ~ "month_range_resolution",
        str_detect(eventDate, "\\d{4}-\\d{4}") ~ "year_range_resolution",
        str_detect(eventDate, "\\d{4}/\\d{4}") ~ "year_range_resolution",
        
        # Mid-priority: Month and single years
        str_detect(eventDate, "\\d{4}-\\d{2}") ~ "month_year",
        str_detect(eventDate, "\\d{4}/\\d{2}") ~ "month_year",
        str_detect(eventDate, "\\d{1,2}/\\d{4}") ~ "month_year", # Corrected rule for MM/YYYY
        str_detect(eventDate, "\\d{4}") ~ "year_only",
        str_detect(eventDate, "\\d{2}-\\d{2}") ~ "month_day",
        
        # Low-priority
        TRUE ~ "missing"
      ),
      temporal_resolution = fct(temporal_resolution, levels = c("day_range_resolution", "full_date", "month_range_resolution", "month_year", "year_range_resolution", "year_only", "month_day", "missing")),
    simple_date = parse_date_time(eventDate, orders = c("%y-%m-%d", "%y/%m/%d", "%m/%d/%y", "%y-%m", "%y-%b", "%y/%m", "%y/%b", "%m-%y", "%b-%y", "%m/%y", "%b/%y", "%y", "%d-%m-%y")))

# day_range_resolution
day_range <- clean_event_dates %>%
  filter(temporal_resolution == "day_range_resolution") %>%
  mutate(
    # Extract the two parts of the date range
    date_parts = str_split(eventDate, "/"),
    # Check for the YYYY-MM-DD/YYYY-MM-DD format
    is_full_range = map_lgl(date_parts, ~ length(.x) == 2 && str_detect(.x[2], "\\d{4}")),
    # Check for the YYYY-MM-DD/DD format
    is_day_range = map_lgl(date_parts, ~ length(.x) == 2 && str_detect(.x[2], "\\d{1,2}") && !str_detect(.x[2], "\\d{4}"))
  ) %>%
  mutate(
    start_date = case_when(
      is_full_range ~ ymd(map_chr(date_parts, ~ .x[1])),
      is_day_range ~ ymd(map_chr(date_parts, ~ .x[1])),
      TRUE ~ NA_Date_
    ),
    end_date = case_when(
      is_full_range ~ ymd(map_chr(date_parts, ~ .x[2])),
      is_day_range ~ ymd(paste0(map_chr(date_parts, ~ str_extract(.x[1], "\\d{4}-\\d{2}-")), map_chr(date_parts, ~ .x[2]))),
      TRUE ~ NA_Date_
    )
  ) %>%
  # Clean up temporary columns
  select(-date_parts, -is_full_range, -is_day_range)

# full_date
full_date <- clean_event_dates %>%
  filter(temporal_resolution == "full_date") %>%
  mutate(temporal_resolution = case_when(str_length(eventDate) >= 20 ~ "day_range_resolution",
                                         TRUE ~ temporal_resolution),
         start_date = case_when(str_detect(temporal_resolution, "day_range_resolution") & str_detect(eventDate, "^\\d{2}/\\d{2}") ~ parse_date_time(str_split(eventDate, "-", simplify = TRUE)[, 1], orders = c("%d/%m/%y", "%m/%d/%y")),
                                str_detect(temporal_resolution, "day_range_resolution") & str_detect(eventDate, "^\\d{4}-\\d{1,2}") ~ parse_date_time(str_split(eventDate, "/", simplify = TRUE)[, 1], orders = c("%y-%m-%d")),
                                TRUE ~ simple_date),
         end_date = case_when(str_detect(temporal_resolution, "day_range_resolution") & str_detect(eventDate, "^\\d{2}/\\d{2}") ~ parse_date_time(str_extract(eventDate, "\\d{2}/\\d{2}/\\d{4}$"), orders = c("%d/%m/%y", "%m/%d/%y")),
                              str_detect(temporal_resolution, "day_range_resolution") & str_detect(eventDate, "^\\d{4}-\\d{1,2}") ~ parse_date_time(str_split(eventDate, "/", simplify = TRUE)[, 2], orders = c("%y-%m-%d")),
                              TRUE ~ simple_date))

# month_range_resolution
month_range <- clean_event_dates %>%
  filter(temporal_resolution == "month_range_resolution") %>%
  mutate(
    # Construct the start date
    start_date = case_when(
      str_detect(eventDate, "^\\d{4}-\\d{2}/\\d{4}-\\d{2}$") ~ ym(str_extract(eventDate, "^\\d{4}-\\d{2}")),
      str_detect(eventDate, "^\\d{1,2}/\\d{4}-\\d{1,2}/\\d{4}$") ~ my(str_extract(eventDate, "^\\d{1,2}/\\d{4}")),
      str_detect(eventDate, "^\\d{4}/\\d{2}-\\d{4}/\\d{2}$") ~ ym(str_extract(eventDate, "^\\d{4}/\\d{2}")),
      # Handle the semicolon case by extracting the first date
      str_detect(eventDate, ";") ~ my(str_extract(eventDate, "^\\d{1,2}/\\d{4}")),
      TRUE ~ NA_Date_
    ),
    # Construct the end date
    end_date = case_when(
      str_detect(eventDate, "^\\d{4}-\\d{2}/\\d{4}-\\d{2}$") ~ ym(str_extract(eventDate, "(?<=[/-])\\d{4}-\\d{2}")) + months(1) - days(1),
      str_detect(eventDate, "^\\d{1,2}/\\d{4}-\\d{1,2}/\\d{4}$") ~ my(str_extract(eventDate, "(?<=[/-])\\d{1,2}/\\d{4}")) + months(1) - days(1),
      str_detect(eventDate, "^\\d{4}/\\d{2}-\\d{4}/\\d{2}$") ~ ym(str_extract(eventDate, "(?<=[/-])\\d{4}/\\d{2}")) + months(1) - days(1),
      # Handle the semicolon case by extracting the last date
      str_detect(eventDate, ";") ~ my(str_extract(eventDate, "\\d{1,2}/\\d{4}$")) + months(1) - days(1),
      TRUE ~ NA_Date_
    )
  )

# month_year
month_date <- clean_event_dates %>%
  filter(temporal_resolution == "month_year") %>%
  mutate(
    temporal_resolution = case_when(str_detect(eventDate, ":") ~ "month_range_resolution",
                                    str_length(eventDate) >= 18 ~ "day_range_resolution",
                                    str_length(eventDate) == 9 ~ "full_date",
                                    str_length(eventDate) >= 8 ~ "month_range_resolution",
                                    TRUE ~ temporal_resolution),
    # Construct start_date and end_date based on the eventDate string
    start_date = case_when(temporal_resolution == "month_year" ~ parse_date_time(eventDate, orders = c("%y-%m", "%m/%y", "%y/%m", "%m/%y")),
                           temporal_resolution == "full_date" ~ parse_date_time(eventDate, orders = c("%y-%m-%d")),
                           temporal_resolution == "day_range_resolution" & str_detect(eventDate, "\\d{1,2}-\\d{1,2}/\\d{4}") ~ parse_date_time(str_split(eventDate, "/", simplify = TRUE)[ ,1], orders = c("%y-%m-%d")),
                           temporal_resolution == "day_range_resolution" & str_detect(eventDate, "\\d{1,2}/\\d{1,2}-\\d{4}") ~ parse_date_time(str_split(eventDate, "-", simplify = TRUE)[ ,1], orders = c("%y-%m-%d")),
                           temporal_resolution == "month_range_resolution" & str_detect(eventDate, "\\d{4}-\\d{1,2}") ~ parse_date_time(str_split(eventDate, "/", simplify = TRUE)[, 1], orders = c("%y-%m")),
                           temporal_resolution == "month_range_resolution" & str_detect(eventDate, ":") ~ parse_date_time(paste0(str_extract_all(eventDate, "\\d{4}", simplify = TRUE)[, 1], "-", str_extract_all(eventDate, "\\d{2}", simplify = TRUE)[, 1]), orders = c("%y-%m")),
                           temporal_resolution == "month_range_resolution" & str_detect(eventDate, "\\d{2}-\\d{2}") ~ parse_date_time(paste0(str_extract_all(eventDate, "\\d{4}", simplify = TRUE)[, 1], "-", str_extract_all(eventDate, "\\d{2}", simplify = TRUE)[, 1]), orders = c("%y-%m")),
                           TRUE ~ NA
                           ),
    end_date = case_when(temporal_resolution == "month_year" ~ ymd(start_date) + months(1) - days(1),
                         temporal_resolution == "full_date" ~ ymd(start_date),
                         temporal_resolution == "day_range_resolution" & str_detect(eventDate, "\\d{1,2}-\\d{1,2}/\\d{4}") ~ parse_date_time(str_split(eventDate, "/", simplify = TRUE)[ ,2], orders = c("%y-%m-%d")),
                         temporal_resolution == "day_range_resolution" & str_detect(eventDate, "\\d{1,2}/\\d{1,2}-\\d{4}") ~ parse_date_time(str_split(eventDate, "-", simplify = TRUE)[ ,2], orders = c("%y-%m-%d")),
                         eventDate == "04-12/2009" ~ parse_date_time("2009-12", orders = c("%y-%m")) + months(1) - days(1),
                         eventDate == "2016-12/2017-6" ~ parse_date_time("2017-06", orders = c("%y-%m")) + months(1) - days(1),
                         eventDate == "07/1995:10/1996:05/1998" ~ parse_date_time("1998-05", orders = c("%y-%m")) + months(1) - days(1),
                         temporal_resolution == "month_range_resolution" & str_length(eventDate) == 10 ~ parse_date_time(paste0(str_extract_all(eventDate, "\\d{4}", simplify = TRUE)[, 1], "-", str_extract_all(eventDate, "\\d{2}$", simplify = TRUE)), orders = c("%y-%m")) + months(1) - days(1),
                         TRUE ~ NA))

# year_range_resolution
year_range <- clean_event_dates %>%
  filter(temporal_resolution == "year_range_resolution") %>%
  mutate(start_date = parse_date_time(str_extract_all(eventDate, "^\\d{4}", simplify = TRUE)[, 1], orders = c("%y")),
         end_date = parse_date_time(str_extract_all(eventDate, "\\d{4}$", simplify = TRUE)[, 1], orders = c("%y")) + years(1) - days(1))

# year_only
year_date <- clean_event_dates %>%
  filter(temporal_resolution == "year_only") %>%
  mutate(temporal_resolution = case_when(str_length(eventDate) == 4 ~ "year_only",
                                         str_detect(eventDate, ";") ~ "year_range_resolution",
                                         str_detect(eventDate, "Spr|Sum|Aut|Win") ~ "season_range_resolution",
                                         str_detect(eventDate, "Aug|May|Nov|Feb") ~ "named_month_range_resolution",
                                         str_detect(eventDate, "\\d{4}-\\d{1,2}$|\\d{4}/\\d{1,2}$") ~ "month_year",
                                         str_detect(eventDate, "2010-1-5–6|2010-1-8–9") ~ "day_range_resolution",
                                         str_detect(eventDate, "\\d{4}-\\d{1,2}-\\d{1,2}") ~ "full_date"),
         start_date = case_when(str_detect(temporal_resolution, "year_only") ~ parse_date_time(eventDate, orders = c("%y")),
                                str_detect(temporal_resolution, "year_range_resolution") ~ parse_date_time(str_extract(eventDate, "^\\d{4}"), orders = c("%y")),
                                str_detect(temporal_resolution, "month_year") ~ parse_date_time(eventDate, orders = c("%y-%m", "%y/%m")),
                                str_detect(temporal_resolution, "full_date") ~ parse_date_time(eventDate, orders = c("%y-%m-%d", "%y/%m/%d")),
                                str_detect(eventDate, "2010-1-5–6") ~ parse_date_time("2010-1-5", orders = c("%y-%m-%d")),
                                str_detect(eventDate, "2010-1-8–9") ~ parse_date_time("2010-1-8", orders = c("%y-%m-%d")),
                                TRUE ~ NA),
         end_date = case_when(str_detect(temporal_resolution, "year_only") ~ start_date + years(1) - days(1),
                              str_detect(temporal_resolution, "year_range_resolution") ~ parse_date_time(str_extract(eventDate, "\\d{4}$"), orders = c("%y")) + years(1) - days(1),
                              str_detect(temporal_resolution, "month_year") ~ start_date + months(1) - days(1),
                              str_detect(temporal_resolution, "full_date") ~ start_date,
                              str_detect(eventDate, "2010-1-5–6") ~ parse_date_time("2010-1-6", orders = c("%y-%m-%d")),
                              str_detect(eventDate, "2010-1-8–9") ~ parse_date_time("2010-1-9", orders = c("%y-%m-%d")),
                              TRUE ~ NA)
         )

# month_day
month_day <- clean_event_dates %>%
  filter(temporal_resolution == "month_day") %>%
  mutate(temporal_resolution = "month_range_resolution",
         start_date = parse_date_time("1996-01-01", orders = c("%y-%m-%d")),
         end_date = parse_date_time("1998-12-31", orders = c("%y-%m-%d")))



# Cleaned Dates -----------------------------------------------------------

combined_date_processing <- bind_rows(day_range,
                                      full_date,
                                      month_range,
                                      month_date,
                                      year_range,
                                      year_date,
                                      month_day) %>%
  mutate(temporal_resolution = fct(temporal_resolution, levels = c("day_range_resolution", "full_date", "month_range_resolution", "named_month_range_resolution", "season_range_resolution", "month_year", "year_only", "year_range_resolution", "undated")))

processed_event_dates <- event_dates %>%
  left_join(combined_date_processing %>%
              bind_rows(event_dates %>%
                          filter(!date_id %in% combined_date_processing$date_id) %>%
                          mutate(temporal_resolution = "undated") %>%
                          distinct(date_id, , .keep_all = TRUE)) %>%
              select(date_id, temporal_resolution, start_date, end_date),
            by = c("date_id")) %>%
  select(-date_id, -eventDate) %>%
  left_join(combined_data$host, by = c("rodent_record_id", "study_id"))

remaining_undated <- processed_event_dates %>%
  filter(is.na(start_date)|is.na(end_date)) %>%
  left_join(combined_data$descriptives_cleaned %>% 
              select(study_id, publication_year), by = "study_id") %>%
  mutate(temporal_resolution = case_when(temporal_resolution == "undated" ~ "publication_derived"),
         end_date = case_when(str_detect(temporal_resolution, "publication_derived") ~ parse_date_time(publication_year, orders = c("%y"))))

# Helper function to parse seasonal dates based on hemisphere
parse_seasonal_date <- function(date_string, latitude) {
  
  # Determine hemisphere
  hemisphere <- if_else(latitude >= 0, "Northern", "Southern")
  
  # Detect and extract year ranges
  if (str_detect(date_string, "\\d{4}-\\d{4}")) {
    years <- as.numeric(str_extract_all(date_string, "\\d{4}")[[1]])
    start_year <- years[1]
    end_year <- years[2]
  } else {
    # If no range, treat it as a single year
    year <- as.numeric(str_extract(date_string, "\\d{4}"))
    start_year <- year
    end_year <- year
  }
  
  # Extract the season
  season_start <- str_extract(date_string, "Spring|Summer|Autumn|Winter")
  season_end <- str_extract(date_string, "(?=-)[A-Za-z]+") %>% str_remove("-")
  
  # Correct logic for a single-season date
  if (is.na(season_end)) {
    season_end <- season_start
  }
  
  if (hemisphere == "Northern") {
    start_date <- case_when(
      season_start == "Spring" ~ ymd(paste0(start_year, "-03-21")),
      season_start == "Summer" ~ ymd(paste0(start_year, "-06-21")),
      season_start == "Autumn" ~ ymd(paste0(start_year, "-09-23")),
      season_start == "Winter" ~ ymd(paste0(start_year, "-12-21")),
      TRUE ~ NA_Date_
    )
    
    end_date <- case_when(
      season_end == "Spring" ~ ymd(paste0(end_year, "-06-20")),
      season_end == "Summer" ~ ymd(paste0(end_year, "-09-22")),
      season_end == "Autumn" ~ ymd(paste0(end_year, "-12-20")),
      season_end == "Winter" ~ ymd(paste0(end_year + 1, "-03-20")),
      TRUE ~ NA_Date_
    )
    
  } else { # Southern Hemisphere
    start_date <- case_when(
      season_start == "Spring" ~ ymd(paste0(start_year, "-09-23")),
      season_start == "Summer" ~ ymd(paste0(start_year, "-12-21")),
      season_start == "Autumn" ~ ymd(paste0(start_year, "-03-21")),
      season_start == "Winter" ~ ymd(paste0(start_year, "-06-21")),
      TRUE ~ NA_Date_
    )
    
    end_date <- case_when(
      season_end == "Spring" ~ ymd(paste0(end_year, "-12-20")),
      season_end == "Summer" ~ ymd(paste0(end_year + 1, "-03-20")),
      season_end == "Autumn" ~ ymd(paste0(end_year, "-06-20")),
      season_end == "Winter" ~ ymd(paste0(end_year, "-09-22")),
      TRUE ~ NA_Date_
    )
  }
  
  return(tibble(start_date, end_date))
}

seasonal_dates_cleaned <- remaining_undated %>%
  filter(str_detect(eventDate, "Spring|Summer|Autumn|Winter")) %>%
  rowwise() %>%
  mutate(
    parsed_seasonal = list(parse_seasonal_date(eventDate, decimalLatitude))
  ) %>%
  ungroup() %>%
  unnest_wider(col = parsed_seasonal, names_sep = "_") %>%
  mutate(
    start_date = coalesce(start_date, parsed_seasonal_start_date),
    end_date = coalesce(end_date, parsed_seasonal_end_date),
    temporal_resolution = "month_range_resolution"
  ) %>%
  select(-parsed_seasonal_start_date, -parsed_seasonal_end_date) %>%
  select(rodent_record_id, study_id, temporal_resolution, start_date, end_date)

publication_derived_dates <- remaining_undated %>%
  filter(temporal_resolution == "publication_derived") %>%
  select(rodent_record_id, study_id, temporal_resolution, start_date, end_date)

final_event_dates <- processed_event_dates %>%
  select(rodent_record_id, study_id, temporal_resolution, start_date, end_date) %>%
  filter(!rodent_record_id %in% c(seasonal_dates_cleaned$rodent_record_id, publication_derived_dates$rodent_record_id)) %>%
  bind_rows(seasonal_dates_cleaned, publication_derived_dates) %>%
  mutate(temporal_resolution = fct(temporal_resolution, levels = c("day_range_resolution", "full_date", "month_range_resolution", "named_month_range_resolution", "season_range_resolution", "month_year", "year_only", "year_range_resolution", "publication_derived", "undated")))

# Manual fixes, may not be needed

combined_data$host <- combined_data$host %>%
  left_join(final_event_dates, by = c("rodent_record_id", "study_id")) %>%
  rename("event_date_raw" = eventDate)

write_rds(combined_data, here("data", "data_cleaning", "03_06_output.rds"))
