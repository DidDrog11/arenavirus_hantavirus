# ArHa_App/R/mod_filters.R
library(shiny)
library(dplyr)
library(tidyr)

# --- UI ---
mod_filters_ui <- function(id) {
  ns <- NS(id)
  tagList(
    accordion(
      accordion_panel(
        "Geography & Study",
        mod_facets_ui(ns("country_counts")),
        selectizeInput(ns("country"), "Country", choices = NULL, multiple = TRUE),
        sliderInput(ns("year"), "Publication Year", min = 1950, max = 2026, value = c(1950, 2026), sep = "")
      ),
      accordion_panel(
        "Host Filters",
        mod_facets_ui(ns("genus_counts")),
        selectizeInput(ns("host_genus"), "Genus", choices = NULL, multiple = TRUE),
        tags$hr(),
        mod_facets_ui(ns("host_counts")),
        selectizeInput(ns("host_species"), "Species", choices = NULL, multiple = TRUE)
      ),
      accordion_panel(
        "Pathogen Filters",
        mod_facets_ui(ns("path_counts")),
        selectizeInput(ns("pathogen_family"), "Family", choices = NULL, multiple = TRUE),
        selectizeInput(ns("pathogen"), "Species", choices = NULL, multiple = TRUE)
      )
    ),
    actionButton(ns("reset"), "Reset Filters", class = "btn-secondary mt-3 w-100")
  )
}

# --- Server ---
mod_filters_server <- function(id, tbl_events, tbl_hosts, tbl_pathogens) {
  moduleServer(id, function(input, output, session) {
    
    # Base join
    base_join <- reactive({
      tbl_events |>
        left_join(tbl_hosts, by = "eventID") |>
        left_join(tbl_pathogens, by = c("eventID", "occurrenceID"))
    })
    
    # Pre-fetch possible values
    all_vals <- list(
      country = tbl_hosts |> filter(!is.na(country)) |> distinct(country) |> pull() |> sort(),
      genus   = tbl_hosts |> filter(!is.na(genus)) |> distinct(genus) |> pull() |> sort(),
      species = tbl_hosts |> filter(!is.na(scientificName)) |> distinct(scientificName) |> pull() |> sort(),
      fam     = tbl_pathogens |> filter(!is.na(family_pathogen)) |> distinct(family_pathogen) |> pull() |> sort(),
      path    = tbl_pathogens |> filter(!is.na(scientificName_pathogen)) |> distinct(scientificName_pathogen) |> pull() |> sort()
    )
    
    # Populate UI choices
    observe({
      updateSelectizeInput(session, "country", choices = all_vals$country)
      updateSelectizeInput(session, "host_genus", choices = all_vals$genus)
      updateSelectizeInput(session, "host_species", choices = all_vals$species)
      updateSelectizeInput(session, "pathogen_family", choices = all_vals$fam)
      updateSelectizeInput(session, "pathogen", choices = all_vals$path)
    })
    
    # Reactive inputs
    curr_inputs <- reactive({
      list(country = input$country, genus = input$host_genus, species = input$host_species, 
           fam = input$pathogen_family, path = input$pathogen)
    })
    
    # Facet data pipeline
    country_facet_data <- reactive({
      get_faceted_counts(base_join(), "country", curr_inputs(), all_vals$country) |> 
        filter(n > 0) |> arrange(desc(n))
    })
    
    genus_facet_data <- reactive({
      get_faceted_counts(base_join(), "genus", curr_inputs(), all_vals$genus) |> 
        filter(n > 0) |> arrange(desc(n))
    })
    
    host_facet_data <- reactive({
      get_faceted_counts(base_join(), "scientificName", curr_inputs(), all_vals$species) |> 
        filter(n > 0) |> arrange(desc(n))
    })
    
    path_facet_data <- reactive({
      get_faceted_counts(base_join(), "scientificName_pathogen", curr_inputs(), all_vals$path) |> 
        filter(n > 0) |> arrange(desc(n))
    })
    
    # Initialise facets
    mod_facets_server("country_counts", country_facet_data, title_text = "Top Countries:")
    mod_facets_server("genus_counts", genus_facet_data, title_text = "Top Genera:")
    mod_facets_server("host_counts", host_facet_data, title_text = "Top Host Species:")
    mod_facets_server("path_counts", path_facet_data, title_text = "Top Pathogens:")
    
    # Final query
    filtered_data <- reactive({
      final_query <- base_join()
      
      if (length(input$country) > 0)        final_query <- final_query |> filter(country %in% !!input$country)
      if (length(input$host_genus) > 0)     final_query <- final_query |> filter(genus %in% !!input$host_genus)
      if (length(input$host_species) > 0)   final_query <- final_query |> filter(scientificName %in% !!input$host_species)
      if (length(input$pathogen_family) > 0) final_query <- final_query |> filter(family_pathogen %in% !!input$pathogen_family)
      if (length(input$pathogen) > 0)       final_query <- final_query |> filter(scientificName_pathogen %in% !!input$pathogen)
      if (!is.null(input$year))             final_query <- final_query |> filter(year >= !!input$year[1] & year <= !!input$year[2])
      
      return(final_query)
    })
    
    # Reset button
    observeEvent(input$reset, {
      updateSelectizeInput(session, "country", selected = character(0))
      updateSelectizeInput(session, "host_genus", selected = character(0))
      updateSelectizeInput(session, "host_species", selected = character(0))
      updateSelectizeInput(session, "pathogen_family", selected = character(0))
      updateSelectizeInput(session, "pathogen", selected = character(0))
      years <- tbl_events |> filter(!is.na(year)) |> distinct(year) |> pull()
      if (length(years) > 0) updateSliderInput(session, "year", value = c(min(years), max(years)))
    })
    
    return(filtered_data)
  })
}