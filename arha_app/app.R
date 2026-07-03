# ArHa_App/app.R
library(shiny)
library(bslib)
library(dplyr)

# Source the Database Engine (DuckDB)
source("R/db_prep.R")

# Source Modules
source("R/mod_home.R")
source("R/utils.R")
source("R/mod_filters.R")
source("R/mod_facets.R")
source("R/mod_dashboard.R")
source("R/mod_map.R")
source("R/mod_table.R")
source("R/mod_genetic_map.R")


# UI ----------------------------------------------------------------------
ui <- page_navbar(
  title = "ArHa Database Explorer",
  id = "main_nav",
  theme = bs_theme(version = 5, preset = "flatly"),
  fillable = TRUE,
  
  # Global Sidebar - Initialised as closed
  sidebar = sidebar(
    id = "global_sidebar", 
    width = 300, 
    open = "closed",
    conditionalPanel(
      condition = "input.main_nav !== 'home'",
      tags$h5("Data Filters", class = "mt-2 mb-3"),
      mod_filters_ui("global_filters"))),
  
  # Home / Landing Page
  nav_panel(
    title = "Home", 
    value = "home",
    icon = icon("house"),
    mod_home_ui("home")
  ),
  
  # Dashboard
  nav_panel(
    title = "Overview Dashboard", 
    value = "dashboard",
    icon = icon("chart-pie"),
    mod_dashboard_ui("dashboard")
  ),
  
  # Spatial Explorer
  nav_panel(
    title = "Spatial Explorer", 
    value = "spatial",
    icon = icon("earth-europe"),
    mod_map_ui("spatial-explorer")
  ),
  
  # Genomic Explorer
  nav_panel(
    title = "Genomic Explorer", 
    value = "genomic", 
    icon = icon("dna"), 
    mod_genetic_map_ui("genetic_map")
  ),
  
  # Data Table & Download
  nav_panel(
    title = "Data & Download", 
    value = "data",
    icon = icon("table"),
    mod_table_ui("data")
  )
)


# Server ------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Static Modules
  mod_home_server("home")
  
  # Reactively open/close sidebar depending on active tab
  observeEvent(input$main_nav, {
    bslib::sidebar_toggle(
      id = "global_sidebar", 
      open = input$main_nav != "home"
    )
  })
  
  # Initialise Core Filter
  filtered_query <- mod_filters_server("global_filters", tbl_events, tbl_hosts, tbl_pathogens)
  
  # Initialise Visualisation Modules
  mod_dashboard_server("dashboard", filtered_query)
  
  mod_map_server("spatial-explorer", filtered_query)
  
  mod_genetic_map_server("genetic_map", filtered_query, tbl_sequences)
  
  # Initialise Data Export Module
  mod_table_server("data", filtered_query)
  
  session$onSessionEnded(function() {
    message("Disconnecting DuckDB...")
    DBI::dbDisconnect(con, shutdown = TRUE)})
}

shinyApp(ui, server)