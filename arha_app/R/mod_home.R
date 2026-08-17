library(shiny)
library(bslib)

# Version constants — update these on each release.
# APP_VERSION:       semantic version of the Explorer application itself.
# DB_SEARCH_DATE:    literature search cutoff date.
# DB_COMPILED_DATE:  date the current Parquet files were produced by the pipeline.
# DB_ZENODO_DOI:     permanent DOI for the archived release on Zenodo.
APP_VERSION      <- "1.0.0"
DB_SEARCH_DATE   <- "24 August 2023"
DB_COMPILED_DATE <- "14 August 2026"
DB_ZENODO_DOI    <- "[PLACEHOLDER — insert DOI]"

# --- UI ---
mod_home_ui <- function(id) {
  ns <- NS(id)
  tagList(
    card(
      class = "bg-primary text-white text-center mb-4", fill = FALSE,
      card_body(
        tags$h1("Welcome to the ArHa Database Explorer", class = "display-5 fw-bold"),
        tags$p(class = "lead mb-0", "An interactive, open-source repository for spatial epidemiology, host ecology, and zoonotic pathogen surveillance.")
      )
    ),
    
    tags$h4("Explore the Database", class = "mb-3"),
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      card(class = "border-info h-100",    fill = FALSE, card_header(icon("chart-pie"),    " Overview Dashboard"), card_body("A high-level ecological summary of the filtered dataset: temporal sampling effort, host taxonomic composition, host-pathogen interaction matrix, and record-level prevalence.")),
      card(class = "border-success h-100", fill = FALSE, card_header(icon("earth-europe"), " Spatial Explorer"),   card_body("An interactive mapping environment for visualising trapping locations, occurrence data, and spatial biases.")),
      card(class = "border-danger h-100",  fill = FALSE, card_header(icon("dna"),          " Genomic Explorer"),   card_body("Geospatial tracking of sequenced genetic material, with polymorphic host/pathogen sequence mapping and NCBI GenBank integration.")),
      card(class = "border-warning h-100", fill = FALSE, card_header(icon("table"),        " Data & Download"),    card_body("Granular data inspection built around FAIR principles. Review quality metrics and export the filtered dataset as a Darwin Core archive (CSV or Parquet)."))
    ),
    tags$br(),
    
    layout_columns(
      col_widths = c(6, 6),
      card(
        fill = FALSE,
        card_header("Filtering Methodology", class = "bg-dark text-white"),
        card_body(
          tags$p("The database uses a non-destructive filtering architecture designed to maintain context during data exploration. Global filters appear on the left-hand sidebar as soon as you navigate to an analytical module."),
          tags$ul(
            tags$li(tags$b("Global State:"), " Filters applied in one tab persist across the entire application."),
            tags$li(tags$b("Non-Destructive Dropdowns:"), " Selecting an option (e.g. a country) filters the dataset but does ", tags$em("not"), " restrict the choices available in other dropdowns — the full taxonomic and geographic dictionary always stays visible."),
            tags$li(tags$b("Dynamic Facets:"), " The 'Top Items' lists above each dropdown update live, showing the number of records remaining in your filtered subset so you can spot data-rich entries quickly.")
          )
        )
      ),
      card(
        fill = FALSE,
        card_header(icon("github"), " Open Source & Local Deployment", class = "bg-dark text-white"),
        card_body(
          tags$p("This application is open-source. If you experience poor performance on the public server, it can be hosted locally on your own machine."),
          tags$ul(
            tags$li(tags$b("Repository: "), tags$a(href = "https://github.com/DidDrog11/arenavirus_hantavirus", target = "_blank", "DidDrog11/arenavirus_hantavirus")),
            tags$li(tags$b("Instructions: "), "Clone the repository, navigate to the ", tags$code("arha_app"), " directory, and launch via ", tags$code("app.R"), "."),
            tags$li(tags$b("Data Payload: "), "All relational data is bundled within the application folder.")
          )
        )
      )
    ),
    
    tags$div(
      class = "mt-4 pt-3 border-top d-flex flex-wrap justify-content-between text-muted",
      style = "font-size: 0.85em;",
      tags$div(tags$b("Literature search current to: "), DB_SEARCH_DATE),
      tags$div(tags$b("Database compiled: "), DB_COMPILED_DATE),
      tags$div(tags$b("Archive: "), tags$a(href = paste0("https://doi.org/", DB_ZENODO_DOI), target = "_blank", DB_ZENODO_DOI)),
      tags$div(tags$b("App version: "), APP_VERSION)
    )
  )
}


# --- Server ---
mod_home_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Static landing page — no server logic required.
  })
}