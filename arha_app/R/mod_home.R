library(shiny)
library(bslib)

# --- UI ---
mod_home_ui <- function(id) {
  ns <- NS(id)
  tagList(card(
    class = "bg-primary text-white text-center mb-4",
    card_body(
      tags$h1("Welcome to the ArHa Database Explorer", class = "display-5 fw-bold"),
      tags$p(class = "lead", "An interactive, open-source repository for spatial epidemiology, host ecology, and zoonotic pathogen surveillance."))
    ),
    
    layout_columns(
      col_widths = c(4, 8),
      card(
        card_header("Filtering Methodology", class = "bg-dark text-white"),
        card_body(
          tags$p("The database utilises a non-destructive filtering architecture designed to maintain context during data exploration. The global filters will become available on the left-hand sidebar as soon as you navigate to an analytical module."),
          tags$h5("How it works:"),
          tags$ul(
            tags$li(tags$b("Global State:"), " Filters applied in one tab persist across the entire application."),
            tags$li(tags$b("Non-Destructive Dropdowns:"), " Selecting an option (e.g., a specific country) filters the underlying dataset, but it ", tags$em("does not"), " restrict the available choices in other dropdown menus. You can always see the full taxonomic and geographic dictionary."),
            tags$li(tags$b("Dynamic Facets:"), " To guide your exploration, the 'Top Items' lists located above each dropdown dynamically update. These facets display the actual number of positive records remaining in your filtered subset, allowing you to quickly identify data-rich entries."))
          )
        ),
      card(
        card_header(icon("github"), " Open Source & Local Deployment", class = "bg-dark text-white"),
        card_body(
          tags$p("This application is open-source. If you experience poor performance on the public server, the application can be hosted locally on your own machine."),
          tags$ul(tags$li(tags$b("Repository: "), tags$a(href="https://github.com/DidDrog11/arenavirus_hantavirus", target="_blank", "DidDrog11/arenavirus_hantavirus")),
                  tags$li(tags$b("Instructions: "), "Clone the repository, navigate to the ", tags$code("arha_app"), " directory, and launch the tool via ", tags$code("app.R"), "."),
                  tags$li(tags$b("Data Payload: "), "All necessary relational data is bundled within the application folder."))
          )
        ),
      card(
        card_header("Application Modules", class = "bg-dark text-white"),
        card_body(
          layout_columns(
            col_widths = c(6, 6),
            card(
              class = "border-info",
              card_header(icon("chart-pie"), " Overview Dashboard"),
              card_body("A high-level ecological summary of the filtered dataset. Features include temporal sampling effort, host taxonomic composition, a dynamic host-pathogen interaction matrix, and record-level prevalence analysis."),
            card(
              class = "border-success",
              card_header(icon("earth-europe"), " Spatial Explorer"),
              card_body("An interactive mapping environment for visualising trapping locations, occurrence data, and spatial biases.")
              ),
            card(
              class = "border-danger",
              card_header(icon("dna"), " Genomic Explorer"),
              card_body("A geospatial interface tracking the distribution of sequenced genetic material. Features polymorphic host/pathogen sequence mapping and integration with NCBI GenBank to bridge the gap between field ecology and molecular epidemiology.")
              ),
            card(
              class = "border-warning",
              card_header(icon("table"), " Data & Download"),
              card_body("A granular data inspection table built around FAIR principles. Review data quality metrics and export the currently filtered dataset as a Darwin Core-compliant relational archive (ZIP) in CSV or Parquet formats."))
            )
            )
          )
        )
      )
  )
}

# --- Server ---
mod_home_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Server logic empty for landing page
  })
}