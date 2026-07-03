# ArHa_App/R/mod_genetic_map.R

library(shiny)
library(bslib)
library(leaflet)
library(dplyr)
library(stringr)

# --- UI ---
mod_genetic_map_ui <- function(id) {
  ns <- NS(id)
  tagList(
    card(
      full_screen = TRUE, class = "mb-3",
      card_header(
        class = "d-flex justify-content-between align-items-center bg-dark text-white p-2",
        tags$span(icon("dna"), " Genomic Spatial Distribution", class = "fs-5"),
        tags$div(
          class = "d-flex align-items-center gap-2",
          tags$span("Filter Type:", class = "text-white-50 small"),
          selectInput(ns("seq_type"), NULL, choices = c("All Sequences", "Host Only", "Pathogen Only"), width = "150px"),
          actionButton(ns("generate_map"), "Render Map", class = "btn-primary", icon = icon("rotate"))
        )
      ),
      card_body(class = "p-0", leafletOutput(ns("seq_map"), height = "650px"))
    )
  )
}

# --- Server ---
mod_genetic_map_server <- function(id, filtered_data, tbl_sequences) {
  moduleServer(id, function(input, output, session) {
    
    map_data <- eventReactive(input$generate_map, {
      id_notify <- showNotification("Isolating genomic records...", duration = NULL, type = "message")
      on.exit(removeNotification(id_notify), add = TRUE)
      
      base_spatial <- filtered_data() |> 
        filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |> 
        select(occurrenceID, measurementID, eventID, scientificName, scientificName_pathogen, year, country, decimalLatitude, decimalLongitude, associatedReferences) |> 
        collect()
      
      seqs <- tbl_sequences |> collect()
      
      df <- bind_rows(inner_join(base_spatial, seqs, by = c("occurrenceID" = "resourceID"), relationship = "many-to-many"),
                      inner_join(base_spatial, seqs, by = c("measurementID" = "resourceID"), relationship = "many-to-many")) |> 
        distinct(accession_primary, occurrenceID, .keep_all = TRUE) |> 
        mutate(display_pathogen = case_when(sequence_type == "Host" ~ "N/A",
                                            sequence_type == "Pathogen" & !is.na(scientificName_pathogen) & scientificName_pathogen != "NA" ~ scientificName_pathogen,
                                            sequence_type == "Pathogen" ~ "Unclassified Pathogen",
                                            TRUE ~ "N/A"))
      
      if (input$seq_type == "Host Only") df <- df |> filter(sequence_type == "Host")
      if (input$seq_type == "Pathogen Only") df <- df |> filter(sequence_type == "Pathogen")
      
      return(df)
    })
    
    output$seq_map <- renderLeaflet({
      df <- map_data()
      m <- leaflet() |> addProviderTiles(providers$CartoDB.Positron) |> setView(lng = 15, lat = 5, zoom = 3)
      if (nrow(df) == 0) return(m)
      
      pal <- colorFactor(palette = c("#2980b9", "#c0392b", "#808080"), domain = c("Host", "Pathogen", "Unknown"), na.color = "#808080")
      
      m |> addCircleMarkers(
        data = df, lng = ~decimalLongitude, lat = ~decimalLatitude, radius = 6,
        fillColor = ~pal(coalesce(sequence_type, "Unknown")), fillOpacity = 0.8, color = "white", weight = 1,
        clusterOptions = markerClusterOptions(maxClusterRadius = 40, spiderfyOnMaxZoom = TRUE),
        popup = ~paste0("<div style='font-family: sans-serif; font-size: 13px; max-width: 300px;'>",
                        "<strong style='color: #2c3e50; font-size: 15px;'>", coalesce(accession_primary, "Unknown Accession"), "</strong><br>",
                        "<span style='color: white; background-color: ", ifelse(sequence_type == "Pathogen", "#c0392b", "#2980b9"), "; padding: 2px 6px; border-radius: 4px; font-size: 11px;'>", sequence_type, " Sequence</span><hr style='margin: 5px 0;'>",
                        "<b>Host:</b> ", coalesce(scientificName, "Unknown Host"), "<br>",
                        "<b>Pathogen Target:</b> ", display_pathogen, "<br>",
                        "<b>Gene Target:</b> ", coalesce(gene_name, "Unspecified"), "<br>",
                        "<b>Location:</b> ", country, "<hr style='margin: 5px 0;'>",
                        "<div style='font-size: 11px; color: #555; line-height: 1.2; margin-bottom: 8px;'><b>Citation:</b> ", coalesce(associatedReferences, "Unpublished / Direct Submission"), "</div>",
                        "<a href='", relatedResourceID, "' target='_blank' style='text-decoration: none; color: #16a085; display: block; text-align: center; background: #ecf0f1; padding: 4px; border-radius: 4px;'><b><i class='fa fa-external-link'></i> View on NCBI GenBank</b></a>",
                        "</div>")) |> 
        addLegend(position = "bottomright", pal = pal, values = coalesce(df$sequence_type, "Unknown"), title = "Sequence Type", opacity = 1)
    })
    
  })
}