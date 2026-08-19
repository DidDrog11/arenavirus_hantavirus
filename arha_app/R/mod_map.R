# ArHa_App/R/mod_map.R
library(shiny)
library(bslib)
library(leaflet)
library(dplyr)
library(RColorBrewer)
library(stringr)

# --- UI ---
mod_map_ui <- function(id) {
  ns <- NS(id)
  tagList(
    card(
      full_screen = TRUE,
      class = "mb-3",
      
      card_header(
        class = "d-flex justify-content-between align-items-center bg-dark text-white p-2",
        tags$span(icon("map-location-dot"), " Spatial Distribution", class = "fs-5"),
        
        tags$div(
          class = "d-flex align-items-center gap-2",
          tags$span("Colour by:", class = "text-white-50 small"),
          selectInput(
            ns("color_var"), 
            label = NULL, 
            choices = c(
              "Positivity Status" = "positivity_status",
              "Pathogen Tested" = "map_colour_pathogen",
              "Host Genus" = "genus",
              "Host Species" = "scientificName"
            ),
            width = "200px"
          ),
          actionButton(ns("generate_map"), "Generate Map", class = "btn-primary", icon = icon("rotate"))
        )
      ),
      
      card_body(
        class = "p-0",
        leafletOutput(ns("map_plot"), height = "650px")
      )
    )
  )
}

# --- Server ---
mod_map_server <- function(id, filtered_data, tbl_sequences) {
  moduleServer(id, function(input, output, session) {
    
    map_data <- eventReactive(input$generate_map, {
      
      id_notify <- showNotification("Querying DuckDB & Rendering Map...", duration = NULL, type = "message")
      on.exit(removeNotification(id_notify), add = TRUE)
      
      base <- filtered_data() |> 
        filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |> 
        select(occurrenceID, measurementID, eventID,
         scientificName, genus, 
         scientificName_pathogen, genus_pathogen, family_pathogen, 
         measurementRemarks, pathogen_species_cleaned, 
         year, eventDate, country, decimalLatitude, decimalLongitude, 
         measurementValue, number_tested,
         associatedReferences)  |> 
        collect() |> 
        mutate(positivity_status = if_else(measurementValue > 0, "Positive Detection", "Negative/No Detection"),
               verbatim_raw = trimws(str_remove(str_extract(measurementRemarks, "Verbatim pathogen: [^|]+"), "Verbatim pathogen: ")),
               verbatim_pathogen = if_else(verbatim_raw %in% c("NA", "", "N/A"), NA_character_, verbatim_raw),
               sci_path_clean = if_else(scientificName_pathogen %in% c("NA", ""), NA_character_, scientificName_pathogen),
               gen_path_clean = if_else(genus_pathogen %in% c("NA", ""), NA_character_, genus_pathogen),
               fam_path_clean = if_else(family_pathogen %in% c("NA", ""), NA_character_, family_pathogen),
               cleaned_path_clean = if_else(pathogen_species_cleaned %in% c("NA", ""), NA_character_, pathogen_species_cleaned),
               display_pathogen = coalesce(sci_path_clean, gen_path_clean, fam_path_clean, verbatim_pathogen, cleaned_path_clean, "Not Specified"),
               map_colour_pathogen = coalesce(sci_path_clean, gen_path_clean, fam_path_clean, "Other / Verbatim"),
               extracted_doi = str_extract(associatedReferences, "(?i)doi:\\s*(10\\.\\S+)"),
               clean_doi     = str_remove(extracted_doi, "(?i)doi:\\s*"),
               event_link    = if_else(
                 !is.na(clean_doi),
                 paste0("(<a href='https://doi.org/", clean_doi, "' target='_blank'>View Source</a>)"),
                 ""
               ))
      
      # Attach linked sequence references -- resourceID is polymorphic
      # (host-keyed OR pathogen-keyed), so check both directions.
      seqs <- tbl_sequences |> collect()
      
      seq_by_host <- seqs |> 
        semi_join(base, by = c("resourceID" = "occurrenceID")) |>
        group_by(resourceID) |>
        summarise(sequence_refs = paste0(
          "<a href='", relatedResourceID, "' target='_blank'>", accession_primary, "</a>",
          collapse = "<br>"), .groups = "drop") |>
        rename(occurrenceID = resourceID)
      
      seq_by_path <- seqs |> 
        semi_join(base, by = c("resourceID" = "measurementID")) |>
        group_by(resourceID) |>
        summarise(sequence_refs = paste0(
          "<a href='", relatedResourceID, "' target='_blank'>", accession_primary, "</a>",
          collapse = "<br>"), .groups = "drop") |>
        rename(measurementID = resourceID)
      
      df <- base |>
        left_join(seq_by_host, by = "occurrenceID") |>
        left_join(seq_by_path, by = "measurementID", suffix = c("_host", "_path")) |>
        mutate(
          host_sequence_refs     = coalesce(sequence_refs_host, "None linked"),
          pathogen_sequence_refs = coalesce(sequence_refs_path, "None linked")
        ) |>
        select(-any_of(c("sequence_refs_host", "sequence_refs_path")))
      
      # Jitter exact-duplicate coordinates so overlapping points stay
      # individually selectable once clustering disables at high zoom.
      df <- df |>
        group_by(decimalLatitude, decimalLongitude) |>
        mutate(
          n_at_point = n(),
          jitter_lat = decimalLatitude + if_else(n_at_point > 1, rnorm(n(), 0, 0.0008), 0),
          jitter_lng = decimalLongitude + if_else(n_at_point > 1, rnorm(n(), 0, 0.0008), 0)
        ) |>
        ungroup()
      
      return(df)
    })
    
    output$map_plot <- renderLeaflet({
      df <- map_data()
      
      m <- leaflet() |> 
        addProviderTiles(providers$CartoDB.Positron) |> 
        setView(lng = 15, lat = 5, zoom = 3)
      
      if (nrow(df) == 0) return(m)
      
      var_to_colour <- input$color_var
      
      mapping_vector <- as.character(df[[var_to_colour]])
      mapping_vector[is.na(mapping_vector) | mapping_vector %in% c("NA", "")] <- "Not Specified"
      df[[var_to_colour]] <- mapping_vector 
      
      if (var_to_colour == "positivity_status") {
        status_levels <- c("Positive Detection", "Negative/No Detection", "Not Specified")
        df[[var_to_colour]] <- factor(df[[var_to_colour]], levels = status_levels)
        pal <- colorFactor(
          palette = c("#e74c3c", "#bdc3c7", "#808080"),
          domain = df[[var_to_colour]]
        )
      }else {
        unique_vals <- unique(mapping_vector)
        n_unique <- length(unique_vals)
        expanding_colours <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(12, max(3, n_unique)), "Paired"))(n_unique)
        pal <- colorFactor(palette = expanding_colours, domain = unique_vals, na.color = "#808080")
      }
      
      m |> 
        addCircleMarkers(
          data = df,
          lng = ~jitter_lng,
          lat = ~jitter_lat,
          radius = 6,
          fillColor = ~pal(get(var_to_colour)), 
          fillOpacity = 0.8,
          color = "white", 
          weight = 1,
          clusterOptions = markerClusterOptions(
            maxClusterRadius = 40, 
            spiderfyOnMaxZoom = TRUE,
            disableClusteringAtZoom = 12,
            spiderfyDistanceMultiplier = 1.5
          ),
          popup = ~paste0(
            "<div style='font-family: sans-serif; font-size: 13px;'>",
            "<strong style='color: #2c3e50; font-size: 15px;'>", coalesce(scientificName, "Not Specified"), "</strong><br>",
            "<em>(", coalesce(genus, "Not Specified"), ")</em><hr style='margin: 5px 0;'>",
            "<b>Event:</b> ", coalesce(eventID, "NA"), " ", event_link, "<br>",
            "<b>Date / Year:</b> ", coalesce(eventDate, as.character(year), "Not Specified"), "<br>",
            "<b>Pathogen Tested:</b> ", display_pathogen, "<br>", 
            "<b>Test Result:</b> <span style='color: ", ifelse(positivity_status == "Positive Detection", "red", "black"), ";'>", positivity_status, "</span><br>",
            "<b>Number Tested:</b> ", coalesce(number_tested, 0), "<br>",
            "<b>Number Positive:</b> ", coalesce(measurementValue, 0),
            "<hr style='margin: 5px 0;'>",
            "<div style='font-size: 11px; color: #7f8c8d; line-height: 1.2;'>",
            "<b>Host ID:</b> ", coalesce(occurrenceID, "NA"), "<br>",
            "<b>Pathogen ID:</b> ", coalesce(measurementID, "NA"), "<br>",
            "<b>Host Sequence(s):</b> ", host_sequence_refs, "<br>",
            "<b>Pathogen Sequence(s):</b> ", pathogen_sequence_refs,
            "</div></div>"
          )
        ) |> 
        addLegend(position = "bottomright", pal = pal, values = df[[var_to_colour]],
                  title = tools::toTitleCase(gsub("_", " ", var_to_colour)), opacity = 1)
    })
  })
}