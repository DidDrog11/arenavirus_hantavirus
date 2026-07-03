# ArHa_App/R/mod_dashboard.R
library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(stringr)
library(DT)
library(ggplot2)

# --- UI ---
mod_dashboard_ui <- function(id) {
  ns <- NS(id)
  tagList(
    
    # Publications Summary Table (Open by default)
    accordion(
      open = "pubs", 
      accordion_panel(title = "Publications Summary", value = "pubs", icon = icon("book"), DTOutput(ns("table_pubs")))
    ),
    tags$br(),
    
    # Temporal Effort & Host Composition (closed by default)
    layout_columns(col_widths = c(6, 6), fill = FALSE,
                   accordion(open = FALSE, accordion_panel(title = "Temporal Sampling Effort", value = "time", icon = icon("chart-bar"), plotlyOutput(ns("plot_time"), height = "350px"))),
                   accordion(open = FALSE, accordion_panel(title = "Host Composition", value = "hosts", icon = icon("paw"), navset_pill(nav_panel("Genus", plotlyOutput(ns("plot_hosts_genus"), height = "350px")), nav_panel("Species", plotlyOutput(ns("plot_hosts_species"), height = "350px")))))
    ),
    tags$br(),
    
    # Pathogen Detections & Matrix (closed by default)
    layout_columns(col_widths = c(5, 7), fill = FALSE,
                   accordion(open = FALSE, accordion_panel(title = "Pathogen Detections", value = "paths", icon = icon("virus"), navset_pill(nav_panel("By Species", plotlyOutput(ns("plot_pathogens"), height = "400px")), nav_panel("By Modality", plotlyOutput(ns("plot_modality"), height = "400px"))))),
                   accordion(open = FALSE, accordion_panel(title = "Host-Pathogen Interaction Matrix", value = "matrix", icon = icon("circle-nodes"), plotlyOutput(ns("plot_matrix"), height = "400px")))
    ),
    tags$br(),
    
    # Record-Level Prevalence Scatter (closed by default)
    accordion(
      open = FALSE, 
      accordion_panel(title = "Record-Level Prevalence by Diagnostic Modality", value = "prev", icon = icon("vial"), plotlyOutput(ns("plot_prevalence"), height = "550px"))
    )
  )
}

# --- Server ---
mod_dashboard_server <- function(id, filtered_data) {
  moduleServer(id, function(input, output, session) {
    
    # Publications Summary Table
    output$table_pubs <- renderDT({
      df <- filtered_data() |> filter(!is.na(associatedReferences)) |> group_by(associatedReferences) |> summarise(n_records = n(), n_host_genera = n_distinct(genus, na.rm = TRUE), n_host_species = n_distinct(scientificName, na.rm = TRUE), n_path_species = n_distinct(scientificName_pathogen, na.rm = TRUE), n_individuals_tested = sum(number_tested, na.rm = TRUE), n_assays = sum(as.integer(measurementValue), na.rm = TRUE)) |> collect()
      req(nrow(df) > 0)
      df <- df |> mutate(extracted_doi = str_extract(associatedReferences, "(?i)doi:\\s*(10\\.\\S+)"), clean_doi = str_remove(extracted_doi, "(?i)doi:\\s*"), safe_ref = str_replace_all(associatedReferences, "\"", "&quot;"), Reference = if_else(!is.na(clean_doi), paste0("<span title=\"", safe_ref, "\" style=\"cursor: help; border-bottom: 1px dotted #999;\">", str_trunc(associatedReferences, 45), "</span> <br><a href='https://doi.org/", clean_doi, "' target='_blank' style='font-size: 0.85em;'><i class='fa fa-external-link'></i> View Source</a>"), paste0("<span title=\"", safe_ref, "\" style=\"cursor: help; border-bottom: 1px dotted #999;\">", str_trunc(associatedReferences, 55), "</span>"))) |> select(Reference, n_records, n_host_genera, n_host_species, n_path_species, n_individuals_tested, n_assays) |> arrange(desc(n_records))
      datatable(df, escape = FALSE, rownames = FALSE, options = list(paging = FALSE, scrollY = "250px", scrollX = TRUE, scrollCollapse = TRUE, dom = 't'), colnames = c("Publication / Reference", "Records", "Host Genera", "Host Species", "Pathogen Species", "Indiv. Tested", "Indiv. Positive"))
    })
    
    # Temporal Sampling Effort
    output$plot_time <- renderPlotly({
      df <- filtered_data() |> filter(!is.na(year), !is.na(scientificName)) |> group_by(year, scientificName) |> summarise(n_records = n(), n_detections = sum(measurementValue, na.rm = TRUE), .groups = 'drop') |> arrange(year, desc(n_records)) |> collect()
      req(nrow(df) > 0)
      n_species <- length(unique(df$scientificName))
      dynamic_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(12, max(3, n_species)), "Paired"))(n_species)
      plot_ly(df, x = ~year, y = ~n_records, color = ~scientificName, colors = dynamic_palette, type = 'bar', hoverinfo = 'text', marker = list(line = list(width = 0)), text = ~paste0("<b>Year:</b> ", year, "<br><b>Species:</b> ", scientificName, "<br><b>Records:</b> ", n_records, "<br><b>Detections:</b> ", n_detections)) |> layout(barmode = 'stack', xaxis = list(title = "Year"), yaxis = list(title = "Records"), legend = list(orientation = "h", x = 0, y = -0.2), margin = list(t = 20, b = 60))
    })
    
    # Host Composition (Genus and Species)
    output$plot_hosts_genus <- renderPlotly({
      df <- filtered_data() |> filter(!is.na(genus)) |> group_by(genus) |> summarise(n_records = n(), n_detections = sum(measurementValue, na.rm = TRUE), .groups = 'drop') |> arrange(desc(n_records)) |> head(15) |> collect()
      req(nrow(df) > 0)
      df$genus <- factor(df$genus, levels = rev(df$genus))
      plot_ly(df, x = ~n_records, y = ~genus, type = 'bar', orientation = 'h', marker = list(color = '#2c3e50'), hoverinfo = 'text', text = ~paste0("<b>Genus:</b> ", genus, "<br><b>Records:</b> ", n_records, "<br><b>Detections:</b> ", n_detections)) |> layout(xaxis = list(title = "Records"), yaxis = list(title = ""), margin = list(t = 10))
    })
    
    output$plot_hosts_species <- renderPlotly({
      df <- filtered_data() |> filter(!is.na(scientificName)) |> group_by(scientificName) |> summarise(n_records = n(), n_detections = sum(measurementValue, na.rm = TRUE), .groups = 'drop') |> arrange(desc(n_records)) |> head(15) |> collect()
      req(nrow(df) > 0)
      df$scientificName <- factor(df$scientificName, levels = rev(df$scientificName))
      plot_ly(df, x = ~n_records, y = ~scientificName, type = 'bar', orientation = 'h', marker = list(color = '#18bc9c'), hoverinfo = 'text', text = ~paste0("<b>Species:</b> ", scientificName, "<br><b>Records:</b> ", n_records, "<br><b>Detections:</b> ", n_detections)) |> layout(xaxis = list(title = "Records"), yaxis = list(title = ""), margin = list(t = 10))
    })
    
    # Pathogen Detections (Species and Modality)
    output$plot_pathogens <- renderPlotly({
      df <- filtered_data() |> filter(!is.na(scientificName_pathogen)) |> group_by(scientificName_pathogen) |> summarise(n_tested = sum(number_tested, na.rm = TRUE), n_pos = sum(measurementValue, na.rm = TRUE), .groups = 'drop') |> mutate(prop_pos = if_else(n_tested > 0, n_pos / n_tested, 0)) |> arrange(desc(n_tested)) |> head(15) |> collect()
      req(nrow(df) > 0)
      df$scientificName_pathogen <- factor(df$scientificName_pathogen, levels = rev(df$scientificName_pathogen))
      plot_ly(df, x = ~n_tested, y = ~scientificName_pathogen, type = 'bar', orientation = 'h', marker = list(color = ~prop_pos, colorscale = 'Reds', cmin = 0, cmax = 1, colorbar = list(title = "Prev."), line = list(width = 0)), hoverinfo = 'text', text = ~paste0("<b>Pathogen:</b> ", scientificName_pathogen, "<br><b>Tested:</b> ", n_tested, "<br><b>Positive:</b> ", n_pos, " (", scales::percent(prop_pos, accuracy = 0.1), ")")) |> layout(xaxis = list(title = "Individuals Tested"), yaxis = list(title = ""), margin = list(t = 10))
    })
    
    output$plot_modality <- renderPlotly({
      df <- filtered_data() |> filter(!is.na(measurementMethod)) |> group_by(measurementMethod) |> summarise(n_tested = sum(number_tested, na.rm = TRUE), n_pos = sum(measurementValue, na.rm = TRUE), .groups = 'drop') |> arrange(desc(n_tested)) |> collect()
      req(nrow(df) > 0)
      df$measurementMethod <- factor(df$measurementMethod, levels = rev(df$measurementMethod))
      plot_ly(df, x = ~n_tested, y = ~measurementMethod, type = 'bar', orientation = 'h', marker = list(color = '#f39c12'), hoverinfo = 'text', text = ~paste0("<b>Modality:</b> ", measurementMethod, "<br><b>Tested:</b> ", n_tested, "<br><b>Positive:</b> ", n_pos)) |> layout(xaxis = list(title = "Individuals Tested"), yaxis = list(title = ""), margin = list(t = 10))
    })
    
    # Host-Pathogen Matrix
    output$plot_matrix <- renderPlotly({
      df <- filtered_data() |> filter(!is.na(scientificName), !is.na(scientificName_pathogen)) |> group_by(scientificName, scientificName_pathogen) |> summarise(n_tested = sum(number_tested, na.rm = TRUE), n_pos = sum(measurementValue, na.rm = TRUE), .groups = 'drop') |> mutate(prop_pos = if_else(n_tested > 0, n_pos / n_tested, 0)) |> filter(n_tested > 0) |> collect()
      req(nrow(df) > 0)
      max_test <- max(df$n_tested, na.rm = TRUE)
      sizeref <- ifelse(max_test > 0, 2.0 * max_test / (40**2), 1)
      plot_ly(data = df, x = ~scientificName_pathogen, y = ~scientificName, type = 'scatter', mode = 'markers', hoverinfo = 'text', text = ~paste0("<b>Host:</b> ", scientificName, "<br><b>Pathogen:</b> ", scientificName_pathogen, "<br><b>Tested:</b> ", formatC(n_tested, format = "d", big.mark = ","), "<br><b>Positive:</b> ", formatC(n_pos, format = "d", big.mark = ","), " (", scales::percent(prop_pos, accuracy = 0.1), ")"), marker = list(size = ~n_tested, color = ~prop_pos, colorscale = 'Reds', showscale = TRUE, colorbar = list(title = "Prev."), sizemode = 'area', sizeref = sizeref, line = list(width = 0, color = 'rgba(0,0,0,0)'))) |> layout(xaxis = list(title = "Pathogen Species", tickangle = 45), yaxis = list(title = "Host Species"), margin = list(b = 100, l = 100))
    })
    
    # Record-level prevalence
    output$plot_prevalence <- renderPlotly({
      id_notify <- showNotification("Calculating diagnostic matrix...", duration = NULL, type = "message")
      on.exit(removeNotification(id_notify), add = TRUE)
      
      df <- filtered_data() |> 
        filter(!is.na(number_tested), number_tested > 0, !is.na(measurementValue), measurementValue >= 0, !is.na(scientificName_pathogen)) |> 
        select(occurrenceID, measurementID, measurementMethod, measurementValue, number_tested, scientificName, genus, scientificName_pathogen) |> 
        collect() |> 
        mutate(modality = case_when(grepl("PCR|Molecular|Sequencing", measurementMethod, ignore.case = TRUE) ~ "Molecular / PCR", 
                                    grepl("Serolog|Antibody|ELISA|IFA|IgG|IgM", measurementMethod, ignore.case = TRUE) ~ "Serology", 
                                    TRUE ~ "Other"), 
               prevalence = (measurementValue / number_tested) * 100, 
               host_display = coalesce(scientificName, genus, "Unknown Host"), 
               hover_text = paste0("<b>Host:</b> ", coalesce(scientificName, genus, "Unknown Host"), "<br><b>Pathogen:</b> ", scientificName_pathogen, "<br><b>Prevalence:</b> ", round(prevalence, 1), "%<br><b>Tested:</b> ", number_tested, " | <b>Positive:</b> ", measurementValue)) |> 
        filter(modality != "Other", prevalence <= 100)
      
      req(nrow(df) > 0)
      
      # Setup Y-axis factor mapping and jitter
      path_levels <- sort(unique(df$scientificName_pathogen))
      df <- df |> mutate(path_num = as.numeric(factor(scientificName_pathogen, levels = path_levels)), y_jitter = path_num + runif(n(), -0.25, 0.25))
      
      # Calculate host-specific pooled prevalence
      df_pooled <- df |> 
        group_by(modality, scientificName_pathogen, path_num, host_display) |> 
        summarise(total_tested = sum(number_tested, na.rm = TRUE), total_pos = sum(measurementValue, na.rm = TRUE), .groups = 'drop') |> 
        mutate(pooled_prev = (total_pos / total_tested) * 100,
               hover_text = paste0("<b>Host:</b> ", host_display, "<br><b>Pathogen:</b> ", scientificName_pathogen, "<br><b>Pooled Prevalence:</b> ", round(pooled_prev, 1), "%<br><b>Total Tested:</b> ", formatC(total_tested, format="d", big.mark=","), " | <b>Total Positive:</b> ", formatC(total_pos, format="d", big.mark=","))) |> 
        filter(total_tested > 0)
      
      # Define global colour palette
      unique_hosts <- sort(unique(df$host_display))
      host_colours <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(12, max(3, length(unique_hosts))), "Paired"))(length(unique_hosts))
      names(host_colours) <- unique_hosts
      
      # Function for building facets
      build_subplot <- function(mod_name, show_leg) {
        d_sub <- df |> filter(modality == mod_name)
        p_sub <- df_pooled |> filter(modality == mod_name)
        
        if (nrow(d_sub) == 0) return(plot_ly(type = 'scattergl', mode = 'markers') |> layout(xaxis = list(title = mod_name, range = c(-5, 105), ticksuffix = "%")))
        
        plot_ly() |> 
          # Raw data points
          add_trace(data = d_sub, x = ~prevalence, y = ~y_jitter, color = ~host_display, colors = host_colours, type = 'scattergl', mode = 'markers', text = ~hover_text, hoverinfo = 'text', legendgroup = ~host_display, showlegend = show_leg, marker = list(size = 6, opacity = 0.4, line = list(width = 0))) |> 
          # Pooled prevalence diamond
          add_trace(data = p_sub, x = ~pooled_prev, y = ~path_num, color = ~host_display, colors = host_colours, type = 'scattergl', mode = 'markers', text = ~hover_text, hoverinfo = 'text', hoverlabel = list(bgcolor = "#1a1a1a", font = list(color = "white"), bordercolor = "white"), legendgroup = ~host_display, showlegend = FALSE, marker = list(symbol = 'diamond', size = 12, line = list(color = 'black', width = 1.5))) |> 
          layout(xaxis = list(title = mod_name, range = c(-5, 105), ticksuffix = "%"))
      }
      
      # Combine
      subplot(build_subplot("Molecular / PCR", TRUE), build_subplot("Serology", FALSE), shareY = TRUE, titleX = TRUE) |> 
        layout(yaxis = list(title = "", tickvals = 1:length(path_levels), ticktext = path_levels, zeroline = FALSE), legend = list(orientation = "h", x = 0, y = -0.2), margin = list(l = 150))
    })
    
  })
}