# ArHa_App/R/mod_table.R
library(shiny)
library(bslib)
library(dplyr)
library(DT)
library(arrow)
library(stringr)
library(zip) 

# --- UI ---
mod_table_ui <- function(id) {
  ns <- NS(id)
  tagList(
    card(
      class = "mb-3",
      card_header(icon("server"), " Dataset Export & Quality Metrics", class = "bg-dark text-white"),
      card_body(
        layout_columns(
          col_widths = c(9, 3), 
          layout_columns(
            fill = FALSE,
            value_box(title = "Total Records", value = textOutput(ns("metric_total")), showcase = icon("database"), theme = "primary"),
            value_box(title = "Precise Spatial", value = textOutput(ns("metric_spatial")), showcase = icon("location-crosshairs"), theme = "success", tags$p("Resolved to 'site' level", class = "mb-0", style = "font-size: 0.8em;")),
            value_box(title = "Precise Temporal", value = textOutput(ns("metric_temporal")), showcase = icon("calendar-day"), theme = "info", tags$p("Resolved to specific day(s)", class = "mb-0", style = "font-size: 0.8em;")),
            value_box(title = "Species-Level ID", value = textOutput(ns("metric_taxa")), showcase = icon("dna"), theme = "warning", tags$p("Host identified to species", class = "mb-0", style = "font-size: 0.8em;"))
          ),
          card(
            class = "border-0 shadow-none bg-light h-100",
            card_body(
              class = "d-flex flex-column justify-content-center align-items-center gap-2",
              tags$h6("Export Filtered Dataset", class = "text-center mb-1"),
              downloadButton(ns("dl_zip_csv"), "Download DwC-A (CSV Zip)", class = "btn-outline-primary w-100"),
              downloadButton(ns("dl_zip_parquet"), "Download DwC-A (Parquet Zip)", class = "btn-outline-success w-100"),
              tags$small("Relational archive containing Events, Hosts, Pathogens & Sequences.", class = "text-muted text-center mt-1")
            )
          )
        )
      )
    ),
    card(
      full_screen = TRUE,
      card_header(icon("table"), " Data Preview (First 1,000 Rows)"),
      card_body(class = "p-0", DTOutput(ns("data_preview")))
    )
  )
}

# --- Server ---
mod_table_server <- function(id, filtered_data) {
  moduleServer(id, function(input, output, session) {
    
    # Calculate summary metrics
    observe({
      metrics <- filtered_data() |> 
        summarise(
          total = n(),
          precise_coords = sum(as.integer(coordinate_resolution_processed == "site"), na.rm = TRUE),
          precise_time = sum(as.integer(temporal_resolution %in% c("full_date", "day_range_resolution")), na.rm = TRUE),
          species_id = sum(as.integer(taxonRank == "species"), na.rm = TRUE)
        ) |> 
        collect()
      
      output$metric_total <- renderText({ formatC(metrics$total, format = "d", big.mark = ",") })
      calc_pct <- function(count, total) { if(total > 0) paste0(formatC(count, format = "d", big.mark = ","), " (", round((count/total)*100, 1), "%)") else "0 (0%)" }
      output$metric_spatial <- renderText({ calc_pct(metrics$precise_coords, metrics$total) })
      output$metric_temporal <- renderText({ calc_pct(metrics$precise_time, metrics$total) })
      output$metric_taxa <- renderText({ calc_pct(metrics$species_id, metrics$total) })
    })
    
    # Render preview
    output$data_preview <- renderDT({
      id_notify <- showNotification("Generating data preview...", duration = NULL, type = "message")
      on.exit(removeNotification(id_notify), add = TRUE)
      
      df_preview <- filtered_data() |> 
        head(1000) |> 
        select(occurrenceID, measurementID, eventID, year, eventDate, country, 
               scientificName, genus, scientificName_pathogen, genus_pathogen, family_pathogen, 
               measurementRemarks, pathogen_species_cleaned, measurementMethod, measurementValue, number_tested,
               decimalLatitude, decimalLongitude, coordinate_resolution_processed, temporal_resolution) |> 
        collect() |> 
        mutate(verbatim_raw = trimws(str_remove(str_extract(measurementRemarks, "Verbatim pathogen: [^|]+"), "Verbatim pathogen: ")),
               verbatim_pathogen = if_else(verbatim_raw %in% c("NA", "", "N/A"), NA_character_, verbatim_raw),
               sci_path_clean = if_else(scientificName_pathogen %in% c("NA", ""), NA_character_, scientificName_pathogen),
               gen_path_clean = if_else(genus_pathogen %in% c("NA", ""), NA_character_, genus_pathogen),
               fam_path_clean = if_else(family_pathogen %in% c("NA", ""), NA_character_, family_pathogen),
               cleaned_path_clean = if_else(pathogen_species_cleaned %in% c("NA", ""), NA_character_, pathogen_species_cleaned),
               Display_Pathogen = coalesce(sci_path_clean, gen_path_clean, fam_path_clean, verbatim_pathogen, cleaned_path_clean, "Not Specified")) |> 
        select(occurrenceID, measurementID, eventID, year, eventDate, country, 
               scientificName, genus, Display_Pathogen, measurementMethod, measurementValue, number_tested,
               decimalLatitude, decimalLongitude, coordinate_resolution_processed, temporal_resolution)
      
      datatable(df_preview, options = list(pageLength = 15, scrollX = TRUE, dom = 'tp'), rownames = FALSE, class = "display compact")
    })
    
    # Zip export
    generate_zip <- function(file, format = "csv") {
      id_notify <- showNotification("Packaging relational database. This may take a moment...", duration = NULL, type = "message")
      on.exit(removeNotification(id_notify))
      
      # Extract valid keys from the filter
      keys <- filtered_data() |> 
        select(eventID, occurrenceID, measurementID) |> 
        distinct() |> 
        collect()
      
      # Subset DuckDB tables 
      exp_events <- tbl_events |> filter(eventID %in% local(keys$eventID)) |> collect()
      exp_hosts <- tbl_hosts |> filter(occurrenceID %in% local(keys$occurrenceID)) |> collect()
      exp_paths <- tbl_pathogens |> filter(measurementID %in% local(keys$measurementID)) |> collect()
      
      # For sequences, we need to handle that resourceID could map to either occurrenceID OR measurementID
      exp_seqs <- tbl_sequences |> 
        filter(resourceID %in% local(keys$occurrenceID) | resourceID %in% local(keys$measurementID)) |> 
        collect()
      
      # Setup temporary directory
      tmp_dir <- tempdir()
      fs <- c("events", "hosts", "pathogens", "sequences")
      ext <- paste0(".", format)
      file_paths <- file.path(tmp_dir, paste0(fs, ext))
      readme_path <- file.path(tmp_dir, "README.txt")
      
      # Write data files
      if (format == "csv") {
        readr::write_csv(exp_events, file_paths[1])
        readr::write_csv(exp_hosts, file_paths[2])
        readr::write_csv(exp_paths, file_paths[3])
        readr::write_csv(exp_seqs, file_paths[4])
      } else {
        arrow::write_parquet(exp_events, file_paths[1])
        arrow::write_parquet(exp_hosts, file_paths[2])
        arrow::write_parquet(exp_paths, file_paths[3])
        arrow::write_parquet(exp_seqs, file_paths[4])
      }
      
      # Generate and write README
      readme_text <- paste0(
        "ArHa Database Explorer - Custom Export\n",
        "Generated: ", Sys.time(), "\n",
        "========================================\n\n",
        
        "CITATION & ATTRIBUTION\n",
        "========================================\n",
        "When utilising this dataset, please cite both the primary ArHa database and the original constituent studies.\n\n",
        "1. PRIMARY CITATION\n",
        "[Placeholder: Authors] (2026). The ArHa Database... DOI: [TBD]\n\n",
        "2. CONSTITUENT CITATIONS\n",
        "The citations for the specific records in this export are preserved in the `associatedReferences` column of the `events` table.\n",
        "Extract them in R using:\n",
        "events |> dplyr::filter(!is.na(associatedReferences)) |> dplyr::distinct(associatedReferences) |> dplyr::pull()\n\n",
        
        "========================================\n",
        "DATA DICTIONARY\n",
        "========================================\n",
        "This Darwin Core-compliant archive contains four relational tables linked by ID fields.\n\n",
        "1. events", ext, "\n",
        "   - eventID: Primary manuscript identifier.\n",
        "   - associatedReferences: Full citation for primary manuscript. \n",
        "   - year: Publication year of primary manuscript. \n",
        "   - samplingProtocol: Sampling effort reported in primary manuscript. Trap nights - The number of trap nights is reported, Time period - The sampling period is reported, Sessions - The number of sampling sessions is reported, Unspecified - An alternative measure of effort is reported, Not reported - No measure of effort is reported.\n",
        "   - sampleSizeValue: Measure of effort associated with samplingProtocol e.g., Trap nights\n",
        "   - samplingEffort: Number of sampling sessions. \n",
        "   - dataAccessLevel: Whether individual level, summary level or an unspecified level of data is available. \n",
        "   - dataResolution: What level of data is included in the records from this study in the dataset. (e.g., study level, site-session). \n\n",
        "2. hosts", ext, "\n",
        "   - occurrenceID: Primary host identifier. Links to resourceID in sequences table.\n",
        "   - eventID: Foreign key linking to the events table.\n",
        "   - scientificName: Resolved species name (if available). \n",
        "   - genus: Resolved genus name (if available). \n",
        "   - family: Resolved family name (if available). \n",
        "   - order: Resolved order name (if available). \n",
        "   - class: Resolved class name (if available). \n",
        "   - taxonRank: Taxon rank of record. \n",
        "   - taxonID: GBIFID for the record at the highest taxonomic resolution. \n",
        "   - individualCount: Number of individual animals detected in the record. \n",
        "   - eventDate: Period in which the record was detected. \n",
        "   - country: Country name(s) associated with the record. \n",
        "   - countryCode: ISO3C code(s) associated with the record. \n",
        "   - stateProvince: Administrative level 1 code(s) associated with the record. \n",
        "   - county: Administrative level 2 code(s) associated with the record. \n",
        "   - municipality: Administrative level 3 code(s) associated with the record. \n",
        "   - decimalLatitude: Latitude of the record, decimal degrees. \n",
        "   - decimalLongitude: Longitude of the record, decimal degrees. \n",
        "   - verbatimIdentification: Reported detection in primary manuscript prior to harmonisation to GBIF. \n",
        "   - verbatimLocality: Reported location of detection in primary manuscript. \n",
        "   - verbatimEventDate: Reported dates of detection in primary manuscript. \n",
        "   - coord_status: QC column reporting whether coordinates are valid or missing for this record. \n",
        "   - dist_from_expected: QC column reporting how far the coordinates are from locations reported in the primary manuscript. \n",
        "   - coordinate_resolution_processed: Level of spatial resolution coordinates are for (e.g., country, adm2, village or site). \n",
        "   - temporal_resolution: Temporal resolution of the record (e.g., full_date, publication_derived, year_only). \n",
        "   - trap_nights_clean: The number of trap nights of effort associated with this record. \n",
        "   - trap_nights_status: Descriptor of the information for trapping effort e.g., nights, trapnight_reported, not_reported. \n\n",
        "3. pathogens", ext, "\n",
        "   - measurementID: Primary assay/pathogen identifier. Links to resourceID in sequences table.\n",
        "   - occurrenceID: Foreign key linking to the hosts table.\n",
        "   - eventID: Foreign key linking to the events table.\n",
        "   - measurementType: Pathogens being detected in the reported assay.\n",
        "   - measurementValue: Number of individuals testing positive in the assay.\n",
        "   - measurementMethod: Method for detection in the assay.\n",
        "   - measurementRemarks: Further information about the assay including the number of individuals tests, the verbatim pathogen reported in the manuscript, the assay type reported in the manuscrip and any associated notes.\n",
        "   - scientificName_pathogen: Harmonised pathogen name at species level.\n",
        "   - genus_pathogen: Harmonised pathogen name at genus level.\n",
        "   - family_pathogen: Harmonised pathogen name at family level.\n",
        "   - taxonRank_pathogen: The taxon rank of the harmonised pathogen name.\n",
        "   - taxonID_pathogen: NCBI ID of the harmonised pathogen name.\n",
        "   - dynamicProperties: Nested column containing additional information on the pathogen taxonomy as applicable.\n",
        "   - number_tested: Number of individuals assayed in this record.\n",
        "   - number_negative: Number of individuals assayed in this record found to be negative for the pathogen.\n",
        "   - number_inconclusive: Number of individuals assayed with inconclusive assay results in this record.\n",
        "   - tested_detected: Binary QC column indicating whether the number reported as tested is greater than the number of reported positive assays.\n",
        "   - positive_tested: Binary QC column indicating whether the number reported as positive is greater than the number of reported individuals assayed.\n",
        "   - pathogen_species_cleaned: Pathogen name prior to harmonisation to NCBI.\n",
        "   - pathogen_family_original: Pathogen family name reported in the primary manuscript.\n\n",
        "4. sequences", ext, "\n",
        "   - resourceRelationshipID: Primary sequence relationship identifier.\n",
        "   - resourceID: Foreign key linking to occurrenceID in the hosts table or measurementID in the pathogen table.\n",
        "   - relatedResourceID: URL associated with the accession of the record on NCBI GenBank.\n",
        "   - relationshipOfResource: Whether this record is a host associated sequence or pathogen associated.\n",
        "   - sequence_type: Whether this record is a host associated sequence or pathogen associated.\n",
        "   - query_accession: Accession extracted from primary manuscript.\n",
        "   - accession_primary: Processed accession related to the primary manuscript.\n"
      )
      writeLines(readme_text, readme_path)
      
      # Bundle into ZIP
      all_files <- c(file_paths, readme_path)
      zip::zip(zipfile = file, files = all_files, mode = "cherry-pick")
    }
    
    # Download handlers
    output$dl_zip_csv <- downloadHandler(
      filename = function() { paste0("ArHa_DwC_CSV_", format(Sys.Date(), "%Y%m%d"), ".zip") },
      content = function(file) { generate_zip(file, format = "csv") }
    )
    
    output$dl_zip_parquet <- downloadHandler(
      filename = function() { paste0("ArHa_DwC_Parquet_", format(Sys.Date(), "%Y%m%d"), ".zip") },
      content = function(file) { generate_zip(file, format = "parquet") }
    )
    
  })
}