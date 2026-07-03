# ArHa_App/R/mod_facets.R
library(shiny)

mod_facets_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("facet_summary"))
}

mod_facets_server <- function(id, facet_data_reactive, title_text = "Items:") {
  moduleServer(id, function(input, output, session) {
    
    output$facet_summary <- renderUI({
      df <- facet_data_reactive()
      
      if (is.null(df) || nrow(df) == 0) return(tags$p(tags$small(tags$em("No matching records."))))
      
      col_name <- names(df)[1]
      
      tagList(
        tags$small(tags$b(title_text)),
        tags$div(
          style = "max-height: 150px; overflow-y: auto; margin-bottom: 10px; padding-right: 5px;",
          tags$ul(class = "list-unstyled ms-2 mb-0 text-muted", 
                  lapply(1:nrow(df), function(i) {
                    tags$li(tags$small(paste0(df[[col_name]][i], " (", df$n[i], ")")))
                  })
          )
        )
      )
    })
  })
}