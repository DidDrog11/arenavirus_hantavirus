# ---
# 03_network_analysis.R
#
# Purpose: To construct, analyze, and visualize host-pathogen bipartite networks
# for Arenaviridae and Hantaviridae. This script performs a stratified analysis,
# creating separate networks for high-confidence "Acute Evidence" (PCR, Culture)
# and for "All Evidence" (including serology).
# ---

# --- 1. Setup ---
# Load necessary packages
library(here)
library(readr)
library(dplyr)
library(stringr)
library(bipartite) # For network analysis
library(igraph)    # For network object creation
library(tidygraph) # For tidy network manipulation
library(ggraph)    # For plotting networks with ggplot2
library(ggrepel)

# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-08-25.rds") # <-- Update this date
arha_db <- read_rds(db_path)
host_data <- arha_db$host
pathogen_data <- arha_db$pathogen

# --- 2. Create Stratified Edge Lists ---

# Create a master edge list that includes the assay type for filtering
master_edge_list <- pathogen_data %>%
  # remove records where the association can not be attributed to single species
  filter(!str_detect(ncbi_id, ",")) %>%
  # remove records that are not resolved at species level for pathogen
  filter(taxonomic_level == "species") %>%
  left_join(host_data %>% select(host_record_id, host_species), by = "host_record_id") %>%
  # remove records where host species is not resolved
  filter(!is.na(host_species)) %>%
  # For this analysis, we only need the presence of a link
  filter(number_positive >= 1) %>%
  distinct(host_species, pathogen_species_cleaned, pathogen_family, assay)

# a) Edge list for high-confidence, acute infections
edge_list_acute <- master_edge_list %>%
  filter(str_detect(assay, "PCR|Culture|Sequencing"))

# b) Edge list for all evidence, including serology
edge_list_all <- master_edge_list


# --- 3. Main Analysis Function ---

# This helper function encapsulates the entire analysis and plotting pipeline
# to be run on different subsets of the data.
run_network_analysis_and_plot <- function(edge_df, family_name, analysis_type) {
  
  message(paste("\n--- Running analysis for:", family_name, "-", analysis_type, "---"))
  
  # --- a) Network Construction ---
  
  # Filter for the specific family
  network_edges <- edge_df %>%
    filter(pathogen_family == family_name) %>%
    select(from = host_species, to = pathogen_species_cleaned)
  
  # Check if there are any edges to plot
  if (nrow(network_edges) == 0) {
    warning(paste("No edges found for", family_name, "-", analysis_type, ". Skipping analysis."))
    return(NULL)
  }
  
  # Create incidence matrix for bipartite analysis
  incidence_table <- table(network_edges$from, network_edges$to)
  incidence_matrix <- matrix(as.numeric(incidence_table), 
                             nrow = nrow(incidence_table), 
                             ncol = ncol(incidence_table),
                             dimnames = dimnames(incidence_table))
  incidence_matrix[incidence_matrix > 1] <- 1
  
  # --- b) Network Analysis ---
  
  # Calculate network-level metrics using the bipartite package
  network_metrics <- networklevel(incidence_matrix)
  
  # Save metrics
  metrics_df <- as_tibble(network_metrics, rownames = "metric")
  write_csv(metrics_df, here("output", "analysis_1", 
                             paste0("network_metrics_", tolower(family_name), "_", analysis_type, ".csv")))
  
  # --- b2) Descriptive Node-level Analysis ---
  
  # Create the tidygraph object to get node degrees
  tidy_net_for_metrics <- graph_from_data_frame(d = network_edges, directed = FALSE) %>%
    as_tbl_graph() %>%
    mutate(degree = centrality_degree())
  
  node_degrees <- tidy_net_for_metrics %>%
    activate(nodes) %>%
    as_tibble()
  
  # 1. Identify top connected hosts
  top_hosts <- node_degrees %>%
    filter(name %in% network_edges$from) %>% # Filter for hosts
    arrange(desc(degree))
  
  write_csv(top_hosts, here("output", "analysis_1",
                            paste0("top_hosts_", tolower(family_name), "_", analysis_type, ".csv")))
  
  # 2. Identify potential recombination hubs
  recombination_hubs <- network_edges %>%
    mutate(pathogen_genus = str_extract(to, "^\\w+")) %>%
    group_by(from, pathogen_genus) %>%
    summarise(n_species_in_genus = n_distinct(to), .groups = "drop") %>%
    filter(n_species_in_genus > 1) %>%
    arrange(desc(n_species_in_genus))
  
  write_csv(recombination_hubs, here("output", "analysis_1",
                                     paste0("recombination_hubs_", tolower(family_name), "_", analysis_type, ".csv")))
  
  # 3. Correlate degree with sampling effort
  # First, calculate total sampling effort per host species
  host_sampling_effort <- pathogen_data %>%
    left_join(host_data %>% select(host_record_id, host_species), by = "host_record_id") %>%
    filter(!is.na(host_species)) %>%
    distinct(host_record_id, host_species, pathogen_species_cleaned, number_tested, assay)
  
  if(analysis_type == "Acute_Evidence") {
    host_sampling_effort <- host_sampling_effort %>%
      filter(str_detect(assay, "PCR|Culture|Sequencing"))
  } else {
    host_sampling_effort <- host_sampling_effort
  }
     
  host_sampling_effort <- host_sampling_effort %>%
    group_by(host_record_id) %>%
    arrange(desc(number_tested)) %>%
    slice(1) %>%
    group_by(host_species) %>%
    summarise(number_tested = sum(number_tested, na.rm = TRUE), .groups = "drop")
    
  degree_vs_effort <- top_hosts %>%
    left_join(host_sampling_effort, by = c("name" = "host_species"))
  
  write_csv(degree_vs_effort, here("output", "analysis_1",
                                   paste0("degree_vs_effort_", tolower(family_name), "_", analysis_type, ".csv")))
  
  # --- c) Network Visualization ---
  host_nodes <- tibble(name = unique(network_edges$from), type = TRUE, node_type = "Host") # TRUE for first mode (hosts)
  pathogen_nodes <- tibble(name = unique(network_edges$to), type = FALSE, node_type = "Pathogen") # FALSE for second mode (pathogens)
  all_nodes <- bind_rows(host_nodes, pathogen_nodes)
  
  igraph_net <- graph_from_data_frame(d = network_edges, vertices = all_nodes, directed = FALSE)
  
  # Convert to tidygraph, create a character version of type for plotting, and calculate degree
  tidy_net <- as_tbl_graph(igraph_net) %>%
    mutate(
      degree = centrality_degree(),
      component = group_components() # Identify distinct network components
    )
  
  # Create a manual layout with concentric circles, ordered to minimize edge crossing
  node_data <- tidy_net %>%
    activate(nodes) %>%
    as_tibble()
  
  # Assign components to proportional wedges of the circle
  component_sizes <- node_data %>%
    count(component, name = "n_nodes") %>%
    arrange(desc(n_nodes))
  
  total_nodes <- sum(component_sizes$n_nodes)
  n_components <- nrow(component_sizes)
  gap_size <- 0.1 # Radians for the gap between wedges
  total_gap_space <- gap_size * n_components
  
  # Calculate the angle range for each component based on its proportional size
  component_angles <- component_sizes %>%
    mutate(
      proportion = n_nodes / total_nodes,
      wedge_size = proportion * (2 * pi - total_gap_space),
      end_angle_raw = cumsum(wedge_size),
      start_angle_raw = lag(end_angle_raw, default = 0),
      start_angle = start_angle_raw + (row_number() - 1) * gap_size,
      end_angle = end_angle_raw + (row_number() - 1) * gap_size
    )
  
  node_data <- node_data %>% left_join(component_angles, by = "component")
  
  # Calculate pathogen coordinates within their assigned wedge
  pathogen_coords <- node_data %>%
    filter(node_type == "Pathogen") %>%
    group_by(component) %>% arrange(name) %>%
    mutate(angle = scales::rescale(row_number(), to = c(start_angle[1], end_angle[1])), x = cos(angle) * 1, y = sin(angle) * 1) %>% ungroup()
  
  # Calculate the target angle for each host based on its connections
  host_target_angles <- network_edges %>%
    left_join(pathogen_coords %>% select(name, pathogen_angle = angle), by = c("to" = "name")) %>%
    group_by(from) %>% summarise(target_angle = atan2(mean(sin(pathogen_angle)), mean(cos(pathogen_angle))), .groups = "drop") %>% rename(name = from)
  
  # --- Conditional Layout Logic ---
  if (analysis_type == "All_Evidence" |
      analysis_type == "Acute_Evidence") {
    # For "All Evidence", use a more flexible "cloud" layout for hosts
    message("Using flexible cloud layout for All_Evidence plot...")
    
    host_coords <- node_data %>%
      filter(node_type == "Host") %>%
      left_join(host_target_angles, by = "name") %>%
      mutate(
        radius = runif(n(), min = 1.8, max = 2.2),
        angle = rnorm(n(), mean = coalesce(target_angle, 0), sd = 0.1),
        x = cos(angle) * radius,
        y = sin(angle) * radius
      )
    
  } else {
    # For "Acute Evidence", use the structured proportional wedge layout
    message("Using structured proportional wedge layout for Acute_Evidence plot...")
    
    host_coords <- node_data %>%
      filter(node_type == "Host") %>%
      left_join(host_target_angles, by = "name") %>%
      mutate(target_angle = coalesce(target_angle, 0)) %>%
      group_by(component) %>% arrange(target_angle) %>%
      mutate(angle = scales::rescale(row_number(), to = c(start_angle[1], end_angle[1])), x = cos(angle) * 2, y = sin(angle) * 2) %>% ungroup()
  }
  
  # Combine coordinates and join back to the original node data to ensure correct order
  all_coords <- bind_rows(host_coords, pathogen_coords)
  
  ordered_node_data <- node_data %>%
    left_join(all_coords %>% select(name, x, y), by = "name")
  
  # Create the layout using the correctly ordered coordinates
  manual_layout <- create_layout(tidy_net, layout = "manual", 
                                 x = ordered_node_data$x, 
                                 y = ordered_node_data$y)
  
  # Create separate data frames for the labels
  host_labels <- manual_layout %>% as_tibble() %>% filter(node_type == 'Host')
  pathogen_labels <- manual_layout %>% 
    filter(node_type == 'Pathogen', degree > 1) %>%
    # Add a line break after the first word (genus) for better spacing
    mutate(name = stringr::str_replace(name, " ", "\n"))
  
  # Create the plot using the manual layout
  p <- ggraph(manual_layout) +
    geom_edge_link(edge_colour = "grey", edge_alpha = 0.8, edge_width = 1) +
    geom_node_point(aes(colour = node_type, size = degree)) +
    # Add separate text layers for hosts (outward) and pathogens (inward)
    geom_text_repel(data = host_labels,
                    aes(x = x, y = y, label = name),
                    size = 2.5,
                    max.overlaps = Inf # Allow all labels to be placed
    ) +
    geom_text_repel(data = pathogen_labels,
                    aes(x = x, y = y, label = name),
                    size = 3, fontface = "bold",
                    max.overlaps = Inf,
                    # Nudge pathogen labels towards the center (0,0)
                    nudge_x = -pathogen_labels$x[pathogen_labels$degree > 1] * 0.4,
                    nudge_y = -pathogen_labels$y[pathogen_labels$degree > 1] * 0.4
    ) +
    scale_color_manual(values = c("Host" = "darkgreen", "Pathogen" = "firebrick"), name = "Node Type") +
    scale_size_continuous(range = c(2, 10), name = "Degree") +
    labs(
      title = paste(family_name, "Host-Pathogen Network"),
      subtitle = paste(analysis_type, "| Inner: Pathogens, Outer: Hosts.")
    ) +
    theme_graph(base_family = 'sans') +
    coord_fixed()
  
  # Save the plot
  ggsave(
    here("output", "analysis_1", paste0("network_plot_", tolower(family_name), "_", analysis_type, ".png")),
    plot = p, width = 20, height = 16, dpi = 300
  )
  
  return(p)
}

# --- 4. Run Analyses ---

# Run for Acute Evidence
arena_acute_plot <- run_network_analysis_and_plot(edge_list_acute, "Arenaviridae", "Acute_Evidence")
hanta_acute_plot <- run_network_analysis_and_plot(edge_list_acute, "Hantaviridae", "Acute_Evidence")

# Run for All Evidence
arena_all_plot <- run_network_analysis_and_plot(edge_list_all, "Arenaviridae", "All_Evidence")
hanta_all_plot <- run_network_analysis_and_plot(edge_list_all, "Hantaviridae", "All_Evidence")

# --- 5. Display Final Plots ---
print(arena_acute_plot)
print(hanta_acute_plot)
print(arena_all_plot)
print(hanta_all_plot)
