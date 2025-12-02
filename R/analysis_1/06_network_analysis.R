# ---
# 05_network_analysis.R
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
db_path <- here("data", "database", "Project_ArHa_database_2025-09-25.rds") # <-- Update this date
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
  filter(number_tested > 0) %>% 
  group_by(host_species, pathogen_species_cleaned, pathogen_family, assay) %>%
  summarise(number_tested = sum(number_tested, na.rm = TRUE),
            number_positive = sum(number_positive, na.rm = TRUE),
            .groups = "drop")

# a) Edge list for high-confidence, acute infections
edge_list_acute <- master_edge_list %>%
  filter(str_detect(assay, "PCR|Culture|Sequencing")) %>%
  group_by(host_species, pathogen_species_cleaned, pathogen_family) %>%
  summarise(number_tested = sum(number_tested, na.rm = TRUE),
            number_positive = sum(number_positive, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(proportion_positive = number_positive / number_tested)

# b) Edge list for all evidence, including serology
edge_list_all <- master_edge_list %>%
  group_by(host_species, pathogen_species_cleaned, pathogen_family) %>%
  summarise(number_tested = sum(number_tested, na.rm = TRUE),
            number_positive = sum(number_positive, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(proportion_positive = number_positive / number_tested)

# --- 3. Main Analysis Function ---

# This helper function encapsulates the entire analysis and plotting pipeline
# to be run on different subsets of the data.
run_network_analysis_and_plot <- function(edge_df, family_name, analysis_type) {
  
  message(paste("\n--- Running analysis for:", family_name, "-", analysis_type, "---"))
  
  # --- a) Network Construction ---
  
  # Filter for the specific family
  network_edges <- edge_df %>%
    filter(pathogen_family == family_name) %>%
    filter(number_positive > 0) %>% 
    select(from = host_species, to = pathogen_species_cleaned, number_tested, proportion_positive)
  
  # Check if there are any edges to plot
  if (nrow(network_edges) == 0) {
    warning(paste("No edges found for", family_name, "-", analysis_type, ". Skipping analysis."))
    return(NULL)
  }
  
  # Create incidence matrix for bipartite analysis
  incidence_table <- table(network_edges$from, network_edges$to)
  # Create a true matrix, ensuring row/col names are preserved
  incidence_matrix <- matrix(
    as.numeric(incidence_table), 
    nrow = nrow(incidence_table), 
    ncol = ncol(incidence_table), 
    dimnames = dimnames(incidence_table)
  )
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
  
  # 3. Correlate degree with sampling effort
  host_sampling_effort <- pathogen_data %>%
    left_join(host_data %>% select(host_record_id, host_species), by = "host_record_id") %>%
    filter(!is.na(host_species)) %>%
    distinct(host_record_id, host_species, pathogen_species_cleaned, number_tested, assay)
  
  if(analysis_type == "Acute_Evidence") {
    host_sampling_effort <- host_sampling_effort %>%
      filter(str_detect(assay, "PCR|Culture|Sequencing"))
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
  
  # Create an explicit node list with the required LOGICAL 'type' column
  host_nodes <- tibble(name = unique(network_edges$from), type = TRUE)
  pathogen_nodes <- tibble(name = unique(network_edges$to), type = FALSE)
  all_nodes <- bind_rows(host_nodes, pathogen_nodes)
  
  igraph_net <- graph_from_data_frame(d = network_edges, vertices = all_nodes, directed = FALSE)
  
  tidy_net <- as_tbl_graph(igraph_net) %>%
    mutate(node_type = if_else(type, "Host", "Pathogen"),
           degree = centrality_degree(),
           component = group_components())
  
  # Create a manual layout
  node_data <- tidy_net %>%
    activate(nodes) %>%
    as_tibble()
  
  # Proportional wedge layout logic
  component_sizes <- node_data %>%
    count(component, name = "n_nodes") %>%
    mutate(proportion = n_nodes / sum(n_nodes)) %>%
    arrange(desc(n_nodes))
  
  total_gap_space <- 0.1 * nrow(component_sizes)
  
  component_angles <- component_sizes %>%
    mutate(wedge_size = proportion * (2 * pi - total_gap_space),
           end_angle_raw = cumsum(wedge_size),
           start_angle_raw = lag(end_angle_raw, default = 0),
           start_angle = start_angle_raw + (row_number() - 1) * 0.1,
           end_angle = end_angle_raw + (row_number() - 1) * 0.1)
  
  node_data <- node_data %>% left_join(component_angles, by = "component")
  
  pathogen_coords <- node_data %>%
    filter(node_type == "Pathogen") %>%
    group_by(component) %>% arrange(name) %>%
    mutate(angle = scales::rescale(row_number(), to = c(start_angle[1], end_angle[1])), x = cos(angle) * 1, y = sin(angle) * 1) %>% ungroup()
  
  host_target_angles <- network_edges %>%
    left_join(pathogen_coords %>% select(name, pathogen_angle = angle), by = c("to" = "name")) %>%
    group_by(from) %>% summarise(target_angle = atan2(mean(sin(pathogen_angle)), mean(cos(pathogen_angle))), .groups = "drop") %>% rename(name = from)
  
  # --- Conditional Layout Logic ---
  if (analysis_type == "All_Evidence") {
    message("... using flexible cloud layout")
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
    message("... using structured proportional wedge layout")
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
  
  manual_layout <- create_layout(tidy_net, layout = "manual", 
                                 x = ordered_node_data$x, 
                                 y = ordered_node_data$y)
  
  # Create separate data frames for the labels
  host_labels <- manual_layout %>% as_tibble() %>% filter(node_type == 'Host')
  pathogen_labels <- manual_layout %>% 
    as_tibble() %>%
    filter(node_type == 'Pathogen', degree > 1) %>%
    mutate(name = stringr::str_replace(name, " ", "\n"))
  
  # Create the plot
  p <- ggraph(manual_layout) +
    geom_edge_link2(aes(edge_width = number_tested, edge_colour = proportion_positive), edge_alpha = 0.9) +
    scale_edge_colour_viridis(name = "Proportion Positive", labels = scales::percent, alpha = 0.9, begin = 0, end = 1, direction = -1,
                             discrete = FALSE, option = "E", aesthetics = "edge_colour", guide = "edge_colourbar") +
    scale_edge_width_continuous(name = "N Tested", range = c(0.25, 2)) + 
    ggnewscale::new_scale_color() +
    geom_node_point(aes(colour = node_type, size = degree)) +
    geom_text_repel(data = host_labels,
                    aes(x = x, y = y, label = name),
                    size = 2.5,
                    max.overlaps = Inf) +
    geom_text_repel(data = pathogen_labels,
                    aes(x = x, y = y, label = name),
                    size = 3, fontface = "bold",
                    max.overlaps = Inf,
                    nudge_x = -pathogen_labels$x * 0.4,
                    nudge_y = -pathogen_labels$y * 0.4) +
    scale_color_manual(values = c("Host" = "darkgreen", "Pathogen" = "firebrick"), name = "Node Type") +
    scale_size_continuous(range = c(2, 10), name = "Degree") +
    labs(
      title = paste(family_name, "Host-Pathogen Network"),
      subtitle = paste(str_to_title(str_replace(analysis_type, "_", " ")), "| Inner: Pathogens, Outer: Hosts.")
    ) +
    theme_graph(base_family = 'sans') +
    coord_fixed()
  
  # Save the plot
  ggsave(
    here("output", "analysis_1", paste0("network_plot_", tolower(family_name), "_", analysis_type, ".png")),
    plot = p, width = 20, height = 16, dpi = 300
  )
  
  # --- d) Return a list of all key outputs ---
  # This allows us to access the network objects later for further analysis (e.g., ERGMs)
  output_list <- list(
    metrics = metrics_df,
    network_matrix = incidence_matrix,
    tidy_network = tidy_net,
    plot = p
  )
  
  return(output_list)
}

# --- 4. Run Analyses ---

network_analysis_results <- list()

network_analysis_results$arena_acute <- run_network_analysis_and_plot(
  edge_df = edge_list_acute, 
  family_name = "Arenaviridae", 
  analysis_type = "Acute_Evidence"
)

network_analysis_results$arena_all <- run_network_analysis_and_plot(
  edge_df = edge_list_all, 
  family_name = "Arenaviridae", 
  analysis_type = "All_Evidence"
)

network_analysis_results$hanta_acute <- run_network_analysis_and_plot(
  edge_df = edge_list_acute, 
  family_name = "Hantaviridae", 
  analysis_type = "Acute_Evidence"
)

network_analysis_results$hanta_all <- run_network_analysis_and_plot(
  edge_df = edge_list_all, 
  family_name = "Hantaviridae", 
  analysis_type = "All_Evidence"
)

write_rds(network_analysis_results, here("output", "analysis_1", "networks.rds"))

# --- 5. Combined analysis ---

combined_matrix_table <- table(master_edge_list$host_species, master_edge_list$pathogen_species_cleaned)

combined_matrix <- matrix(
  as.numeric(combined_matrix_table), 
  nrow = nrow(combined_matrix_table), 
  dimnames = dimnames(combined_matrix_table)
)

combined_matrix[combined_matrix > 1] <- 1

combined_modules <- bipartite::computeModules(combined_matrix)

modularity_score <- combined_modules@likelihood
module_assignments_vector <- combined_modules@modules[1, ]

n_hosts <- nrow(combined_matrix)
n_pathogens <- ncol(combined_matrix)

host_module_ids_ordered <- module_assignments_vector[1:n_hosts]
pathogen_module_ids_ordered <- module_assignments_vector[(n_hosts + 1):(n_hosts + n_pathogens)]

host_names_original <- rownames(combined_matrix)
pathogen_names_original <- colnames(combined_matrix)

host_order_indices <- combined_modules@orderA
pathogen_order_indices <- combined_modules@orderB

host_modules <- tibble(host_species = host_names_original,
                       order = host_order_indices) %>%
  arrange(order) %>%
  mutate(module_id = host_module_ids_ordered) %>%
  select(host_species, module_id)

pathogen_modules <- tibble(pathogen_species = pathogen_names_original,
                           order = pathogen_order_indices) %>%
  arrange(order) %>%
  mutate(module_id = pathogen_module_ids_ordered) %>%
  select(pathogen_species, module_id)

pathogen_family_lookup <- master_edge_list %>%
  distinct(pathogen_species_cleaned, pathogen_family)

pathogen_module_composition <- pathogen_modules %>%
  left_join(pathogen_family_lookup, by = c("pathogen_species" = "pathogen_species_cleaned"))

module_summary <- pathogen_module_composition %>%
  group_by(module_id) %>%
  count(pathogen_family, name = "n_species") %>%
  tidyr::pivot_wider(
    names_from = "pathogen_family",
    values_from = "n_species",
    values_fill = 0
  ) %>%
  mutate(total_species = Arenaviridae + Hantaviridae) %>%
  filter(total_species >= 1) %>%
  arrange(module_id)

nrow(module_summary)
nrow(module_summary %>% filter(Hantaviridae > 0 & Arenaviridae == 0))
nrow(module_summary %>% filter(Hantaviridae == 0 & Arenaviridae > 0))
nrow(module_summary %>% filter(Hantaviridae > 0 & Arenaviridae > 0))
