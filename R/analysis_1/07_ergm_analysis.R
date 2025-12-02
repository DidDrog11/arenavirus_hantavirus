# ---
# 07_ergm_analysis.R
#
# Purpose: To fit Exponential Random Graph Models (ERGMs) to test hypotheses
# about the drivers of host-pathogen interactions (positive links), while
# controlling for sampling effort bias.
# ---

# --- 1. Setup ---
# Load necessary packages
library(here)
library(readr)
library(dplyr)
library(tidyr) 
library(stringr)
library(statnet)
library(tibble)
library(brms)
library(tidybayes)
library(modelr)

# Create output directory for ERGM results
output_dir_ergm <- here("output", "analysis_1")

# --- 2. Load Data ---
# Load the final database object
db_path <- here("data", "database", "Project_ArHa_database_2025-09-25.rds") # <-- Check date
arha_db <- read_rds(db_path)
host_data <- arha_db$host
pathogen_data <- arha_db$pathogen

# Load harmonized traits data
host_traits <- read_rds(here("data", "external", "host_traits_cleaned.rds")) %>%
  select(host_species = scientific_name,
         species_gbifid,
         range_size_km2,
         synanthropy_status = synanthropy,
         s_index) %>%
  distinct(host_species, .keep_all = TRUE)


# --- 3. Calculate Sampling Effort Covariate ---

# Calculate total number tested per host species across ALL assays
host_sampling_effort <- pathogen_data %>%
  left_join(host_data %>% select(host_record_id, host_species), by = "host_record_id") %>%
  filter(!is.na(host_species)) %>%
  # Sum tests across all pathogens/assays for each host record
  group_by(host_record_id, host_species) %>%
  summarise(record_tests = sum(number_tested, na.rm = TRUE), .groups = "drop") %>%
  # Now sum across all records for each host species
  group_by(host_species) %>%
  summarise(total_individuals_tested = sum(record_tests, na.rm = TRUE), .groups = "drop") %>%
  # Add log(x+1) version for the model
  mutate(log_sampling_effort = log1p(total_individuals_tested)) %>%
  # Ensure we only have hosts that exist in host_traits (for merging later)
  filter(host_species %in% host_traits$host_species)

# --- 4. Data Preparation Function for ERGM ---

prepare_ergm_data <- function(family_name, pathogen_data, host_sampling_effort, host_traits) {
  
  message(paste("--- Preparing data for:", family_name, "---"))
  
  # --- a) Identify Hosts and Pathogens for this family ---
  family_pathogens <- pathogen_data %>%
    filter(pathogen_family == family_name) %>%
    filter(!is.na(pathogen_species_cleaned)) %>%
    distinct(pathogen_species_cleaned) %>%
    pull(pathogen_species_cleaned)
  
  family_positives <- pathogen_data %>%
    filter(pathogen_family == family_name) %>%
    filter(number_positive > 0) %>%
    left_join(host_data %>% select(host_record_id, host_species), by = "host_record_id") %>%
    filter(!is.na(host_species)) %>%
    distinct(host_species, pathogen_species_cleaned)
  
  # Get all hosts that have *any* sampling effort recorded and are in traits
  all_sampled_hosts <- host_sampling_effort$host_species
  
  # --- b) Create the Binary Incidence Matrix ---
  
  # Create all possible dyads between sampled hosts and family pathogens
  all_dyads <- expand_grid(host_species = all_sampled_hosts,
                           pathogen_species_cleaned = family_pathogens)
  
  # Mark positive interactions
  interaction_data <- all_dyads %>%
    left_join(family_positives %>% mutate(positive = 1),
              by = c("host_species", "pathogen_species_cleaned")) %>%
    mutate(interaction = ifelse(is.na(positive), 0, 1)) %>%
    select(host_species, pathogen_species_cleaned, interaction)
  
  # Pivot to wide matrix format
  incidence_matrix <- interaction_data %>%
    pivot_wider(names_from = pathogen_species_cleaned,
                values_from = interaction,
                values_fill = 0) %>% # Fill missing dyads explicitly as 0
    column_to_rownames("host_species") %>%
    as.matrix()
  
  # --- c) Prepare Host Covariates ---
  
  # Filter traits and effort data for hosts present in the matrix rows
  matrix_hosts <- rownames(incidence_matrix)
  
  host_covariates <- host_traits %>%
    filter(host_species %in% matrix_hosts) %>%
    left_join(host_sampling_effort %>% select(host_species, log_sampling_effort), by = "host_species") %>%
    
    # Prepare range size (log, handle zeros/NA)
    mutate(log_range_size = log1p(range_size_km2)) %>%
    select(host_species, log_sampling_effort, log_range_size, synanthropy_status)
  
  # Handle NAs in covariates: ERGM cannot handle NAs. We must exclude hosts with NA.
  hosts_with_complete_range <- host_covariates %>%
    filter(!is.na(log_range_size)) %>%
    pull(host_species)
  
  hosts_with_complete_synanthropy <- host_covariates %>%
    filter(!is.na(synanthropy_status), synanthropy_status != "") %>% # Also exclude empty strings
    pull(host_species)
  
  # Align covariates with matrix order
  # For the "full" model (using range size)
  covariates_full <- host_covariates %>%
    filter(host_species %in% hosts_with_complete_range) %>%
    arrange(match(host_species, matrix_hosts))
  
  # For the "synanthropy" subset model
  covariates_synanthropy <- host_covariates %>%
    filter(host_species %in% hosts_with_complete_synanthropy) %>%
    filter(host_species %in% hosts_with_complete_range) %>% # Ensure range is also complete
    mutate(synanthropy_status = factor(synanthropy_status, levels = c("Not Synanthropic", "Occasionally Synanthropic", "Totally Synanthropic"))) %>% # Convert to factor for nodefactor
    arrange(match(host_species, matrix_hosts))
  
  # --- d) Subset Matrix for Models ---
  
  # Matrix for the "full" model
  matrix_full <- incidence_matrix[rownames(incidence_matrix) %in% covariates_full$host_species, ]
  
  # Matrix for the "synanthropy" subset model
  matrix_synanthropy <- incidence_matrix[rownames(incidence_matrix) %in% covariates_synanthropy$host_species, ]
  
  # --- e) Return objects ---
  return(list(
    matrix_full = matrix_full,
    covariates_full = covariates_full,
    matrix_synanthropy = matrix_synanthropy,
    covariates_synanthropy = covariates_synanthropy,
    family = family_name
  ))
}


# --- 5. Prepare Data for Both Families ---

arena_data <- prepare_ergm_data("Arenaviridae", pathogen_data, host_sampling_effort, host_traits)
hanta_data <- prepare_ergm_data("Hantaviridae", pathogen_data, host_sampling_effort, host_traits)


# --- 6. Create Network Objects ---

# --- Arenaviridae ---
net_arena_full_base <- network(arena_data$matrix_full,
                               matrix.type = "adjacency",
                               bipartite = TRUE,
                               ignore.eval = FALSE,
                               names.eval = "interaction",
                               directed = FALSE)

aligned_covariates_full <- arena_data$covariates_full %>%
  arrange(match(host_species, rownames(arena_data$matrix_full)))

net_arena_full <- set.vertex.attribute(net_arena_full_base,
                                       attrname = colnames(aligned_covariates_full)[-1], # Exclude host_species column
                                       value = aligned_covariates_full[, -1],           # Pass the data frame of attributes
                                       v = 1:nrow(arena_data$matrix_full))

net_arena_synanthropy_base <- network(arena_data$matrix_synanthropy,
                                      matrix.type = "adjacency",
                                      bipartite = TRUE,
                                      ignore.eval = FALSE,
                                      names.eval = "interaction",
                                      directed = FALSE)

aligned_covariates_synanthropy <- arena_data$covariates_synanthropy %>%
  arrange(match(host_species, rownames(arena_data$matrix_synanthropy)))

net_arena_synanthropy <- set.vertex.attribute(net_arena_synanthropy_base,
                                              attrname = colnames(aligned_covariates_synanthropy)[-1],
                                              value = aligned_covariates_synanthropy[, -1],
                                              v = 1:nrow(arena_data$matrix_synanthropy))

# --- Hantaviridae ---
net_hanta_full_base <- network(hanta_data$matrix_full,
                               matrix.type = "adjacency",
                               bipartite = TRUE,
                               ignore.eval = FALSE,
                               names.eval = "interaction",
                               directed = FALSE)

aligned_hanta_covariates_full <- hanta_data$covariates_full %>%
  arrange(match(host_species, rownames(hanta_data$matrix_full)))

net_hanta_full <- set.vertex.attribute(net_hanta_full_base,
                                       attrname = colnames(aligned_hanta_covariates_full)[-1],
                                       value = aligned_hanta_covariates_full[, -1],
                                       v = 1:nrow(hanta_data$matrix_full))

net_hanta_synanthropy_base <- network(hanta_data$matrix_synanthropy,
                                      matrix.type = "adjacency",
                                      bipartite = TRUE,
                                      ignore.eval = FALSE,
                                      names.eval = "interaction",
                                      directed = FALSE)

aligned_hanta_covariates_synanthropy <- hanta_data$covariates_synanthropy %>%
  arrange(match(host_species, rownames(hanta_data$matrix_synanthropy)))

net_hanta_synanthropy <- set.vertex.attribute(net_hanta_synanthropy_base,
                                              attrname = colnames(aligned_hanta_covariates_synanthropy)[-1],
                                              value = aligned_hanta_covariates_synanthropy[, -1],
                                              v = 1:nrow(hanta_data$matrix_synanthropy))

# --- 7. Define and Fit ERGM Models ---

# --- Arenaviridae Models ---

# Model 1: Full dataset (Effort + Range Size + edge density)
formula_arena_full<- net_arena_full ~ b1cov("log_sampling_effort") + b1cov("log_range_size") + edges

ergm_arena_full <- ergm(formula_arena_full, control = control.ergm(seed = 123))

# Model 2: Synanthropy subset (Effort + Range Size + Synanthropy)
formula_arena_synanthropy <- net_arena_synanthropy ~ b1cov("log_sampling_effort") + b1cov("log_range_size") + b1factor("synanthropy_status") + edges

ergm_arena_synanthropy <- ergm(formula_arena_synanthropy, control = control.ergm(seed = 123))

# --- Hantaviridae Models ---

# Model 3: Full dataset (Effort + Range Size)
formula_hanta_full <- net_hanta_full ~ b1cov("log_sampling_effort") + b1cov("log_range_size") + edges

ergm_hanta_full <- ergm(formula_hanta_full, control = control.ergm(seed = 123))

# Model 4: Synanthropy subset (Effort + Range Size + Synanthropy)
formula_hanta_synanthropy <- net_hanta_synanthropy ~ b1cov("log_sampling_effort") + b1cov("log_range_size") + b1factor("synanthropy_status") + edges

ergm_hanta_synanthropy <- ergm(formula_hanta_synanthropy, control = control.ergm(seed = 123))


# --- 8. Summarize Results ---

summary_arena_full <- summary(ergm_arena_full)
summary_arena_synanthropy <- summary(ergm_arena_synanthropy)
summary_hanta_full <- summary(ergm_hanta_full)
summary_hanta_synanthropy <- summary(ergm_hanta_synanthropy)

# Save summaries
all_summaries <- list(arena_full = summary_arena_full,
                      arena_synanthropy = summary_arena_synanthropy,
                      hanta_full = summary_hanta_full,
                      hanta_synanthropy = summary_hanta_synanthropy)

write_rds(all_summaries, here(output_dir_ergm, "ergm_summaries.rds"))


# --- 9. Goodness-of-Fit Diagnostics ---
# Example for one model:
gof_arena_full <- gof(ergm_arena_full)
gof_arena_synanthropy <- gof(ergm_arena_synanthropy)
gof_hanta_full <- gof(ergm_hanta_full)
gof_hanta_synanthropy <- gof(ergm_hanta_synanthropy)
#plot(gof_arena_full)

# --- 10. Prepare Data Frame for Dyadic GLMM ---

prepare_dyadic_df <- function(matrix_data, covariates_data) {
  
  host_names <- rownames(matrix_data)
  pathogen_names <- colnames(matrix_data)
  
  # Create all possible dyads
  dyad_df <- expand_grid(host_species = host_names, 
                         pathogen_species = pathogen_names)
  
  # Add the observed interaction status (0 or 1)
  interaction_long <- as.data.frame.table(matrix_data, responseName = "interaction")
  colnames(interaction_long) <- c("host_species", "pathogen_species", "interaction")
  
  dyad_df <- dyad_df %>%
    left_join(interaction_long, by = c("host_species", "pathogen_species"))
  
  # Add host covariates
  # Ensure covariate data only contains columns to join + actual covariates
  covs_to_join <- covariates_data %>% 
    select(host_species, log_sampling_effort, log_range_size, any_of(contains("synanthropy_status")))
  
  dyad_df <- dyad_df %>%
    left_join(covs_to_join, by = "host_species")
  
  # Ensure factors are correctly set
  dyad_df <- dyad_df %>%
    mutate(host_species = factor(host_species),
           pathogen_species = factor(pathogen_species))
  
  # Re-level synanthropy if it exists, setting "Not Synanthropic" as reference
  if("synanthropy_status" %in% names(dyad_df)){
    dyad_df <- dyad_df %>% 
      mutate(synanthropy_status = factor(synanthropy_status, levels = c("Not Synanthropic", "Occasionally Synanthropic", "Totally Synanthropic")))
  }
  
  return(dyad_df)
}

dyadic_df_arena_full <- prepare_dyadic_df(arena_data$matrix_full, arena_data$covariates_full)
dyadic_df_arena_synanthropy <- prepare_dyadic_df(arena_data$matrix_synanthropy, arena_data$covariates_synanthropy)
dyadic_df_hanta_full <- prepare_dyadic_df(hanta_data$matrix_full, hanta_data$covariates_full)
dyadic_df_hanta_synanthropy <- prepare_dyadic_df(hanta_data$matrix_synanthropy, hanta_data$covariates_synanthropy)

# --- 11. Define and Fit brms Models ---
# Arenaviridae
formula_glmm_arena_full <- bf(interaction ~ scale(log_sampling_effort) + scale(log_range_size) + (1 | host_species) + (1 | pathogen_species))

brm_arena_full <- brm(formula = formula_glmm_arena_full,
                      data = dyadic_df_arena_full,
                      family = bernoulli(link = "logit"),
                      cores = 4, chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.95),
                      seed = 123)

formula_glmm_arena_synanthropy <- bf(interaction ~ scale(log_sampling_effort) + scale(log_range_size) + synanthropy_status + (1 | host_species) + (1 | pathogen_species))

brm_arena_synanthropy <- brm(formula = formula_glmm_arena_synanthropy,
                             data = dyadic_df_arena_synanthropy,
                             family = bernoulli(link = "logit"),
                             cores = 4, chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.95),
                             seed = 123)

# Hantaviridae
formula_glmm_hanta_full <- bf(interaction ~ scale(log_sampling_effort) + scale(log_range_size) + (1 | host_species) + (1 | pathogen_species))

brm_hanta_full <- brm(formula = formula_glmm_hanta_full,
                      data = dyadic_df_hanta_full,
                      family = bernoulli(link = "logit"),
                      cores = 4, chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.95),
                      seed = 123)

formula_glmm_hanta_synanthropy <- bf(interaction ~ scale(log_sampling_effort) + scale(log_range_size) + synanthropy_status + (1 | host_species) + (1 | pathogen_species))

brm_hanta_synanthropy <- brm(formula = formula_glmm_hanta_synanthropy,
                             data = dyadic_df_hanta_synanthropy,
                             family = bernoulli(link = "logit"),
                             cores = 4, chains = 4, iter = 2000, warmup = 1000, control = list(adapt_delta = 0.95),
                             seed = 123)

summary_glmm_arena_full <- summary(brm_arena_full)
summary_glmm_arena_synanthropy <- summary(brm_arena_synanthropy)
summary_glmm_hanta_full <- summary(brm_hanta_full)
summary_glmm_hanta_synanthropy <- summary(brm_hanta_synanthropy)

all_summaries_glmm <- list(arena_full_glmm = summary_glmm_arena_full,
                           arena_synanthropy_glmm = summary_glmm_arena_synanthropy,
                           hanta_full_glmm = summary_glmm_hanta_full,
                           hanta_synanthropy_glmm = summary_glmm_hanta_synanthropy)

write_rds(all_summaries_glmm, here(output_dir_ergm, "glmm_summaries.rds"))

# Visualise
draws_arena <- brm_arena_synanthropy %>%
  gather_draws(`b_synanthropy_status.*`, regex = TRUE) %>%
  mutate(virus_family = "Arenaviridae")

draws_hanta <- brm_hanta_synanthropy %>%
  gather_draws(`b_synanthropy_status.*`, regex = TRUE) %>%
  mutate(virus_family = "Hantaviridae")

all_draws <- bind_rows(draws_arena, draws_hanta) %>%
  mutate(OR = exp(.value)) %>%
  mutate(Comparison = str_remove(.variable, "b_synanthropy_status") %>%
           str_replace_all("([A-Z])", " \\1") %>%
           str_trim())

p_odds_ratios_halfeye <- ggplot(all_draws, aes(x = OR, y = Comparison, fill = virus_family)) +
  stat_halfeye(adjust = 0.5, .width = c(0.66, 0.95), point_interval = "median_qi") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32)) + 
  scale_fill_manual(values = c("Arenaviridae" = "#E69F00", "Hantaviridae" = "#56B4E9"),
                    name = "Virus Family") +
  facet_wrap(~ virus_family, ncol = 1) + 
  labs(y = element_blank(),
       x = "OR with 66% & 95% CrI (log scale)") +
  theme_minimal() +
  theme(legend.position = "none")

# Save the plot
ggsave(here("output", "analysis_1", "synanthropy_odds_ratio_plot.png"), 
       plot = p_odds_ratios_halfeye, width = 7, height = 5, dpi = 300)
