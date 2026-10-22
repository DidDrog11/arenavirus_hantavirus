# Project ArHa: Figures for the data release paper
# 09_figures.R
# Presentation only: reads the released tables (data/data_release/parquet/) and draws the
# main figures into output/figures/. Supplementary Figure S1 (PRISMA) is 99_prisma_fig.R.
#   Figure 1: sampling localities, coloured by the resolution of their coordinates
#   Figure 2: what the database holds (studies by year, host records by order, assay records)
#   Figure 3: the four relational tables and how they link

if (!require("pacman")) install.packages("pacman")
pacman::p_load(c("here", "tidyverse", "arrow", "terra", "tidyterra", "geodata", "cowplot",
                 "DiagrammeR", "DiagrammeRsvg", "rsvg"), character.only = TRUE)

parquet_dir <- here("data", "data_release", "parquet")
sampling_events        <- read_parquet(file.path(parquet_dir, "sampling_events.parquet"))
host_occurrences       <- read_parquet(file.path(parquet_dir, "host_occurrences.parquet"))
pathogen_mof           <- read_parquet(file.path(parquet_dir, "pathogen_mof.parquet"))
resource_relationships <- read_parquet(file.path(parquet_dir, "resource_relationships.parquet"))

figure_dir <- here("output", "figures")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Colours: the first three categorical slots of the reference palette (they stay distinct
# for colour-blind readers when all pairs are compared, as on a map), a single hue for
# bars, and grey ink for text and the basemap
col_series <- c("#2a78d6", "#eb6834", "#1baf7a")
col_bar    <- "#2a78d6"
col_text   <- "#52514e"
col_land   <- "#e4e3df"

theme_arha <- theme_minimal(base_size = 9) +
  theme(text = element_text(colour = col_text),
        axis.text = element_text(colour = col_text),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 9, colour = "black"),
        plot.title.position = "plot")

fmt <- function(n) format(n, big.mark = ",", trim = TRUE)

# --- Figure 3: relational schema ----------------------------------------------
# One box per table with its Darwin Core role, key and row count; arrows run from the
# referenced table to the table holding the foreign key
table_box <- function(name, role, key, n) {
  paste0('<<table border="0" cellborder="0" cellpadding="2">',
         '<tr><td><b>', name, '</b></td></tr>',
         '<tr><td><font point-size="9">', role, '</font></td></tr>',
         '<tr><td><font point-size="9">key: ', key, '</font></td></tr>',
         '<tr><td><font point-size="9">', fmt(n), ' records</font></td></tr></table>>')
}

schema_dot <- paste0('
  digraph schema {
    graph [rankdir = LR, nodesep = 0.5, ranksep = 0.9]
    node [shape = box, style = "rounded", fontname = Arial, fontsize = 11, margin = "0.2,0.1"]
    edge [fontname = Arial, fontsize = 9, color = "#52514e", fontcolor = "#52514e"]

    events [label = ', table_box("sampling_events", "Event core", "eventID", nrow(sampling_events)), ']
    hosts [label = ', table_box("host_occurrences", "Occurrence extension", "occurrenceID", nrow(host_occurrences)), ']
    assays [label = ', table_box("pathogen_mof", "Extended MeasurementOrFact", "measurementID", nrow(pathogen_mof)), ']
    sequences [label = ', table_box("resource_relationships", "ResourceRelationship extension", "resourceRelationshipID", nrow(resource_relationships)), ']

    events -> hosts [label = "eventID"]
    hosts -> assays [label = "occurrenceID"]
    events -> assays [label = "eventID", style = dashed]
    hosts -> sequences [label = "resourceID (host)"]
    assays -> sequences [label = "resourceID (assay)"]
  }
')

grViz(schema_dot) |>
  export_svg() |>
  charToRaw() |>
  rsvg_png(file.path(figure_dir, "figure_3_schema.png"), width = 2400)

# --- Figure 1: sampling localities --------------------------------------------
# One point per distinct coordinate (locality), grouped into three levels of precision; a
# coordinate used at more than one resolution takes its most precise. Points are drawn
# coarsest first, so the most precise localities sit on top.
resolution_levels <- c("Site or settlement", "Study area or administrative unit", "Country or not specified")

localities <- host_occurrences |>
  mutate(resolution = case_when(
    coordinate_resolution_processed %in% c("site", "village", "town", "city") ~ 1L,
    coordinate_resolution_processed %in% c("study_area", "adm3", "adm2", "adm1") ~ 2L,
    TRUE ~ 3L)) |>
  group_by(decimalLatitude, decimalLongitude) |>
  summarise(resolution = factor(resolution_levels[min(resolution)], levels = resolution_levels), .groups = "drop")

resolution_counts <- count(localities, resolution)
resolution_labels <- setNames(paste0(resolution_counts$resolution, " (", fmt(resolution_counts$n), ")"),
                              resolution_counts$resolution)

# Equal Earth projection; the coarse GADM world outline is enough for a basemap
equal_earth <- "+proj=eqearth"
world <- world(path = here("data")) |> crop(ext(-180, 180, -60, 90)) |> project(equal_earth)
locality_points <- vect(localities, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326") |>
  project(equal_earth)
localities <- localities |> bind_cols(as_tibble(crds(locality_points))) |>
  arrange(desc(resolution))

map_plot <- ggplot() +
  geom_spatvector(data = world, fill = col_land, colour = "white", linewidth = 0.1) +
  geom_point(data = localities, aes(x, y, fill = resolution), shape = 21, size = 1.1, stroke = 0.15, colour = "white") +
  scale_fill_manual(values = setNames(col_series, resolution_levels), labels = resolution_labels,
                    name = "Resolution of the coordinate (localities)") +
  coord_sf(crs = equal_earth, datum = NA) +
  guides(fill = guide_legend(override.aes = list(size = 2.5), nrow = 1)) +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", legend.title.position = "top",
        legend.text = element_text(colour = col_text), legend.title = element_text(colour = col_text))

ggsave(file.path(figure_dir, "figure_1_map.png"), map_plot, width = 180, height = 105, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure_1_map.pdf"), map_plot, width = 180, height = 105, units = "mm", bg = "white")

# --- Figure 2: contents ------------------------------------------------------
# (a) Studies by the year sampling began. Dates inferred from the publication are left out.
first_year <- host_occurrences |>
  filter(temporal_resolution != "publication_derived") |>
  mutate(start_year = as.integer(substr(eventDate, 1, 4))) |>
  group_by(eventID) |>
  summarise(start_year = min(start_year), .groups = "drop")
n_studies_dated <- nrow(first_year)

panel_years <- ggplot(count(first_year, start_year), aes(start_year, n)) +
  geom_col(fill = col_bar, width = 0.8) +
  scale_x_continuous(breaks = seq(1960, 2020, by = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = paste0("a  Studies by the year sampling began (n = ", fmt(n_studies_dated), ")"), x = NULL, y = "Studies") +
  theme_arha + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "#ecebe7"))

# (b) Host records by order, with the smallest orders folded into "Other"
orders <- host_occurrences |>
  mutate(order = coalesce(order, "Not resolved to order")) |>
  count(order) |>
  mutate(order = if_else(rank(-n, ties.method = "first") <= 6, order, "Other orders")) |>
  group_by(order) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(order = fct_reorder(order, n) |> fct_relevel("Other orders"))

panel_orders <- ggplot(orders, aes(n, order)) +
  geom_col(fill = col_bar, width = 0.7) +
  geom_text(aes(x = n + 0.015 * max(n), label = fmt(n)), hjust = 0, size = 2.6, colour = col_text) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)), breaks = c(0, 25000, 50000), labels = fmt) +
  labs(title = "b  Host records by order", x = "Host records", y = NULL) +
  theme_arha

# (c) Assay records by assay type, for the two target families
assays <- pathogen_mof |>
  filter(family_pathogen %in% c("Arenaviridae", "Hantaviridae")) |>
  mutate(assay = case_when(measurementMethod %in% c("PCR", "Serology", "Antigen") ~ measurementMethod,
                           TRUE ~ "Other or not reported")) |>
  count(family_pathogen, assay) |>
  mutate(assay = factor(assay, levels = rev(c("PCR", "Serology", "Antigen", "Other or not reported"))))
n_other_families <- sum(!pathogen_mof$family_pathogen %in% c("Arenaviridae", "Hantaviridae"))

panel_assays <- ggplot(assays, aes(n, assay)) +
  geom_col(fill = col_bar, width = 0.7) +
  geom_text(aes(x = n + 0.015 * max(n), label = fmt(n)), hjust = 0, size = 2.6, colour = col_text) +
  facet_wrap(~ family_pathogen, ncol = 1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2)), labels = fmt) +
  labs(title = "c  Assay records by virus family and assay", x = "Assay records", y = NULL) +
  theme_arha + theme(strip.text = element_text(hjust = 0, face = "italic", colour = col_text))

contents_plot <- plot_grid(panel_years, plot_grid(panel_orders, panel_assays, nrow = 1, align = "h", axis = "b"),
                           ncol = 1, rel_heights = c(0.8, 1))

ggsave(file.path(figure_dir, "figure_2_contents.png"), contents_plot, width = 180, height = 150, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "figure_2_contents.pdf"), contents_plot, width = 180, height = 150, units = "mm", bg = "white")

cat("Figures written to ", figure_dir, "\n",
    "Figure 2a: ", n_studies_dated, " of ", nrow(sampling_events), " studies have a sampling date\n",
    "Figure 2c: ", n_other_families, " assay records for other or unrecorded virus families left out\n", sep = "")
