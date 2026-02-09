suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript prune_upham_trees.R <upham_tree> <hanta_matrix_csv> <arena_matrix_csv> <out_dir> <prefix>")
}

upham_tree_path <- args[1]
hanta_matrix_path <- args[2]
arena_matrix_path <- args[3]
out_dir <- args[4]
prefix <- args[5]

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

hanta <- read.csv(hanta_matrix_path, row.names = 1, check.names = FALSE)
arena <- read.csv(arena_matrix_path, row.names = 1, check.names = FALSE)

hanta_hosts <- trimws(rownames(hanta))
arena_hosts <- trimws(rownames(arena))

all_hosts <- sort(unique(c(hanta_hosts, arena_hosts)))

upham <- read.nexus(upham_tree_path)
clean_tip <- function(label) {
  parts <- unlist(strsplit(label, "_"))
  parts <- parts[parts != ""]
  if (length(parts) >= 2) {
    return(paste(parts[1:2], collapse = " "))
  }
  return(label)
}
upham$tip.label <- vapply(upham$tip.label, clean_tip, character(1))
upham$tip.label <- trimws(upham$tip.label)

dup_labels <- duplicated(upham$tip.label)
if (any(dup_labels)) {
  upham <- drop.tip(upham, upham$tip.label[dup_labels])
}

in_all <- intersect(all_hosts, upham$tip.label)
if (length(in_all) == 0) {
  stop("No overlap between Upham tree tips and host lists.")
}

in_arena <- intersect(arena_hosts, upham$tip.label)

all_tree <- keep.tip(upham, in_all)
arena_tree <- keep.tip(upham, in_arena)

all_path <- file.path(out_dir, paste0(prefix, "_upham_hosts_hanta_arena.nwk"))
arena_path <- file.path(out_dir, paste0(prefix, "_upham_hosts_arena.nwk"))

write.tree(all_tree, all_path)
write.tree(arena_tree, arena_path)

summary_path <- file.path(out_dir, paste0(prefix, "_upham_prune_summary.txt"))
summary_lines <- c(
  paste("Total tree tips:", length(upham$tip.label)),
  paste("Hanta hosts:", length(hanta_hosts)),
  paste("Arena hosts:", length(arena_hosts)),
  paste("Unique hosts (Hanta+Arena):", length(all_hosts)),
  paste("Hosts kept (Hanta+Arena):", length(in_all)),
  paste("Hosts kept (Arena only):", length(in_arena))
)
writeLines(summary_lines, con = summary_path)
