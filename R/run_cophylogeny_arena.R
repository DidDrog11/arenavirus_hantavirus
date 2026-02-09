suppressPackageStartupMessages({
  library(ape)
  library(paco)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_cophylogeny_arena.R <assoc_matrix_csv> <host_tree_nwk> <virus_tree_mcc> <ncbi_map_csv> <out_dir>")
}

assoc_path <- args[1]
host_tree_path <- args[2]
virus_tree_path <- args[3]
map_path <- args[4]
out_dir <- args[5]

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

assoc <- read.csv(assoc_path, row.names = 1, check.names = FALSE)
assoc <- as.matrix(assoc)
assoc[assoc > 0] <- 1
assoc <- assoc[rowSums(assoc) > 0, colSums(assoc) > 0, drop = FALSE]

htree <- read.tree(host_tree_path)
htree$tip.label <- gsub("_", " ", htree$tip.label)
htree$tip.label <- trimws(htree$tip.label)

common_hosts <- intersect(rownames(assoc), htree$tip.label)
if (length(common_hosts) == 0) {
  stop("No overlap between host tree tips and association matrix hosts.")
}
assoc <- assoc[common_hosts, , drop = FALSE]
htree <- keep.tip(htree, common_hosts)

ptree <- read.nexus(virus_tree_path)
ptree$tip.label <- gsub("'", "", ptree$tip.label)

ptree$tip.label <- sapply(ptree$tip.label, function(x) {
  if (grepl("\\|", x)) {
    parts <- trimws(strsplit(x, "\\|", fixed = FALSE)[[1]])
    hit <- parts[grepl("arenavirus", tolower(parts))]
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  x
})

if (any(grepl("arenavirus", tolower(ptree$tip.label)))) {
  ptree$tip.label <- gsub("_", " ", ptree$tip.label)
} else {
  ncbi <- read.csv(map_path)
  map <- setNames(ncbi$Species, ncbi$Accession)
  ptree$tip.label <- sapply(ptree$tip.label, function(x) if (x %in% names(map)) map[[x]] else x)
  ptree$tip.label <- gsub("_", " ", ptree$tip.label)
}
ptree$tip.label <- trimws(ptree$tip.label)

common_viruses <- intersect(colnames(assoc), ptree$tip.label)
if (length(common_viruses) == 0) {
  stop("No overlap between virus tree tips and association matrix pathogens.")
}
assoc <- assoc[, common_viruses, drop = FALSE]
ptree <- keep.tip(ptree, common_viruses)

assoc <- assoc[rowSums(assoc) > 0, colSums(assoc) > 0, drop = FALSE]
htree <- keep.tip(htree, rownames(assoc))
ptree <- keep.tip(ptree, colnames(assoc))

write.csv(assoc, file = file.path(out_dir, "arena_assoc_matrix_pruned.csv"))
write.tree(htree, file = file.path(out_dir, "arena_host_tree_pruned.nwk"))
write.tree(ptree, file = file.path(out_dir, "arena_virus_tree_pruned.nwk"))

hdist <- cophenetic.phylo(htree)
pdist <- cophenetic.phylo(ptree)

hdist <- hdist[rownames(assoc), rownames(assoc), drop = FALSE]
pdist <- pdist[colnames(assoc), colnames(assoc), drop = FALSE]

set.seed(1)
pfit <- parafit(host.D = hdist, para.D = pdist, HP = assoc,
                correction = "cailliez", nperm = 999, test.links = TRUE)

D <- prepare_paco_data(H = hdist, P = pdist, HP = assoc)
D <- add_pcoord(D, correction = "cailliez")
set.seed(1)
pac <- PACo(D, nperm = 999, seed = 1, method = "r0", symmetric = FALSE)
pac_links <- paco_links(pac)
res.l <- residuals_paco(pac_links$proc)

imat <- which(assoc > 0, arr.ind = TRUE)
assoc_df <- data.frame(
  host = rownames(assoc)[imat[, 1]],
  virus = colnames(assoc)[imat[, 2]],
  link = 1
)

min_alpha <- 0.2
max_alpha <- 0.8
min_lwd <- 0.2
max_lwd <- 0.8

inv_res <- 1 / (1 + res.l)
inv_res_scaled <- (inv_res - min(inv_res)) / (max(inv_res) - min(inv_res))
alpha_vals <- min_alpha + inv_res_scaled * (max_alpha - min_alpha)
lwd_vals <- min_lwd + inv_res_scaled * (max_lwd - min_lwd)

jpeg(file.path(out_dir, "cophylo_arena.jpeg"), width = 2500, height = 1600, res = 300)
par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
co <- cophylo(htree, ptree, assoc = assoc_df, rotate = TRUE)
plot(co, link.type = "curved", link.lty = "solid",
     fsize = c(0.4, 0.5),
     link.lwd = lwd_vals,
     link.col = grDevices::rgb(0, 0, 0, alpha = alpha_vals))
dev.off()

summary_path <- file.path(out_dir, "cophylo_summary.txt")
summary_lines <- c(
  paste("Hosts:", length(htree$tip.label)),
  paste("Viruses:", length(ptree$tip.label)),
  paste("Links:", nrow(assoc_df)),
  paste("PACo_gof:", paste(capture.output(print(pac_links$gof)), collapse = " "))
)
writeLines(summary_lines, con = summary_path)

saveRDS(list(parafit = pfit, paco = pac, paco_links = pac_links, residuals = res.l),
        file = file.path(out_dir, "cophylo_results.rds"))

pair_names <- names(res.l)
parts <- strsplit(pair_names, "-", fixed = TRUE)
host <- vapply(parts, `[`, character(1), 1)
parasite <- vapply(parts, `[`, character(1), 2)
res_df <- data.frame(host = host, virus = parasite, residual = as.numeric(res.l))
write.csv(res_df, file = file.path(out_dir, "cophylo_pairs_residuals.csv"), row.names = FALSE)
