#Cophylogeny analysis Hanta - Project ArHa, March, 2025

# Co-phylogenetic analysis of Hantavirus and hosts
# Clean workspace and load libraries
rm(list = ls())
graphics.off()
gc()

# Load necessary libraries
library(tidyverse)
library(janitor)
library(ape)
library(paco)
library(ggtree)
library(reshape2)
library(treedataverse)
library(taxonomizr)
library(camcorder)
library(extrafont)
library(extrafontdb)
library(seriation)
library(plotrix)
library(phytools)
library(MoMAColors)
library(ragg)
library(ggplot2)
library(hrbrthemes)
library(tidyr)
library(dplyr)

#---------------------------
# FUNCTIONS
#---------------------------

# Function to build and binarize the host-virus association matrix
build_assoc_matrix <- function(data, host_col = "host", virus_col = "organism_name") {
  assoc <- with(data, table(get(host_col), get(virus_col)))
  assoc_bin <- ifelse(assoc > 0, 1, 0)
  return(assoc_bin)
}

# Function to extract species name from a tip label for viruses (between | symbols)
extract_and_clean <- function(vector) {
  vector %>% 
    str_extract("(?<=\\|)[^|]+(?=\\|)") %>%  # Extract string between first and second "|"
    str_replace_all("_", " ")                # Replace "_" with " "
}

# Function to extract clean host species names from tip labels (before the second underscore)
split_before_second_underscore <- function(vector) {
  vector %>%
    str_extract("^[^_]+_[^_]+") %>%            # Extract everything before the second underscore
    str_replace_all("_", " ")                  # Replace "_" with " "
}

# Function to filter association matrix based on tree tip labels
filter_assoc_matrix <- function(am, host_tips, virus_tips) {
  am_filtered <- am[rownames(am) %in% host_tips, ]
  am_filtered <- am_filtered[, colnames(am_filtered) %in% virus_tips]
  return(am_filtered)
}

#---------------------------
# DATA LOADING & PREPROCESSING
#---------------------------

#Read ArHa data to build host-virus association data, filtering the sequences without info
arha <- read_csv('data/ncbi_virus/arha_enriched.csv')
#Separate into arena and Hanta
hanta <- arha %>% 
  filter(virus_type == 'hantavirus',
         species != 'Homo sapiens')

arena <- arha %>% 
  filter(virus_type == 'arenavirus')

#Create host-virus pairs

hv_hanta_matrix <- build_assoc_matrix(hanta, 
                                      host_col = 'species', 
                                      virus_col = 'virus_clean')

hnames <- rownames(hv_hanta_matrix)
pnames <- colnames(hv_hanta_matrix)

#Phylogenetic trees, for Host we'll use the Upham tree
tree <- read.nexus('~/Documents/Hanta/Data/Trees/Upham/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

#Fix taxa name in the Upham tree
taxa=read.csv('~/Documents/Hanta/Data/Trees/Upham/taxonomy_mamPhy_5911species.csv',header=T)
taxa$tip=taxa$Species_Name

#Remove underscores
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

#Filter data not available
species <- trimws(tree$tip.label, 'both') #Remove leading and trailing whitespaces

host_species_upham <- species[species %in% hnames]

tree=keep.tip(tree,host_species_upham)
plotTree(tree)

#save tree
ape::write.tree(tree, 'data/ncbi_virus/uphamTree_hantaHosts.nwk')
#Save taxa list to get genomes (continuing in Python)
write(pnames, file='data/ncbi_virus/associatedViralSpecies_cophylo.txt', sep = '\t')

#Load BEAST MCC tree for Hanta
hanta_tree <- read.beast('data/ncbi_virus/Hanta_S-hanta_S_aligned.mcc.tree') %>% 
  as.phylo()

hanta_tree$tip.label <- sapply(hanta_tree$tip.label, function(x) str_split(x, "\\.")[[1]][1])


#Set the tipnames to species using full NCBI data
ncbi_hanta <- read_csv('data/ncbi_virus/hantaviridae_ncbi_03062025.csv')

mapping <- setNames(ncbi_hanta$Species, ncbi_hanta$Accession)
hanta_tree$tip.label <- sapply(hanta_tree$tip.label, function(x) {
  if (x %in% names(mapping)) {
    mapping[x]
  } else {
    x
  }
})

plotTree(hanta_tree)

longdata <- melt(hv_hanta_matrix)

hlist <- tree$tip.label
plist <- hanta_tree$tip.label

# Reorder association matrix using seriation
o <- seriate(hv_hanta_matrix, method = "BEA_TSP")
longdata$host <- factor(longdata$Var1, levels = names(unlist(o[[1]][])))
longdata$organism_name <- factor(longdata$Var2, levels = names(unlist(o[[2]][])))

# Filter longdata to include only rows with tree tip labels and cast to matrix
ld <- longdata %>% 
  filter(host %in% hlist, organism_name %in% plist) %>% 
  droplevels()
amtab <- acast(ld, host ~ organism_name)

amtab_filtered <- filter_assoc_matrix(amtab, hnames, pnames)
htree <- keep.tip(tree, rownames(amtab_filtered))
ptree <- keep.tip(hanta_tree, colnames(amtab_filtered))

cat("Dimensions of filtered association matrix:", dim(amtab_filtered), "\n")
cat("Number of host tips:", length(htree$tip.label), "\n")
cat("Number of virus tips:", length(ptree$tip.label), "\n")


#---------------------------
# PHYLOGENETIC ANALYSES
#---------------------------

# Compute cophenetic distances
hdist <- cophenetic.phylo(htree)
pdist <- cophenetic.phylo(ptree)


# Parafit test
set.seed(1)
pfit <- parafit(host.D = hdist, para.D = pdist, HP = amtab_filtered,
                correction = "cailliez", nperm = 999, test.links = TRUE)
print(pfit)

# Prepare PACo data and perform PACo analysis
D <- prepare_paco_data(H = hdist, P = pdist, HP = amtab_filtered)
D <- add_pcoord(D, correction = "cailliez")
set.seed(1)
pac <- PACo(D, nperm = 999, seed = 1, method = "r0", symmetric = FALSE)
pac_links <- paco_links(pac)
res.l <- residuals_paco(pac_links$proc)
cat("PACo goodness-of-fit:\n")
print(pac_links$gof)

# Calculate weights for links based on residuals
wei <- ((res.l^-2) / 2)

# Create interaction matrix for cophyloplot format
imat.l <- melt(amtab)
names(imat.l) <- c("host", "virus", "link")
imat.l <- imat.l[which(imat.l$link > 0), ]

# Plot cophyloplot
par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
cophyloplot(htree, ptree, assoc = imat.l, show.tip.label = TRUE,
            use.edge.length = FALSE, lwd = 1 / wei, space = 200, gap = 5, length.line = 20)

# Generate color palette for virus tips
pal <- MoMAColors::moma.colors("Ohchi", length(unique(ptree$tip.label)))

# Create cophylogenetic object and plot using phytools
co.phylo.l <- cophylo(htree, ptree, assoc = imat.l, rotate = TRUE)
link <- co.phylo.l$assoc

ragg::agg_jpeg(filename = "output/cophylo_iter1.jpeg",
               height = 25 * 10/16, width = 25, units = "in", res = 750)

plot(co.phylo.l, link.type = "curved", link.lty = "solid",
     fsize = c(0.7, 0.5), link.lwd = ifelse(res.l < 50, 1, 0.1),
     link.col = pal)

dev.off()
