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
