#!/usr/bin/env Rscript

# Author: Henry Ward, Xiang Zhang
# Clustering to get layer 0 and layer 1 on fullFF 

# Loads all packages in a way that allows exporting to child environments
source("/project/csbio/henry/Documents/libraries/lib_util.R")
.libPaths("/project/csbio/henry/Documents/libraries/R")
packages <- c("ggplot2", "ggthemes")
for (p in packages) {
  library(p, character.only = TRUE)
}

######
# DATA PREP
######

# Loads qGI scores and sets output folder
setwd("/project/GIN/clustering/")
source(file.path("../src", "clusterv2.R"))
scores <- read.csv("/project/GIN/input/qGI_20211111_fullFF.txt", sep = "\t", row.names = 1)
output_folder <- file.path("output_fullFF", "clusters_v2_expanded", "2021_11_11_qGI_exp_redundant_fullFF_snr")
if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }

# Sets parameters that determine the number of clusters for each run
n_layers <- 1
metric_thresholds <- c(0.01, 0.01)
#metric_thresholds <- c(0.01, 0.01, 0.01, 0.01)  # One more than the number of layers if we don't load existing clusters
metric <- "snr"
sort_metric <- "size"
sort_decreasing <- FALSE
member_threshold <- 3
#filter_file <- file.path(output_folder, "orig_layer1_clusters_to_filter.txt")
filter_file <- NULL

# Sets whether or not to load existing clusters by giving a path to that file or generate new ones,
# by setting existing_clusters to NULL
existing_clusters <- NULL
#existing_clusters <- file.path(output_folder, "layer0", "best_clusters.txt")

# Replaces NAs with row means
row_means <- rowMeans(scores, na.rm = TRUE)
for (i in 1:nrow(scores)) {
  na_ind <- which(is.na(scores[i,]))
  scores[i, na_ind] <- row_means[i]
}

######
# FIRST LAYER PREP
######

# Either loads in existing clusters or clusters qGI scores directly
first_clust <- NULL
layer_output_folder <- file.path(output_folder, "layer0")
if (!dir.exists(layer_output_folder)) { dir.create(layer_output_folder) }
if (is.null(existing_clusters)) {
  gi_cluster(scores, layer_output_folder,
             metric_threshold = metric_thresholds[1],
             metric = metric,
             sort_metric = sort_metric,
             sort_decreasing = sort_decreasing,
             pick_num = 1,
             member_threshold = member_threshold,
             verbose = TRUE)
  first_clust <- readLines(file.path(layer_output_folder, "best_clusters.txt"))
} else {
  first_clust <- readLines(existing_clusters)
}

# Manually filters the first layer of clusters if file specified
if (!is.null(filter_file)) {
  to_remove <- as.numeric(readLines(filter_file))
  to_keep <- 1:length(first_clust)
  to_keep <- to_keep[!(to_keep %in% to_remove)]
  first_clust <- first_clust[to_keep]
}

######
# MANUAL CHANGES
######

# Loads in layer0 manually to map layer1 back to gene names
#writeLines(first_clust, file.path(layer_output_folder, "best_clusters_subset.txt"))
##layer0_clust <- readLines(file.path(output_folder, "orig_layer0", "best_clusters.txt"))
#layer0_clust <- readLines(file.path(output_folder, "layer0", "best_clusters.txt"))
#gene_list <- c()
#for (i in 1:length(first_clust)) {
#  clust <- strsplit(first_clust[i], ";")[[1]]
#  all_genes <- NULL
#  for (id in clust) {
#    clust_id <- as.numeric(strsplit(id, "_")[[1]][2])
#    clust_genes <- strsplit(layer0_clust[clust_id], ";")[[1]]
#    if(is.null(all_genes)) {
#      all_genes <- clust_genes
#    } else {
#      all_genes <- c(all_genes, clust_genes)
#    }
#  }
#  all_genes <- paste(sort(all_genes), collapse = ";")
#  gene_list <- c(gene_list, all_genes)
#  first_clust[i] <- all_genes
#}
#writeLines(gene_list, file.path(layer_output_folder, "clust_genes.txt"))

######
# MAIN SCRIPT
######

# Creates metagene dataframe for the first layer
layer_df <- data.frame(matrix(nrow = length(first_clust), ncol = ncol(scores)))
for (i in 1:length(first_clust)) {
  clust <- strsplit(first_clust[i], ";")[[1]]
  layer_df[i,] <- colMeans(scores[clust,])
}
rownames(layer_df) <- paste0("Clust0_", 1:nrow(layer_df))
colnames(layer_df) <- colnames(scores)
write.table(layer_df, file.path(layer_output_folder, "merged_scores.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Clusters each layer
previous_layer_clust <- first_clust
previous_layer_genes <- first_clust
previous_layer_df <- layer_df
for (layer in 1:n_layers) {

  # Makes output folder for the current layer
  layer_output_folder <- file.path(output_folder, paste0("layer", layer))
  if (!dir.exists(layer_output_folder)) { dir.create(layer_output_folder) }

  # Clusters the current layer
  gi_cluster(layer_df, layer_output_folder,
             metric_threshold = metric_thresholds[2],
             metric = metric,
             sort_metric = sort_metric,
             sort_decreasing = sort_decreasing,
             pick_num = 1,
             member_threshold = member_threshold,
             verbose = TRUE)

  # Loads clusters
  current_clust <- readLines(file.path(layer_output_folder, "best_clusters.txt"))

  # Creates metagene dataframe for the current layer
  layer_df <- data.frame(matrix(nrow = length(current_clust), ncol = ncol(scores)))
  for (i in 1:length(current_clust)) {
    clust <- strsplit(current_clust[i], ";")[[1]]
    layer_df[i,] <- colMeans(previous_layer_df[clust,])
  }
  clust_str <- paste0("Clust", layer, "_")
  rownames(layer_df) <- paste0(clust_str, 1:nrow(layer_df))
  colnames(layer_df) <- colnames(scores)
  write.table(layer_df, file.path(layer_output_folder, "merged_scores.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Maps clusters to gene names
  gene_list <- c()
  for (i in 1:length(current_clust)) {
    clust <- strsplit(current_clust[i], ";")[[1]]
    all_genes <- NULL
    for (id in clust) {
      clust_id <- as.numeric(strsplit(id, "_")[[1]][2])
      clust_genes <- strsplit(previous_layer_genes[clust_id], ";")[[1]]
      if(is.null(all_genes)) {
        all_genes <- clust_genes
      } else {
        all_genes <- c(all_genes, clust_genes)
      }
    }
    all_genes <- paste(sort(all_genes), collapse = ";")
    gene_list <- c(gene_list, all_genes)
  }
  writeLines(gene_list, file.path(layer_output_folder, "clust_genes.txt"))

  # Saves variables for next layer of clustering
  previous_layer_df <- layer_df
  previous_layer_clust <- current_clust
  previous_layer_genes <- gene_list
}

# Writes clustering data to zip
#fname <- file.path(output_folder, paste0(basename(output_folder), ".zip"))
#zip(fname, dir(output_folder, full.names = TRUE))
#stop("Script finished")
