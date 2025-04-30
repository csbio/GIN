
library(ggplot2)
library(ggthemes)

# Clusters scores by ranking valid hierarchical clusters by the mean of their
# Pearson correlation SDs across many iterations.
# metric: one of "sd" or "snr" (default "sd")
# sort_metric: one of "size", "sd" or "snr" (default "size")
# sort_decreasing: whether to sort by the metric in decreasing order (default TRUE)
gi_cluster <- function(scores, output_folder, metric_threshold, 
                       pick_num = 5, member_threshold = 3,
                       metric = "sd", sort_metric = "size",
                       sort_decreasing = TRUE, verbose = TRUE,
                       save_intermediate = TRUE, load_intermediate = FALSE) {

  # Checks arguments
  if (!(metric %in% c("sd", "snr", "pcc"))) {
    stop("ERROR: metric must be one of 'sd', 'pcc' or 'snr'")
  }
  if (!(sort_metric %in% c("size", "sd", "snr"))) {
    stop("ERROR: sort_metric must be one of 'size,' 'sd' or 'snr'")
  }

  # Makes output folder if necessary
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  # Sets path to optional intermediate output file
  intermediate_file <- file.path(output_folder, "intermediate_output.rda")
  
  # Computes the PCC distance matrix and gets initial clusters
  cor_mat <- cor(t(scores))

  # Main algorithm
  remaining_genes <- rownames(scores)
  remaining_genes_list <- -1
  prev_remaining_genes <- -1
  final_clusters <- rep(NA, nrow(scores) + 1)
  count <- 1
  iter <- 0
  dist_mat <- NA
  while (length(remaining_genes) > 0) {

    # Saves to file if specified or loads intermediate file if specified
    if (save_intermediate & !load_intermediate) {
      if (verbose) {
        cat(paste0("Saving intermediate output to ", intermediate_file, "...\n"))
      }
      save(remaining_genes, remaining_genes_list, prev_remaining_genes, final_clusters, count,
           iter, file = intermediate_file)
    } else if (load_intermediate) {
      load(intermediate_file)
    }

    # Gets number of remaining genes
    remaining_genes_list <- c(remaining_genes_list, length(remaining_genes))
    if (verbose) {
      cat(paste("Picking", pick_num, "clusters from", length(remaining_genes), "remaining genes\n"))
    }
    
    # Gets dendrogram
    dist_mat <- as.dist(1 - cor_mat[rownames(cor_mat) %in% remaining_genes,
                                    colnames(cor_mat) %in% remaining_genes])
    clusts <- hclust(dist_mat, method = "average")
    dend <- as.dendrogram(clusts)
    
    # Traverses tree to get all clusters
    clust_df <- get_clusters_outer(dend)
    
    # Subsets clusters to those with at least 2 members
    larger_clust <- clust_df[unlist(lapply(clust_df, function(x) length(x) > 1))]
    
    # Converts list of clusters to dataframe
    if (length(larger_clust) > 1) {
      clust_df <- data.frame(genes = I(larger_clust), 
                             PCC = NA,
                             sd = NA,
                             members = NA)
    } else if (length(larger_clust) == 1) {
      clust_df <- data.frame(genes = larger_clust[[1]], 
                             PCC = NA,
                             sd = NA,
                             members = NA)
    } else {
      new_clust_df <- data.frame(genes = NA, PCC = NA, sd = NA, members = NA)
      new_clust_df$genes[1] <- list(clust_df)
      clust_df <- new_clust_df
    }
    
    # Gets the mean PCC, mean SD of PCCs, SNR and # of members for all clusters
    for (i in 1:nrow(clust_df)) {
      temp_cor <- cor_mat[rownames(cor_mat) %in% clust_df$genes[i][[1]],
                          colnames(cor_mat) %in% clust_df$genes[i][[1]]]
      temp_cor[upper.tri(temp_cor, diag = TRUE)] <- NA
      clust_df$PCC[i] <- mean(temp_cor, na.rm = TRUE)
      clust_df$sd[i] <- sd(temp_cor, na.rm = TRUE)
      clust_df$members[i] <- length(clust_df$genes[i][[1]])
      clust_df$snr[i] <- clust_df$PCC[i] / clust_df$sd[i]
    }
    
    # For the typical case, we filter clusters to those with a set number of members
    # (typically 3 or more) and those which have SDs less than the specified quantile
    if (nrow(clust_df) > pick_num) {
      clust_df <- clust_df[clust_df$members >= member_threshold,]
      if (metric == "sd") {
        thresh <- quantile(clust_df$sd, metric_threshold, na.rm = TRUE)
        if (sum(clust_df$sd <= thresh, na.rm = TRUE) > 0) {
          clust_df <- clust_df[clust_df$sd <= thresh & complete.cases(clust_df$sd),]
        }
      } else if (metric == "snr") {
        thresh <- quantile(clust_df$snr, 1 - metric_threshold, na.rm = TRUE)
        if (sum(clust_df$snr >= thresh, na.rm = TRUE) > 0) {
          clust_df <- clust_df[clust_df$snr >= thresh & complete.cases(clust_df$snr),]
        }
      } else if (metric == "pcc") {
        thresh <- quantile(clust_df$PCC, 1 - metric_threshold, na.rm = TRUE)
        if (sum(clust_df$PCC >= thresh, na.rm = TRUE) > 0) {
          clust_df <- clust_df[clust_df$PCC >= thresh & complete.cases(clust_df$PCC),]
        }
      }
    }
    
    # Sorts clusters by # of members, SD or SNR
    if (sort_metric == "size") {
      clust_df <- clust_df[order(clust_df$members, decreasing = sort_decreasing),]
    } else if (sort_metric == "sd") {
      clust_df <- clust_df[order(clust_df$sd, decreasing = sort_decreasing),]
    } else if (sort_metric == "snr") {
      clust_df <- clust_df[order(clust_df$snr, decreasing = sort_decreasing),]
    }
    
    # Takes out the highest-ranked clusters and adds to list of final clusters unless
    # we have no clusters left, in which case we break
    if (nrow(clust_df) > 0) {
      best_clusts <- clust_df$genes[1:min(pick_num, nrow(clust_df))]
      for (j in 1:length(best_clusts)) {
        current_clust <- best_clusts[j][[1]]
        final_clusters[count] <- paste(current_clust, collapse = ";") 
        remaining_genes <- remaining_genes[!(remaining_genes %in% current_clust)]
        count <- count + 1
      }
    }

    # Additional check for the end of the loop - if we end up unable to pick any clusters, 
    # we return the remaining genes as a single cluster and finish the algorithm
    if (length(remaining_genes) == prev_remaining_genes | length(remaining_genes) == 1 |
        nrow(clust_df) == 0) {
      final_clusters[count] <- paste(remaining_genes, collapse = ";") 
      remaining_genes <- c()
      count <- count + 1
      if (verbose) {
        cat(paste("Unable to sort clusters - returning leftover genes and finishing\n"))
      }
      break
    }

    # Iterates counters
    iter <- iter + 1
    prev_remaining_genes <- length(remaining_genes)
  }
  
  # Formats resulting dataframe and writes to file
  final_clusters <- final_clusters[!is.na(final_clusters)]
  empty_ind <- unlist(lapply(final_clusters, function(x) x == ""))
  final_clusters <- final_clusters[!empty_ind]
  sink(file.path(output_folder, "best_clusters.txt"))
  cat(paste(final_clusters, collapse = "\n"))
  cat("\n")
  sink()
  
  # Makes plot of iteration vs. remaining genes and saves to file
  remaining_genes_list <- remaining_genes_list[!(remaining_genes_list == -1)]
  remaining_genes_list <- remaining_genes_list[1:iter]
  temp <- data.frame(iter = 1:iter, remaining_genes = remaining_genes_list)
  ggplot(temp, aes(x = iter, y = remaining_genes)) + 
    geom_line() + 
    xlab("Iteration") +
    ylab("Remaining genes") +
    theme_tufte(base_size = 18)
  ggsave(file.path(output_folder, "remaining_genes.png"))
  
  # Explicitly returns clusters
  return(final_clusters)
}

# Recursive function to get all clusters in a dendrogram
get_clusters_outer <- function(dend) {
  type1 <- typeof(dend[[1]])
  type2 <- typeof(dend[[2]])
  if (type1 == "integer" & type2 == "integer") {
    gene1 <- get_clusters_inner(dend[[1]])
    gene2 <- get_clusters_inner(dend[[2]])
    return(c(list(gene1),
             list(gene2)))
  } else if (type1 == "integer") {
    return(c(list(get_clusters_inner(dend[[1]])),
             list(get_clusters_inner(dend[[2]])),
             get_clusters_outer(dend[[2]])))
  } else if (type2 == "integer") {
    return(c(list(get_clusters_inner(dend[[1]])),
             list(get_clusters_inner(dend[[2]])),
             get_clusters_outer(dend[[1]])))
  } else {
    return(c(list(get_clusters_inner(dend[[1]])),
             list(get_clusters_inner(dend[[2]])),
             get_clusters_outer(dend[[1]]),
             get_clusters_outer(dend[[2]]))) 
  }
}

# Inner recursive function for the above
get_clusters_inner <- function(dend) {
  if (!is.null(attr(dend, "leaf"))) {
    return(c(attr(dend, "label")))
  } else {
    return(c(get_clusters_inner(dend[[1]]), 
             get_clusters_inner(dend[[2]])))
  }
}
