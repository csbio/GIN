# Authors: Henry Ward, Xiang Zhang
# Last update: 12/13/2022

# Loads packages
packages <- c("reshape2", "dplyr", "rjson", "igraph", "RCy3")#, "reticulate", "rjson", "igraph", "NetPathMiner", "RCy3", "tools")
for (p in packages) {
  library(p, character.only = TRUE)
}

######
# UTILITY FUNCTIONS
######

# Removes edges for nodes that appear X or fewer times
filter_small_nodes <- function(adj_list, x) {
  adj_list[,1] <- as.character(adj_list[,1])
  adj_list[,2] <- as.character(adj_list[,2])
  unique_nodes <- unique(c(adj_list[,1], adj_list[,2]))
  node_count <- c()
  for (node in unique_nodes) {
    node_count <- c(node_count, sum(adj_list[,1] == node) + sum(adj_list[,2] == node))
  }
  to_remove <- unique_nodes[node_count <= x]
  adj_list <- adj_list[!(adj_list[,1] %in% to_remove) | !(adj_list[,2] %in% to_remove),]
  return(adj_list)
}

# Removes the bottom connections for each node if they have X or more nodes
filter_weak_edges <- function(adj_list, x) {
  adj_list[,1] <- as.character(adj_list[,1])
  adj_list[,2] <- as.character(adj_list[,2])
  unique_nodes <- unique(c(adj_list[,1], adj_list[,2]))
  to_remove <- c()
  for (node in unique_nodes) {
    ind <- adj_list[,1] == node | adj_list[,2] == node
    temp <- adj_list[ind,]
    temp <- temp[order(temp$value, decreasing = TRUE),]
    if (nrow(temp) > x) {
      to_remove <- c(to_remove, rownames(temp)[(x+1):nrow(temp)])
    }
  }
  to_remove <- unique(to_remove)
  adj_list <- adj_list[!(rownames(adj_list) %in% to_remove),]
  return(adj_list)
}

# Expands network by adding edges from mutual KNN
mutual_knn <- function(cor_mat, adj_list, k = 10, pcc = 0.4, min_pcc = 0.15, 
                       max_additions = NULL, max_edges_new_nodes = 1) {
  
  # Sets max connections to k if unspecified
  if (is.null(max_additions)) {
    max_additions <- k
  }
  
  # Gets correlation matrix and unique edges
  for (i in 1:nrow(adj_list)) {
    key <- paste(sort(c(adj_list[i,1], adj_list[i,2])), collapse = "_")
    adj_list$key[i] <- key
  }
  
  # Expands network by adding mutual KNN
  knn_adj_list <- data.frame(node1 = NA, node2 = NA, value = NA, connections = NA)
  unique_nodes <- unique(c(adj_list[,1], adj_list[,2]))
  counter <- 1
  for (node in unique_nodes) {
    knn <- sort(cor_mat[node,], decreasing = TRUE)[2:(k+1)]
    knn <- knn[knn > min_pcc]
    additions <- 0
    if (length(knn) > 0) {
      for (neighbor in names(knn)) {
        neighbor_knn <- sort(cor_mat[neighbor,], decreasing = TRUE)[2:(k+1)]
        neighbor_knn <- neighbor_knn[neighbor_knn > min_pcc]
        if (length(neighbor_knn) > 0) {
          if (node %in% names(neighbor_knn)) {
            key <- paste(sort(c(node, neighbor)), collapse = "_")
            if (!(key %in% adj_list$key)) {
              if (additions < max_additions) {
                connections <- names(neighbor_knn)[names(neighbor_knn) %in% unique_nodes]
                connections <- paste(connections, collapse = ";")
                knn_adj_list[counter,] <- c(node, neighbor, knn[neighbor][[1]], connections) 
                counter <- counter + 1
                additions <- additions + 1
              }
            }
          }
        }
      } 
    }
  }
  
  # Subsets to 1 edge per new node
  unique_new_nodes <- unique(knn_adj_list[,2])
  to_keep <- c()
  for (node in unique_new_nodes) {
    temp <- knn_adj_list[knn_adj_list[,2] == node,]
    temp <- temp[order(temp[,3], decreasing = TRUE),]
    max_edges <- min(nrow(temp), max_edges_new_nodes)
    to_keep <- c(to_keep, rownames(temp)[1:max_edges])
  }
  knn_adj_list <- knn_adj_list[to_keep,]
  
  # Removes anchor-anchor edges that don't meet the original PCC threshold
  to_keep <- rep(TRUE, nrow(knn_adj_list))
  for (i in 1:nrow(knn_adj_list)) {
    if (knn_adj_list[i,2] %in% unique_nodes & knn_adj_list[i,3] < pcc) {
      to_keep[i] <- FALSE
    }
  }
  knn_adj_list <- knn_adj_list[to_keep,]
  return(knn_adj_list)
}

# Expands network by adding edges from KNN run only once
single_knn <- function(cor_mat, adj_list, k = 10, pcc = 0.4, min_pcc = 0.15, 
                       max_additions = NULL, max_edges_new_nodes = 1) {
  
  # Sets max connections to k if unspecified
  if (is.null(max_additions)) {
    max_additions <- k
  }
  
  # Gets unique edges
  for (i in 1:nrow(adj_list)) {
    key <- paste(sort(c(adj_list[i,1], adj_list[i,2])), collapse = "_")
    adj_list$key[i] <- key
  }
  
  # Expands network by adding single KNN
  knn_adj_list <- data.frame(node1 = NA, node2 = NA, value = NA)
  unique_nodes <- unique(c(adj_list[,1], adj_list[,2]))
  counter <- 1
  for (node in unique_nodes) {
    knn <- sort(cor_mat[node,], decreasing = TRUE)[2:(k+1)]
    knn <- knn[knn > min_pcc]
    additions <- 0
    if (length(knn) > 0) {
      for (neighbor in names(knn)) {
        key <- paste(sort(c(node, neighbor)), collapse = "_")
        if (!(key %in% adj_list$key)) {
          if (additions < max_additions) {
            knn_adj_list[counter,] <- c(node, neighbor, knn[neighbor][[1]]) 
            counter <- counter + 1
            additions <- additions + 1
          }
        } 
      }
    }
  }
  
  # Removes anchor-anchor edges that don't meet the original PCC threshold
  to_keep <- rep(TRUE, nrow(knn_adj_list))
  for (i in 1:nrow(knn_adj_list)) {
    if (knn_adj_list[i,2] %in% unique_nodes & knn_adj_list[i,3] < pcc) {
      to_keep[i] <- FALSE
    }
  }
  knn_adj_list <- knn_adj_list[to_keep,]
  
  # Subsets to a single edge per new node
  unique_new_nodes <- unique(knn_adj_list[,2])
  to_keep <- c()
  for (node in unique_new_nodes) {
    temp <- knn_adj_list[knn_adj_list[,2] == node,]
    temp <- temp[order(temp[,3], decreasing = TRUE),]
    max_edges <- min(nrow(temp), max_edges_new_nodes)
    to_keep <- c(to_keep, rownames(temp)[1:max_edges])
  }
  knn_adj_list <- knn_adj_list[to_keep,]
  
  return(knn_adj_list)
}

# Adds additional nodes and edges to skeleton network
add_edges_to_network <- function(nodes, edges, new_edges, scale_dist = 0.1,
                                 jitter = 0.33) {
  
  # Gets distance at which to place new vertices
  rand_nodes <- sample(1:nrow(nodes), 50)
  pairs <- combn(rand_nodes, 2)
  rand_dist <- c()
  for (i in 1:ncol(pairs)) {
    point1 <- c(nodes$x[pairs[1,i]], nodes$y[pairs[1,i]])
    point2 <- c(nodes$x[pairs[2,i]], nodes$y[pairs[2,i]])
    rand_dist <- c(rand_dist, sqrt(sum((point1 - point2)^2)))
  }
  rand_dist <- mean(rand_dist) * scale_dist
  
  # Places new vertices
  current_vertices <- unique(c(edges$from, edges$to))
  to_place <- new_edges$node2[!(new_edges$node2 %in% current_vertices)]
  counter <- nrow(nodes) + 1
  new_nodes <- nodes
  for (gene in to_place) {
    partners <- new_edges$node1[new_edges$node2 == gene]
    if (length(partners) == 1) {
      partner_x <- as.numeric(nodes$x[nodes$id == partners[1]])
      partner_y <- as.numeric(nodes$y[nodes$id == partners[1]])
      angle <- runif(1) * 2 * pi
      jitter_x <- rnorm(1, 0, jitter * rand_dist)
      jitter_y <- rnorm(1, 0, jitter * rand_dist)
      new_x <- partner_x + (cos(angle) * (rand_dist + jitter_x))
      new_y <- partner_y + (sin(angle) * (rand_dist + jitter_y))
      new_nodes[counter,] <- c(gene, new_x, new_y)
    } else if (length(partners) > 1) {
      print(gene)
      edge_weights <- as.numeric(new_edges$value[new_edges$node2 == gene])
      max_weight <- max(edge_weights)
      min_weight <- min(edge_weights)
      norm_weights <- lapply(edge_weights, function(x) (x - min_weight) / (max_weight - min_weight))
    }
    counter <- counter + 1
  }
  new_nodes$x <- as.double(new_nodes$x)
  new_nodes$y <- as.double(new_nodes$y)
  
  # Adds new edges to existing edges
  colnames(new_edges) <- c("from", "to", "weight")
  new_edges <- rbind(edges, new_edges)
  new_edges$weight <- as.numeric(new_edges$weight)
  
  # Returns new nodes and edges
  results <- list()
  results[["new_nodes"]] <- new_nodes
  results[["new_edges"]] <- new_edges
  return(results)
}

######
# CORE NETWORK GENERATION
######

# Begin script
setwd("C:\\project")

output_folder <- file.path("output", "GI_network")
if (!dir.exists(output_folder)) { dir.create(output_folder) }

# Sets important parameters
filter_method <- "pcc"

# Reads in data (original scores or pcc matrix)
df <- read.csv(file.path("input", "qGI_centroid12_rm_no507Var2.txt"), sep = "\t", row.names = 1)
#cor_mat <- data.matrix(df)

# Writes node labels to file
node_file <- file.path("output", "node_labels.txt")
node_labels <- data.frame(name = rownames(df), gene = rownames(df))
write.table(node_labels, node_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Generates adjacency list from PCC matrix for all edges
cor_mat <- cor(t(df), use = "complete.obs")
cor_mat[lower.tri(cor_mat, diag = TRUE)] <- 50
all_edges <- melt(cor_mat)
all_edges <- all_edges %>%
  filter(value != 50)
rm(cor_mat)
cor_mat <- cor(t(df), use = "complete.obs")

# Sets PCC thresholds and filtering parameters
thresh_list <- c(0.41)
edge_filter <- c(30)

# Generates PCC networks with different edge cutoffs and filtering parameters
for (i in 1:length(thresh_list)) {
  for (edge_thresh in edge_filter) {

    # Generates PCC network
    thresh <- thresh_list[i]
    edges <- all_edges %>%
      filter(value > thresh) %>%
      filter_weak_edges(edge_thresh)
    colnames(edges) <- c("node1", "node2", "value")
    edges <- edges[complete.cases(edges),]

    # Expands network with KNN approach if specified
    unique_nodes <- unique(c(edges[,1], edges[,2]))
    print(paste("Edges:", nrow(edges)))
    print(paste("Nodes:", length(unique_nodes)))
    file_str <- paste0("filter_", edge_thresh, "_", gsub("\\.", "", thresh), ".txt")
    edge_file <- file.path(output_folder, paste0("edges_", file_str))
    node_file <- file.path(output_folder, paste0("nodes_", file_str))
    write.table(edges, edge_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.table(unique_nodes, node_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
stop()

######
# MANUAL EDITS
######

# Before the expanded network is generated, we load the network up in cytoscape,
# remove all connected components save for the largest, and write to file again.

######
# EXPANDED NETWORK GENERATION
######

# Loads largest connected component from the network as well as GIN data
setwd("C:\\project")

output_folder <- file.path("output", "GI_network")
if (!dir.exists(output_folder)) { dir.create(output_folder) }
json_network <- fromJSON(file = file.path(output_folder, "mc_core_c12rm_no507Var2_041_30.cyjs"))

# Generates network from skeleton network
json_nodes <- json_network$elements$nodes
json_edges <- json_network$elements$edges
nodes <- data.frame(id = NA, x = NA, y = NA)
edges <- data.frame(from = NA, to = NA, weight = NA)
for (i in 1:length(json_nodes)) {
  nodes[i,1] <- json_nodes[[i]][["data"]][["name"]]
  nodes[i,2] <- json_nodes[[i]][["position"]][["x"]]
  nodes[i,3] <- json_nodes[[i]][["position"]][["y"]]
}
for (i in 1:length(json_edges)) {
  shared_name <- strsplit(json_edges[[i]][["data"]][["shared_name"]], " ")[[1]]
  edges[i,1] <- shared_name[1]
  edges[i,2] <- shared_name[4]
  edges[i,3] <- json_edges[[i]][["data"]][["value"]]
}
nodes$x <- as.numeric(nodes$x)
nodes$y <- as.numeric(nodes$y)
nodes$x <- nodes$x + abs(min(nodes$x))
nodes$y <- nodes$y + abs(min(nodes$y))
network <- graph_from_data_frame(vertices = nodes, d = edges, directed = FALSE)
plot(network, vertex.size = 1, vertex.label = NA)

# Runs KNN on network
temp_edges <- edges
colnames(temp_edges) <- c("node1", "node2", "value")
knn_edges <- single_knn(cor_mat, temp_edges, k = 35, pcc = 0.41, min_pcc = 0.10)
results <- add_edges_to_network(nodes, edges, knn_edges, scale_dist = 0.1, jitter = 0)
results$new_nodes$x <- results$new_nodes$x + abs(min(results$new_nodes$x))
results$new_nodes$y <- results$new_nodes$y + abs(min(results$new_nodes$y))
cat(paste("Number of genes post-KNN:", nrow(results$new_nodes), "\n"))

# Generates new layout within open cytoscape connection and writes to file
knn_network <- graph_from_data_frame(vertices = results$new_nodes, d = results$new_edges, directed = FALSE)
plot(knn_network, vertex.size = 1, vertex.label = NA)
name_str <- "c12rm_no507Var2_041new_knn35_min010"
createNetworkFromIgraph(knn_network, name_str, collection = "expanded")
write_graph(knn_network, file.path(output_folder, paste0(name_str, ".gml")), format = "gml")