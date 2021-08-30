

does_node_precede <- function(nodes, positions, node_i, node_j) {
  precedes = positions[which(nodes == node_i)] < positions[which(nodes == node_j)]
  return(precedes)
}

read_networks <- function(path, edge_delim='-', pair_delim = ',') {
  networks <- read_lines(path)
  res <- data.frame(from = c(), to = c(), t = c())
  for (i in 1:length(networks)) {
    network = networks[i]
    pairs = str_split(network, pattern = pair_delim, simplify = TRUE)
    split_pairs = str_split(pairs, pattern = edge_delim, simplify = TRUE)
    from = split_pairs[,1]
    to = split_pairs[,2]
    res = rbind(res, data.frame(from = from, to = to, t = i))
  }
  return(res)
}

summarize_network_connectivity <- function(network) {
  connectivity <- network %>% group_by(from, to) %>% summarize(count = n())
  return(connectivity)
}

generate_precedents_matrix <- function(orderings) {
  long_ord <- orderings %>%
    rowid_to_column() %>%
    pivot_longer(cols = -rowid, names_to = "position", values_to = "node") %>%
    mutate(position = as.numeric(str_remove(position, pattern = "X")))

  nodes <- sort(unique(long_ord$node))

  precedents_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes))

  for (i in 1:(length(nodes) - 1)) {
    node_i = nodes[i]
    print(i)
    for (j in (i + 1):length(nodes)) {
      node_j = nodes[j]
      tmp <- long_ord %>% group_by(rowid) %>% summarise(i_before_j = does_node_precede(node, position, node_i, node_j))
      precedents_matrix[i, j] = mean(tmp$i_before_j)
      precedents_matrix[j, i] = 1 - precedents_matrix[i, j]
    }
  }
  return(precedents_matrix)
}

generate_clustering_matrix <- function(cluster_membership) {
  node_labels <- names(cluster_membership)
  mat <- matrix(nrow = length(node_labels), ncol = length(node_labels))

  for (i in 1:length(node_labels)) {
    node_i_label <- node_labels[i]
    for (j in i:length(node_labels)) {
      node_j_label <- node_labels[j]
      mat[i,j] <- mat[j,i] <- as.numeric(cluster_membership[node_i_label] == cluster_membership[node_j_label])
    }
  }
  colnames(mat) <- rownames(mat) <- node_labels
  return(mat)
}

hamming_distance <- function(genotype1, genotype2) {
  g1 <- as.logical(as.numeric(str_split(genotype1, "", simplify = T)))
  g2 <- as.logical(as.numeric(str_split(genotype2, "", simplify = T)))
  return(sum(xor(g1, g2)))
}

calculate_genetic_dist_matrix <- function(nodes, metric) {
  total_nodes <- length(nodes)
  latent_distance_matrix <- matrix(nrow = total_nodes, ncol = total_nodes)
  observed_distance_matrix <- matrix(nrow = total_nodes, ncol = total_nodes)
  node_names <- c()

  for (i in 1:total_nodes) {

    node_i = nodes[[i]]
    total_genotypes <- length(node_i$observed_genotype)
    node_names <- c(node_names, node_i$id)

    for (j in i:total_nodes) {
      node_j = nodes[[j]]
      latent_distance = 0
      observed_distance = 0
      if (i != j) {
        for (k in 1:total_genotypes) {
          latent_distance = latent_distance + metric(
            node_i$latent_genotype[[k]]$genotype,
            node_j$latent_genotype[[k]]$genotype
          )
          observed_distance = observed_distance + metric(
            node_i$observed_genotype[[k]]$genotype,
            node_j$observed_genotype[[k]]$genotype
          )
        }
      }
      latent_distance_matrix[i, j] = latent_distance
      latent_distance_matrix[j, i] = latent_distance
      observed_distance_matrix[i, j] = observed_distance
      observed_distance_matrix[j, i] = observed_distance
    }
  }

  rownames(latent_distance_matrix) <- node_names
  colnames(latent_distance_matrix) <- node_names
  rownames(observed_distance_matrix) <- node_names
  colnames(observed_distance_matrix) <- node_names
  return(list(observed_distance_matrix = observed_distance_matrix,
              latent_distance_matrix = latent_distance_matrix))
}

extract_parent_sets <- function(dir_path) {
  ps_prob_list = list()
  for (file in dir(dir_path)) {
    sets <- read_csv(paste0(dir_path, file),
                     col_names = c("set", "llik", "iter"),
                     col_types = cols(set = col_character(), llik = col_double(), iter = col_double()), skip = 1)
    ps_probs <- sets %>%
      mutate(llik = if_else(is.na(llik), 0, llik)) %>%
      pivot_wider(names_from = set, values_from = llik, values_fill = 0)
    node = str_split(file, "_")[[1]][1]
    prop = ps_probs %>% select(-iter) %>% colSums() / nrow(ps_probs)
    ps_prob_list[[node]] = prop
  }
  return(ps_prob_list)
}

network_from_parent_sets <- function(parent_set_list) {
  network <- data.frame(from=c(), to=c(), weight=c())
  for(node in names(parent_set_list)) {
    parent_set_lliks = parent_set_list[[node]]
    # Doesn't handle multiple parents yet
    parent_sets = unlist(str_split(str_remove(str_remove(names(parent_set_lliks), "\\{"), "\\}"), ","))
    network <- rbind(network, data.frame(from=parent_sets, to=node, weight=parent_set_lliks))
  }
  return(network)
}

rc_from_weighted_network <- function(network) {
 rc <- network %>%
  filter(from != "S") %>%
  group_by(from) %>%
  summarize(rc = sum(weight))
 return(rc)
}

rc_from_network <- function(network) {
  rc <- network %>%
    group_by(from) %>%
    summarize(rc = sum(n()))
  return(rc)
}

generate_matched_layouts <- function(target_network, other_network, ...) {
  target_layout <- create_layout(target_network, ...)
  other_layout <- create_layout(other_network, ...)
  target_layout_points <- target_layout %>% select(x, y, name)
  other_layout_points <- data.frame(name = other_layout$name) %>% left_join(target_layout_points)
  other_layout$x = other_layout_points$x
  other_layout$y = other_layout_points$y
 return(list(target = target_layout, other = other_layout))
}

match_layout <- function(network, target_layout) {
  other_layout <- create_layout(network, layout = "nicely")
  target_layout_points <- target_layout %>% select(x, y, name)
  other_layout_points <- data.frame(name = other_layout$name) %>% left_join(target_layout_points)
  other_layout$x = other_layout_points$x
  other_layout$y = other_layout_points$y
  return(other_layout)
}

extract_true_network_tibble <- function(path) {
  input_data <- read_json(path)
  true_network_edges <- data.frame(from=c(), to=c())
  for (el in input_data$network) {
    true_network_edges <- rbind(true_network_edges, data.frame(from=el$from, to=el$to))
  }

  true_network_vertices <- tibble(name=c(), infection_duration=c(), infection_time=c(), observation_time=c())
  for (el in input_data$nodes) {
    true_network_vertices <- rbind(true_network_vertices,
                                   tibble(name=el$id, infection_duration=el$infection_duration,
                                          infection_time=el$infection_time, observation_time=el$observation_time,
                                          observed_genotype=list(el$observed_genotype), latent_genotype=list(el$latent_genotype)))
  }

  return(list(edges = true_network_edges, vertices = true_network_vertices))
}
