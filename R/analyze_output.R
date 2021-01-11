library(tidyverse)
library(ggplot2)
library(superheat)
library(igraph)
library(jsonlite)
library(ggraph)
library(tidygraph)

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

summarize_network_connectivity <- function(network_df) {
  network_df %>% group_by(from, to) %>% summarize(count = n())
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



setwd('~/Workspace/transmission_nets/test/outputs/ModelFiveTests/CoreTest/')
# setwd('~/Workspace/transmission_nets/test/outputs/ModelThreeTests/CoreTest/')
llik <- read_csv('params/likelihood.csv', col_names = c("llik"))
geo_gen_prob <- read_csv('params/geo_gen_prob.csv', col_names = c("geo_gen"))
mean_coi <- read_csv('params/mean_coi.csv', col_names = c("coi"))
loss_prob <- read_csv('params/loss_prob.csv', col_names = c("prob"))
mutation_prob <- read_csv('params/mutation_prob.csv', col_names = c("prob"))
orderings <- read_csv('params/infection_order.csv', col_names = F)
sets <- read_csv('stats/1_ps.csv', col_names = c("set", "llik", "iter"))
networks <- read_networks('network.csv')
fpr <- read_csv('params/fpr.csv', col_names = c("fpr"))
fnr <- read_csv('params/fnr.csv', col_names = c("fnr"))


plot(1:length(llik$llik), llik$llik, type = "l")
plot(1:length(mean_coi$coi), mean_coi$coi, type = "l")
qplot(geo_gen_prob$geo_gen)
qplot(mean_coi$coi)
qplot(mutation_prob$prob)
qplot(loss_prob$prob)
qplot(fpr$fpr)
qplot(fnr$fnr)



precedes_matrix <- generate_precedents_matrix(orderings)

superheat(precedes_matrix, heat.pal = c("red", "white", "blue"), X.text = round(precedes_matrix, 1))
superheat(precedes_matrix, heat.pal = c("red", "white", "blue"))

res = data.frame(before = c(), after = c(), prop = c())
for (i in 1:dim(precedes_matrix)[1]) {
  for (j in i:dim(precedes_matrix)[1]) {
      res <- rbind(res, data.frame(before = i, after = j, prop = precedes_matrix[i, j]))
  }
}

res %>% filter(prop >= 1) %>% View

edge_data_frame %>% filter

network_summary <- summarize_network_connectivity(networks)
inferred_tnet <- graph_from_data_frame(network_summary)
inferred_tnet <- delete_vertices(inferred_tnet, "S")
inferred_tnet <- set_vertex_attr(inferred_tnet, "size", value = 8)

plot(inferred_tnet)


inferred_clusters <- clusters(inferred_tnet)
true_clusters <- clusters(tnet)

inferred_cluster_matrix <- generate_clustering_matrix(inferred_clusters$membership)
true_cluster_matrix <- generate_clustering_matrix(true_clusters$membership)

inferred_cluster_matrix <- inferred_cluster_matrix[order(rownames(inferred_cluster_matrix)), order(colnames(inferred_cluster_matrix))]
true_cluster_matrix <- true_cluster_matrix[order(rownames(true_cluster_matrix)), order(colnames(true_cluster_matrix))]

inferred_true_diff <- inferred_cluster_matrix - true_cluster_matrix

superheat(inferred_true_diff, pretty.order.rows = T, pretty.order.cols = T)

superheat(inferred_cluster_matrix, pretty.order.rows = T, pretty.order.cols = T)
superheat(true_cluster_matrix, pretty.order.rows = T, pretty.order.cols = T)

input_data <- read_json("~/Workspace/transmission_nets/test/resources/JSON/nodes3.json")
input_data2 <- read_json("~/Workspace/transmission_nets/test/resources/JSON/nodes2.json")


# Genotype is a string of zeros and ones
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

distance_matrices_100 <- calculate_genetic_dist_matrix(input_data$nodes, hamming_distance)
distance_matrices_20 <- calculate_genetic_dist_matrix(input_data2$nodes, hamming_distance)

superheat(distance_matrices_20$latent_distance_matrix, X.text = distance_matrices_20$latent_distance_matrix)
superheat(distance_matrices_20$observed_distance_matrix, X.text = distance_matrices_20$observed_distance_matrix)


ps_prob_list = list()
for (file in dir('stats/')) {
  sets <- read_csv(paste0('stats/', file), col_names = c("set", "llik", "iter"))
  ps_probs <- sets %>%
    # filter(iter > 3500) %>%
    mutate(llik = if_else(is.na(llik), 0, llik)) %>%
    pivot_wider(names_from = set, values_from = llik, values_fill = 0)
  node = str_split(file, "_")[[1]][1]
  prop = ps_probs %>% select(-iter) %>% colSums() / nrow(ps_probs)
  ps_prob_list[[node]] = prop
}


network <- data.frame(from=c(), to=c(), weight=c())
for(node in names(ps_prob_list)) {
  parent_set_lliks = ps_prob_list[[node]]
  # Doesn't handle multiple parents yet
  parent_sets = unlist(str_split(str_remove(str_remove(names(parent_set_lliks), "\\{"), "\\}"), ","))
  names(parent_set_lliks) <- NULL
  print(parent_sets)
  network <- rbind(network, data.frame(from=parent_sets, to=node, weight=parent_set_lliks))
}

r_c <- network %>%
  filter(from != "S") %>%
  group_by(from) %>%
  summarize(bf = sum(weight))

pos_nets_max_llik <- network %>% group_by(to) %>% filter(weight == max(weight)) %>% arrange(to)
max_llik_graph <- graph_from_data_frame(pos_nets_max_llik)
inf_graph <- graph_from_data_frame(network)

tmax_llik_graph <- tidygraph::as_tbl_graph(max_llik_graph)
tinf_graph <- tidygraph::as_tbl_graph(inf_graph)

ggraph(tmax_llik_graph %>% activate(nodes) %>% filter(name != "S")) +
  geom_edge_fan(show.legend = FALSE) +
  geom_node_point()

ggraph(tinf_graph %>% activate(nodes) %>% filter(name != "S") %>% activate(edges) %>% filter(weight >= .25)) +
  geom_edge_fan(show.legend = FALSE, aes(edge_alpha = weight)) +
  geom_node_point()


plot(delete_vertices(max_llik_graph, "S"))
plot(delete_vertices(inf_graph, "S"))

input_data2 <- read_json("~/Workspace/transmission_nets/test/resources/JSON/nodes2.json")

true_dat2 <- data.frame(from=c(), to=c())
for (el in input_data2$network) {
  print(el)
  true_dat2 <- rbind(true_dat2, data.frame(from=el$from, to=el$to))
}

true_graph <- graph_from_data_frame(true_dat2)
plot(delete_vertices(true_graph, "source"))

