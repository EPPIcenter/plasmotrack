library(tidyverse)
library(ggplot2)
library(superheat)
library(igraph)
library(jsonlite)
library(ggraph)
library(tidygraph)

setwd('~/Workspace/transmission_nets/test/outputs/ModelFiveTests/CoreTest/')
# setwd('~/Workspace/transmission_nets/test/outputs/ModelThreeTests/CoreTest/')
llik <- read_csv('parameters/likelihood.csv', col_names = c("llik"))
geo_gen_prob <- read_csv('parameters/geo_gen_prob.csv', col_names = c("geo_gen"))
mean_coi <- read_csv('parameters/mean_coi.csv', col_names = c("coi"))
loss_prob <- read_csv('parameters/loss_prob.csv', col_names = c("prob"))
mutation_prob <- read_csv('parameters/mutation_prob.csv', col_names = c("prob"))
orderings <- read_csv('parameters/infection_order.csv', col_names = F)
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



distance_matrices_100 <- calculate_genetic_dist_matrix(input_data$nodes, hamming_distance)
distance_matrices_20 <- calculate_genetic_dist_matrix(input_data2$nodes, hamming_distance)

superheat(distance_matrices_20$latent_distance_matrix, X.text = distance_matrices_20$latent_distance_matrix)
superheat(distance_matrices_20$observed_distance_matrix, X.text = distance_matrices_20$observed_distance_matrix)








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

input_data2 <- read_json("~/Workspace/transmission_nets/test/resources/JSON/nodes7.json")

true_dat2 <- data.frame(from=c(), to=c())
for (el in input_data2$network) {
  print(el)
  true_dat2 <- rbind(true_dat2, data.frame(from=el$from, to=el$to))
}

true_graph <- graph_from_data_frame(true_dat2)
tinf_graph_sub <- tinf_graph %>%
  activate(nodes) %>%
  filter(name != "S") %>%
  activate(edges) %>%
  filter(weight >= .25)

layouts <- generate_matched_layouts(true_graph, tinf_graph_sub)
true_layout <- layouts$target
tinf_layout <- layouts$other


ggraph(true_layout) +
  geom_edge_fan() +
  geom_node_point(size = 8) +
  geom_node_text(aes(label = name), color = "white")

ggraph(tinf_layout) +
  geom_edge_fan() +
  geom_node_point(size = 8) +
  geom_node_text(aes(label = name), color = "white")




r_c_cmp <- r_c %>% left_join(true_r_c, by="from") %>% replace_na(list(bf.y = 0))
mean(r_c_cmp$bf.x)
mean(r_c_cmp$bf.y)
