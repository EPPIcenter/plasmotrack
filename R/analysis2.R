library(tidyverse)
library(ggplot2)
library(superheat)
library(igraph)
library(jsonlite)
library(ggraph)
library(tidygraph)


setwd('~/Workspace/transmission_nets/test/outputs/ModelFiveTests_small2/CoreTest/')
llik <- read_csv('parameters/likelihood.csv', col_names = c("llik"))
geo_gen_prob <- read_csv('parameters/geo_gen_prob.csv', col_names = c("geo_gen"))
mean_coi <- read_csv('parameters/mean_coi.csv', col_names = c("coi"))
loss_prob <- read_csv('parameters/loss_prob.csv', col_names = c("prob"))
mutation_prob <- read_csv('parameters/mutation_prob.csv', col_names = c("prob"))
orderings <- read_csv('parameters/infection_order.csv', col_names = F)
infection_duration <- cbind(read_csv('parameters/infection_duration_shape.csv', col_names = "shape"),
                            read_csv('parameters/infection_duration_scale.csv', col_names = "scale"))

input_data <- read_json("~/Workspace/transmission_nets/test/resources/JSON/nodes8.json")
true_network_edges <- data.frame(from=c(), to=c())
for (el in input_data$network) {
  true_network_edges <- rbind(true_network_edges, data.frame(from=el$from, to=el$to))
}

true_network_vertices <- data.frame(name=c(), infection_duration=c(), infection_time=c(), observation_time=c())
for (el in input_data$nodes) {
  true_network_vertices <- rbind(true_network_vertices,
                                 data.frame(name=el$id, infection_duration=el$infection_duration,
                                            infection_time=el$infection_time, observation_time=el$observation_time))
}

from_nodes <- true_network_edges %>% pull(from) %>% unique()
to_nodes <- true_network_edges %>% pull(to) %>% unique()
all_nodes <- true_network_vertices %>% pull(name) %>% unique()
is_source_node <- !(all_nodes %in% to_nodes)
true_network_vertices <- true_network_vertices %>% mutate(is_source = is_source_node)

tidy_network <- tidygraph::as_tbl_graph(graph_from_data_frame(true_network_edges, vertices = true_network_vertices))

parent_sets <- extract_parent_sets("stats/")
inf_network_df <- network_from_parent_sets(parent_sets)

is_source <- inf_network_df %>% filter(from == "S") %>% mutate(is_source = weight > .5) %>% select(to, is_source) %>% rename(name = to)
inf_network_vertices <- true_network_vertices %>% select(-is_source) %>% left_join(is_source, by = "name")
inf_network <- graph_from_data_frame(inf_network_df %>% filter(from != "S"), vertices = inf_network_vertices)
tinf_network <- tidygraph::as_tbl_graph(inf_network)


tinf_graph_sub <- tinf_network %>%
  activate(edges) %>%
  filter(weight > 0)


layouts <- generate_matched_layouts(tidy_network, tinf_graph_sub, layout = "kk")
# layouts <- generate_matched_layouts(tidy_network, tidy_network, layout = "kk")

true_layout <- layouts$target
tinf_layout <- layouts$other

google_pres_theme <- theme(panel.background = element_rect(fill = "#1A212B", color = "#1A212B"),
                           plot.background = element_rect(fill = "#1A212B", color = "#1A212B"),
                           plot.title = element_text(color = "white"),
                           legend.background = element_rect(fill = "#1A212B", color = "#1A212B"),
                           legend.key = element_rect(fill = "#1A212B", color = "#1A212B"),
                           legend.text = element_text(color = "white"),
                           legend.title = element_text(color = "white"))


white_theme <-  theme(panel.background = element_rect(fill = "white", color = "white"),
                      plot.background = element_rect(fill = "white", color = "white"),
                      plot.title = element_text(color = "black"),
                      legend.background = element_rect(fill = "white", color = "white"),
                      legend.key = element_rect(fill = "white", color = "white"),
                      legend.text = element_text(color = "black"),
                      legend.title = element_text(color = "black"))

ggraph(true_layout) +
  geom_edge_fan(color = "black") +
  geom_node_point(size = 8, aes(color = is_source)) +
  geom_node_text(aes(label = name), color = "black") +
  ggtitle("True Network - 10 Loci") +
  white_theme +
  # google_pres_theme +
  ggsave("~/Desktop/lab meeting presentation/small_true.pdf", dpi = "retina")

ggraph(tinf_layout) +
  geom_edge_fan(aes(alpha = weight), color = "black") +
  geom_node_point(size = 8, aes(color = is_source)) +
  geom_node_text(aes(label = name), color = "black") +
  ggtitle("Inferred Edge Probabilities - 10 Loci") +
  white_theme +
  # google_pres_theme +
  ggsave("~/Desktop/lab meeting presentation/small_inf.pdf", dpi = "retina")

#1A212B




inf_rc <- rc_from_weighted_network(inf_network_df)
true_rc <- rc_from_network(true_network_edges)
merged_rc <- inf_rc %>% left_join(true_rc, by = "from") %>% replace_na(list(rc.y = 0))
