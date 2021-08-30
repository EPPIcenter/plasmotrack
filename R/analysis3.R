library(tidyverse)
library(ggplot2)
# library(superheat)
library(igraph)
library(jsonlite)
library(ggraph)
library(tidygraph)
library(ggalt)

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


# setwd('~/Workspace/transmission_nets/test/outputs/ModelFiveTests_large2/CoreTest/')
# setwd('~/Workspace/transmission_nets/test/outputs/ModelFiveTests/CoreTest/')
# setwd('~/Workspace/transmission_nets/test/outputs/obs_exp/CoreTest90/')
setwd("~/Workspace/transmission_nets/test")
true_network <- extract_true_network_tibble("resources/JSON/locus_subset_exp/full_nodes.json")
true_network_vertices <- true_network$vertices
true_network_edges <- true_network$edges
all_nodes <- true_network_vertices %>% pull(name) %>% unique()
from_nodes <- true_network_edges %>% pull(from) %>% unique()
to_nodes <- true_network_edges %>% pull(to) %>% unique()
is_source_node <- !(all_nodes %in% to_nodes)
true_network_vertices <- true_network_vertices %>% mutate(is_source = is_source_node)
tidy_true_network <- tidygraph::as_tbl_graph(graph_from_data_frame(true_network_edges, vertices = true_network_vertices))
tidy_true_network <- tidy_true_network %>% activate(what = "nodes") %>% mutate(cluster = factor(igraph::components(tidy_true_network)$membership))
true_layout <- create_layout(tidy_true_network, layout = "nicely")

true_rc <- rc_from_network(true_network_edges)

# to_analyze = c("CoreTestFull", "CoreTest90", "CoreTest75", "CoreTest50", "CoreTest25", "CoreTest125")
to_analyze = c("CoreTestFull")

cluster_convex_hulls <- tibble(true_layout) %>%
  group_by(cluster) %>%
  filter(n() > 2) %>%
  summarize(cvx_hull = list(chull(x, y)), cvx_x=x[unlist(cvx_hull)], cvx_y=y[unlist(cvx_hull)]) %>%
  select(-cvx_hull)

rc <- data.frame(from=c(), rc=c(), obs_label=c())
for (path in to_analyze) {
  dir_path = paste0("outputs/locus_subset_exp/", path)
  llik <- read_csv(paste0(dir_path, '/parameters/likelihood.csv'))
  obs_nodes <- extract_true_network_tibble(paste0("resources/JSON/locus_subset_exp/nodes_", str_remove(path, "CoreTest"), ".json"))$vertices$name
  # geo_gen_prob <- read_csv(paste0(dir_path, '/parameters/geo_gen_prob.csv'))
  # mean_coi <- read_csv(paste0(dir_path, '/parameters/mean_coi.csv'))
  # loss_prob <- read_csv(paste0(dir_path, '/parameters/loss_prob.csv'))
  # mutation_prob <- read_csv(paste0(dir_path, '/parameters/mutation_prob.csv'))
  # orderings <- read_csv(paste0(dir_path, '/parameters/infection_order.csv', col_names = F, skip = 1))
  # infection_duration <- cbind(read_csv(paste0(dir_path, '/parameters/infection_duration_shape.csv')),
  #                             read_csv(paste0(dir_path, '/parameters/infection_duration_scale.csv')))

  parent_sets <- extract_parent_sets(paste0(dir_path, "/stats/"))

  inf_network_df <- network_from_parent_sets(parent_sets) %>%
    filter(weight > 0)

  is_source <- inf_network_df %>%
    filter(from == "S") %>%
    mutate(is_source = weight, is_source_dich = weight > .5) %>%
    select(to, is_source, is_source_dich) %>%
    rename(name = to)

  inf_network_vertices <- true_network_vertices %>%
    select(-is_source) %>%
    left_join(is_source, by = "name") %>%
    mutate(is_observed = name %in% obs_nodes)
  inf_network <- graph_from_data_frame(inf_network_df %>% filter(from != "S"), vertices = inf_network_vertices)
  tidy_inf_network <- tidygraph::as_tbl_graph(inf_network)

  # layouts <- generate_matched_layouts(tidy_true_network, tinf_graph_sub, layout = "nicely")
  # layouts <- generate_matched_layouts(tidy_network, tinf_graph_sub, layout = "auto")
  # true_layout <- layouts$target
  # tinf_layout <- layouts$other
  tinf_layout <- match_layout(tidy_inf_network, true_layout)

  ggraph(true_layout) +
    geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .2, spread = 0.005) +
    geom_edge_fan(color = "black", width = 1) +
    geom_node_point(aes(color = is_source), size = 6) +
    geom_node_text(aes(label = name), size = 3, color = "black") +
    ggtitle("True Network - 25 Loci") +
    guides(color = FALSE, fill = FALSE) +
    white_theme +
    # google_pres_theme +
    # ggsave("~/Desktop/lab meeting presentation/large_true.png", dpi = "retina")
    ggsave("~/Desktop/tnet_tests/25loci/true_network.png", dpi = "retina")

  if (!all(inf_network_vertices$is_observed == 1)) {
    ggraph(tinf_layout) +
      geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .2, expand = 0.005) +
      geom_edge_fan(aes(alpha = weight), color = "black", width = 1) +
      geom_node_point(size = 6, aes(color = is_source,  alpha = if_else(is_observed, 1, .25))) +
      geom_node_text(aes(label = name), size = 3, color = "grey") +
      ggtitle(paste0("Inferred Edge Probabilities - 25 Loci - ", path)) +
      white_theme +
      guides(color = FALSE, alpha = FALSE, fill = FALSE)
    # google_pres_theme +
    # ggsave("~/Desktop/lab meeting presentation/large_inf.png", dpi = "retina")
    ggsave(paste0("~/Desktop/tnet_tests/25loci/", path, ".png"), dpi = "retina")

  } else {
    ggraph(tinf_layout) +
      geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .2, expand = 0.005) +
      geom_edge_fan(aes(alpha = weight), color = "black", width = 1) +
      geom_node_point(size = 6, aes(color = is_source)) +
      geom_node_text(aes(label = name), size = 3, color = "grey") +
      ggtitle(paste0("Inferred Edge Probabilities - 25 Loci - ", path)) +
      white_theme +
      guides(color = FALSE, alpha = FALSE, fill = FALSE)
    # google_pres_theme +
    # ggsave("~/Desktop/lab meeting presentation/large_inf.png", dpi = "retina")
    ggsave(paste0("~/Desktop/tnet_tests/25loci/", path, ".png"), dpi = "retina")
  }
  inf_rc <- rc_from_weighted_network(inf_network_df)
  inf_rc$obs_label <- path
  rc <- rbind(rc, inf_rc)
}


ggraph(true_layout) +
  geom_edge_fan(color = "black", width = 1) +
  geom_node_point(size = 6, aes(color = is_source)) +
  geom_node_text(aes(label = name), size = 3, color = "black") +
  ggtitle("True Network - 25 Loci") +
  white_theme
# google_pres_theme +
# ggsave("~/Desktop/lab meeting presentation/large_true.png", dpi = "retina")
# ggsave("~/Desktop/tnet_tests/large_true.png", dpi = "retina")


ggraph(true_layout) +
  geom_edge_fan(color = "black", width = 1) +
  geom_node_point(size = 3, aes(color = is_source)) +
  geom_node_text(aes(label = name), size = 1, color = "black") +
  ggtitle("True Network - 25 Loci") +
  white_theme
# google_pres_theme +
# ggsave("~/Desktop/lab meeting presentation/large_true_dendro.png", dpi = "retina")
# ggsave("~/Desktop/tnet_tests/large_true.png", dpi = "retina")

ggraph(tinf_layout) +
  geom_edge_fan(aes(alpha = weight), color = "black", width = 1) +
  geom_node_point(size = 6, aes(color = is_source)) +
  geom_node_text(aes(label = name), size = 3, color = "grey") +
  ggtitle("Inferred Edge Probabilities - 25 Loci") +
  white_theme
# google_pres_theme +
# ggsave("~/Desktop/lab meeting presentation/large_inf.png", dpi = "retina")
# ggsave("~/Desktop/tnet_tests/large_inf.png", dpi = "retina")

ggraph(tinf_layout) +
  geom_edge_fan(aes(alpha = weight), color = "black", width = 1) +
  geom_node_point(size = 3, aes(color = is_source)) +
  geom_node_text(aes(label = name), size = 1, color = "grey") +
  ggtitle("Inferred Edge Probabilities - 25 Loci") +
  white_theme
# google_pres_theme +
# ggsave("~/Desktop/lab meeting presentation/large_inf_dendro.png", dpi = "retina")
# ggsave("~/Desktop/tnet_tests/large_true.png", dpi = "retina")

inf_rc <- rc_from_weighted_network(inf_network_df)
true_rc <- rc_from_network(true_network_edges)
merged_rc <- inf_rc %>% left_join(true_rc, by = "from") %>% replace_na(list(rc.y = 0))

genetic_dist_matrix <- calculate_genetic_dist_matrix(input_data$nodes, hamming_distance)

infection_durations_full_dat <- read_csv('outputs/locus_subset_exp/CoreTestFull/parameters/infection_duration.csv')

ggplot(infection_durations_full_dat, aes(x = iter, y = infection_duration, group = id, color = id)) + geom_line()

qplot(x = iter, y = infection_duration, data = infection_durations_full_dat %>% filter(id == 2), geom = "line",)

infection_durations <- infection_durations_full_dat %>%
  mutate(id = as.character(id)) %>%
  group_by(id) %>%
  summarize(mean = mean(infection_duration),
            upper = quantile(infection_duration, .975), lower = quantile(infection_duration, .025),
            min = min(infection_duration), max = max(infection_duration)) %>%
  left_join(true_network$vertices %>% select(name, infection_duration, infection_time, observation_time), by = c("id" = "name")) %>%
  mutate(inferred_infection_time = observation_time - mean) %>%
  mutate(ci_width = abs(lower - upper)) %>%
  arrange(infection_time) %>%
  mutate(true_ord = row_number()) %>%
  arrange(inferred_infection_time) %>%
  mutate(inferred_ord = row_number()) %>%
  arrange(observation_time) %>%
  mutate(naive_ord = row_number())

infection_ord <- infection_durations_full_dat %>%
  mutate(id = as.character(id)) %>%
  left_join(true_network$vertices %>% select(name, infection_time, observation_time), by = c("id" = "name")) %>%
  mutate(inferred_infection_time = observation_time - infection_duration) %>%
  group_by(iter) %>%
  arrange(inferred_infection_time) %>%
  mutate(inferred_ord = row_number()) %>%
  arrange(infection_time) %>%
  mutate(true_ord = row_number()) %>%
  arrange(observation_time) %>%
  mutate(naive_ord = row_number()) %>%
  ungroup() %>%
  group_by(id) %>%
  summarize(mean_ord = mean(inferred_ord), true_ord = mean(true_ord), naive_ord = mean(naive_ord))

ggplot(infection_ord, aes(x = true_ord)) +
  geom_point(aes(y = mean_ord), color = "blue") +
  geom_point(aes(y = naive_ord), color = "red") +
  geom_abline(slope = 1)


