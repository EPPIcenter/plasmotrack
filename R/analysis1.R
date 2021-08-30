library(tidyverse)
library(ggplot2)
# library(superheat)
library(igraph)
library(jsonlite)
library(ggraph)
library(tidygraph)
library(ggalt)

google_pres_theme <- theme(panel.background = element_rect(fill = "#1A212C", color = "#1A212C"),
                           plot.background = element_rect(fill = "#1A212C", color = "#1A212C"),
                           plot.title = element_text(color = "white", size = 24),
                           # axis.title = element_text(color = "white"),
                           # axis.text = element_text(color = "white"),
                           # axis.line = element_line(color = "white"),
                           panel.border = element_blank(),
                           panel.grid = element_blank(),
                           legend.background = element_rect(fill = "#1A212C", color = "#1A212C"),
                           legend.key = element_rect(fill = "#1A212C", color = "#1A212C"),
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
true_network <- extract_true_network_tibble("resources/JSON/obs_exp/nodes_full.json")
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

to_analyze = c("CoreTestFull", "CoreTest90", "CoreTest75", "CoreTest50", "CoreTest25", "CoreTest125")
titles <- c("Fully Observed", "90% Observed", "75% Observed", "50% Observed", "25% Observed", "12.5% Observed")

cluster_convex_hulls <- tibble(true_layout) %>%
  group_by(cluster) %>%
  filter(n() > 2) %>%
  summarize(cvx_hull = list(chull(x, y)), cvx_x=x[unlist(cvx_hull)], cvx_y=y[unlist(cvx_hull)]) %>%
  select(-cvx_hull)

ggraph(true_layout) +
  geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .4, expand = 0.005, show.legend = FALSE) +
  geom_edge_fan(color = "white", width = 1, show.legend = FALSE) +
  geom_node_point(color = "white", size = 6) +
  geom_node_point(aes(color = is_source), size = 5, show.legend = FALSE) +
  geom_node_text(aes(label = name), size = 2.5, color = "black") +
  ggtitle("True Network") +
  # white_theme +
  google_pres_theme +
  # ggsave("~/Desktop/lab meeting presentation/large_true.png", dpi = "retina")
  ggsave("~/Desktop/tnet_tests/true_network.png", dpi = "retina", scale = 1.25)

rc <- data.frame(from=c(), rc=c(), obs_label=c())
for (i in 1:length(to_analyze)) {
  path = to_analyze[i]
  title = titles[i]
  dir_path = paste0("outputs/obs_exp/", path)
  llik <- read_csv(paste0(dir_path, '/parameters/likelihood.csv')) %>% mutate(iter = row_number())
  # ggplot(llik, aes(x = iter, y = llik)) +
  #   geom_line(color = "white") +
  #   xlab("Iter.") +
  #   ylab("LLik.") +
  #   theme_minimal() +
  #   google_pres_theme +
  #   ggsave(paste0("~/Desktop/tnet_tests/100loci/llik_", path, ".png"), dpi = "retina", scale = 1.25)
  obs_nodes <- extract_true_network_tibble(paste0("resources/JSON/obs_exp/nodes_", str_remove(path, "CoreTest"), ".json"))$vertices$name
  # geo_gen_prob <- read_csv(paste0(dir_path, '/parameters/geo_gen_prob.csv'))
  # mean_coi <- read_csv(paste0(dir_path, '/parameters/mean_coi.csv'))
  # loss_prob <- read_csv(paste0(dir_path, '/parameters/loss_prob.csv'))
  # mutation_prob <- read_csv(paste0(dir_path, '/parameters/mutation_prob.csv'))
  # orderings <- read_csv(paste0(dir_path, '/parameters/infection_order.csv', col_names = F, skip = 1))
  infection_duration <- cbind(read_csv(paste0(dir_path, '/parameters/infection_duration_shape.csv')),
                              read_csv(paste0(dir_path, '/parameters/infection_duration_scale.csv')))

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

  tinf_layout <- match_layout(tidy_inf_network, true_layout)

  if (!all(inf_network_vertices$is_observed == 1)) {
    ggraph(tinf_layout) +
      geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .4, expand = 0.005, show.legend = FALSE) +
      geom_edge_fan(aes(alpha = weight), color = "white", width = 1, show.legend = FALSE) +
      geom_node_point(aes(alpha = if_else(is_observed, 1, .5)), size = 6, color = "white", show.legend = FALSE) +
      geom_node_point(aes(color = is_source, alpha = if_else(is_observed, 1, .5)), size = 5, show.legend = FALSE) +
      geom_node_text(aes(label = name), size = 2.5, color = "black") +
      ggtitle(paste0("Inferred Edge Probabilities - 100 Loci - ", title)) +
      # white_theme +
      guides(color = FALSE, alpha = FALSE, fill = FALSE) +
      google_pres_theme +
      ggsave(paste0("~/Desktop/tnet_tests/100loci/", path, ".png"), dpi = "retina", scale = 1.25)
  } else {
   ggraph(tinf_layout) +
    geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .4, expand = 0.005, show.legend = FALSE) +
    geom_edge_fan(aes(alpha = weight), color = "white", width = 1, show.legend = FALSE) +
    geom_node_point(size = 6, color = "white") +
    geom_node_point(size = 5, aes(color = is_source), show.legend = FALSE) +
    geom_node_text(aes(label = name), size = 2.5, color = "black") +
    ggtitle(paste0("Inferred Edge Probabilities - 100 Loci - ", title)) +
    # white_theme +
    google_pres_theme +
    ggsave(paste0("~/Desktop/tnet_tests/100loci/", path, ".png"), dpi = "retina", scale = 1.25)
  }
  inf_rc <- rc_from_weighted_network(inf_network_df)
  inf_rc$obs_label <- path
  rc <- rbind(rc, inf_rc)
}


rc2 <- data.frame(from=c(), rc=c(), obs_label=c())
for (i in 1:length(to_analyze)) {
  path = to_analyze[i]
  title = titles[i]
  dir_path = paste0("outputs/locus_subset_exp/", path)
  llik <- read_csv(paste0(dir_path, '/parameters/likelihood.csv')) %>% mutate(iter = row_number())
  obs_nodes <- extract_true_network_tibble(paste0("resources/JSON/locus_subset_exp/nodes_", str_remove(path, "CoreTest"), ".json"))$vertices$name
  # geo_gen_prob <- read_csv(paste0(dir_path, '/parameters/geo_gen_prob.csv'))
  # mean_coi <- read_csv(paste0(dir_path, '/parameters/mean_coi.csv'))
  # loss_prob <- read_csv(paste0(dir_path, '/parameters/loss_prob.csv'))
  # mutation_prob <- read_csv(paste0(dir_path, '/parameters/mutation_prob.csv'))
  # orderings <- read_csv(paste0(dir_path, '/parameters/infection_order.csv', col_names = F, skip = 1))
  # infection_duration <- cbind(read_csv(paste0(dir_path, '/parameters/infection_duration_shape.csv')),
  #                             read_csv(paste0(dir_path, '/parameters/infection_duration_scale.csv')))

  # ggplot(llik, aes(x = iter, y = llik)) +
  #   geom_line(color = "white") +
  #   xlab("Iter.") +
  #   ylab("LLik.") +
  #   theme_minimal() +
  #   google_pres_theme +
  #   ggsave(paste0("~/Desktop/tnet_tests/25loci/llik_", path, ".png"), dpi = "retina", scale = 1.25)
  parent_sets <- extract_parent_sets(paste0(dir_path, "/stats/"))
  #
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

  tinf_layout <- match_layout(tidy_inf_network, true_layout)

  if (!all(inf_network_vertices$is_observed == 1)) {
    ggraph(tinf_layout) +
      geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .4, expand = 0.005, show.legend = FALSE) +
      geom_edge_fan(aes(alpha = weight), color = "white", width = 1, show.legend = FALSE) +
      geom_node_point(aes(alpha = if_else(is_observed, 1, .5)), size = 6, color = "white", show.legend = FALSE) +
      geom_node_point(aes(color = is_source, alpha = if_else(is_observed, 1, .5)), size = 5, show.legend = FALSE) +
      geom_node_text(aes(label = name), size = 2.5, color = "black") +
      ggtitle(paste0("Inferred Edge Probabilities - 25 Loci - ", title)) +
      # white_theme +
      guides(color = FALSE, alpha = FALSE, fill = FALSE) +
      google_pres_theme +
      ggsave(paste0("~/Desktop/tnet_tests/25loci/", path, ".png"), dpi = "retina", scale = 1.25)
  } else {
    ggraph(tinf_layout) +
      geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .4, expand = 0.005, show.legend = FALSE) +
      geom_edge_fan(aes(alpha = weight), color = "white", width = 1, show.legend = FALSE) +
      geom_node_point(size = 6, color = "white") +
      geom_node_point(size = 5, aes(color = is_source), show.legend = FALSE) +
      geom_node_text(aes(label = name), size = 2.5, color = "black") +
      ggtitle(paste0("Inferred Edge Probabilities - 25 Loci - ", title)) +
      # white_theme +
      google_pres_theme +
      ggsave(paste0("~/Desktop/tnet_tests/25loci/", path, ".png"), dpi = "retina", scale = 1.25)
  }
  inf_rc <- rc_from_weighted_network(inf_network_df)
  inf_rc$obs_label <- path
  rc2 <- rbind(rc2, inf_rc)
}




ggraph(true_layout) +
  geom_edge_fan(color = "black", width = 1) +
  geom_node_point(size = 6, aes(color = is_source)) +
  geom_node_text(aes(label = name), size = 3, color = "black") +
  ggtitle("True Network - 100 Loci") +
  white_theme
  # google_pres_theme +
  # ggsave("~/Desktop/lab meeting presentation/large_true.png", dpi = "retina")
  # ggsave("~/Desktop/tnet_tests/large_true.png", dpi = "retina")


ggraph(true_layout) +
  geom_edge_fan(color = "black", width = 1) +
  geom_node_point(size = 3, aes(color = is_source)) +
  geom_node_text(aes(label = name), size = 1, color = "black") +
  ggtitle("True Network - 100 Loci") +
  white_theme
  # google_pres_theme +
  # ggsave("~/Desktop/lab meeting presentation/large_true_dendro.png", dpi = "retina")
  # ggsave("~/Desktop/tnet_tests/large_true.png", dpi = "retina")

ggraph(tinf_layout) +
  geom_edge_fan(aes(alpha = weight), color = "black", width = 1) +
  geom_node_point(size = 6, aes(color = is_source)) +
  geom_node_text(aes(label = name), size = 3, color = "grey") +
  ggtitle("Inferred Edge Probabilities - 100 Loci") +
  white_theme
  # google_pres_theme +
  # ggsave("~/Desktop/lab meeting presentation/large_inf.png", dpi = "retina")
  # ggsave("~/Desktop/tnet_tests/large_inf.png", dpi = "retina")

ggraph(tinf_layout) +
  geom_edge_fan(aes(alpha = weight), color = "black", width = 1) +
  geom_node_point(size = 3, aes(color = is_source)) +
  geom_node_text(aes(label = name), size = 1, color = "grey") +
  ggtitle("Inferred Edge Probabilities - 100 Loci") +
  white_theme
  # google_pres_theme +
  # ggsave("~/Desktop/lab meeting presentation/large_inf_dendro.png", dpi = "retina")
  # ggsave("~/Desktop/tnet_tests/large_true.png", dpi = "retina")

inf_rc <- rc_from_weighted_network(inf_network_df)
true_rc <- rc_from_network(true_network_edges)
merged_rc <- inf_rc %>% left_join(true_rc, by = "from") %>% replace_na(list(rc.y = 0))


true_network_json <- read_json("resources/JSON/obs_exp/nodes_full.json")
genetic_dist_matrix <- calculate_genetic_dist_matrix(true_network_json$nodes, hamming_distance)

subset_network_json <- read_json("resources/JSON/locus_subset_exp/nodes_full.json")
sub_genetic_dist_matrix <- calculate_genetic_dist_matrix(subset_network_json$nodes, hamming_distance)

obs_dist <- as.data.frame(as.table(genetic_dist_matrix$observed_distance_matrix)) %>%
  rename(to = Var1, from = Var2, dist = Freq)

sub_obs_dist <- as.data.frame(as.table(sub_genetic_dist_matrix$observed_distance_matrix)) %>%
  rename(to = Var1, from = Var2, dist = Freq)

obs_dist_network <- obs_dist %>%
  filter(dist > 0) %>%
  group_by(to) %>%
  summarize(from = from[which.min(dist)])

sub_obs_dist_network <- sub_obs_dist %>%
  filter(dist > 0) %>%
  group_by(to) %>%
  summarize(from = from[which.min(dist)])

tidy_obs_dist_network <- tidygraph::as_tbl_graph(graph_from_data_frame(obs_dist_network))
tidy_obs_layout <- match_layout(tidy_obs_dist_network, true_layout)

tidy_sub_obs_dist_network <- tidygraph::as_tbl_graph(graph_from_data_frame(sub_obs_dist_network))
tidy_sub_obs_layout <- match_layout(tidy_sub_obs_dist_network, true_layout)

ggraph(tidy_obs_layout) +
  geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .4, expand = 0.005, show.legend = FALSE) +
  geom_edge_fan(color = "white", width = 1, show.legend = FALSE) +
  geom_node_point(color = "white", size = 6) +
  geom_node_text(aes(label = name), size = 2.5, color = "black") +
  ggtitle("Naive Network - 100 Loci") +
  # white_theme +
  google_pres_theme +
  ggsave("~/Desktop/tnet_tests/100loci/naive.png", dpi = "retina", scale = 1.25)


ggraph(tidy_sub_obs_layout) +
  geom_encircle(data = cluster_convex_hulls, aes(x = cvx_x, y = cvx_y, fill = cluster), alpha = .4, expand = 0.005, show.legend = FALSE) +
  geom_edge_fan(color = "white", width = 1, show.legend = FALSE) +
  geom_node_point(color = "white", size = 6) +
  geom_node_text(aes(label = name), size = 2.5, color = "black") +
  ggtitle("Naive Network - 25 Loci") +
  # white_theme +
  google_pres_theme +
  ggsave("~/Desktop/tnet_tests/25loci/naive.png", dpi = "retina", scale = 1.25)

