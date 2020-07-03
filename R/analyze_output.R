library(tidyverse)
library(ggplot2)
library(superheat)

does_node_precede <- function(nodes, positions, node_i, node_j) {
  precedes = positions[which(nodes == node_i)] < positions[which(nodes == node_j)]
  return(precedes)
}

setwd('~/Workspace/transmission_nets/test/outputs/ModelTwoTests/CoreTest/')
llik <- read_csv('likelihood.csv', col_names = c("llik"))
geo_gen_prob <- read_csv('geo_gen_prob.csv', col_names = c("geo_gen"))
mean_coi <- read_csv('mean_coi.csv', col_names = c("coi"))
loss_prob <- read_csv('loss_prob.csv', col_names = c("prob"))
mutation_prob <- read_csv('mutation_prob.csv', col_names = c("prob"))
orderings <- read_csv('infection_order.csv', col_names = F)
fpr <- read_csv('fpr.csv', col_names = c("fpr"))
fnr <- read_csv('fnr.csv', col_names = c("fnr"))

plot(1:length(llik$llik), llik$llik, type = "l")
qplot(geo_gen_prob$geo_gen)
qplot(mean_coi$coi)
qplot(mutation_prob$prob)
qplot(loss_prob$prob)
qplot(fpr$fpr)
qplot(fnr$fnr)

long_ord <- orderings %>%
  rowid_to_column() %>%
  pivot_longer(cols = -rowid, names_to = "position", values_to = "node") %>%
  mutate(position = as.numeric(str_remove(position, pattern = "X")))

nodes <- sort(unique(long_ord$node))

precedes_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes))

for (i in 1:(length(nodes) - 1)) {
  node_i = nodes[i]
  print(i)
  for (j in (i + 1):length(nodes)) {
    node_j = nodes[j]
    tmp <- long_ord %>% group_by(rowid) %>% summarise(i_before_j = does_node_precede(node, position, node_i, node_j))
    precedes_matrix[i, j] = mean(tmp$i_before_j)
    precedes_matrix[j, i] = 1 - precedes_matrix[i, j]
  }
}

superheat(precedes_matrix, heat.pal = c("red", "white", "blue"))



qplot(seq(0, 8, .1), dgamma(seq(0, 8, .1), 100, scale = .01, log = T))

