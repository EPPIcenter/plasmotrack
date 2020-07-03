library(gtools)
library(jsonlite)
library(tidyverse)
library(rmutil)
library(igraph)

## Transmission Network Simulator
##


simulate_source <- function(num_loci, min_alleles, max_alleles) {
  num_alleles <- runif(num_loci, min_alleles, max_alleles)
  locus_allele_frequencies <- lapply(num_alleles, function(n) rdirichlet(1, rep(1,n)))
  names(locus_allele_frequencies) <- paste0("L", 1:num_loci)
  return(locus_allele_frequencies)
}

simulate_founder_infection <- function(source, coi, fpr, fnr) {
  total_nodes_ <<- total_nodes_ + 1 # global
  infection = list(id = as.character(total_nodes_), parent = "source")
  infection$latent_genotype <- lapply(source, function(s) as.numeric(rmultinom(1, coi, s) > 0))
  names(infection$latent_genotype) <- names(source)

  infection$obs_genotype <- lapply(infection$latent_genotype, function(l) simulate_observed_infection(l, fpr, fnr))
  names(infection$obs_genotype) <- names(infection$latent_genotype)

  return(infection)
}

simulate_child_infection <- function(parent, pr, assoc, fpr, fnr) {
  total_nodes_ <<- total_nodes_ + 1 # global
  infection = list(id = as.character(total_nodes_), parent = parent$id)
  infection$latent_genotype <- lapply(parent$latent_genotype, function(p) {
    total_retained <- rmultbinom(1, sum(p) - 1, pr, assoc) + 1
    retained <- sample(which(as.logical(p)), total_retained)
    g <- rep(0, length(p))
    g[retained] = 1
    return(g)
  })
  names(infection$latent_genotype) <- names(parent$latent_genotype)

  infection$obs_genotype <- lapply(infection$latent_genotype, function(l) simulate_observed_infection(l, fpr, fnr))
  names(infection$obs_genotype) <- names(infection$latent_genotype)

  return(infection)
}

simulate_observed_infection <- function(latent_genotype, fpr, fnr) {
  positives <- which(as.logical(latent_genotype))
  negatives <- which(!as.logical(latent_genotype))
  false_positives <- sample(negatives, rbinom(1, length(negatives), fpr))
  false_negatives <- sample(positives, rbinom(1, length(positives), fnr))
  obs_genotype <- latent_genotype
  obs_genotype[false_positives] = 1
  obs_genotype[false_negatives] = 0
  return(obs_genotype)
}

total_founders <- 20
num_loci <- 10
min_alleles <- 5
max_alleles <- 16
fpr = .05
fnr = .05
pr = .8
assoc = 1
total_nodes_ = 0

s1 <- simulate_source(num_loci, min_alleles, max_alleles)
nodes <- lapply(1:total_founders, function(x) simulate_founder_infection(s1, runif(1, 1, 5), fpr = fpr, fnr = fnr))

for(i in 1:20) {
  parents <- sample(nodes, 1)
  children <- lapply(parents, function(p) simulate_child_infection(p, pr, assoc, fpr, fnr))
  nodes <- append(nodes, children)
}

nodes_input_data <- lapply(nodes, function(n) {
  observed_genotype <- lapply(names(n$obs_genotype), function(loc) {
    list(locus=loc, genotype=paste0(n$obs_genotype[[loc]], collapse = ''))
  })
  id <- n$id
  res = list(observed_genotype = observed_genotype, id = id)
  return(res)
})

loci_data <- lapply(names(s1), function(x){
 res <- list(locus = x, num_alleles=length(s1[[x]]))
})

input_data <- list(loci=loci_data, nodes=nodes_input_data)

edge_list <- t(sapply(nodes, function(n) c(n$parent, n$id)))
tnet <- graph_from_edgelist(edge_list)
plot(delete_vertices(tnet, "source"), vertex.size=5)

write_json(input_data, "~/Workspace/transmission_nets/test/resources/JSON/nodes.json", auto_unbox=T)



