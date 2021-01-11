library(gtools)
library(jsonlite)
library(tidyverse)
library(rmutil)
library(igraph)

## Transmission Network Simulator
##

simulate_source <- function(num_loci, min_alleles, max_alleles) {
  num_alleles <- runif(num_loci, min_alleles, max_alleles)
  locus_allele_frequencies <-
    map(num_alleles, function(n)
      rdirichlet(1, rep(1, n)))
  names(locus_allele_frequencies) <- paste0("L", 1:num_loci)
  return(locus_allele_frequencies)
}

simulate_founder_infection <- function(id, source, coi, fpr, fnr) {
  infection = list(id = as.character(id), parent = "source")
  infection$latent_genotype <-
    map(source, function(s)
      as.numeric(rmultinom(1, coi, s) > 0))
  names(infection$latent_genotype) <- names(source)

  infection$observed_genotype <-
    map(infection$latent_genotype, function(l)
      simulate_observed_infection(l, fpr, fnr))
  names(infection$observed_genotype) <-
    names(infection$latent_genotype)

  return(infection)
}

simulate_child_infection <-
  function(id, parent, pr, assoc, fpr, fnr) {
    infection = list(id = as.character(id), parent = parent$id)
    infection$latent_genotype <-
      map(parent$latent_genotype, function(p) {
        total_retained <- rmultbinom(1, sum(p) - 1, pr, assoc) + 1
        retained <- sample(which(as.logical(p)), total_retained)
        g <- rep(0, length(p))
        g[retained] = 1
        return(g)
      })
    names(infection$latent_genotype) <- names(parent$latent_genotype)

    infection$observed_genotype <-
      map(infection$latent_genotype, function(l)
        simulate_observed_infection(l, fpr, fnr))
    names(infection$observed_genotype) <-
      names(infection$latent_genotype)

    return(infection)
  }

simulate_observed_infection <- function(latent_genotype, fpr, fnr) {
  positives <- which(as.logical(latent_genotype))
  negatives <- which(!as.logical(latent_genotype))
  false_positives <-
    sample(negatives, rbinom(1, length(negatives), fpr))
  false_negatives <-
    sample(positives, rbinom(1, length(positives), fnr))
  obs_genotype <- latent_genotype
  obs_genotype[false_positives] = 1
  obs_genotype[false_negatives] = 0
  return(obs_genotype)
}

# fully observed network
simulate_network <- function(source_allele_freqs, mean_founder_coi, fpr, fnr, pr, assoc, num_founders, total_nodes) {
    founder_node_args_list <- list(coi = rpois(num_founders, mean_founder_coi), id = 1:num_founders)
    founder_nodes <- pmap(founder_node_args_list, function(id, coi) {
      simulate_founder_infection(id, source_allele_freqs, coi, fpr, fnr)
    })
    nodes = append(list(), founder_nodes)
    for (i in (num_founders + 1):total_nodes) {
      parent <- sample(nodes, 1)[[1]]
      child <- list(simulate_child_infection(i, parent, pr, assoc, fpr, fnr))
      nodes <- append(nodes, child)
    }
    return(nodes)
  }

genotype = paste0(n$obs_genotype[[loc]], collapse = '')

generate_data <-
  function(source_allele_freqs,
           mean_founder_coi,
           fpr,
           fnr,
           pr,
           assoc,
           num_founders,
           total_nodes) {
    simulated_nodes <- simulate_network(
      source_allele_freqs,
      mean_founder_coi,
      fpr,
      fnr,
      pr,
      assoc,
      num_founders,
      total_nodes
    )
    loci <- pmap(list(
      locus = names(source_allele_freqs),
      allele_freqs = source_allele_freqs
    ),
    function(locus, allele_freqs) {
      list(
        locus = locus,
        allele_freqs = allele_freqs[1,],
        num_alleles = length(allele_freqs)
      )
    })
    nodes <- map(simulated_nodes, function(node) {
      list(
        id = node$id,
        latent_genotype = pmap(list(
          locus = names(node$latent_genotype),
          genotype = node$latent_genotype
        ),
        function(locus, genotype)
          list(
            locus = locus,
            genotype = paste0(genotype, collapse = '')
          )),
        observed_genotype = pmap(list(
          locus = names(node$observed_genotype),
          genotype = node$observed_genotype
        ),
        function(locus, genotype)
          list(
            locus = locus,
            genotype = paste0(genotype, collapse = '')
          ))
      )
    })

    network <- map(simulated_nodes, function(node) {
      list(from = node$parent, to = node$id)
    })

    return(
      list(
        loci = loci,
        nodes = nodes,
        network = network,
        mean_founder_coi = mean_founder_coi,
        fpr = fpr,
        fnr = fnr,
        pr = pr,
        assoc = assoc
      )
    )
  }


# Schema:
# {
#   loci: [
#     {
#       locus: locus_name,
#       allele_freqs: [allele1, allele2...],
#       num_alleles: int
#     }
#   ],
#   nodes: [
#     {
#       id: node_name,
#       latent_genotype: [
#         {
#           locus: locus_name,
#           genotype: genotype_string
#         }
#       ],
#       observed_genotype: [
#         {
#           locus: locus_name,
#           genotype: genotype_string
#         }
#       ]
#     }
#   ],
#   network: [
#     {
#       from: node_id,
#       to: node_id
#     }
#   ]
# }



num_loci <- 20
min_alleles <- 6
max_alleles <- 24
s1 <- simulate_source(num_loci, min_alleles, max_alleles)

mean_founder_coi <- 4
fpr = .005
fnr = .05
pr = .7
assoc = 1
num_founders <- 50
total_nodes = 100

dat <-
  generate_data(s1,
                mean_founder_coi,
                fpr,
                fnr,
                pr,
                assoc,
                num_founders,
                total_nodes)
write_json(
  dat,
  "~/Workspace/transmission_nets/test/resources/JSON/nodes3.json",
  auto_unbox = T
)

edge_data_frame <-
  map_dfr(dat$network, function(edge)
    data.frame(from = edge$from, to = edge$to))
tnet <- graph_from_data_frame(edge_data_frame)
tnet <- delete_vertices(tnet, "source")
tnet <- set_vertex_attr(tnet, "size", value = 8)
plot(tnet)


nodes_input_data <- lapply(nodes, function(n) {
  observed_genotype <- lapply(names(n$obs_genotype), function(loc) {
    list(locus = loc,
         genotype = paste0(n$obs_genotype[[loc]], collapse = ''))
  })
  id <- n$id
  res = list(observed_genotype = observed_genotype, id = id)
  return(res)
})

loci_data <- lapply(names(s1), function(x) {
  res <- list(locus = x, num_alleles = length(s1[[x]]))
})

input_data <- list(loci = loci_data, nodes = nodes_input_data)

transmission_edge_list <-
  t(sapply(nodes, function(n)
    c(n$parent, n$id)))
obs_edge_list <-
  t(sapply(nodes, function(n)
    c(n$id, paste0(n$id, "obs"))))
edge_list = rbind(transmission_edge_list, obs_edge_list)

# tnet <- graph_from_edgelist(edge_list)
tnet <- graph_from_edgelist(transmission_edge_list)
tnet <- delete_vertices(tnet, "source")
plot(tnet)

vertices <- V(tnet)
edges <- E(tnet)
latent_vertices <- vertices[1:105]
obs_vertices <- vertices[106:length(vertices)]

random_latent_vertices <-
  sample(latent_vertices, length(latent_vertices) * .5)

colored_tnet <-
  set_vertex_attr(tnet, "color", index = latent_vertices, "#AFD2E9")
colored_tnet <-
  set_vertex_attr(colored_tnet, "color", index = obs_vertices, "#C98686")

png(
  filename = "~/Google Drive//Gradschool/Presentations/model-it meeting/latent_network.png",
  width = 2000,
  height = 2000,
  bg = NA
)
plot(
  delete_vertices(colored_tnet, obs_vertices),
  vertex.size = 5,
  vertex.label = NA,
  edge.arrow.size = 2
)
dev.off()

png(
  filename = "~/Google Drive//Gradschool/Presentations/model-it meeting/obs_network.png",
  width = 2000,
  height = 2000,
  bg = NA
)
plot(
  colored_tnet,
  vertex.size = 5,
  vertex.label = NA,
  edge.arrow.size = 2
)
dev.off()

write_json(
  input_data,
  "~/Workspace/transmission_nets/test/resources/JSON/nodes.json",
  auto_unbox = T
)

input_data2 <-
  read_json("~/Workspace/transmission_nets/test/resources/JSON/nodes2.json")
