# Hand-written simulator truth JSON, shaped like the thesis simulations' sim_*_s*.json.gz (keys network, nodes,
# full_nodes, loci).
#
#   truth_node(id, infection_time = 0, observation_time = 10, duration = 10, symptomatic = TRUE, genotype = list(),
#              observed = list())
#     One node. `genotype` is a list of strains, each a named list Locus -> 0/1 vector over the simulator's Alleles.
#     `observed` is the node's observed_genotype: a named list Locus -> 0/1 vector, written as a bitstring like the
#     simulator's ("0101"), in the simulator's numbering. Omitted when empty.
#   write_truth(edges, nodes, loci = list(), path = NULL)
#     `edges` is a list of c(src, dest); `loci` a named list Locus -> Allele count. Writes plain JSON to `path`
#     (default: a temp file removed when the caller exits) and returns the path.
#   gzip_copy(path)  a gzipped copy next to `path`, named `<path>.gz`.

truth_node <- function(id, infection_time = 0, observation_time = 10, duration = 10, symptomatic = TRUE,
                       genotype = list(), observed = list()) {
  strains <- lapply(genotype, function(strain) {
    lapply(names(strain), function(locus) list(locus = locus, genotype = I(as.integer(strain[[locus]]))))
  })
  node <- list(id = as.character(id), latent_genotype = strains, infection_duration = duration,
               infection_time = infection_time, observation_time = observation_time, symptomatic = symptomatic)
  if (length(observed)) {
    node$observed_genotype <- lapply(names(observed), function(locus) {
      list(locus = locus, genotype = paste(as.integer(observed[[locus]]), collapse = ""))
    })
  }
  node
}

write_truth <- function(edges, nodes, loci = list(), path = NULL, .local_envir = parent.frame()) {
  if (is.null(path)) path <- withr::local_tempfile(fileext = ".json", .local_envir = .local_envir)
  json <- list(
    network = lapply(edges, function(e) list(src = as.character(e[[1]]), dest = as.character(e[[2]]))),
    nodes = nodes,
    full_nodes = nodes,
    loci = lapply(names(loci), function(l) {
      list(locus = l, allele_freqs = I(rep(1 / loci[[l]], loci[[l]])), num_alleles = loci[[l]])
    })
  )
  jsonlite::write_json(json, path, auto_unbox = TRUE, digits = NA)
  path
}

gzip_copy <- function(path) {
  out <- paste0(path, ".gz")
  con <- gzfile(out, "wb")
  writeLines(readLines(path), con)
  close(con)
  out
}
