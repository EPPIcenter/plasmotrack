# Network export: the posterior network as an igraph or a tidygraph, vertices from pt_infection_summary(), edges
# from pt_edge_marginals() or pt_map_graph().

#' Export the posterior network
#'
#' Builds a directed graph of the posterior transmission network, for drawing or graph algorithms in igraph,
#' tidygraph or ggraph. The package exports networks; it does not draw them.
#'
#' * **Vertices** are every Infection, edge or not, in id string order. Their attributes are the columns of
#'   [pt_infection_summary()], with `infection` as the vertex `name` (`map_members` stays a list).
#' * **Edges** run parent to child (tidygraph's `from` and `to`) and carry the pooled Edge `marginal` and the `lag`
#'   (child's observation time minus the parent's, `NA` where either is missing), as in [pt_edge_marginals()].
#'   * `graph = "marginal"`: every edge whose pooled Edge marginal is at least `floor`.
#'   * `graph = "map"`: the MAP graph of [pt_map_graph()], one edge from each observed member of an Infection's MAP
#'     Parent set. `floor` does not apply.
#'
#' The Latent parent is never an edge: an Infection whose MAP Parent set is the Latent parent alone has no incoming
#' MAP edge. Values pool the sampling Draws of every Run; Runs with no sampling Draws are skipped with a message.
#'
#' igraph (and, for `pt_as_tbl_graph()`, tidygraph) must be installed; they are suggested, not imported.
#'
#' @inheritParams pt_edge_marginals
#' @param floor For `graph = "marginal"`, keep edges whose pooled marginal is at least `floor`.
#' @param graph `"marginal"` for the Edge-marginal graph, `"map"` for the MAP graph.
#' @return `pt_as_igraph()`: a directed `igraph`. `pt_as_tbl_graph()`: the same graph as a `tbl_graph`.
#' @seealso [pt_infection_summary()], [pt_edge_marginals()], [pt_map_graph()].
#' @export
#' @examplesIf rlang::is_installed("igraph")
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' g <- pt_as_igraph(a, floor = 0.05)
#' g
#' pt_as_igraph(a, graph = "map")
#' pt_close(a)
pt_as_igraph <- function(a, floor = 0.05, graph = c("marginal", "map"), runs = NULL) {
  rlang::check_installed("igraph", reason = "to export the network.")
  network <- network_tables(a, floor, rlang::arg_match(graph), runs)
  igraph::graph_from_data_frame(network$edges, directed = TRUE, vertices = network$vertices)
}

#' @rdname pt_as_igraph
#' @export
#' @examplesIf rlang::is_installed(c("igraph", "tidygraph"))
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_as_tbl_graph(a, floor = 0.05)
#' pt_close(a)
pt_as_tbl_graph <- function(a, floor = 0.05, graph = c("marginal", "map"), runs = NULL) {
  rlang::check_installed(c("igraph", "tidygraph"), reason = "to export the network.")
  tidygraph::as_tbl_graph(pt_as_igraph(a, floor, rlang::arg_match(graph), runs))
}

# Vertex and edge data frames of the network: vertices are pt_infection_summary() with `infection` renamed `name`;
# edges are from (parent), to (child), marginal, lag, sorted by child then parent.
network_tables <- function(a, floor, graph, runs, call = rlang::caller_env()) {
  if (!is.numeric(floor) || length(floor) != 1 || is.na(floor) || floor < 0 || floor > 1) {
    cli::cli_abort("{.arg floor} must be a number between 0 and 1.", call = call)
  }
  runs <- runs_with_draws(a, runs)
  vertices <- pt_infection_summary(a, runs = runs)
  names(vertices)[names(vertices) == "infection"] <- "name"

  marginals <- pt_edge_marginals(a, floor = if (graph == "marginal") floor else 0, runs = runs)
  if (graph == "map") {
    map <- pt_map_graph(a, runs = runs)
    child <- rep(map$infection, lengths(map$members))
    parent <- as.character(unlist(map$members))
    o <- order(child, parent, method = "radix")
    at <- match(paste(child[o], parent[o], sep = "\r"), paste(marginals$child, marginals$parent, sep = "\r"))
    marginals <- marginals[at, ]
  }
  edges <- data.frame(
    from = marginals$parent, to = marginals$child, marginal = marginals$marginal, lag = marginals$lag
  )
  list(vertices = as.data.frame(vertices), edges = edges)
}
