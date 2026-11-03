# pt_as_igraph() and pt_as_tbl_graph(): vertices are pt_infection_summary() rows, edges are hand-computed Edge
# marginals (marginal graph) or MAP Parent-set members (MAP graph), parent to child.

skip_if_not_installed("igraph")
skip_if_not_installed("tidygraph")

# Infections "10", "2", "3", "4" over 2 sampling Draws, observed at 1, 5, 3, 7:
#   "2":  Draw 1 {10} 0.5, {3} 0.5;  Draw 2 {10} 1   -> 10->2 0.75, 3->2 0.25; MAP {10}
#   "3":  Draw 1 {10} 0.6, {L} 0.4;  Draw 2 {L} 1    -> 10->3 0.3;             MAP {L}
#   "4":  {3, 10} 1 at both Draws                    -> 10->4 1, 3->4 1;       MAP {10, 3}
#   "10": {L} 1 at both Draws                        -> no edge;               MAP {L}
local_network_store <- function(.local_envir = parent.frame()) {
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  path <- fake_store(
    infections = c("10", "2", "3", "4"),
    parent_sets = rbind(
      sets(1, "2", "10", FALSE, 0.5), sets(1, "2", "3", FALSE, 0.5), sets(2, "2", "10", FALSE, 1),
      sets(1, "3", "10", FALSE, 0.6), sets(1, "3", "", TRUE, 0.4), sets(2, "3", "", TRUE, 1),
      sets(1:2, "4", "3,10", FALSE, 1),
      sets(1:2, "10", "", TRUE, 1)
    ),
    observation_time = c("10" = 1, "2" = 5, "3" = 3, "4" = 7),
    .local_envir = .local_envir
  )
  a <- pt_read(c(r = path))
  withr::defer(pt_close(a), envir = .local_envir)
  a
}

# The edges of an igraph as a data frame of from, to, marginal, lag, sorted by child then parent.
igraph_edges <- function(g) {
  e <- igraph::as_data_frame(g, what = "edges")
  e <- e[order(e$to, e$from, method = "radix"), c("from", "to", "marginal", "lag")]
  rownames(e) <- NULL
  e
}

test_that("pt_as_igraph() keeps marginal edges at or above floor, parent to child, with marginal and lag", {
  a <- local_network_store()
  g <- pt_as_igraph(a)
  expect_s3_class(g, "igraph")
  expect_true(igraph::is_directed(g))
  expect_equal(igraph_edges(g), data.frame(
    from = c("10", "3", "10", "10", "3"), to = c("2", "2", "3", "4", "4"),
    marginal = c(0.75, 0.25, 0.3, 1, 1), lag = c(4, 2, 2, 6, 4)
  ), tolerance = 1e-12)

  # 0.3 is kept (at the floor), 0.25 is dropped.
  g <- pt_as_igraph(a, floor = 0.3)
  expect_equal(igraph_edges(g), data.frame(
    from = c("10", "10", "10", "3"), to = c("2", "3", "4", "4"),
    marginal = c(0.75, 0.3, 1, 1), lag = c(4, 2, 6, 4)
  ), tolerance = 1e-12)
  # Every Infection stays a vertex, edge or not.
  expect_equal(igraph::V(g)$name, c("10", "2", "3", "4"))
})

test_that("pt_as_igraph() vertices are every Infection with pt_infection_summary() columns", {
  a <- local_network_store()
  g <- pt_as_igraph(a, floor = 0.9)
  s <- pt_infection_summary(a)
  v <- igraph::as_data_frame(g, what = "vertices")
  expect_equal(names(v), c("name", setdiff(names(s), "infection")))
  expect_equal(v$name, s$infection)
  for (column in setdiff(names(s), "infection")) {
    expect_equal(v[[column]], s[[column]], info = column)
  }
  expect_equal(igraph::V(g)$map_members, list(character(), "10", character(), c("10", "3")))
})

test_that("pt_as_igraph(graph = 'map') exports the MAP graph's observed-parent edges", {
  a <- local_network_store()
  g <- pt_as_igraph(a, graph = "map")
  expect_equal(igraph::V(g)$name, c("10", "2", "3", "4"))
  # The Latent parent is never an edge; floor does not apply to the MAP graph.
  expect_equal(igraph_edges(g), data.frame(
    from = c("10", "10", "3"), to = c("2", "4", "4"), marginal = c(0.75, 1, 1), lag = c(4, 6, 4)
  ), tolerance = 1e-12)
  expect_equal(igraph_edges(pt_as_igraph(a, floor = 1, graph = "map")), igraph_edges(g))
})

test_that("pt_as_tbl_graph() builds the same graph with from and to as parent and child", {
  a <- local_network_store()
  for (graph in c("marginal", "map")) {
    tg <- pt_as_tbl_graph(a, floor = 0.3, graph = graph)
    expect_s3_class(tg, "tbl_graph")
    nodes <- tibble::as_tibble(tg, active = "nodes")
    edges <- tibble::as_tibble(tg, active = "edges")
    expect_equal(nodes$name, c("10", "2", "3", "4"))
    expect_equal(nodes$map_members, list(character(), "10", character(), c("10", "3")))
    g <- pt_as_igraph(a, floor = 0.3, graph = graph)
    expect_equal(
      data.frame(from = nodes$name[edges$from], to = nodes$name[edges$to], marginal = edges$marginal, lag = edges$lag),
      igraph::as_data_frame(g, what = "edges")
    )
  }
})

test_that("network export pools Runs and checks its arguments", {
  a <- read_golden(c(one = golden_store("single-chain"), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  g <- pt_as_igraph(a, floor = 0.05)
  m <- pt_edge_marginals(a, floor = 0.05)
  expect_equal(igraph::vcount(g), nrow(pt_infections(a)))
  e <- igraph::as_data_frame(g, what = "edges")
  expect_equal(nrow(e), nrow(m))
  expect_equal(e$marginal, m$marginal[match(paste(e$from, e$to), paste(m$parent, m$child))])

  map <- pt_map_graph(a)
  expect_equal(igraph::ecount(pt_as_igraph(a, graph = "map")), sum(lengths(map$members)))

  expect_error(pt_as_igraph(a, graph = "tree"), "graph")
  expect_error(pt_as_igraph(a, floor = 2), "floor")
})
