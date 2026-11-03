# Flows between strata (#127): attribution, enrichment, net flow, long-lag reassignment.

# One Run, Draws 1 and 2, Infections 1-4. Strata: x = {1, 2}, y = {3}, NA = {4} (missing from epi).
# Times 1 = 0, 2 = 2, 3 = 10, 4 = NA, so lags: 2<-1 is 2, 3<-1 is 10, 3<-2 is 8, 4<-3 is NA.
#   Draw 1: 1 {L} 1; 2 {1} 0.6, {L} 0.4; 3 {1, L} 0.7, {2} 0.3; 4 {3} 1.
#   Draw 2: 1 {L} 1; 2 {1} 1;            3 {L} 0.5, {1, 2} 0.5; 4 {L} 1.
local_flow_store <- function(path = NULL, .local_envir = parent.frame()) {
  if (is.null(path)) path <- file.path(withr::local_tempdir(.local_envir = .local_envir), "flow")
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  fake_store(
    path,
    infections = c("1", "2", "3", "4"),
    parent_sets = rbind(
      sets(1, "1", "", TRUE, 1), sets(1, "2", "1", FALSE, 0.6), sets(1, "2", "", TRUE, 0.4),
      sets(1, "3", "1", TRUE, 0.7), sets(1, "3", "2", FALSE, 0.3), sets(1, "4", "3", FALSE, 1),
      sets(2, "1", "", TRUE, 1), sets(2, "2", "1", FALSE, 1), sets(2, "3", "", TRUE, 0.5),
      sets(2, "3", "1,2", FALSE, 0.5), sets(2, "4", "", TRUE, 1)
    ),
    observation_time = c("1" = 0, "2" = 2, "3" = 10),
    allowed_parents = list("2" = "1", "3" = c("1", "2"), "4" = c("3", "2"))
  )
}

flow_epi <- function() data.frame(infection = c("1", "2", "3"), district = c("x", "x", "y"))

read_flow_store <- function(path = local_flow_store(.local_envir = parent.frame()), epi = flow_epi()) {
  pt_read(path, epi = epi)
}

sorted <- function(x, keys) {
  x <- as.data.frame(dplyr::collect(x))
  x <- x[do.call(order, c(unname(as.list(x[keys])), method = "radix")), , drop = FALSE]
  rownames(x) <- NULL
  x
}

test_that("pt_attribution splits each recipient stratum's parent mass by source, Latent parent last", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  x <- pt_attribution(a, by = "district")
  expect_s3_class(x, "tbl_lazy")
  x <- sorted(x, c("run", "iteration", "recipient", "latent", "source"))
  expect_equal(names(x), c("run", "iteration", "recipient", "source", "latent", "share"))
  one <- x[x$iteration == 1, ]
  # Draw 1. x = {1, 2}: 2<-1 0.6 from x; Latent parent 1 + 0.4. y = {3}: 0.7 + 0.3 from x; Latent parent 0.7.
  # NA = {4}: 1 from y.
  expect_equal(one$recipient, rep(c("x", "y", NA), each = 4))
  expect_equal(one$source, rep(c("x", "y", NA, NA), 3))
  expect_equal(one$latent, rep(c(FALSE, FALSE, FALSE, TRUE), 3))
  expect_equal(one$share, c(0.3, 0, 0, 0.7, 1 / 1.7, 0, 0, 0.7 / 1.7, 0, 1, 0, 0), tolerance = 1e-12)
  # Draw 2. x: 2<-1 1, Latent parent 1. y: 0.5 + 0.5 from x, Latent parent 0.5. NA: Latent parent 1.
  two <- x[x$iteration == 2, ]
  expect_equal(two$share, c(0.5, 0, 0, 0.5, 1 / 1.5, 0, 0, 0.5 / 1.5, 0, 0, 0, 1), tolerance = 1e-12)
  sums <- tapply(x$share, paste(x$iteration, x$recipient), sum)
  expect_equal(as.vector(sums), rep(1, 6), tolerance = 1e-12)
})

test_that("pt_enrichment divides each cell's edge share by its candidate share, NA without candidates", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  x <- pt_enrichment(a, by = "district")
  expect_s3_class(x, "tbl_lazy")
  x <- sorted(x, c("run", "iteration", "source", "recipient"))
  expect_equal(names(x), c("run", "iteration", "source", "recipient", "enrichment"))
  # Candidates 2<-1 (x, x), 3<-1 and 3<-2 (x, y), 4<-3 (y, NA), 4<-2 (x, NA): shares 0.2, 0.4, 0.2, 0.2.
  # Draw 1 edge mass: (x, x) 0.6, (x, y) 1, (y, NA) 1, of 2.6. Draw 2: (x, x) 1, (x, y) 1, of 2.
  expect_equal(x$source, rep(rep(c("x", "y", NA), each = 3), 2))
  expect_equal(x$recipient, rep(c("x", "y", NA), 6))
  expect_equal(x$enrichment, c(0.6 / 2.6 / 0.2, 1 / 2.6 / 0.4, 0, NA, NA, 1 / 2.6 / 0.2, NA, NA, NA,
                               2.5, 1.25, 0, NA, NA, 0, NA, NA, NA), tolerance = 1e-12)
  # Candidates given by hand; unknown ids and repeats are dropped.
  candidates <- data.frame(parent = c("1", "1", "9", "3"), child = c("2", "2", "3", "4"))
  y <- sorted(pt_enrichment(a, by = "district", candidates = candidates), c("run", "iteration", "source", "recipient"))
  expect_equal(y$enrichment[y$iteration == 1], c(0.6 / 2.6 / 0.5, NA, NA, NA, NA, 1 / 2.6 / 0.5, NA, NA, NA),
               tolerance = 1e-12)
  expect_error(pt_enrichment(a, by = "district", candidates = data.frame(from = "1")), "parent")
})

test_that("pt_net_flow counts expected out-edges minus in-edges between different strata", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  x <- pt_net_flow(a, by = "district")
  expect_s3_class(x, "tbl_lazy")
  x <- sorted(x, c("run", "iteration", "stratum"))
  expect_equal(names(x), c("run", "iteration", "stratum", "net_flow"))
  # Draw 1: x -> y 0.7 + 0.3, y -> NA 1 (2<-1 stays within x). Draw 2: x -> y 1.
  expect_equal(x$stratum, rep(c("x", "y", NA), 2))
  expect_equal(x$net_flow, c(1, 0, -1, 1, -1, 0), tolerance = 1e-12)
})

test_that("pt_reassign_long_lags moves edges beyond the lag to the Latent parent, clamping at 1", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  x <- pt_reassign_long_lags(a, max_lag = 5)
  expect_s3_class(x, "tbl_lazy")
  x <- sorted(x, c("run", "iteration"))
  expect_equal(names(x), c("run", "iteration", "mean_out_degree", "mean_p_latent_parent"))
  # Lags 3<-1 (10) and 3<-2 (8) move; 2<-1 (2) stays; 4<-3 (NA) never moves.
  # Draw 1 out-degrees (1, 2, 3, 4): 0.6, 0, 1, 0. Latent parent: 1, 0.4, min(0.7 + 1, 1) (set {1, L} holds a
  # moved parent and the Latent parent), 0. Draw 2: out 1, 0, 0, 0; Latent parent 1, 0, min(0.5 + 1, 1), 1.
  expect_equal(x$mean_out_degree, c(1.6 / 4, 0.25), tolerance = 1e-12)
  expect_equal(x$mean_p_latent_parent, c(2.4 / 4, 0.75), tolerance = 1e-12)

  # Per stratum: out-degree over the stratum's sources, Latent parent over its recipients.
  y <- sorted(pt_reassign_long_lags(a, max_lag = 5, by = "district"), c("run", "iteration", "stratum"))
  expect_equal(names(y), c("run", "iteration", "stratum", "mean_out_degree", "mean_p_latent_parent"))
  expect_equal(y$stratum, rep(c("x", "y", NA), 2))
  expect_equal(y$mean_out_degree, c(0.3, 1, 0, 0.5, 0, 0), tolerance = 1e-12)
  expect_equal(y$mean_p_latent_parent, c(0.7, 1, 0, 0.5, 1, 1), tolerance = 1e-12)
})

test_that("pt_reassign_long_lags never moves an edge with a missing lag", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  # Every edge with a lag moves; 4<-3 (no time for 4) stays. Draw 1 out: 3 has 1. Latent parent: 1, 0.4 + 0.6,
  # clamped 1, 0.
  x <- sorted(pt_reassign_long_lags(a, max_lag = -100), c("run", "iteration"))
  expect_equal(x$mean_out_degree[[1]], 0.25, tolerance = 1e-12)
  expect_equal(x$mean_p_latent_parent[[1]], 0.75, tolerance = 1e-12)
  # Nothing beyond Inf: the plain functionals.
  y <- sorted(pt_reassign_long_lags(a, max_lag = Inf), c("run", "iteration"))
  expect_equal(y$mean_out_degree, sorted(pt_mean_out_degree(a), c("run", "iteration"))$mean_out_degree,
               tolerance = 1e-12)
  expect_equal(y$mean_p_latent_parent,
               sorted(pt_mean_p_latent_parent(a), c("run", "iteration"))$mean_p_latent_parent, tolerance = 1e-12)
  expect_error(pt_reassign_long_lags(a, max_lag = NA), "max_lag")
  expect_error(pt_reassign_long_lags(a, max_lag = c(1, 2)), "max_lag")
})

# Two Runs a, b of 2 Draws on the flow Store's input: b's Draws differ from a's.
local_flow_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  local_flow_store(file.path(parent, "a"))
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  fake_store(
    file.path(parent, "b"),
    infections = c("1", "2", "3", "4"),
    parent_sets = rbind(
      sets(1, "1", "4", TRUE, 0.25), sets(1, "1", "", TRUE, 0.75), sets(1, "2", "3,4", TRUE, 1),
      sets(1, "3", "2", FALSE, 1), sets(1, "4", "1", TRUE, 0.5), sets(1, "4", "2", FALSE, 0.5),
      sets(2, "1", "", TRUE, 1), sets(2, "2", "1", TRUE, 1), sets(2, "3", "1,2", TRUE, 0.9),
      sets(2, "3", "", TRUE, 0.1), sets(2, "4", "2,3", FALSE, 1)
    ),
    observation_time = c("1" = 0, "2" = 2, "3" = 10),
    allowed_parents = list("2" = "1", "3" = c("1", "2"), "4" = c("3", "2"))
  )
  file.path(parent, c("a", "b"))
}

flow_golden_epi <- function() data.frame(infection = c("1", "2", "7"), district = c("north", "south", "east"),
                                         travel = c(TRUE, NA, FALSE), zone = factor(c("b", "b", "a")))

golden_candidates <- function() {
  nodes <- jsonlite::read_json(golden_input())$nodes
  parents <- lapply(nodes, function(n) as.character(unlist(n$allowed_parents)))
  data.frame(parent = unlist(parents), child = rep(vapply(nodes, function(n) n$id, ""), lengths(parents)))
}

test_that("flows match the oracle per Draw on fake and golden Stores", {
  flow_times <- c("1" = 0, "2" = 2, "3" = 10, "4" = NA)
  flow_candidates <- data.frame(parent = c("1", "1", "2", "3", "2"), child = c("2", "3", "3", "4", "4"))
  fixtures <- list(
    list(paths = local_flow_runs(), epi = flow_epi(), by = "district", times = flow_times,
         candidates = flow_candidates, lags = c(-Inf, 5, 9, Inf)),
    list(paths = c(one = golden_store(), four = golden_store("replica-exchange")), epi = flow_golden_epi(),
         by = "district", times = golden_observation_times(), candidates = golden_candidates(), lags = c(50, 100)),
    list(paths = golden_store("replica-exchange"), epi = flow_golden_epi(), by = "travel",
         times = golden_observation_times(), candidates = golden_candidates(), lags = 100),
    list(paths = golden_store(), epi = flow_golden_epi(), by = "zone",
         times = golden_observation_times(), candidates = golden_candidates(), lags = 60)
  )
  for (f in fixtures) {
    a <- pt_read(f$paths, epi = f$epi)
    o <- oracle_read(f$paths)
    strata <- oracle_strata(o, f$epi, f$by)
    keys <- c("run", "iteration")
    expect_oracle_equal(pt_attribution(a, f$by), oracle_attribution(o, strata),
                        c(keys, "recipient", "latent", "source"))
    expect_oracle_equal(pt_enrichment(a, f$by), oracle_enrichment(o, strata, f$candidates),
                        c(keys, "source", "recipient"))
    expect_oracle_equal(pt_net_flow(a, f$by), oracle_net_flow(o, strata), c(keys, "stratum"))
    for (lag in f$lags) {
      expect_oracle_equal(pt_reassign_long_lags(a, lag), oracle_reassign_long_lags(o, f$times, lag), keys)
      expect_oracle_equal(pt_reassign_long_lags(a, lag, by = f$by),
                          oracle_reassign_long_lags(o, f$times, lag, strata), c(keys, "stratum"))
    }
    pt_close(a)
  }
})

test_that("flows take runs and phase", {
  paths <- c(one = golden_store(), four = golden_store("replica-exchange"))
  a <- pt_read(paths, epi = flow_golden_epi())
  on.exit(pt_close(a))
  o <- oracle_read(paths)
  strata <- oracle_strata(o, flow_golden_epi(), "district")
  both <- c("burn_in", "sampling")
  keys <- c("run", "iteration")
  expect_oracle_equal(pt_attribution(a, "district", runs = "four", phase = both),
                      oracle_attribution(o, strata, runs = "four", phase = both),
                      c(keys, "recipient", "latent", "source"))
  expect_oracle_equal(pt_enrichment(a, "district", runs = "one", phase = "burn_in"),
                      oracle_enrichment(o, strata, golden_candidates(), runs = "one", phase = "burn_in"),
                      c(keys, "source", "recipient"))
  expect_oracle_equal(pt_net_flow(a, "district", phase = both), oracle_net_flow(o, strata, phase = both),
                      c(keys, "stratum"))
  expect_oracle_equal(pt_reassign_long_lags(a, 100, "district", runs = "one", phase = both),
                      oracle_reassign_long_lags(o, golden_observation_times(), 100, strata, runs = "one",
                                                phase = both), c(keys, "stratum"))
})

test_that("net flow sums to 0 and attribution to 1 at every Draw", {
  a <- pt_read(local_flow_runs(), epi = flow_epi())
  on.exit(pt_close(a))
  net <- dplyr::collect(pt_net_flow(a, "district"))
  expect_equal(as.vector(tapply(net$net_flow, paste(net$run, net$iteration), sum)), rep(0, 4), tolerance = 1e-12)
  share <- dplyr::collect(pt_attribution(a, "district"))
  expect_equal(as.vector(tapply(share$share, paste(share$run, share$iteration, share$recipient), sum)),
               rep(1, 12), tolerance = 1e-12)
  latent <- share[share$latent, ]
  expect_true(all(is.na(latent$source)))
  expect_equal(nrow(latent), 12)
})

test_that("flows need one epi column; reassignment warns without observation times", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  expect_error(pt_attribution(a, by = NULL), "between strata")
  expect_error(pt_net_flow(a, by = "region"), "region")
  expect_error(pt_enrichment(a, by = c("district", "district")), "one epi column")
  b <- suppressMessages(pt_read(local_golden_copy_without_input(), epi = flow_golden_epi()))
  on.exit(pt_close(b), add = TRUE)
  expect_error(pt_enrichment(b, "district"), "input JSON")
  expect_warning(x <- pt_reassign_long_lags(b, 10), "observation time")
  expect_equal(dplyr::collect(x)$mean_out_degree, dplyr::collect(pt_mean_out_degree(b))$mean_out_degree)
})
