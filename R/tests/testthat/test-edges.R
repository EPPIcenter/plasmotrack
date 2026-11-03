# Per-Draw Edges and pooled Edge marginals, on fake Stores (hand-computed) and golden Stores (oracle).

test_that("fake_store() writes a Store pt_read() accepts", {
  sets <- data.frame(iteration = c(1, 1, 2), infection = "2", parents = c("10", "", "10"),
                     latent = c(FALSE, TRUE, FALSE), prob = c(0.7, 0.3, 1))
  path <- fake_store(parent_sets = sets, infections = c("10", "2"), burn_in = 0, thin = 1)
  a <- pt_read(c(r = path))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$draws, 2)
  expect_equal(pt_infections(a)$infection, c("10", "2"))
  expect_equal(nrow(dplyr::collect(pt_table(a, "parent_set_members"))), 3)
})

test_that("pt_edge_draws() gives observed parents only, run and iteration first", {
  sets <- data.frame(
    iteration = c(1, 1, 1, 2, 2),
    infection = c("2", "2", "2", "2", "3"),
    parents = c("10", "10,3", "", "3", ""),
    latent = c(FALSE, TRUE, TRUE, FALSE, TRUE),
    prob = c(0.5, 0.2, 0.3, 1, 1)
  )
  path <- fake_store(parent_sets = sets, infections = c("10", "2", "3"))
  a <- pt_read(c(r = path))
  on.exit(pt_close(a))
  draws <- pt_edge_draws(a)
  expect_s3_class(draws, "tbl_lazy")
  expect_equal(colnames(draws), c("run", "iteration", "child", "parent", "prob"))
  got <- dplyr::arrange(dplyr::collect(draws), iteration, child, parent)
  expect_equal(as.data.frame(got), data.frame(
    run = "r", iteration = c(1, 1, 2), child = "2", parent = c("10", "3", "3"), prob = c(0.7, 0.2, 1)
  ))
})

golden_paths <- function() c(one = golden_store(), four = golden_store("replica-exchange"))

test_that("per-Draw tables match the oracle on both golden Stores", {
  a <- read_golden(golden_paths())
  on.exit(pt_close(a))
  o <- oracle_read(golden_paths())
  expect_oracle_equal(pt_edge_draws(a), oracle_edge_draws(o), c("run", "iteration", "child", "parent"))
  expect_oracle_equal(pt_infection_draws(a), oracle_infection_draws(o), c("run", "iteration", "infection"))
  expect_oracle_equal(pt_class_draws(a), oracle_class_draws(o), c("run", "iteration", "infection", "class"))
  expect_oracle_equal(pt_out_degree_draws(a), oracle_out_degree_draws(o), c("run", "iteration", "infection"))
})

test_that("per-Draw tables take runs and phase", {
  a <- read_golden(golden_paths())
  on.exit(pt_close(a))
  o <- oracle_read(golden_paths())
  both <- c("burn_in", "sampling")
  expect_oracle_equal(pt_edge_draws(a, runs = "four", phase = both), oracle_edge_draws(o, "four", both),
                      c("run", "iteration", "child", "parent"))
  expect_oracle_equal(pt_infection_draws(a, phase = "burn_in"), oracle_infection_draws(o, phase = "burn_in"),
                      c("run", "iteration", "infection"))
  expect_oracle_equal(pt_class_draws(a, runs = "one", phase = both), oracle_class_draws(o, "one", both),
                      c("run", "iteration", "infection", "class"))
  burn_in <- dplyr::collect(pt_out_degree_draws(a, phase = "burn_in"))
  expect_true(0 %in% burn_in$iteration)
  expect_oracle_equal(burn_in, oracle_out_degree_draws(o, phase = "burn_in"), c("run", "iteration", "infection"))
  expect_error(pt_edge_draws(a, runs = "five"), "five")
  expect_error(pt_class_draws(a, phase = "warmup"), "phase")
})

# Three Runs of unequal length over Infections "10", "2", "3" (string order: "10" < "2" < "3"). Weights are as
# stored: some Draws lost a set to the floor, so an Infection's weights sum to less than 1 there.
local_three_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  ids <- c("10", "2", "3")
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  # a: 2 Draws. Draw 1 dropped a set of weight 1e-10 under the floor.
  fake_store(file.path(parent, "a"), infections = ids, parent_sets = rbind(
    sets(1, "2", "10", FALSE, 0.6), sets(1, "2", "", TRUE, 0.4 - 1e-10),
    sets(2, "2", "10,3", FALSE, 0.5), sets(2, "2", "3", FALSE, 0.5),
    sets(1:2, "3", "", TRUE, 1), sets(1:2, "10", "", TRUE, 1)
  ))
  # b: 3 Draws.
  fake_store(file.path(parent, "b"), infections = ids, parent_sets = rbind(
    sets(1:3, "2", "10", FALSE, 1),
    sets(1, "10", "3", FALSE, 0.25), sets(1, "10", "", TRUE, 0.75), sets(2:3, "10", "", TRUE, 1),
    sets(1:3, "3", "", TRUE, 1)
  ))
  # c: 1 Draw, no observed edge; Infection 2's only set falls short of 1.
  fake_store(file.path(parent, "c"), infections = ids, parent_sets = rbind(
    sets(1, "2", "", TRUE, 0.999999), sets(1, "3", "", TRUE, 1), sets(1, "10", "", TRUE, 1)
  ))
  parent
}

# The three Runs differ in Draw count, which pt_read() warns about.
read_unequal <- function(paths) {
  expect_warning(a <- pt_read(paths), "Draw counts differ")
  a
}

test_that("pt_edge_marginals() pools every Draw of every Run, weights as stored", {
  a <- read_unequal(local_three_runs())
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$draws, c(2, 3, 1))
  m <- pt_edge_marginals(a)
  expect_s3_class(m, "tbl_df")
  expect_equal(as.data.frame(m), data.frame(
    child = c("10", "2", "2"), parent = c("3", "10", "3"),
    marginal = c(0.25 / 6, (0.6 + 0.5 + 3) / 6, (0.5 + 0.5) / 6), lag = NA_real_
  ), tolerance = 1e-12)
})

test_that("pt_edge_marginals(by_run = TRUE) has a 0 row in each Run lacking an edge", {
  a <- read_unequal(local_three_runs())
  on.exit(pt_close(a))
  m <- pt_edge_marginals(a, by_run = TRUE)
  expect_equal(as.data.frame(m), data.frame(
    run = rep(c("a", "b", "c"), each = 3),
    child = rep(c("10", "2", "2"), 3), parent = rep(c("3", "10", "3"), 3),
    marginal = c(0, 1.1 / 2, 1 / 2, 0.25 / 3, 1, 0, 0, 0, 0), lag = NA_real_
  ), tolerance = 1e-12)
})

test_that("floor drops edges below it: pooled, or below it in every Run", {
  a <- read_unequal(local_three_runs())
  on.exit(pt_close(a))
  expect_equal(pt_edge_marginals(a, floor = 0.1)$parent, c("10", "3"))
  by_run <- pt_edge_marginals(a, floor = 0.1, by_run = TRUE)
  expect_equal(unique(paste(by_run$child, by_run$parent)), c("2 10", "2 3"))
  expect_equal(nrow(by_run), 6)
  expect_error(pt_edge_marginals(a, floor = -1), "floor")
})

test_that("pt_edge_marginals() matches the oracle on fake and golden Stores", {
  paths <- local_three_runs()
  runs <- file.path(paths, c("a", "b", "c"))
  a <- read_unequal(runs)
  on.exit(pt_close(a))
  o <- oracle_read(runs)
  expect_oracle_equal(pt_edge_marginals(a), oracle_edge_marginals(o), c("child", "parent"))
  expect_oracle_equal(pt_edge_marginals(a, by_run = TRUE, floor = 0.05),
                      oracle_edge_marginals(o, by_run = TRUE, floor = 0.05), c("run", "child", "parent"))
  expect_oracle_equal(pt_edge_draws(a), oracle_edge_draws(o), c("run", "iteration", "child", "parent"))
  expect_oracle_equal(pt_infection_draws(a), oracle_infection_draws(o), c("run", "iteration", "infection"))
  expect_oracle_equal(pt_class_draws(a), oracle_class_draws(o), c("run", "iteration", "infection", "class"))
  expect_oracle_equal(pt_out_degree_draws(a), oracle_out_degree_draws(o), c("run", "iteration", "infection"))

  g <- read_golden(golden_paths())
  on.exit(pt_close(g), add = TRUE)
  og <- oracle_read(golden_paths())
  times <- golden_observation_times()
  expect_oracle_equal(pt_edge_marginals(g), oracle_edge_marginals(og, observation_time = times), c("child", "parent"))
  expect_oracle_equal(pt_edge_marginals(g, by_run = TRUE),
                      oracle_edge_marginals(og, by_run = TRUE, observation_time = times), c("run", "child", "parent"))
  expect_oracle_equal(pt_edge_marginals(g, runs = "four"),
                      oracle_edge_marginals(og, runs = "four", observation_time = times), c("child", "parent"))
  expect_false(anyNA(dplyr::collect(pt_edge_marginals(g))$lag))
})

test_that("pt_edge_marginals() skips a Run with no sampling Draws", {
  parent <- local_three_runs()
  # A fourth Run whose only sampling Window is unclosed: it has 0 Draws.
  fake_store(file.path(parent, "d"), infections = c("10", "2", "3"), unclosed = 1,
             parent_sets = data.frame(iteration = 1, infection = "2", parents = "3", prob = 1))
  a <- suppressWarnings(pt_read(parent))
  on.exit(pt_close(a))
  expect_message(m <- pt_edge_marginals(a, by_run = TRUE), "d")
  expect_equal(unique(m$run), c("a", "b", "c"))
})

test_that("the Latent parent never appears in any output", {
  for (paths in list(local_three_runs(), golden_paths())) {
    a <- suppressWarnings(suppressMessages(pt_read(paths)))
    ids <- pt_infections(a)$infection
    outputs <- list(
      dplyr::collect(pt_edge_draws(a, phase = c("burn_in", "sampling"))), pt_edge_marginals(a, by_run = TRUE),
      dplyr::collect(pt_infection_draws(a)), dplyr::collect(pt_class_draws(a)),
      dplyr::collect(pt_out_degree_draws(a))
    )
    pt_close(a)
    for (x in outputs) {
      for (column in intersect(c("child", "parent", "infection"), names(x))) expect_true(all(x[[column]] %in% ids))
      text <- unlist(lapply(x, function(v) if (is.character(v)) v))
      expect_false(any(grepl("_copy", text, fixed = TRUE)))
    }
  }
})

test_that("pt_class_draws() and pt_out_degree_draws() are full zero-filled grids", {
  a <- read_unequal(local_three_runs())
  on.exit(pt_close(a))
  classes <- dplyr::collect(pt_class_draws(a))
  # 6 Draws x 3 Infections x the 3 classes with weight anywhere.
  expect_equal(sort(unique(classes$class)), c("latent_only", "observed_1", "observed_2"))
  expect_equal(nrow(classes), 6 * 3 * 3)
  draw <- classes[classes$run == "a" & classes$iteration == 2 & classes$infection == "2", ]
  draw <- draw[order(draw$class), ]
  expect_equal(draw$share, c(0, 0.5, 0.5))
  expect_equal(draw$n_observed, c(0, 1, 2))
  expect_equal(draw$has_latent_parent, c(TRUE, FALSE, FALSE))
  # Run c's Infection 2 keeps its stored weight, short of 1.
  expect_equal(classes$share[classes$run == "c" & classes$infection == "2" & classes$class == "latent_only"], 0.999999)

  degrees <- dplyr::collect(pt_out_degree_draws(a))
  expect_equal(nrow(degrees), 6 * 3)
  expect_equal(degrees$out_degree[degrees$run == "a" & degrees$iteration == 2 & degrees$infection == "3"], 1)
  expect_equal(degrees$out_degree[degrees$run == "a" & degrees$iteration == 1 & degrees$infection == "10"], 0.6)
  expect_equal(sum(degrees$out_degree[degrees$run == "c"]), 0)
})

test_that("pt_infection_draws() follows the spec's formulas", {
  a <- read_unequal(local_three_runs())
  on.exit(pt_close(a))
  x <- dplyr::collect(pt_infection_draws(a, runs = "a"))
  one <- x[x$iteration == 1 & x$infection == "2", ]
  w <- c(0.6, 0.4 - 1e-10)
  expect_equal(one$p_latent_parent, 0.4 - 1e-10)
  expect_equal(one$expected_observed_parents, 0.6)
  expect_equal(one$entropy, -sum(w * log(w)))
  expect_equal(one$effective_size, exp(-sum(w * log(w))))
  expect_equal(one$multi_source_augmented, 0)
  two <- x[x$iteration == 2 & x$infection == "2", ]
  expect_equal(two$expected_observed_parents, 1.5)
  expect_equal(two$multi_source_observed, 0.5)
  expect_equal(two$multi_source_augmented, 0.5)
})
