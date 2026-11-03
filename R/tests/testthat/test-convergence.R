# Convergence report (#124): the Draw set, per-quantity status, R-hat and ESS, pt_as_draws().

golden_pair <- function() golden_store(c("single-chain", "replica-exchange"))

test_that("pt_convergence() quantities equal posterior on the oracle's matrices (golden Runs)", {
  paths <- golden_pair()
  a <- read_golden(paths)
  on.exit(pt_close(a))
  conv <- pt_convergence(a)
  expect_s3_class(conv, "pt_convergence")
  o <- oracle_read(paths)
  expected <- oracle_convergence_quantities(oracle_convergence_matrices(o, oracle_convergence_draws(o)))
  expect_s3_class(conv$quantities, "tbl_df")
  expect_identical(conv$quantities$name, expected$name)
  expect_oracle_equal(conv$quantities, expected, keys = c("group", "name"))
})

# Run scalars over iterations 0 (burn-in) to 4 (sampling 1-4), with the given sampling values.
scalar_rows <- function(mean_coi = 1.5, mean_strains_tx = 1.2, parent_set_size_prob = 0.5, prior = -1,
                        likelihood = -2) {
  sampling <- data.frame(iteration = 1:4, mean_coi = mean_coi, mean_strains_tx = mean_strains_tx,
                         parent_set_size_prob = parent_set_size_prob, prior = prior, likelihood = likelihood)
  rows <- rbind(transform(sampling[1, ], iteration = 0), sampling)
  cbind(rows, observation = -0.5, transmission = -1.5, source = 0)
}

# Two Runs of 4 Draws: mean_strains_tx constant, mean_coi stuck (1 in a, 2 in b), parent_set_size_prob mixing,
# and a log posterior exactly half 0 and half 1 (the fold case).
local_status_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  fake_store(file.path(parent, "a"), infections = "1", samples = 4, tables = list(run_scalars = scalar_rows(
    mean_coi = 1, parent_set_size_prob = c(0.1, 0.4, 0.2, 0.3), prior = c(0, 1, 0, 1), likelihood = 0
  )))
  fake_store(file.path(parent, "b"), infections = "1", samples = 4, tables = list(run_scalars = scalar_rows(
    mean_coi = 2, parent_set_size_prob = c(0.3, 0.1, 0.4, 0.2), prior = c(1, 0, 1, 0), likelihood = 0
  )))
  file.path(parent, c("a", "b"))
}

quantity <- function(conv, name) as.list(conv$quantities[conv$quantities$name == name, ])

test_that("a quantity with one value everywhere is constant: NA diagnostics, never flagged", {
  a <- pt_read(local_status_runs())
  on.exit(pt_close(a))
  q <- quantity(pt_convergence(a), "mean_strains_tx")
  expect_identical(q$status, "constant")
  expect_identical(c(q$rhat, q$ess_bulk, q$ess_tail), rep(NA_real_, 3))
  expect_false(q$flagged)
})

test_that("Runs constant at different values are stuck: R-hat Inf, always flagged", {
  a <- pt_read(local_status_runs())
  on.exit(pt_close(a))
  q <- quantity(pt_convergence(a, rhat = 1e6, ess = 1e-6), "mean_coi")
  expect_identical(q$status, "stuck")
  expect_identical(q$rhat, Inf)
  expect_true(q$flagged)
})

test_that("a mixing quantity gets posterior's rhat and ESS on the Runs' matrix", {
  a <- pt_read(local_status_runs())
  on.exit(pt_close(a))
  conv <- pt_convergence(a)
  m <- cbind(c(0.1, 0.4, 0.2, 0.3), c(0.3, 0.1, 0.4, 0.2))
  q <- quantity(conv, "parent_set_size_prob")
  expect_identical(q$status, "mixing")
  expect_equal(q$rhat, posterior::rhat(m), tolerance = 1e-12)
  expect_equal(q$ess_bulk, posterior::ess_bulk(m), tolerance = 1e-12)
  expect_equal(q$ess_tail, posterior::ess_tail(m), tolerance = 1e-12)
  expect_identical(q$flagged, q$rhat > 1.01 || q$ess_bulk < 200)
})

test_that("R-hat falls back to rhat_basic when rhat() is NA (half 0, half 1)", {
  a <- pt_read(local_status_runs())
  on.exit(pt_close(a))
  m <- cbind(c(0, 1, 0, 1), c(1, 0, 1, 0))
  expect_true(is.na(posterior::rhat(m)))
  q <- quantity(pt_convergence(a), "log_posterior")
  expect_identical(q$status, "mixing")
  expect_equal(q$rhat, posterior::rhat_basic(m), tolerance = 1e-12)
})

test_that("Runs with different thinning strides error", {
  parent <- withr::local_tempdir()
  fake_store(file.path(parent, "a"), infections = "1", samples = 4, thin = 1)
  fake_store(file.path(parent, "b"), infections = "1", samples = 8, thin = 2)
  a <- pt_read(parent)
  on.exit(pt_close(a))
  expect_error(pt_convergence(a), "different thinning strides")
  expect_error(pt_as_draws(a), "different thinning strides")
  expect_s3_class(pt_convergence(a, runs = "a"), "pt_convergence")
})

# Run a has 4 Draws, b has 6: every diagnostic uses each Run's first 4.
local_unequal_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  fake_store(file.path(parent, "a"), infections = "1", samples = 4, tables = list(run_scalars = scalar_rows(
    mean_coi = c(1, 3, 2, 4), prior = c(-5, -4, -6, -5)
  )))
  b <- scalar_rows(mean_coi = c(2, 1, 4, 3), prior = c(-7, -8, -6, -7))
  b <- rbind(b, transform(b[2:3, ], iteration = 5:6, mean_coi = 100, prior = 100))
  fake_store(file.path(parent, "b"), infections = "1", samples = 6, tables = list(run_scalars = b))
  parent
}

test_that("unequal Runs are cut to the shortest Run's first Draws", {
  expect_warning(a <- pt_read(local_unequal_runs()), "Draw counts differ")
  on.exit(pt_close(a))
  conv <- pt_convergence(a)
  expect_equal(as.data.frame(conv$runs), data.frame(
    run = c("a", "b"), draws = 4L, first_iteration = 1, last_iteration = 4,
    log_posterior = c(-5, -7) - 2, gap = c(0, 2)
  ))
  m <- cbind(c(1, 3, 2, 4), c(2, 1, 4, 3))
  q <- quantity(conv, "mean_coi")
  expect_equal(q$rhat, posterior::rhat(m), tolerance = 1e-12)
  expect_equal(q$ess_bulk, posterior::ess_bulk(m), tolerance = 1e-12)
})

test_that("pt_as_draws() cuts to the shortest Run with a message", {
  expect_warning(a <- pt_read(local_unequal_runs()), "Draw counts differ")
  on.exit(pt_close(a))
  expect_message(d <- pt_as_draws(a, groups = c("log_posterior", "global")), "shortest Run's 4 Draws")
  expect_s3_class(d, "draws_array")
  expect_identical(posterior::variables(d), c("log_posterior", "mean_coi", "mean_strains_tx",
                                              "parent_set_size_prob"))
  expect_identical(posterior::niterations(d), 4L)
  expect_equal(unname(posterior::extract_variable_matrix(d, "mean_coi")), cbind(c(1, 3, 2, 4), c(2, 1, 4, 3)))
})

test_that("pt_as_draws() variables match pt_convergence() names (golden Runs)", {
  a <- read_golden(golden_pair())
  on.exit(pt_close(a))
  expect_silent(d <- pt_as_draws(a))
  q <- pt_convergence(a)$quantities
  expect_identical(posterior::variables(d), q$name)
  expect_identical(posterior::nchains(d), 2L)
  expect_equal(posterior::rhat(posterior::extract_variable_matrix(d, "mean_coi")),
               q$rhat[q$name == "mean_coi"], tolerance = 1e-12)
})

# Infections 1, 2, 3 over two Runs of 4 Draws. Edge 2<-1: 0.5 at a's Draws 1-2, never in b (Run means 0.25, 0).
# Edge 3<-1: 0.02 at a's Draw 1 and 0.03 at b's Draw 1 (Run means 0.005, 0.0075). Edge 3<-2: 1 at every Draw of
# a, never in b.
local_edge_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  sets <- function(iteration, infection, parents, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = FALSE, prob = prob)
  }
  fake_store(file.path(parent, "a"), infections = c("1", "2", "3"), samples = 4, parent_sets = rbind(
    sets(1:2, "2", "1", 0.5), sets(1, "3", "1", 0.02), sets(1:4, "3", "2", 1)
  ))
  fake_store(file.path(parent, "b"), infections = c("1", "2", "3"), samples = 4, parent_sets = rbind(
    sets(1, "3", "1", 0.03), sets(1:4, "2", "", 0)
  ))
  file.path(parent, c("a", "b"))
}

test_that("edges below edge_floor in every Run are left out; edges are zero-filled", {
  a <- pt_read(local_edge_runs())
  on.exit(pt_close(a))
  q <- pt_convergence(a)$quantities
  expect_identical(q$name[q$group == "edge"], c("edge[2,1]", "edge[3,2]"))
  m <- cbind(c(0.5, 0.5, 0, 0), c(0, 0, 0, 0))
  expect_equal(quantity(pt_convergence(a), "edge[2,1]")$rhat, posterior::rhat(m), tolerance = 1e-12)
  stuck <- quantity(pt_convergence(a), "edge[3,2]")
  expect_identical(c(stuck$status, format(stuck$rhat)), c("stuck", "Inf"))
  low <- pt_convergence(a, edge_floor = 0.0075)$quantities
  expect_identical(low$name[low$group == "edge"], c("edge[2,1]", "edge[3,1]", "edge[3,2]"))
  expect_identical(pt_convergence(a, edge_floor = 0.00751)$quantities$name[low$group == "edge"][1:2],
                   c("edge[2,1]", "edge[3,2]"))
})

test_that("Runs with 0 sampling Draws are skipped with a message", {
  parent <- withr::local_tempdir()
  fake_store(file.path(parent, "a"), infections = "1", samples = 4)
  fake_store(file.path(parent, "b"), infections = "1", samples = 4)
  drop_windows(file.path(parent, "b"), 1)
  a <- suppressWarnings(pt_read(parent))
  on.exit(pt_close(a))
  expect_message(conv <- pt_convergence(a), "Skipping Run .*b.*: no sampling Draws")
  expect_identical(conv$runs$run, "a")
  expect_message(d <- pt_as_draws(a, groups = "global"), "Skipping Run")
  expect_identical(posterior::nchains(d), 1L)
})

test_that("iterations restricts every diagnostic to a range of iterations", {
  expect_warning(a <- pt_read(local_unequal_runs()), "Draw counts differ")
  on.exit(pt_close(a))
  conv <- pt_convergence(a, iterations = c(2, 3))
  expect_equal(conv$runs$draws, c(2L, 2L))
  expect_equal(conv$runs$first_iteration, c(2, 2))
  expect_equal(conv$runs$last_iteration, c(3, 3))
  expect_equal(quantity(conv, "mean_coi")$rhat, posterior::rhat(cbind(c(3, 2), c(1, 4))), tolerance = 1e-12)
  expect_message(conv <- pt_convergence(a, iterations = c(5, 6)), "Skipping Run .*a")
  expect_identical(conv$runs$run, "b")
  expect_error(pt_convergence(a, iterations = c(50, 60)), "No sampling Draws")
  expect_error(pt_convergence(a, iterations = c(3, 2)), "first <= last")
})

test_that("pt_convergence() and pt_as_draws() check their arguments", {
  a <- pt_read(local_status_runs())
  on.exit(pt_close(a))
  expect_error(pt_convergence(a, edge_floor = 2), "edge_floor")
  expect_error(pt_convergence(a, rhat = 0), "rhat")
  expect_error(pt_convergence(a, ess = -1), "ess")
  expect_error(pt_convergence(a, spread = NA), "spread")
  expect_error(pt_as_draws(a, groups = "edges"), "groups")
  expect_error(pt_convergence(a, runs = "z"), "Unknown Run")
})

test_that("print() shows the Runs, the Draws used, the gap and flagged counts per group", {
  expect_warning(a <- pt_read(local_unequal_runs()), "Draw counts differ")
  on.exit(pt_close(a))
  expect_snapshot(print(pt_convergence(a)))
  expect_snapshot(print(pt_convergence(a, iterations = c(2, 4))))
})

test_that("print() on the golden Runs", {
  a <- read_golden(golden_pair())
  on.exit(pt_close(a))
  expect_snapshot(print(pt_convergence(a)))
})
