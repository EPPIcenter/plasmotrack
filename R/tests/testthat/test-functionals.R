# Functionals (a mean over Infections at each Draw), on a fake Store (hand-computed) and golden Stores (oracle).

# Two Runs over Infections "10", "2", "3" (string order "10" < "2" < "3"): a has 2 Draws, b has 1.
# At Run a's Draw 1, Infection 10 has no Parent set row at all: it counts 0 (Convention 2).
# At Run a's Draw 2, a set of weight 1e-10 was dropped by the floor: weights as stored.
local_functional_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  ids <- c("10", "2", "3")
  fake_store(file.path(parent, "a"), infections = ids, parent_sets = rbind(
    sets(1, "2", "10", FALSE, 0.6), sets(1, "2", "", TRUE, 0.4),
    sets(1, "3", "10,2", FALSE, 0.5), sets(1, "3", "2", TRUE, 0.5),
    sets(2, "2", "", TRUE, 1), sets(2, "3", "2", FALSE, 0.25), sets(2, "3", "", TRUE, 0.75 - 1e-10),
    sets(2, "10", "", TRUE, 1)
  ))
  fake_store(file.path(parent, "b"), infections = ids, parent_sets = rbind(
    sets(1, "2", "3", FALSE, 1), sets(1, "3", "", TRUE, 1), sets(1, "10", "2,3", TRUE, 1)
  ))
  file.path(parent, c("a", "b"))
}

read_functional_runs <- function(paths) {
  expect_warning(a <- pt_read(paths), "Draw counts differ")
  a
}

by_draw <- function(x) {
  x <- as.data.frame(dplyr::collect(x))
  x <- x[order(x$run, x$iteration, method = "radix"), , drop = FALSE]
  rownames(x) <- NULL
  x
}

golden_functional_paths <- function() c(one = golden_store(), four = golden_store("replica-exchange"))

test_that("pt_mean_out_degree() is the mean over every Infection of its out-degree, per Draw", {
  a <- read_functional_runs(local_functional_runs())
  on.exit(pt_close(a))
  x <- pt_mean_out_degree(a)
  expect_s3_class(x, "tbl_lazy")
  expect_equal(by_draw(x), data.frame(
    run = c("a", "a", "b"), iteration = c(1, 2, 1),
    mean_out_degree = c((0.6 + 0.5 + 0.5 + 0.5) / 3, 0.25 / 3, 3 / 3)
  ), tolerance = 1e-12)
})

test_that("pt_mean_p_latent_parent() averages the Latent parent probability over every Infection", {
  a <- read_functional_runs(local_functional_runs())
  on.exit(pt_close(a))
  expect_equal(by_draw(pt_mean_p_latent_parent(a)), data.frame(
    run = c("a", "a", "b"), iteration = c(1, 2, 1),
    mean_p_latent_parent = c((0.4 + 0.5) / 3, (1 + 0.75 - 1e-10 + 1) / 3, 2 / 3)
  ), tolerance = 1e-12)
})

test_that("pt_multi_source() gives augmented and observed-only multi-source probabilities", {
  a <- read_functional_runs(local_functional_runs())
  on.exit(pt_close(a))
  expect_equal(by_draw(pt_multi_source(a)), data.frame(
    run = c("a", "a", "b"), iteration = c(1, 2, 1),
    augmented = c(1 / 3, 0, 1 / 3), observed_only = c(0.5 / 3, 0, 1 / 3)
  ), tolerance = 1e-12)
})

test_that("pt_entropy() gives the mean entropy and the mean of exp(entropy)", {
  a <- read_functional_runs(local_functional_runs())
  on.exit(pt_close(a))
  h2 <- -(0.6 * log(0.6) + 0.4 * log(0.4))
  w <- c(0.25, 0.75 - 1e-10)
  h3 <- -sum(w * log(w))
  expect_equal(by_draw(pt_entropy(a)), data.frame(
    run = c("a", "a", "b"), iteration = c(1, 2, 1),
    # Run a, Draw 1: Infection 10 has no row, so entropy 0 and effective size 1.
    mean_entropy = c((h2 + log(2)) / 3, h3 / 3, 0),
    mean_effective_size = c((exp(h2) + 2 + 1) / 3, (1 + exp(h3) + 1) / 3, 1)
  ), tolerance = 1e-12)
})

test_that("pt_composition() is long: every class with weight anywhere, at every Draw", {
  a <- read_functional_runs(local_functional_runs())
  on.exit(pt_close(a))
  x <- dplyr::collect(pt_composition(a))
  expect_equal(names(x), c("run", "iteration", "class", "share"))
  expect_equal(nrow(x), 3 * 5)
  one <- x[x$run == "a" & x$iteration == 1, ]
  one <- one[order(one$class, method = "radix"), ]
  expect_equal(one$class, c("latent_only", "observed_1", "observed_1_latent", "observed_2", "observed_2_latent"))
  expect_equal(one$share, c(0.4, 0.6, 0.5, 0.5, 0) / 3, tolerance = 1e-12)
})

test_that("composition shares sum to each Draw's mean total weight", {
  a <- read_functional_runs(local_functional_runs())
  on.exit(pt_close(a))
  x <- dplyr::collect(pt_composition(a))
  total <- tapply(x$share, paste(x$run, x$iteration), sum)
  # Run a Draw 1: Infection 10 holds no weight; Draw 2: the floor dropped 1e-10.
  expect_equal(as.vector(total[c("a 1", "a 2", "b 1")]), c(2 / 3, 1 - 1e-10 / 3, 1), tolerance = 1e-12)

  g <- pt_read(golden_functional_paths())
  on.exit(pt_close(g), add = TRUE)
  shares <- dplyr::collect(pt_composition(g))
  weights <- dplyr::collect(pt_table(g, "parent_sets"))
  n <- nrow(pt_infections(g))
  expect_equal(nrow(shares) / length(unique(shares$class)), sum(pt_runs(g)$draws))
  key <- function(x) paste(x$run, x$iteration)
  sums <- tapply(shares$share, key(shares), sum)
  mean_total <- tapply(weights$prob, key(weights), sum) / n
  expect_equal(sums[names(mean_total)], mean_total, tolerance = 1e-12)
})

test_that("mean out-degree equals the mean over Infections of pt_out_degree_draws() at every Draw", {
  for (paths in list(local_functional_runs(), golden_functional_paths())) {
    a <- suppressWarnings(pt_read(paths))
    degrees <- dplyr::collect(pt_out_degree_draws(a))
    means <- by_draw(pt_mean_out_degree(a))
    pt_close(a)
    expected <- tapply(degrees$out_degree, paste(degrees$run, degrees$iteration, sep = "\r"), mean)
    expect_equal(means$mean_out_degree, as.vector(expected[paste(means$run, means$iteration, sep = "\r")]),
                 tolerance = 1e-12)
  }
})

test_that("functionals match the oracle per Draw on fake and golden Stores", {
  for (paths in list(local_functional_runs(), golden_functional_paths())) {
    a <- suppressWarnings(pt_read(paths))
    o <- oracle_read(paths)
    draw <- c("run", "iteration")
    expect_oracle_equal(pt_mean_out_degree(a), oracle_mean_out_degree(o), draw)
    expect_oracle_equal(pt_mean_p_latent_parent(a), oracle_mean_p_latent_parent(o), draw)
    expect_oracle_equal(pt_multi_source(a), oracle_multi_source(o), draw)
    expect_oracle_equal(pt_entropy(a), oracle_entropy(o), draw)
    expect_oracle_equal(pt_composition(a), oracle_composition(o), c(draw, "class"))
    pt_close(a)
  }
})

test_that("functionals take runs and phase", {
  a <- pt_read(golden_functional_paths())
  on.exit(pt_close(a))
  o <- oracle_read(golden_functional_paths())
  both <- c("burn_in", "sampling")
  draw <- c("run", "iteration")
  expect_oracle_equal(pt_mean_out_degree(a, runs = "four", phase = both), oracle_mean_out_degree(o, "four", both),
                      draw)
  expect_oracle_equal(pt_entropy(a, phase = "burn_in"), oracle_entropy(o, phase = "burn_in"), draw)
  expect_oracle_equal(pt_composition(a, runs = "one", phase = both), oracle_composition(o, "one", both),
                      c(draw, "class"))
  expect_error(pt_mean_p_latent_parent(a, runs = "five"), "five")
  expect_error(pt_multi_source(a, phase = "warmup"), "phase")
})
