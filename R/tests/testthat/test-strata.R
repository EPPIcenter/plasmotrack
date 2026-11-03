# Strata (#126): `by =` on the functionals, ψ_MAP and the size range. A stratum is one epi column's value; Infections
# with an NA value or missing from epi form the NA stratum.

# Two Runs over Infections "10", "2", "3": a has 2 Draws, b has 1 (as in test-functionals.R).
local_strata_runs <- function(.local_envir = parent.frame()) {
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

# Infection 3 is missing from epi: the NA stratum.
strata_epi <- function() data.frame(infection = c("10", "2", "99"), district = c("x", "x", "y"))

read_strata_runs <- function(paths, epi = strata_epi()) {
  expect_warning(a <- pt_read(paths, epi = epi), "Draw counts differ")
  a
}

by_draw_stratum <- function(x, keys = "stratum") {
  x <- as.data.frame(dplyr::collect(x))
  x <- x[do.call(order, c(unname(as.list(x[c("run", "iteration", keys)])), method = "radix")), , drop = FALSE]
  rownames(x) <- NULL
  x
}

test_that("pt_mean_out_degree(by =) averages out-degree over each stratum's Infections as sources", {
  a <- read_strata_runs(local_strata_runs())
  on.exit(pt_close(a))
  x <- pt_mean_out_degree(a, by = "district")
  expect_s3_class(x, "tbl_lazy")
  # Out-degrees (10, 2, 3): a1 (1.1, 1, 0), a2 (0, 0.25, 0), b1 (0, 1, 2). Stratum x = {10, 2}, NA = {3}.
  expect_equal(by_draw_stratum(x), data.frame(
    run = rep(c("a", "a", "b"), each = 2), iteration = rep(c(1, 2, 1), each = 2),
    stratum = rep(c("x", NA), 3),
    mean_out_degree = c(2.1 / 2, 0, 0.25 / 2, 0, 1 / 2, 2)
  ), tolerance = 1e-12)
})

test_that("recipient functionals average over each stratum's Infections as recipients", {
  a <- read_strata_runs(local_strata_runs())
  on.exit(pt_close(a))
  # Latent parent probability (10, 2, 3): a1 (0, 0.4, 0.5), a2 (1, 1, 0.75 - 1e-10), b1 (1, 0, 1).
  expect_equal(by_draw_stratum(pt_mean_p_latent_parent(a, by = "district")), data.frame(
    run = rep(c("a", "a", "b"), each = 2), iteration = rep(c(1, 2, 1), each = 2),
    stratum = rep(c("x", NA), 3),
    mean_p_latent_parent = c(0.2, 0.5, 1, 0.75 - 1e-10, 0.5, 1)
  ), tolerance = 1e-12)
  # Composition at a1, stratum x: 10 holds nothing, 2 has {10} 0.6 and {L} 0.4.
  x <- by_draw_stratum(pt_composition(a, by = "district"), c("stratum", "class"))
  expect_equal(names(x), c("run", "iteration", "stratum", "class", "share"))
  one <- x[x$run == "a" & x$iteration == 1 & x$stratum %in% "x", ]
  expect_equal(one$class, c("latent_only", "observed_1", "observed_1_latent", "observed_2", "observed_2_latent"))
  expect_equal(one$share, c(0.2, 0.3, 0, 0, 0), tolerance = 1e-12)
  na <- x[x$run == "a" & x$iteration == 1 & is.na(x$stratum), ]
  expect_equal(na$share, c(0, 0, 0.5, 0.5, 0), tolerance = 1e-12)
})

test_that("Infections with an NA epi value join those missing from epi in the NA stratum", {
  epi <- data.frame(infection = c("10", "2"), district = c("x", NA))
  a <- read_strata_runs(local_strata_runs(), epi)
  on.exit(pt_close(a))
  # NA stratum = {2 (NA value), 3 (missing)}; out-degrees as above.
  expect_equal(by_draw_stratum(pt_mean_out_degree(a, by = "district")), data.frame(
    run = rep(c("a", "a", "b"), each = 2), iteration = rep(c(1, 2, 1), each = 2),
    stratum = rep(c("x", NA), 3),
    mean_out_degree = c(1.1, 0.5, 0, 0.125, 0, 1.5)
  ), tolerance = 1e-12)
})

test_that("by = must name one epi column", {
  a <- read_strata_runs(local_strata_runs())
  on.exit(pt_close(a))
  expect_error(pt_mean_out_degree(a, by = "region"), "region")
  expect_error(pt_entropy(a, by = c("district", "district")), "one epi column")
  expect_error(pt_multi_source(a, by = "infection"), "epi column")
  b <- read_strata_runs(local_strata_runs(), epi = NULL)
  on.exit(pt_close(b), add = TRUE)
  expect_error(pt_composition(b, by = "district"), "epi")
})

golden_strata_paths <- function() c(one = golden_store(), four = golden_store("replica-exchange"))

# Golden Infections are "1", "2", "3": 3 is missing, so three strata including NA.
golden_epi <- function() {
  data.frame(infection = c("1", "2", "7"), district = c("north", "south", "east"), travel = c(TRUE, NA, FALSE),
             zone = factor(c("b", "a", "a")), size = c(2.5, 2.5, 1))
}

test_that("stratified functionals match the oracle per Draw on fake and golden Stores", {
  fixtures <- list(
    list(paths = local_strata_runs(), epi = strata_epi(), by = "district"),
    list(paths = golden_strata_paths(), epi = golden_epi(), by = "district"),
    list(paths = golden_strata_paths(), epi = golden_epi(), by = "travel"),
    list(paths = golden_strata_paths(), epi = golden_epi(), by = "zone"),
    list(paths = golden_strata_paths(), epi = golden_epi(), by = "size")
  )
  for (f in fixtures) {
    a <- suppressWarnings(pt_read(f$paths, epi = f$epi))
    o <- oracle_read(f$paths)
    strata <- oracle_strata(o, f$epi, f$by)
    keys <- c("run", "iteration", "stratum")
    by <- f$by
    expect_oracle_equal(pt_mean_out_degree(a, by = by),
                        oracle_stratified_mean(o, oracle_out_degree_draws(o), c(mean_out_degree = "out_degree"),
                                               c(mean_out_degree = 0), strata), keys)
    infection <- oracle_infection_draws(o)
    expect_oracle_equal(pt_mean_p_latent_parent(a, by = by),
                        oracle_stratified_mean(o, infection, c(mean_p_latent_parent = "p_latent_parent"),
                                               c(mean_p_latent_parent = 0), strata), keys)
    expect_oracle_equal(pt_multi_source(a, by = by),
                        oracle_stratified_mean(o, infection, c(augmented = "multi_source_augmented",
                                                               observed_only = "multi_source_observed"),
                                               c(augmented = 0, observed_only = 0), strata), keys)
    expect_oracle_equal(pt_entropy(a, by = by),
                        oracle_stratified_mean(o, infection, c(mean_entropy = "entropy",
                                                               mean_effective_size = "effective_size"),
                                               c(mean_entropy = 0, mean_effective_size = 1), strata), keys)
    expect_oracle_equal(pt_composition(a, by = by), oracle_stratified_composition(o, strata), c(keys, "class"))
    pt_close(a)
  }
})

test_that("stratified functionals take runs and phase", {
  a <- pt_read(golden_strata_paths(), epi = golden_epi())
  on.exit(pt_close(a))
  o <- oracle_read(golden_strata_paths())
  strata <- oracle_strata(o, golden_epi(), "district")
  both <- c("burn_in", "sampling")
  expect_oracle_equal(pt_mean_out_degree(a, by = "district", runs = "four", phase = both),
                      oracle_stratified_mean(o, oracle_out_degree_draws(o, "four", both),
                                             c(mean_out_degree = "out_degree"), c(mean_out_degree = 0), strata,
                                             "four", both),
                      c("run", "iteration", "stratum"))
})

# Σ_S n_S ψ(S) / n = ψ at every Draw.
expect_strata_average_back <- function(a, by, f, columns, keys = character()) {
  n <- nrow(pt_infections(a))
  sizes <- table(infection_strata_sizes(a, by), useNA = "ifany")
  overall <- as.data.frame(dplyr::collect(f(a)))
  strat <- as.data.frame(dplyr::collect(f(a, by = by)))
  weight <- as.vector(sizes[match(as.character(strat$stratum), names(sizes))])
  weight[is.na(strat$stratum)] <- sizes[is.na(names(sizes))]
  draw <- do.call(paste, c(unname(as.list(strat[c("run", "iteration", keys)])), sep = "\r"))
  for (column in columns) {
    back <- tapply(strat[[column]] * weight / n, draw, sum)
    want <- do.call(paste, c(unname(as.list(overall[c("run", "iteration", keys)])), sep = "\r"))
    expect_equal(length(back), nrow(overall))
    expect_true(all(abs(back[want] - overall[[column]]) <= 1e-12), label = column)
  }
}

infection_strata_sizes <- function(a, by) {
  value <- a$epi[[by]]
  if (is.factor(value)) value <- as.character(value)
  as.character(value[match(pt_infections(a)$infection, a$epi$infection)])
}

test_that("per-stratum values weighted by n_S average back to the overall value, Draw by Draw", {
  for (f in list(list(paths = local_strata_runs(), epi = strata_epi()),
                 list(paths = golden_strata_paths(), epi = golden_epi()))) {
    a <- suppressWarnings(pt_read(f$paths, epi = f$epi))
    for (by in c("district", if (!is.null(f$epi$travel)) "travel")) {
      expect_strata_average_back(a, by, pt_mean_out_degree, "mean_out_degree")
      expect_strata_average_back(a, by, pt_mean_p_latent_parent, "mean_p_latent_parent")
      expect_strata_average_back(a, by, pt_multi_source, c("augmented", "observed_only"))
      expect_strata_average_back(a, by, pt_entropy, c("mean_entropy", "mean_effective_size"))
      expect_strata_average_back(a, by, pt_composition, "share", keys = "class")
    }
    pt_close(a)
  }
})
