# pt_interval(): pooled mean and equal-tailed type-7 interval of any per-Draw table.

local_interval_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  ids <- c("10", "2", "3")
  # Run a, 2 Draws: mean out-degree 2.1 / 3 = 0.7, then 0.25 / 3. Run b, 1 Draw: 3 / 3 = 1.
  fake_store(file.path(parent, "a"), infections = ids, parent_sets = rbind(
    sets(1, "2", "10", FALSE, 0.6), sets(1, "2", "", TRUE, 0.4),
    sets(1, "3", "10,2", FALSE, 0.5), sets(1, "3", "2", TRUE, 0.5), sets(1, "10", "", TRUE, 1),
    sets(2, "2", "", TRUE, 1), sets(2, "3", "2", FALSE, 0.25), sets(2, "3", "", TRUE, 0.75),
    sets(2, "10", "", TRUE, 1)
  ))
  fake_store(file.path(parent, "b"), infections = ids, parent_sets = rbind(
    sets(1, "2", "3", FALSE, 1), sets(1, "3", "", TRUE, 1), sets(1, "10", "2,3", TRUE, 1)
  ))
  file.path(parent, c("a", "b"))
}

# Quantiles must match the oracle exactly, not just within tolerance.
expect_same_quantiles <- function(actual, expected, keys) {
  actual <- expect_oracle_equal(actual, expected, keys)
  expected <- as.data.frame(expected)
  expected <- expected[do.call(order, c(unname(as.list(expected[keys])), method = "radix")), , drop = FALSE]
  expect_identical(actual$lower, expected$lower)
  expect_identical(actual$upper, expected$upper)
  expect_identical(as.integer(actual$draws), as.integer(expected$draws))
}

test_that("pt_interval() pools every Draw of every Run: mean and type-7 interval", {
  expect_warning(a <- pt_read(local_interval_runs()), "Draw counts differ")
  on.exit(pt_close(a))
  x <- pt_interval(pt_mean_out_degree(a), level = 0.5)
  expect_s3_class(x, "tbl_df")
  values <- c(0.7, 0.25 / 3, 1)
  expect_equal(as.data.frame(x), data.frame(
    estimate = sum(values) / 3, lower = (0.25 / 3 + 0.7) / 2, upper = (0.7 + 1) / 2, level = 0.5, draws = 3L
  ), tolerance = 1e-12)
  expect_equal(pt_interval(pt_mean_out_degree(a))$level, 0.95)
})

test_that("pt_interval() matches the oracle on functionals, with and without by", {
  paths <- c(one = golden_store(), four = golden_store("replica-exchange"))
  a <- read_golden(paths)
  on.exit(pt_close(a))
  degree <- dplyr::collect(pt_mean_out_degree(a))
  expect_same_quantiles(pt_interval(pt_mean_out_degree(a)), oracle_interval(degree, value = "mean_out_degree"),
                        character())
  composition <- dplyr::collect(pt_composition(a))
  expect_same_quantiles(pt_interval(pt_composition(a), by = class, level = 0.9),
                        oracle_interval(composition, "class", "share", 0.9), "class")
  expect_same_quantiles(pt_interval(pt_composition(a), by = "class", level = 0.9),
                        oracle_interval(composition, "class", "share", 0.9), "class")
  infections <- dplyr::collect(pt_infection_draws(a))
  expect_same_quantiles(pt_interval(pt_infection_draws(a), by = infection, value = entropy, level = 0.8),
                        oracle_interval(infections, "infection", "entropy", 0.8), "infection")
})

test_that("several value columns give a value column naming each", {
  expect_warning(a <- pt_read(local_interval_runs()), "Draw counts differ")
  on.exit(pt_close(a))
  x <- pt_interval(pt_entropy(a))
  expect_equal(names(x), c("value", "estimate", "lower", "upper", "level", "draws"))
  expect_same_quantiles(x, oracle_interval(dplyr::collect(pt_entropy(a)), character(),
                                           c("mean_entropy", "mean_effective_size")), "value")
  one <- pt_interval(pt_entropy(a), value = mean_entropy)
  expect_equal(names(one), c("estimate", "lower", "upper", "level", "draws"))
})

test_that("pt_interval() takes a per-Draw table the user builds, lazy or not", {
  paths <- c(one = golden_store(), four = golden_store("replica-exchange"))
  a <- read_golden(paths)
  on.exit(pt_close(a))
  # A raw Trace: window is bookkeeping, not a value.
  freqs <- pt_table(a, "allele_frequencies")
  x <- pt_interval(freqs, by = c(locus, allele))
  expect_equal(names(x), c("locus", "allele", "estimate", "lower", "upper", "level", "draws"))
  expect_same_quantiles(x, oracle_interval(dplyr::collect(freqs), c("locus", "allele"), "frequency"),
                        c("locus", "allele"))
  expect_equal(unique(x$draws), sum(pt_runs(a)$draws))

  # Built in R: a derived quantity over two Runs of unequal length.
  user <- data.frame(run = rep(c("x", "y"), c(4, 7)), iteration = c(1:4, 1:7),
                     group = rep(c("p", "q"), length.out = 11), score = c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5),
                     other = seq(0.5, 5.5, by = 0.5))
  expect_same_quantiles(pt_interval(user, by = group, level = 0.6),
                        oracle_interval(user, "group", c("score", "other"), 0.6), c("group", "value"))
  lazy <- dplyr::copy_to(dbplyr::remote_con(freqs), user, "user_draws")
  expect_same_quantiles(pt_interval(lazy, by = group, value = score),
                        oracle_interval(user, "group", "score"), "group")
})

test_that("pt_interval() sorts by keys as strings and handles a single Draw", {
  user <- data.frame(run = "r", iteration = c(1, 1, 1, 2, 2, 2), id = c("2", "10", "3"), v = c(1, 2, 3, 4, 5, 6))
  x <- pt_interval(user, by = id)
  expect_equal(x$id, c("10", "2", "3"))
  expect_equal(x$estimate, c(3.5, 2.5, 4.5))
  one <- pt_interval(user[1, ])
  expect_equal(c(one$estimate, one$lower, one$upper, one$draws), c(1, 1, 1, 1))
})

test_that("pt_interval() rejects tables that are not per Draw", {
  user <- data.frame(run = "r", iteration = c(1, 1, 2), id = c("a", "b", "a"), v = c(1, 2, 3), label = "x")
  expect_error(pt_interval(user), "more than one row per Draw")
  expect_error(pt_interval(user[c("iteration", "v")]), "run")
  expect_error(pt_interval(user, by = id, value = label), "numeric")
  expect_error(pt_interval(user[c("run", "iteration", "id")], by = id), "value column")
  expect_error(pt_interval(user, by = id, level = 1), "level")
  expect_error(pt_interval(user, by = id, level = c(0.5, 0.9)), "level")
  user$v[1] <- NA
  expect_error(pt_interval(user, by = id), "missing")
})

test_that("pt_interval() names the keys with missing values and how to drop them", {
  user <- data.frame(run = "r", iteration = rep(1:2, each = 3), id = c("a", "b", "c"), v = c(NA, 1, 2, 3, NA, 4),
                     w = 1)
  err <- expect_error(pt_interval(user, by = id), class = "rlang_error")
  msg <- conditionMessage(err)
  expect_match(msg, "id = \"a\"", fixed = FALSE)
  expect_match(msg, "id = \"b\"")
  expect_no_match(msg, "id = \"c\"")
  expect_match(msg, "filter(!is.na(v))", fixed = TRUE)
  expect_no_match(msg, "is.na(w)", fixed = TRUE)
  # Without `by`, the Draws are counted.
  expect_match(conditionMessage(expect_error(pt_interval(user[user$id == "a", ]))), "1 Draw")
})
