# Infection times and the Ordering, on a fake Store (hand-computed) and golden Stores (oracle).

# Infections 10, 2, 3; sampling Draws 1-3 (iteration 0 is burn-in). Ranks per Draw (0 = acquired first):
#   iteration 0: 10, 2, 3    iteration 1: 10, 2, 3    iteration 2: 2, 10, 3    iteration 3: 3, 2, 10
ordering_rows <- function() {
  data.frame(
    iteration = rep(0:3, each = 3),
    infection = rep(c("10", "2", "3"), 4),
    rank = c(0, 1, 2, 0, 1, 2, 1, 0, 2, 2, 1, 0)
  )
}

# Durations: 10 -> 1, 4, 2 at Draws 1-3; 2 -> 5 throughout; 3 -> 3, 3, 6. Iteration 0 has 100 (burn-in only).
duration_rows <- function() {
  data.frame(
    iteration = rep(0:3, each = 3),
    infection = rep(c("10", "2", "3"), 4),
    expected_false_positives = 0, expected_false_negatives = 0,
    infection_duration = c(100, 100, 100, 1, 5, 3, 4, 5, 3, 2, 5, 6)
  )
}

local_ordering_store <- function(path = NULL, observation_time = c("10" = 20, "2" = 30), .local_envir = parent.frame()) {
  fake_store(path = path, infections = c("10", "2", "3"), samples = 3, observation_time = observation_time,
             tables = list(ordering = ordering_rows(), infections = duration_rows()), .local_envir = .local_envir)
}

local_ordering_analysis <- function(.local_envir = parent.frame()) {
  a <- pt_read(c(r = local_ordering_store(.local_envir = .local_envir)))
  withr::defer(pt_close(a), envir = .local_envir)
  a
}

test_that("pt_order_pair_draws() says whether i, the smaller id as a string, was acquired before j", {
  a <- local_ordering_analysis()
  pairs <- pt_order_pair_draws(a)
  expect_s3_class(pairs, "tbl_lazy")
  expect_equal(colnames(pairs), c("run", "iteration", "i", "j", "before"))
  got <- as.data.frame(dplyr::arrange(dplyr::collect(pairs), iteration, i, j))
  expect_equal(got, data.frame(
    run = "r", iteration = rep(1:3, each = 3),
    i = c("10", "10", "2"), j = c("2", "3", "3"),
    before = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)
  ))
  burn_in <- dplyr::collect(pt_order_pair_draws(a, phase = "burn_in"))
  expect_equal(burn_in$iteration, c(0, 0, 0))
  expect_true(all(burn_in$before))
  # Filtering before collect() is how the O(n^2) table is meant to be used.
  one <- dplyr::collect(dplyr::arrange(dplyr::filter(pairs, i == "2"), iteration))
  expect_equal(one$before, c(TRUE, TRUE, FALSE))
})

# Runs r (the flipping Store above) and s (3 Draws, Ordering always 10, 2, 3) under one parent directory.
local_two_ordering_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  local_ordering_store(file.path(parent, "r"), .local_envir = .local_envir)
  fake_store(file.path(parent, "s"), infections = c("10", "2", "3"), samples = 3,
             observation_time = c("10" = 20, "2" = 30), .local_envir = .local_envir)
  a <- pt_read(parent)
  withr::defer(pt_close(a), envir = .local_envir)
  a
}

test_that("pt_order_pairs() gives the share of Draws at which i preceded j", {
  a <- local_ordering_analysis()
  got <- pt_order_pairs(a)
  expect_s3_class(got, "tbl_df")
  expect_equal(as.data.frame(got), data.frame(
    i = c("10", "10", "2"), j = c("2", "3", "3"), p_before = c(1, 2, 2) / 3, draws = 3
  ))
})

test_that("pt_order_pairs() pools Draws over Runs, or keeps each Run with by_run", {
  a <- local_two_ordering_runs()
  expect_equal(as.data.frame(pt_order_pairs(a)), data.frame(
    i = c("10", "10", "2"), j = c("2", "3", "3"), p_before = c(4, 5, 5) / 6, draws = 6
  ))
  expect_equal(as.data.frame(pt_order_pairs(a, by_run = TRUE)), data.frame(
    run = rep(c("r", "s"), each = 3), i = c("10", "10", "2"), j = c("2", "3", "3"),
    p_before = c(1 / 3, 2 / 3, 2 / 3, 1, 1, 1), draws = 3
  ))
  expect_equal(pt_order_pairs(a, runs = "s")$p_before, c(1, 1, 1))
})

test_that("pt_order_pairs(infections =) keeps only pairs touching those Infections", {
  a <- local_ordering_analysis()
  expect_equal(as.data.frame(pt_order_pairs(a, infections = "3")), data.frame(
    i = c("10", "2"), j = c("3", "3"), p_before = c(2, 2) / 3, draws = 3
  ))
  expect_equal(as.data.frame(pt_order_pairs(a, infections = c("2", "10"))), data.frame(
    i = c("10", "10", "2"), j = c("2", "3", "3"), p_before = c(1, 2, 2) / 3, draws = 3
  ))
  expect_equal(pt_order_pairs(a, infections = "2")[c("i", "j")], tibble::tibble(i = c("10", "2"), j = c("2", "3")))
  expect_error(pt_order_pairs(a, infections = c("2", "99")), "99")
})
