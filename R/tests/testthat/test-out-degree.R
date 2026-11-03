# Out-degree posteriors (#128): pt_out_degree_dist(), pt_realized_out_degree(), pt_sample_networks().

# Two Runs over Infections "1" to "4": a has 2 Draws, b has 1. At a's Draw 1, Infection 2's {Latent} set lost 1e-10
# to the floor. Infection 4 is no Infection's parent at any Draw.
local_out_degree_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  ids <- c("1", "2", "3", "4")
  fake_store(file.path(parent, "a"), infections = ids, parent_sets = rbind(
    sets(1, "2", "1", FALSE, 0.6), sets(1, "2", "", TRUE, 0.4 - 1e-10),
    sets(1, "3", "1,2", FALSE, 0.5), sets(1, "3", "2", TRUE, 0.3), sets(1, "3", "", TRUE, 0.2),
    sets(1, "4", "1", FALSE, 0.5), sets(1, "4", "3", FALSE, 0.5),
    sets(2, "2", "1", FALSE, 1), sets(2, "3", "", TRUE, 1), sets(2, "4", "2", FALSE, 1),
    sets(1:2, "1", "", TRUE, 1)
  ))
  fake_store(file.path(parent, "b"), infections = ids, parent_sets = rbind(
    sets(1, "2", "", TRUE, 1), sets(1, "3", "1", FALSE, 0.9), sets(1, "3", "", TRUE, 0.1),
    sets(1, "1", "", TRUE, 1), sets(1, "4", "", TRUE, 1)
  ))
  file.path(parent, c("a", "b"))
}

read_out_degree_runs <- function(paths = NULL, epi = NULL, .local_envir = parent.frame()) {
  if (is.null(paths)) paths <- local_out_degree_runs(.local_envir)
  expect_warning(a <- pt_read(paths, epi = epi), "Draw counts differ")
  a
}

test_that("pt_out_degree_dist() mixes each Infection's Poisson-binomial out-degree over Draws", {
  a <- read_out_degree_runs()
  on.exit(pt_close(a))
  x <- pt_out_degree_dist(a, infections = c("4", "1"))
  expect_s3_class(x, "tbl_df")
  # Infection 1's children: a1 2 (0.6), 3 (0.5), 4 (0.5), pmf (0.1, 0.35, 0.4, 0.15), mean 1.6, var 0.74;
  # a2 2 (1), mean 1, var 0; b1 3 (0.9), mean 0.9, var 0.09. Mixture (0.2, 2.25, 0.4, 0.15) / 3.
  expect_equal(as.data.frame(x), data.frame(
    infection = c("1", "4"), estimate = c(3.5 / 3, 0), lower = c(0, 0), upper = c(3, 0), level = 0.95,
    var_within = c(0.83 / 3, 0), var_between = c(0.86 / 9, 0), p_zero = c(0.2 / 3, 1), draws = 3L
  ), tolerance = 1e-12)

  pmf <- pt_out_degree_dist(a, infections = "1", pmf = TRUE)
  expect_equal(as.data.frame(pmf), data.frame(
    infection = "1", k = 0:3, probability = c(0.2, 2.25, 0.4, 0.15) / 3
  ), tolerance = 1e-12)
  # A narrower level moves the bounds: cdf (0.067, 0.817, 0.95, 1), quantiles at 0.25 and 0.75.
  expect_equal(pt_out_degree_dist(a, infections = "1", level = 0.5)[c("lower", "upper")],
               tibble::tibble(lower = 1, upper = 1))
})
