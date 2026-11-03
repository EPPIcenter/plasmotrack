# Plots: each returns a ggplot. Tests check the layers (geoms) and the data they map, never pixels.

skip_if_not_installed("ggplot2")

geoms <- function(p) unname(vapply(p$layers, function(l) class(l$geom)[[1]], character(1)))

# Two Runs over iterations 0-3; iterations 0-1 are burn-in.
trace_rows <- function() {
  tibble::tibble(
    run = rep(c("a", "b"), each = 4), iteration = rep(0:3, 2), phase = rep(c("burn_in", "burn_in", "sampling",
                                                                             "sampling"), 2),
    value = c(1, 2, 3, 4, 5, 6, 7, 8), other = 0
  )
}

test_that("pt_plot_trace() draws one line per Run over iteration", {
  rows <- trace_rows()[trace_rows()$phase == "sampling", c("run", "iteration", "value", "other")]
  p <- pt_plot_trace(rows, value = value)
  expect_s3_class(p, "ggplot")
  expect_identical(geoms(p), "GeomLine")
  d <- ggplot2::layer_data(p, 1)
  expect_identical(d$x, c(2, 3, 2, 3))
  expect_identical(d$y, c(3, 4, 7, 8))
  expect_identical(length(unique(d$colour)), 2L)
  expect_identical(length(unique(d$group)), 2L)
  # Only the columns it needs.
  expect_named(p$data, c("run", "iteration", "value"), ignore.order = TRUE)
})

test_that("pt_plot_trace() shades burn-in when the table has a phase column", {
  p <- pt_plot_trace(trace_rows(), value = "value")
  expect_identical(geoms(p), c("GeomRect", "GeomLine"))
  shade <- ggplot2::layer_data(p, 1)
  expect_identical(nrow(shade), 1L)
  expect_identical(c(shade$xmin, shade$xmax), c(0, 1))
  expect_identical(nrow(ggplot2::layer_data(p, 2)), 8L)
})

test_that("pt_plot_trace() facets by `by` and refuses more than max_series keys", {
  rows <- tibble::tibble(run = "a", iteration = rep(1:2, each = 3), key = rep(c("x", "y", "z"), 2), v = 1:6)
  p <- pt_plot_trace(rows, value = v, by = key)
  expect_s3_class(p$facet, "FacetWrap")
  expect_identical(sort(unique(ggplot2::layer_data(p, 1)$PANEL)), factor(1:3))
  expect_error(pt_plot_trace(rows, value = v, by = key, max_series = 2), "3 series")
  expect_error(pt_plot_trace(rows, value = v, by = key, max_series = 2), "[Ff]ilter")
  expect_error(pt_plot_trace(rows, value = v), "more than one row per Draw")
})

test_that("pt_plot_trace() takes a lazy table and collects only the columns it needs", {
  a <- read_golden(golden_store())
  on.exit(pt_close(a))
  x <- pt_table(a, "run_scalars", phase = c("burn_in", "sampling"))
  p <- pt_plot_trace(x, value = likelihood)
  expect_named(p$data, c("run", "iteration", "phase", "likelihood"), ignore.order = TRUE)
  expect_identical(geoms(p), c("GeomRect", "GeomLine"))
  burn <- dplyr::collect(dplyr::filter(x, phase == "burn_in"))
  shade <- ggplot2::layer_data(p, 1)
  expect_identical(c(shade$xmin, shade$xmax), as.numeric(range(burn$iteration)))
})

test_that("a plot without ggplot2 asks to install it", {
  local_mocked_bindings(plot_package = function() "ggplot2notinstalled")
  expect_error(pt_plot_trace(trace_rows(), value = value), class = "rlib_error_package_not_found")
  expect_error(pt_plot_trace(trace_rows(), value = value), "ggplot2notinstalled.*to draw plots")
})

test_that("pt_plot_rank() bins the pooled ranks per Run, zero bins kept", {
  rows <- tibble::tibble(run = rep(c("a", "b"), each = 3), iteration = rep(1:3, 2), v = c(3, 1, 2, 6, 4, 5))
  p <- pt_plot_rank(rows, value = v, bins = 2)
  expect_s3_class(p, "ggplot")
  expect_identical(geoms(p), "GeomCol")
  expect_identical(p$data$run, c("a", "a", "b", "b"))
  expect_identical(p$data$bin, c(1L, 2L, 1L, 2L))
  expect_identical(p$data$count, c(3L, 0L, 0L, 3L))
  d <- ggplot2::layer_data(p, 1)
  expect_identical(length(unique(d$fill)), 2L)
  expect_identical(length(unique(d$PANEL)), 2L)
})

test_that("pt_plot_rank() averages tied ranks and ranks within each `by` key", {
  rows <- tibble::tibble(
    run = rep(c("a", "a", "b", "b"), 2), iteration = rep(1:2, 4), key = rep(c("x", "y"), each = 4),
    v = c(1, 1, 1, 2, 9, 8, 7, 6)
  )
  # Key x: three 1s share rank 2 (bin 1 of 2), the 2 has rank 4 (bin 2).
  # Key y: a has ranks 4, 3 (bin 2), b has 2, 1 (bin 1).
  p <- pt_plot_rank(rows, value = v, by = key, bins = 2)
  expect_identical(p$data$key, rep(c("x", "y"), each = 4))
  expect_identical(p$data$run, rep(c("a", "a", "b", "b"), 2))
  expect_identical(p$data$count, c(2L, 0L, 1L, 1L, 0L, 2L, 2L, 0L))
  expect_s3_class(p$facet, "FacetGrid")
  expect_error(pt_plot_rank(rows, value = v, by = key, max_series = 1), "2 series")
  expect_error(pt_plot_rank(rows, value = v, bins = 0), "bins")
})

interval_rows <- function() {
  tibble::tibble(class = c("latent_only", "observed_1"), estimate = c(0.25, 0.75), lower = c(0.1, 0.5),
                 upper = c(0.4, 0.9), level = 0.95, draws = 10L)
}

test_that("pt_plot_interval() draws one pointrange per key", {
  p <- pt_plot_interval(interval_rows())
  expect_s3_class(p, "ggplot")
  expect_identical(geoms(p), "GeomPointrange")
  d <- ggplot2::layer_data(p, 1)
  expect_identical(d$x, c(0.25, 0.75))
  expect_identical(d$xmin, c(0.1, 0.5))
  expect_identical(d$xmax, c(0.4, 0.9))
  expect_identical(as.character(p$data$label), c("latent_only", "observed_1"))
  # The first row is at the top.
  expect_true(d$y[[1]] > d$y[[2]])
})

test_that("pt_plot_interval() adds a truth line or truth points", {
  p <- pt_plot_interval(interval_rows(), truth = 0.5)
  expect_identical(geoms(p), c("GeomPointrange", "GeomVline"))
  expect_identical(ggplot2::layer_data(p, 2)$xintercept, 0.5)

  truth <- tibble::tibble(class = c("observed_1", "latent_only"), truth = c(0.8, 0.2))
  p <- pt_plot_interval(interval_rows(), truth = truth)
  expect_identical(geoms(p), c("GeomPointrange", "GeomPoint"))
  d <- ggplot2::layer_data(p, 2)
  ranges <- ggplot2::layer_data(p, 1)
  expect_identical(d$x[order(-d$y)], c(0.2, 0.8))
  expect_identical(sort(d$y), sort(ranges$y))

  expect_error(pt_plot_interval(interval_rows(), truth = tibble::tibble(cls = "x", truth = 1)), "cls")
  expect_error(pt_plot_interval(interval_rows(), truth = "a"), "truth")
  expect_error(pt_plot_interval(interval_rows()[c("class", "estimate")]), "lower")
})

test_that("pt_plot_interval() colours by run and facets several value columns", {
  rows <- tibble::tibble(run = c("a", "a", "b", "b"), value = c("u", "v", "u", "v"), estimate = 1:4,
                         lower = 0:3, upper = 2:5, level = 0.9, draws = 3L)
  p <- pt_plot_interval(rows)
  expect_s3_class(p$facet, "FacetWrap")
  d <- ggplot2::layer_data(p, 1)
  expect_identical(length(unique(d$colour)), 2L)
  expect_identical(length(unique(d$PANEL)), 2L)
  expect_identical(sort(d$x), c(1, 2, 3, 4))
  # Runs at the same key are dodged apart.
  expect_identical(length(unique(d$y[d$PANEL == 1])), 2L)
})

test_that("pt_plot_interval() takes the out-degree posteriors and pt_interval() output", {
  a <- read_golden(golden_store())
  on.exit(pt_close(a))
  dist <- pt_out_degree_dist(a)
  p <- pt_plot_interval(dist)
  expect_identical(as.character(p$data$label), dist$infection)
  expect_identical(ggplot2::layer_data(p, 1)$x, dist$estimate)

  one <- pt_interval(pt_mean_out_degree(a))
  p <- pt_plot_interval(one, truth = 1)
  expect_identical(nrow(ggplot2::layer_data(p, 1)), 1L)
})

test_that("pt_plot_run_pairs() scatters each pair of Runs' values on the y = x line", {
  rows <- tibble::tibble(run = c("a", "a", "b", "b", "c"), infection = c("1", "2", "1", "2", "1"),
                         note = "n", v = c(1, 2, 3, 4, 5), other = 0)
  p <- pt_plot_run_pairs(rows, value = v)
  expect_s3_class(p, "ggplot")
  expect_identical(geoms(p), c("GeomAbline", "GeomPoint"))
  expect_identical(unique(ggplot2::layer_data(p, 1)$slope), 1)
  expect_identical(unique(ggplot2::layer_data(p, 1)$intercept), 0)
  # Rows match on the key columns before `value`; a key missing in one Run of a pair is dropped.
  expect_identical(as.character(p$data$pair), c("a vs b", "a vs b", "a vs c", "b vs c"))
  expect_identical(p$data$x, c(1, 2, 1, 3))
  expect_identical(p$data$y, c(3, 4, 5, 5))
  expect_s3_class(p$facet, "FacetWrap")
  expect_identical(length(unique(ggplot2::layer_data(p, 2)$PANEL)), 3L)
})

test_that("pt_plot_run_pairs() takes explicit keys and needs two Runs", {
  rows <- tibble::tibble(run = rep(c("a", "b"), each = 2), locus = "L", allele = c(0, 1, 0, 1), v = 1:4)
  expect_error(pt_plot_run_pairs(rows, value = v), "more than one row per Run")
  p <- pt_plot_run_pairs(rows, value = v, by = c(locus, allele))
  expect_identical(p$data$x, c(1L, 2L))
  expect_identical(p$data$y, c(3L, 4L))
  expect_error(pt_plot_run_pairs(rows[rows$run == "a", ], value = v, by = allele), "two Runs")
})

test_that("pt_plot_run_pairs() takes a by_run summary", {
  a <- read_golden(c(golden_store("single-chain"), golden_store("replica-exchange")))
  on.exit(pt_close(a))
  marginals <- pt_edge_marginals(a, by_run = TRUE)
  p <- pt_plot_run_pairs(marginals, value = marginal)
  wide <- merge(marginals[marginals$run == "single-chain", ], marginals[marginals$run == "replica-exchange", ],
                by = c("child", "parent"))
  expect_identical(nrow(p$data), nrow(wide))
  expect_identical(sort(p$data$x), sort(wide$marginal.x))
})
