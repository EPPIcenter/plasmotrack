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
