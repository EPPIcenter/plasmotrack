# ggplot2 plots of per-Draw tables, intervals, Run agreement and a pt_convergence. ggplot2 is in Suggests: every
# plot calls check_ggplot2() first. Plots return the ggplot object, never print it.
#
#   check_ggplot2(call)               rlang::check_installed() on plot_package().
#   plot_column(x, value, arg, call) the one numeric column `value` selects (tidy-select quosure).
#   plot_keys(x, by, max_series)     the `by` columns; errors past max_series distinct keys (counted lazily).
#   plot_rows(x, columns)            collects only those columns.

# The package plots need; a function so tests can stand in a package that isn't installed.
plot_package <- function() "ggplot2"

check_ggplot2 <- function(call = rlang::caller_env()) {
  rlang::check_installed(plot_package(), reason = "to draw plots.", call = call)
}

check_draw_table <- function(x, call = rlang::caller_env()) {
  missing <- setdiff(c("run", "iteration"), colnames(x))
  if (length(missing)) {
    cli::cli_abort("{.arg x} must be a per-Draw table: it has no {.field {missing}} column{?s}.", call = call)
  }
}

plot_column <- function(x, value, arg = "value", call = rlang::caller_env()) {
  if (rlang::quo_is_missing(value) || rlang::quo_is_null(value)) {
    cli::cli_abort("{.arg {arg}} must name a numeric column of {.arg x}.", call = call)
  }
  column <- colnames(dplyr::select(x, !!value))
  if (length(column) != 1) {
    cli::cli_abort("{.arg {arg}} must select one column, not {length(column)}.", call = call)
  }
  types <- vapply(dplyr::collect(utils::head(dplyr::select(x, dplyr::all_of(column)), 0)), is.numeric, logical(1))
  if (!types[[1]]) cli::cli_abort("Column {.field {column}} must be numeric.", call = call)
  column
}

plot_keys <- function(x, by, max_series, call = rlang::caller_env()) {
  if (!is.numeric(max_series) || length(max_series) != 1 || is.na(max_series) || max_series < 1) {
    cli::cli_abort("{.arg max_series} must be a positive number.", call = call)
  }
  by <- colnames(dplyr::select(x, !!by))
  if (length(by) == 0) return(by)
  n <- nrow(dplyr::collect(dplyr::distinct(dplyr::select(x, dplyr::all_of(by)))))
  if (n > max_series) {
    cli::cli_abort(c(
      "{n} series to plot, more than {.arg max_series} = {max_series}.",
      i = "Filter {.arg x} to the keys you want first, such as with {.fn dplyr::filter}."
    ), call = call)
  }
  by
}

plot_rows <- function(x, columns) {
  dplyr::collect(dplyr::select(x, dplyr::all_of(unique(columns))))
}

check_one_per_draw <- function(rows, by, call = rlang::caller_env()) {
  if (anyDuplicated(rows[c("run", "iteration", by)])) {
    cli::cli_abort(c("{.arg x} has more than one row per Draw per {.arg by} key.",
                     i = "Give every key column in {.arg by}."), call = call)
  }
}

#' Trace plot of a per-Draw table
#'
#' One line per Run over iteration, for one value column of any per-Draw table: a functional such as
#' [pt_mean_out_degree()], a per-Draw table such as [pt_infection_draws()] with `by = infection`, or a raw Trace from
#' [pt_table()]. When the table has a `phase` column (a raw Trace read with `phase = c("burn_in", "sampling")`), the
#' burn-in iterations are shaded.
#'
#' A lazy table is collected with only the columns the plot needs. Filter it first: the plot refuses more than
#' `max_series` keys.
#'
#' @param x A per-Draw table (lazy or a data frame) with `run` and `iteration` columns.
#' @param value The numeric column to plot (tidy-select, bare or a string).
#' @param by Key columns (tidy-select); one facet per key. `x` must have one row per Draw per key.
#' @param max_series Error when `by` has more distinct keys than this.
#' @return A ggplot object, colour mapped to `run`.
#' @seealso [pt_plot_rank()] for rank histograms of the same tables.
#' @export
#' @examplesIf rlang::is_installed("ggplot2")
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_plot_trace(pt_mean_out_degree(a), value = mean_out_degree)
#' pt_table(a, "run_scalars", phase = c("burn_in", "sampling")) |>
#'   pt_plot_trace(value = likelihood)
#' pt_close(a)
pt_plot_trace <- function(x, value, by = NULL, max_series = 25) {
  check_ggplot2()
  check_draw_table(x)
  value <- plot_column(x, rlang::enquo(value))
  by <- plot_keys(x, rlang::enquo(by), max_series)
  phased <- "phase" %in% colnames(x)
  rows <- plot_rows(x, c("run", "iteration", if (phased) "phase", by, value))
  check_one_per_draw(rows, by)

  p <- ggplot2::ggplot(rows, ggplot2::aes(.data$iteration, .data[[value]], colour = .data$run, group = .data$run))
  if (phased && any(rows$phase == "burn_in")) {
    burn <- rows[rows$phase == "burn_in", c("run", "iteration")]
    shade <- dplyr::distinct(dplyr::summarise(dplyr::group_by(burn, .data$run), xmin = min(.data$iteration),
                                              xmax = max(.data$iteration)), .data$xmin, .data$xmax)
    p <- p + ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax), data = shade, ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE, fill = "grey50", alpha = 0.2
    )
  }
  p <- p + ggplot2::geom_line() + ggplot2::labs(x = "iteration", y = value, colour = "run")
  if (length(by)) p <- p + ggplot2::facet_wrap(ggplot2::vars(!!!rlang::syms(by)), scales = "free_y")
  p
}
