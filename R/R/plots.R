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

#' Rank histograms of a per-Draw table
#'
#' For one value column, pools the values of every Run (within each `by` key), ranks them (ties get their average
#' rank), splits the ranks 1 to *N* into `bins` equal bins, and draws each Run's histogram of bins. When the Runs mix
#' well, every Run's histogram is flat; a Run stuck high or low piles up at one end. See Vehtari et al. (2021).
#'
#' Takes the same tables as [pt_plot_trace()], and collects only the columns it needs.
#'
#' @inheritParams pt_plot_trace
#' @param by Key columns (tidy-select); one row of facets per key, ranked separately. `x` must have one row per Draw
#'   per key.
#' @param bins Number of rank bins.
#' @return A ggplot object of columns, one facet per Run (and key), fill mapped to `run`. Its data has the `by`
#'   columns, `run`, `bin` (1 to `bins`) and `count` (empty bins included).
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C. (2021). Rank-normalization,
#'   folding, and localization: an improved R-hat for assessing convergence of MCMC. *Bayesian Analysis* 16(2).
#' @seealso [pt_plot_trace()].
#' @export
#' @examplesIf rlang::is_installed("ggplot2")
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c("single-chain", "replica-exchange")), threads = 1)
#' pt_plot_rank(pt_mean_p_latent_parent(a), value = mean_p_latent_parent, bins = 5)
#' pt_close(a)
pt_plot_rank <- function(x, value, by = NULL, bins = 20, max_series = 25) {
  check_ggplot2()
  check_draw_table(x)
  if (!is.numeric(bins) || length(bins) != 1 || is.na(bins) || bins < 1 || bins != round(bins)) {
    cli::cli_abort("{.arg bins} must be a whole number of at least 1.")
  }
  bins <- as.integer(bins)
  value <- plot_column(x, rlang::enquo(value))
  by <- plot_keys(x, rlang::enquo(by), max_series)
  rows <- plot_rows(x, c("run", "iteration", by, value))
  check_one_per_draw(rows, by)
  if (anyNA(rows[[value]])) cli::cli_abort("Column {.field {value}} has missing values; filter them first.")

  runs <- unique(rows$run)
  counts <- dplyr::group_by(rows, dplyr::across(dplyr::all_of(by))) |>
    dplyr::group_modify(function(group, key) {
      ranks <- rank(group[[value]], ties.method = "average")
      bin <- pmin(pmax(as.integer(ceiling(ranks * bins / length(ranks))), 1L), bins)
      grid <- expand.grid(bin = seq_len(bins), run = runs, stringsAsFactors = FALSE)
      tally <- table(factor(group$run, runs), factor(bin, seq_len(bins)))
      tibble::tibble(run = grid$run, bin = grid$bin, count = as.integer(tally[cbind(match(grid$run, runs), grid$bin)]))
    }) |>
    dplyr::ungroup()

  facets <- if (length(by)) ggplot2::vars(!!!rlang::syms(by)) else NULL
  ggplot2::ggplot(counts, ggplot2::aes(.data$bin, .data$count, fill = .data$run)) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(rows = facets, cols = ggplot2::vars(.data$run)) +
    ggplot2::labs(x = paste0("rank bin of ", value), y = "Draws", fill = "run")
}

# Columns of an interval table that are statistics, not keys.
interval_statistics <- c("estimate", "lower", "upper", "level", "draws", "n", "var_within", "var_between", "p_zero")

#' Interval plot
#'
#' Pointranges from any table with `estimate`, `lower` and `upper` columns: the output of [pt_interval()], of
#' [pt_out_degree_dist()] and [pt_realized_out_degree()], or a table the user builds. One row per key, the first row
#' at the top.
#'
#' Every column other than the statistics (`estimate`, `lower`, `upper`, `level`, `draws`, `n`, `var_within`,
#' `var_between`, `p_zero`) is a key. A `run` key (from `pt_interval(by = run)`) is mapped to colour; a `value` key
#' (from `pt_interval()` on several value columns) gives one facet per value column, with free scales; the other keys
#' label the rows.
#'
#' @param x A data frame with `estimate`, `lower` and `upper`.
#' @param truth A reference: `NULL` for none, a number for one vertical line, or a data frame with a `truth` column and
#'   key columns of `x` to match rows on (such as `class` and `truth`), for one point per matched row.
#' @return A ggplot object.
#' @seealso [pt_interval()].
#' @export
#' @examplesIf rlang::is_installed("ggplot2")
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_composition(a) |> pt_interval(by = class) |> pt_plot_interval()
#' pt_mean_out_degree(a) |> pt_interval(by = run) |> pt_plot_interval(truth = 1)
#' pt_plot_interval(pt_out_degree_dist(a))
#' pt_close(a)
pt_plot_interval <- function(x, truth = NULL) {
  check_ggplot2()
  x <- dplyr::collect(x)
  missing <- setdiff(c("estimate", "lower", "upper"), names(x))
  if (length(missing)) cli::cli_abort("{.arg x} has no {.field {missing}} column{?s}.")
  keys <- setdiff(names(x), interval_statistics)
  facet <- intersect("value", keys)
  colour <- intersect("run", keys)
  labels <- setdiff(keys, c(facet, colour))

  label <- if (length(labels)) do.call(paste, c(lapply(x[labels], as.character), sep = ", ")) else rep("all", nrow(x))
  rows <- x[c(keys, "estimate", "lower", "upper")]
  rows$label <- factor(label, levels = rev(unique(label)))

  aesthetics <- ggplot2::aes(x = .data$estimate, xmin = .data$lower, xmax = .data$upper, y = .data$label)
  if (length(colour)) aesthetics$colour <- rlang::quo(.data$run)
  dodge <- if (length(colour)) ggplot2::position_dodge(width = 0.5) else "identity"
  p <- ggplot2::ggplot(rows, aesthetics) +
    ggplot2::geom_pointrange(position = dodge) +
    ggplot2::labs(x = "estimate", y = NULL)

  if (is.numeric(truth) && length(truth) == 1 && !is.na(truth)) {
    p <- p + ggplot2::geom_vline(xintercept = truth, linetype = "dashed")
  } else if (is.data.frame(truth)) {
    on <- setdiff(names(truth), "truth")
    unknown <- setdiff(on, keys)
    if (!"truth" %in% names(truth) || length(unknown)) {
      cli::cli_abort(c("{.arg truth} must have a {.field truth} column and key columns of {.arg x}.",
                       x = if (length(unknown)) "Not keys of {.arg x}: {.field {unknown}}.",
                       i = "Keys of {.arg x}: {.field {keys}}."))
    }
    points <- dplyr::inner_join(rows[unique(c(on, facet, "label"))], truth, by = on)
    p <- p + ggplot2::geom_point(ggplot2::aes(x = .data$truth, y = .data$label), data = points,
                                 inherit.aes = FALSE, shape = 4, size = 3)
  } else if (!is.null(truth)) {
    cli::cli_abort("{.arg truth} must be {.code NULL}, a number, or a data frame with a {.field truth} column.")
  }
  if (length(facet)) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$value), scales = "free")
  p
}
