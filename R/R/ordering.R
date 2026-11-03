# Infection times and the Ordering: from the Store's `infections` (infection_duration) and `ordering` Traces.

ordering_rows <- function(a, runs, phase) {
  pt_table(a, "ordering", runs = runs, phase = phase) |>
    dplyr::select("run", "iteration", "infection", "rank")
}

# Ordered pairs of a lazy Ordering table: one row per Draw and pair i < j (as strings) with `before`. With
# `infections`, only pairs with at least one of them, built by joining those Infections' rows to every row, so the
# O(n^2) join is never formed in full.
order_pairs_lazy <- function(a, runs, phase, infections = NULL) {
  rows <- ordering_rows(a, runs, phase)
  x <- rows
  if (!is.null(infections)) x <- dplyr::filter(x, .data$infection %in% !!infections)
  x <- dplyr::rename(x, x_infection = "infection", x_rank = "rank")
  y <- dplyr::rename(rows, y_infection = "infection", y_rank = "rank")
  pairs <- dplyr::inner_join(x, y, by = c("run", "iteration"))
  pairs <- if (is.null(infections)) {
    dplyr::filter(pairs, .data$x_infection < .data$y_infection)
  } else {
    # Each pair once: from the smaller id's side, or from x's side when y is not among `infections`.
    dplyr::filter(pairs, .data$x_infection != .data$y_infection &
                    (.data$x_infection < .data$y_infection | !(.data$y_infection %in% !!infections)))
  }
  pairs |>
    dplyr::mutate(
      i = dplyr::if_else(.data$x_infection < .data$y_infection, .data$x_infection, .data$y_infection),
      j = dplyr::if_else(.data$x_infection < .data$y_infection, .data$y_infection, .data$x_infection),
      before = (.data$x_infection < .data$y_infection) == (.data$x_rank < .data$y_rank)
    ) |>
    dplyr::select("run", "iteration", "i", "j", "before")
}

#' Per-Draw Ordering of Infection pairs
#'
#' A lazy table with one row per Draw (`run`, `iteration`) per pair of Infections, from the Store's `ordering`
#' Trace: the order in which the model has the Infections acquired. `i` is the smaller id and `j` the larger, as
#' strings (`"10"` before `"2"`), and `before` is `TRUE` when `i` was acquired before `j` at that Draw.
#'
#' The table has n(n - 1) / 2 rows per Draw for n Infections: filter it (by `i`, `j` or `iteration`) before
#' [dplyr::collect()]. For the probability that one Infection preceded another, use [pt_order_pairs()].
#'
#' @inheritParams pt_table
#' @return A lazy table (`tbl_lazy`) of `run`, `iteration`, `i`, `j` and `before`.
#' @seealso [pt_order_pairs()] for the pooled probabilities, [pt_infection_times()] for Infection times.
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_order_pair_draws(a) |> dplyr::filter(i == "1") |> dplyr::collect()
#' pt_close(a)
pt_order_pair_draws <- function(a, runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  order_pairs_lazy(a, runs, phase)
}

#' Pooled Ordering of Infection pairs
#'
#' For each pair of Infections, the posterior probability that one was acquired before the other: the share of
#' sampling Draws whose Ordering (the Store's `ordering` Trace) puts `i` before `j`. `i` is the smaller id and `j` the
#' larger, as strings (`"10"` before `"2"`); the probability that `j` came first is `1 - p_before`.
#'
#' Draws pool as everywhere in the package: counts sum over Runs and divide by the total Draw count, so each Run
#' weighs by its Draws. Runs with no sampling Draws are skipped with a message. With n Infections there are
#' n(n - 1) / 2 pairs; `infections` restricts the work, not just the output, to pairs with at least one of them.
#'
#' @inheritParams pt_coi
#' @param infections Infection ids (character). Keep only the pairs with at least one of them. `NULL`, the default,
#'   keeps every pair.
#' @return A tibble of `i`, `j`, `p_before` and `draws` (the Draws summarised), after `run` when `by_run = TRUE`,
#'   sorted by `run` (in [pt_runs()] order), then `i` and `j` as strings.
#' @seealso [pt_order_pair_draws()] for the per-Draw pairs, [pt_infection_times()] for Infection times.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_order_pairs(a)
#' pt_order_pairs(a, infections = "2", by_run = TRUE)
#' pt_close(a)
pt_order_pairs <- function(a, infections = NULL, by_run = FALSE, runs = NULL) {
  check_flag(by_run)
  infections <- check_infection_subset(a, infections)
  runs <- runs_with_draws(a, runs)
  groups <- c(if (by_run) "run", "i", "j")
  counts <- order_pairs_lazy(a, runs, "sampling", infections) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
    dplyr::summarise(n_before = sum(dplyr::if_else(.data$before, 1, 0), na.rm = TRUE), .groups = "drop") |>
    dplyr::collect()
  draws <- stats::setNames(as.numeric(a$runs$draws), a$runs$run)
  counts$draws <- if (by_run) unname(draws[counts$run]) else rep(sum(draws[runs]), nrow(counts))
  counts$p_before <- as.numeric(counts$n_before) / counts$draws
  order_by <- lapply(groups, function(g) if (g == "run") match(counts$run, a$runs$run) else counts[[g]])
  counts <- counts[do.call(order, c(unname(order_by), method = "radix")), ]
  tibble::as_tibble(counts[c(groups, "p_before", "draws")])
}

# NULL, or `infections` as unique character ids of the Analysis; unknown ids are an error.
check_infection_subset <- function(a, infections, arg = rlang::caller_arg(infections), call = rlang::caller_env()) {
  if (is.null(infections)) return(NULL)
  if (!is.character(infections) || anyNA(infections)) {
    cli::cli_abort("{.arg {arg}} must be a character vector of Infection ids.", call = call)
  }
  unknown <- setdiff(infections, a$infections$infection)
  if (length(unknown)) {
    cli::cli_abort("Unknown Infection{?s}: {.val {unknown}}.", call = call)
  }
  unique(infections)
}

#' Infection times
#'
#' When each Infection was acquired: its observation time minus the infection duration (the Store's
#' `infection_duration`, the time from infection to observation in the time unit of the IDP files), summarised over
#' the sampling Draws. `median` is the observation time minus the median duration, and `lower` and `upper` bound the
#' equal-tailed interval: the observation time minus the type-7 quantiles (R's default, see [stats::quantile()]) of
#' the pooled durations at \eqn{(1 + level) / 2} and \eqn{(1 - level) / 2}.
#'
#' Observation times come from the input JSON (see [pt_infections()]); an Infection without one, or every Infection
#' when [pt_read()] found no input, gets `NA` times. Draws pool as everywhere in the package: the quantiles are over
#' every Draw of every Run, so each Run weighs by its Draws. Runs with no sampling Draws are skipped with a message.
#'
#' @inheritParams pt_coi
#' @param level Interval probability, strictly between 0 and 1.
#' @return A tibble of `infection`, `lower`, `median` and `upper`, after `run` when `by_run = TRUE`, sorted by `run`
#'   (in [pt_runs()] order), then `infection` as strings.
#' @seealso [pt_order_pairs()] for the probability that one Infection was acquired before another;
#'   `pt_table(a, "infections") |> pt_interval(by = infection, value = infection_duration)` for the durations.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' input <- system.file("extdata", "example_input.json", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), input = input, threads = 1)
#' pt_infection_times(a)
#' pt_infection_times(a, level = 0.5, by_run = TRUE)
#' pt_close(a)
pt_infection_times <- function(a, level = 0.95, by_run = FALSE, runs = NULL) {
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be a number strictly between 0 and 1.")
  }
  check_flag(by_run)
  runs <- runs_with_draws(a, runs)
  durations <- pt_table(a, "infections", runs = runs) |>
    dplyr::select("run", "iteration", "infection", duration = "infection_duration")
  # Distinct durations and their Draw counts: a duration often stays put between Draws.
  counts <- tally_draws(a, durations, c("infection", "duration"), by_run, runs)
  keys <- c(if (by_run) "run", "infection")
  group <- do.call(paste, c(unname(as.list(counts[keys])), sep = "\r"))
  rows <- split(seq_len(nrow(counts)), factor(group, levels = unique(group)))
  probs <- c((1 + level) / 2, 0.5, (1 - level) / 2)
  quantiles <- vapply(rows, function(i) {
    stats::quantile(rep(counts$duration[i], counts$n[i]), probs, names = FALSE, type = 7)
  }, numeric(3), USE.NAMES = FALSE)
  out <- counts[!duplicated(group), keys]
  time <- a$infections$observation_time[match(out$infection, a$infections$infection)]
  out$lower <- time - quantiles[1, ]
  out$median <- time - quantiles[2, ]
  out$upper <- time - quantiles[3, ]
  tibble::as_tibble(out)
}
