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
