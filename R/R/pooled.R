# Pooled summaries: eager tibbles averaged over the sampling Draws of every Run (Conventions 1-3).

#' Edge marginals
#'
#' The posterior probability that one Infection infected another (the manuscript's edge probability, whose mean
#' over Infections is ψ_edge): the total weight of the child's Parent sets that hold the parent, averaged over
#' sampling Draws. The Latent parent is never an edge; its share is the Latent parent probability (see
#' [pt_infection_draws()]).
#'
#' A Run's Draws are its sampling iterations (`draws` in [pt_runs()]). An edge with no weight at a Draw counts 0
#' there. Pooled marginals sum over every Draw of every Run and divide by the total Draw count, so each Run weighs by
#' its Draws. Runs with no sampling Draws are skipped with a message.
#'
#' Weights are used as stored: the Store drops Parent sets below the Run's `parent_set_floor` (default 1e-9) without
#' renormalising, so a marginal can be low by at most the floor times the number of the child's dropped sets.
#'
#' @inheritParams pt_table
#' @param floor Keep edges whose marginal is at least `floor`. With `by_run = TRUE`, keep edges whose marginal in
#'   some Run is at least `floor`, with every Run's row. The default 0 keeps every edge with weight at any Draw.
#' @param by_run `FALSE` for marginals pooled over Runs; `TRUE` for one row per Run per edge, where an edge seen in
#'   any Run has a row (0 if need be) in every Run.
#' @return A tibble of `child`, `parent`, `marginal` and `lag`, after `run` when `by_run = TRUE`, sorted by `run`
#'   (in [pt_runs()] order), then `child` and `parent` as strings. `lag` is the child's observation time minus the
#'   parent's, from [pt_infections()]; it is `NA` where either time is missing.
#' @seealso [pt_edge_draws()] for the per-Draw edge weights.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_edge_marginals(a)
#' pt_edge_marginals(a, by_run = TRUE)
#' pt_close(a)
pt_edge_marginals <- function(a, floor = 0, by_run = FALSE, runs = NULL) {
  if (!is.numeric(floor) || length(floor) != 1 || is.na(floor) || floor < 0 || floor > 1) {
    cli::cli_abort("{.arg floor} must be a number between 0 and 1.")
  }
  check_flag(by_run)
  runs <- runs_with_draws(a, runs)
  totals <- observed_edges(a, runs, "sampling") |>
    dplyr::group_by(.data$run, .data$child, .data$parent) |>
    dplyr::summarise(total = sum(.data$prob, na.rm = TRUE), .groups = "drop") |>
    dplyr::collect()
  draws <- stats::setNames(as.numeric(a$runs$draws), a$runs$run)[runs]

  edges <- unique(totals[c("child", "parent")])
  edges <- edges[order(edges$child, edges$parent, method = "radix"), ]
  # Edge x Run totals, 0 where a Run lacks the edge.
  total <- matrix(0, nrow(edges), length(runs), dimnames = list(NULL, runs))
  total[cbind(match(paste(totals$child, totals$parent, sep = "\r"), paste(edges$child, edges$parent, sep = "\r")),
              match(totals$run, runs))] <- totals$total

  lag <- edge_lag(a, edges$child, edges$parent)
  if (by_run) {
    marginal <- sweep(total, 2, draws, "/")
    keep <- if (length(runs)) apply(marginal, 1, max) >= floor else logical()
    keep[is.na(keep)] <- FALSE
    k <- sum(keep)
    tibble::tibble(
      run = rep(runs, each = k),
      child = rep(edges$child[keep], length(runs)),
      parent = rep(edges$parent[keep], length(runs)),
      marginal = as.vector(marginal[keep, , drop = FALSE]),
      lag = rep(lag[keep], length(runs))
    )
  } else {
    marginal <- rowSums(total) / sum(draws)
    keep <- marginal >= floor
    tibble::tibble(child = edges$child[keep], parent = edges$parent[keep], marginal = marginal[keep], lag = lag[keep])
  }
}

# Child's observation time minus the parent's; NA where either is missing.
edge_lag <- function(a, child, parent) {
  time <- stats::setNames(a$infections$observation_time, a$infections$infection)
  unname(as.numeric(time[child]) - as.numeric(time[parent]))
}
