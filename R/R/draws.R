# Per-Draw Parent-set quantities: lazy tables built on pt_table(), never collected here.
#
# Building blocks shared by the per-Draw tables, the functionals and the pooled summaries:
#   structure_infections(a)   lazy: infection, latent_parent (structure/, the same in every Run)
#   set_composition(a, runs)  lazy: run, infection, set_id, n_observed, has_latent_parent, class
#   weighted_sets(a, runs, phase)  lazy: run, iteration, infection, set_id, prob + the composition columns
#   observed_edges(a, runs, phase) lazy: run, iteration, child, parent, prob (summed over the child's sets)
#   draw_grid(a, runs, phase) lazy: run, iteration, one row per Draw (from run_scalars: Convention 1)

structure_infections <- function(a) {
  con <- analysis_con(a)
  file <- file.path(a$runs$path[[1]], "structure", "infections.parquet")
  dplyr::tbl(con, dbplyr::sql(sprintf(
    "SELECT infection, latent_parent FROM read_parquet(%s)", sql_string(con, file)
  )))
}

# Each Parent set's composition: observed members are those other than the Infection's Latent parent. A set's
# members are read over every closed Window of both Phases (pt_table() does that for parent_set_members).
set_composition <- function(a, runs = NULL) {
  members <- pt_table(a, "parent_set_members", runs = runs)
  members |>
    dplyr::inner_join(structure_infections(a), by = "infection") |>
    dplyr::group_by(.data$run, .data$infection, .data$set_id) |>
    dplyr::summarise(
      n_observed = sum(as.integer(.data$member != .data$latent_parent), na.rm = TRUE),
      has_latent_parent = sum(as.integer(.data$member == .data$latent_parent), na.rm = TRUE) > 0L,
      .groups = "drop"
    ) |>
    dplyr::mutate(class = composition_class_sql)
}

# Composition class names: latent_only, observed_1, observed_1_latent, observed_2, ... (open-ended). A set with no
# member at all cannot occur in a Store; it would be no_parent. Raw SQL: dbplyr cannot translate an R helper.
composition_class_sql <- dbplyr::sql(paste(
  "CASE WHEN n_observed = 0 AND has_latent_parent THEN 'latent_only'",
  "WHEN n_observed = 0 THEN 'no_parent'",
  "WHEN has_latent_parent THEN 'observed_' || CAST(n_observed AS VARCHAR) || '_latent'",
  "ELSE 'observed_' || CAST(n_observed AS VARCHAR) END"
))

weighted_sets <- function(a, runs = NULL, phase = "sampling") {
  pt_table(a, "parent_sets", runs = runs, phase = phase) |>
    dplyr::select("run", "iteration", "infection", "set_id", "prob") |>
    dplyr::inner_join(set_composition(a, runs), by = c("run", "infection", "set_id"))
}

observed_edges <- function(a, runs = NULL, phase = "sampling") {
  members <- pt_table(a, "parent_set_members", runs = runs) |>
    dplyr::inner_join(structure_infections(a), by = "infection") |>
    dplyr::filter(.data$member != .data$latent_parent) |>
    dplyr::select("run", "infection", "set_id", parent = "member")
  pt_table(a, "parent_sets", runs = runs, phase = phase) |>
    dplyr::select("run", "iteration", "infection", "set_id", "prob") |>
    dplyr::inner_join(members, by = c("run", "infection", "set_id")) |>
    dplyr::group_by(.data$run, .data$iteration, child = .data$infection, .data$parent) |>
    dplyr::summarise(prob = sum(.data$prob, na.rm = TRUE), .groups = "drop")
}

draw_grid <- function(a, runs = NULL, phase = "sampling") {
  pt_table(a, "run_scalars", runs = runs, phase = phase) |>
    dplyr::distinct(.data$run, .data$iteration)
}

#' Per-Draw Parent-set tables
#'
#' Lazy tables with one row per Draw (`run`, `iteration`) per key, computed from the Store's `parent_sets` and
#' `parent_set_members`. Filter them before [dplyr::collect()]; reuse one with [dplyr::compute()].
#'
#' At each Draw, an Infection's Parent sets carry the stored weights `w` (the `prob` column of `parent_sets`). The
#' Store drops sets below the Run's `parent_set_floor` (default 1e-9) without renormalising, and these tables use
#' the weights as stored, so an Infection's weights may sum to slightly less than 1: by at most the floor times the
#' number of its dropped sets. A set with no row at a Draw has weight 0 there.
#'
#' The Latent parent is never an edge: it never appears as a `parent`, and its `<id>_copy` id never appears in any
#' output. The manuscript calls it the unobserved source: its share of an Infection's Parent sets, q_i there, is
#' `p_latent_parent` here.
#'
#' * `pt_edge_draws()`: `child`, `parent`, `prob`, where `prob` is the total weight of the child's Parent sets that
#'   hold `parent`. Observed parents only, and only edges with weight at the Draw (no zero rows).
#' * `pt_infection_draws()`, one row per Infection with Parent sets at the Draw:
#'   * `p_latent_parent` \eqn{= \sum w [\text{has Latent parent}]}, the Latent parent probability (q_i).
#'   * `expected_observed_parents` \eqn{= \sum w \, n_{obs}}.
#'   * `entropy` \eqn{= -\sum w \ln w} over \eqn{w > 0}, and `effective_size` \eqn{= \exp(\text{entropy})}.
#'   * `multi_source_augmented` \eqn{= \sum w [n_{obs} + \text{has Latent parent} > 1]}: more than one parent,
#'     counting the Latent parent.
#'   * `multi_source_observed` \eqn{= \sum w [n_{obs} > 1]}: more than one observed parent.
#' * `pt_class_draws()`: `class`, `n_observed`, `has_latent_parent`, `share` \eqn{= \sum w [\text{class}]}. Composition
#'   classes are named by their observed parents and the Latent parent: `latent_only` (the manuscript's "unobserved
#'   source only"), `observed_1`, `observed_1_latent`, `observed_2`, and so on. The classes are those of any Parent
#'   set with weight in the Runs and Phases read, and every one of them has a row for every Infection at every Draw,
#'   0 when the Infection has no set of that class.
#' * `pt_out_degree_draws()`: `infection`, `out_degree`, the sum over every child of the probability, at the Draw, that
#'   the Infection is its parent (the per-Draw Edge marginals). A full grid: every Infection at every Draw, 0 when it is no one's parent.
#'
#' The Draws of a Run are the iterations of its `run_scalars` rows in the Phases read.
#'
#' @inheritParams pt_table
#' @return A lazy table (`tbl_lazy`) with `run` and `iteration` first, then the columns above.
#' @seealso [pt_edge_marginals()] for Edge marginals pooled over Draws; [pt_table()] for the raw tables.
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_edge_draws(a) |> dplyr::filter(child == "3") |> dplyr::collect()
#' pt_infection_draws(a) |> dplyr::collect()
#' pt_class_draws(a) |> dplyr::filter(share > 0) |> dplyr::collect()
#' pt_out_degree_draws(a) |> dplyr::collect()
#' pt_close(a)
pt_edge_draws <- function(a, runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  observed_edges(a, runs, phase) |>
    dplyr::select("run", "iteration", "child", "parent", "prob")
}

#' @rdname pt_edge_draws
#' @export
pt_infection_draws <- function(a, runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  weighted_sets(a, runs, phase) |>
    dplyr::mutate(
      latent = as.numeric(.data$has_latent_parent),
      w_log_w = dplyr::if_else(.data$prob > 0, .data$prob * log(.data$prob), 0)
    ) |>
    dplyr::group_by(.data$run, .data$iteration, .data$infection) |>
    dplyr::summarise(
      p_latent_parent = sum(.data$prob * .data$latent, na.rm = TRUE),
      expected_observed_parents = sum(.data$prob * .data$n_observed, na.rm = TRUE),
      entropy = -sum(.data$w_log_w, na.rm = TRUE),
      multi_source_augmented = sum(dplyr::if_else(.data$n_observed + .data$latent > 1, .data$prob, 0), na.rm = TRUE),
      multi_source_observed = sum(dplyr::if_else(.data$n_observed > 1, .data$prob, 0), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(effective_size = exp(.data$entropy)) |>
    dplyr::select("run", "iteration", "infection", "p_latent_parent", "expected_observed_parents", "entropy",
                  "effective_size", "multi_source_augmented", "multi_source_observed")
}

#' @rdname pt_edge_draws
#' @export
pt_class_draws <- function(a, runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  sets <- weighted_sets(a, runs, phase)
  shares <- sets |>
    dplyr::group_by(.data$run, .data$iteration, .data$infection, .data$class) |>
    dplyr::summarise(share = sum(.data$prob, na.rm = TRUE), .groups = "drop")
  classes <- sets |> dplyr::distinct(.data$class, .data$n_observed, .data$has_latent_parent)
  draw_grid(a, runs, phase) |>
    dplyr::cross_join(dplyr::select(structure_infections(a), "infection")) |>
    dplyr::cross_join(classes) |>
    dplyr::left_join(shares, by = c("run", "iteration", "infection", "class")) |>
    dplyr::mutate(share = dplyr::coalesce(.data$share, 0)) |>
    dplyr::select("run", "iteration", "infection", "class", "n_observed", "has_latent_parent", "share")
}

#' @rdname pt_edge_draws
#' @export
pt_out_degree_draws <- function(a, runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  degrees <- observed_edges(a, runs, phase) |>
    dplyr::group_by(.data$run, .data$iteration, infection = .data$parent) |>
    dplyr::summarise(out_degree = sum(.data$prob, na.rm = TRUE), .groups = "drop")
  draw_grid(a, runs, phase) |>
    dplyr::cross_join(dplyr::select(structure_infections(a), "infection")) |>
    dplyr::left_join(degrees, by = c("run", "iteration", "infection")) |>
    dplyr::mutate(out_degree = dplyr::coalesce(.data$out_degree, 0)) |>
    dplyr::select("run", "iteration", "infection", "out_degree")
}
