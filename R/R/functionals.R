# Functionals: one value per Draw, a mean over Infections of a per-Draw value. Lazy, like the per-Draw tables.
#
#   per_draw_mean(values, columns, keys)  lazy: mean of `columns` over the rows of each Draw (and `keys`). Every
#                                         caller hands it a full Infection grid, so the mean is over every Infection.
#   infection_grid(a, runs, phase)        lazy: run, iteration, infection: every Infection at every Draw.

per_draw_mean <- function(values, columns, keys = character()) {
  values |>
    dplyr::group_by(.data$run, .data$iteration, dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(columns), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

infection_grid <- function(a, runs = NULL, phase = "sampling") {
  draw_grid(a, runs, phase) |>
    dplyr::cross_join(dplyr::select(structure_infections(a), "infection"))
}

#' Functionals
#'
#' A Functional is one number per Draw computed from the whole network, here a mean over Infections. Each function
#' returns a lazy table with one row per Draw (`run`, `iteration`, then the value columns); [pt_interval()] turns it
#' into a posterior mean with an interval.
#'
#' Every Infection counts at every Draw, and the mean divides by the number of Infections. An Infection with no
#' Parent set row at a Draw (which an engine Store never writes) counts 0, with entropy 0 and effective size 1.
#' Weights are used as stored: see [pt_infection_draws()] for the per-Infection values and the floor.
#'
#' * `pt_mean_out_degree()`: `mean_out_degree`, the mean of [pt_out_degree_draws()] over Infections, which equals
#'   the mean number of observed parents per Infection. The manuscript's ψ_edge. Every Infection counts as a source.
#' * `pt_mean_p_latent_parent()`: `mean_p_latent_parent`, the mean Latent parent probability (the manuscript's mean
#'   q, the probability of an unobserved source).
#' * `pt_composition()`: long, one row per Draw per composition class: `class`, `share`, the mean over Infections of
#'   the class's weight (see [pt_class_draws()] for the classes). At each Draw the shares sum to the mean total weight
#'   of an Infection's Parent sets, 1 less what the floor dropped.
#' * `pt_multi_source()`: `augmented`, the mean probability of more than one parent counting the Latent parent, and
#'   `observed_only`, of more than one observed parent. The manuscript's multi-source probability.
#' * `pt_entropy()`: `mean_entropy`, the mean Parent-set entropy \eqn{-\sum w \ln w}, and `mean_effective_size`, the
#'   mean of its exponential (not the exponential of the mean).
#'
#' @inheritParams pt_table
#' @return A lazy table (`tbl_lazy`) with `run` and `iteration` first, then the value columns above.
#' @seealso [pt_interval()] for the posterior mean and interval; [pt_infection_draws()] for per-Infection values.
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_mean_out_degree(a) |> dplyr::collect()
#' pt_mean_out_degree(a) |> pt_interval()
#' pt_composition(a) |> pt_interval(by = class, level = 0.9)
#' pt_entropy(a) |> pt_interval()
#' pt_close(a)
pt_mean_out_degree <- function(a, runs = NULL, phase = "sampling") {
  pt_out_degree_draws(a, runs, phase) |>
    dplyr::rename(mean_out_degree = "out_degree") |>
    per_draw_mean("mean_out_degree")
}

# pt_infection_draws() over the full Infection grid: an Infection with no row at a Draw counts 0, entropy 0 (so
# effective size 1).
infection_values <- function(a, runs, phase) {
  infection_grid(a, runs, phase) |>
    dplyr::left_join(pt_infection_draws(a, runs, phase), by = c("run", "iteration", "infection")) |>
    dplyr::mutate(
      p_latent_parent = dplyr::coalesce(.data$p_latent_parent, 0),
      entropy = dplyr::coalesce(.data$entropy, 0),
      effective_size = dplyr::coalesce(.data$effective_size, 1),
      multi_source_augmented = dplyr::coalesce(.data$multi_source_augmented, 0),
      multi_source_observed = dplyr::coalesce(.data$multi_source_observed, 0)
    )
}

#' @rdname pt_mean_out_degree
#' @export
pt_mean_p_latent_parent <- function(a, runs = NULL, phase = "sampling") {
  infection_values(a, runs, phase) |>
    dplyr::select("run", "iteration", "infection", mean_p_latent_parent = "p_latent_parent") |>
    per_draw_mean("mean_p_latent_parent")
}

#' @rdname pt_mean_out_degree
#' @export
pt_composition <- function(a, runs = NULL, phase = "sampling") {
  pt_class_draws(a, runs, phase) |>
    per_draw_mean("share", keys = "class")
}

#' @rdname pt_mean_out_degree
#' @export
pt_multi_source <- function(a, runs = NULL, phase = "sampling") {
  infection_values(a, runs, phase) |>
    dplyr::select("run", "iteration", "infection", augmented = "multi_source_augmented",
                  observed_only = "multi_source_observed") |>
    per_draw_mean(c("augmented", "observed_only"))
}

#' @rdname pt_mean_out_degree
#' @export
pt_entropy <- function(a, runs = NULL, phase = "sampling") {
  infection_values(a, runs, phase) |>
    dplyr::select("run", "iteration", "infection", mean_entropy = "entropy",
                  mean_effective_size = "effective_size") |>
    per_draw_mean(c("mean_entropy", "mean_effective_size"))
}
