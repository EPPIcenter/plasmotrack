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
#' **Strata.** With `by`, each value is a mean over one stratum's Infections instead of all of them, with a
#' `stratum` column after `iteration` (and before `class`). A stratum is a value of the epi column `by`; Infections
#' with an NA value, or missing from epi, form the NA stratum, which is kept. `pt_mean_out_degree()` counts a stratum's
#' Infections as sources (children anywhere); the others count them as recipients. So the per-stratum values weighted
#' by the stratum sizes n_S average back to the overall value at every Draw. Cross-strata come from a combined column
#' added to epi.
#'
#' @inheritParams pt_table
#' @param by `NULL`, or the name of one epi column (a string; see [pt_read()]) to stratify by.
#' @return A lazy table (`tbl_lazy`) with `run` and `iteration` first, then `stratum` with `by`, then the value
#'   columns above.
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
#'
#' epi <- data.frame(infection = c("1", "2"), district = c("north", "south"))
#' a <- pt_read(store, epi = epi, threads = 1)
#' pt_mean_out_degree(a, by = "district") |> pt_interval(by = stratum)
#' pt_composition(a, by = "district") |> pt_interval(by = c(stratum, class))
#' pt_close(a)
pt_mean_out_degree <- function(a, by = NULL, runs = NULL, phase = "sampling") {
  check_by(a, by)
  pt_out_degree_draws(a, runs, phase) |>
    dplyr::rename(mean_out_degree = "out_degree") |>
    with_strata(a, by) |>
    per_draw_mean("mean_out_degree", keys = strata_keys(by))
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
pt_mean_p_latent_parent <- function(a, by = NULL, runs = NULL, phase = "sampling") {
  check_by(a, by)
  infection_values(a, runs, phase) |>
    dplyr::select("run", "iteration", "infection", mean_p_latent_parent = "p_latent_parent") |>
    with_strata(a, by) |>
    per_draw_mean("mean_p_latent_parent", keys = strata_keys(by))
}

#' @rdname pt_mean_out_degree
#' @export
pt_composition <- function(a, by = NULL, runs = NULL, phase = "sampling") {
  check_by(a, by)
  pt_class_draws(a, runs, phase) |>
    with_strata(a, by) |>
    per_draw_mean("share", keys = c(strata_keys(by), "class"))
}

#' @rdname pt_mean_out_degree
#' @export
pt_multi_source <- function(a, by = NULL, runs = NULL, phase = "sampling") {
  check_by(a, by)
  infection_values(a, runs, phase) |>
    dplyr::select("run", "iteration", "infection", augmented = "multi_source_augmented",
                  observed_only = "multi_source_observed") |>
    with_strata(a, by) |>
    per_draw_mean(c("augmented", "observed_only"), keys = strata_keys(by))
}

#' @rdname pt_mean_out_degree
#' @export
pt_entropy <- function(a, by = NULL, runs = NULL, phase = "sampling") {
  check_by(a, by)
  infection_values(a, runs, phase) |>
    dplyr::select("run", "iteration", "infection", mean_entropy = "entropy",
                  mean_effective_size = "effective_size") |>
    with_strata(a, by) |>
    per_draw_mean(c("mean_entropy", "mean_effective_size"), keys = strata_keys(by))
}

#' Posterior mean and interval of a per-Draw table
#'
#' Turns any per-Draw table into estimates with intervals: the step from lazy per-Draw values to an eager summary.
#' The input needs `run` and `iteration` columns (the Draw) and numeric value columns. It can be a functional such
#' as [pt_mean_out_degree()], a per-Draw table such as [pt_infection_draws()] with `by = infection`, a raw Trace from
#' [pt_table()], or a table the user builds, lazy or in memory.
#'
#' Within each `by` group, the estimate is the mean over every row, pooled over the Draws of every Run, so each Run
#' weighs by its Draws. The interval is equal-tailed: the type-7 quantiles ([stats::quantile()]'s default) at
#' \eqn{(1 - level) / 2} and \eqn{(1 + level) / 2} of the same pooled values. It covers the variation between Draws
#' only (see the spec's out-degree posteriors for the variance within a Draw).
#'
#' Each group must have at most one row per Draw: give every key column in `by`. Value columns may not hold missing
#' values; filter them first.
#'
#' @param x A per-Draw table (lazy or a data frame) with `run`, `iteration` and numeric value columns.
#' @param by Key columns (tidy-select, such as `class` or `c(locus, allele)`); one row per key in the result.
#'   `NULL` summarises the whole table. `by = run` gives one interval per Run.
#' @param value Value columns (tidy-select). `NULL` takes every numeric column other than `run`, `iteration`,
#'   `window` and the `by` columns.
#' @param level Interval probability, between 0 and 1.
#' @return A tibble with the `by` columns (sorted, strings in C order), then `value` (the value column's name, only
#'   when there are several), `estimate`, `lower`, `upper`, `level` and `draws` (the number of Draws summarised).
#' @seealso [pt_mean_out_degree()] and the other functionals.
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_mean_out_degree(a) |> pt_interval()
#' pt_composition(a) |> pt_interval(by = class, level = 0.9)
#' pt_table(a, "allele_frequencies") |> pt_interval(by = c(locus, allele))
#' pt_infection_draws(a) |> pt_interval(by = infection, value = p_latent_parent)
#' pt_close(a)
pt_interval <- function(x, by = NULL, value = NULL, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be a number strictly between 0 and 1.")
  }
  missing <- setdiff(c("run", "iteration"), colnames(x))
  if (length(missing)) {
    cli::cli_abort("{.arg x} must be a per-Draw table: it has no {.field {missing}} column{?s}.")
  }
  by <- colnames(dplyr::select(x, {{ by }}))
  types <- vapply(dplyr::collect(utils::head(x, 0)), is.numeric, logical(1))
  value_quo <- rlang::enquo(value)
  if (rlang::quo_is_null(value_quo)) {
    value <- setdiff(names(types)[types], c("run", "iteration", "window", by))
    if (length(value) == 0) {
      cli::cli_abort("{.arg x} has no numeric value column besides {.field {c('run', 'iteration', by)}}.")
    }
  } else {
    value <- colnames(dplyr::select(x, !!value_quo))
    if (!all(types[value])) {
      cli::cli_abort("Value column{?s} {.field {value[!types[value]]}} must be numeric.")
    }
  }

  rows <- x |>
    dplyr::select(dplyr::all_of(unique(c("run", "iteration", by, value)))) |>
    dplyr::collect()
  bad <- value[vapply(value, function(v) anyNA(rows[[v]]), logical(1))]
  if (length(bad)) {
    cli::cli_abort("Value column{?s} {.field {bad}} {?has/have} missing values; filter them first.")
  }
  if (anyDuplicated(rows[unique(c("run", "iteration", by))])) {
    cli::cli_abort(c("{.arg x} has more than one row per Draw per {.arg by} key.",
                     i = "Give every key column in {.arg by}."))
  }

  # Locals are injected with !!: a value or key column may share their names.
  probs <- c((1 - level) / 2, (1 + level) / 2)
  grouped <- dplyr::group_by(rows, dplyr::across(dplyr::all_of(by)))
  out <- lapply(value, function(v) {
    dplyr::summarise(
      grouped,
      value = !!v,
      estimate = mean(.data[[!!v]]),
      lower = stats::quantile(.data[[!!v]], !!probs[[1]], type = 7, names = FALSE),
      upper = stats::quantile(.data[[!!v]], !!probs[[2]], type = 7, names = FALSE),
      level = !!level,
      draws = dplyr::n(),
      .groups = "drop"
    )
  })
  out <- dplyr::arrange(dplyr::bind_rows(out), dplyr::across(dplyr::all_of(by)))
  if (length(value) == 1) out$value <- NULL
  out
}
