# The fit summary: one row per fit, in the column layout of the thesis simulations' fit summary.

# The thesis factor columns `factors` replaces, typed, NA when `factors` is NULL.
fit_summary_factors <- tibble::tibble(
  run_key = NA_character_, simulation_model = NA_character_, total_nodes = NA_real_,
  founder_probability = NA_real_, superinfection_probability = NA_real_, secondary_case_rate = NA_real_,
  founder_coi_rate = NA_real_, asymptomatic_probability = NA_real_, panel_name = NA_character_,
  panel_alleles_per_locus = NA_real_, simulation_index = NA_integer_, seed = NA_integer_
)

# Functional columns: name stem -> the per-Draw value it summarises (see fit_summary_functionals()).
fit_summary_functional_names <- c(
  estimated_average_outdegree = "mean_out_degree",
  estimated_unobserved_source_proportion = "mean_p_latent_parent",
  mean_observed_only_probability = "observed_only",
  mean_mixed_observed_unobserved_probability = "mixed",
  mean_unobserved_only_probability = "unobserved_only",
  mean_no_parent_probability = "no_parent",
  mean_multi_source_probability_augmented = "augmented",
  mean_multi_source_probability_observed_only = "multi_observed_only",
  mean_parent_set_entropy = "mean_entropy",
  mean_effective_parent_sets = "mean_effective_size"
)

#' Summarise a simulation fit in one row
#'
#' Scores one fit of a simulation against its truth and returns one row in the column layout of the thesis
#' simulations' fit summary (70 columns, names, order and types unchanged), so the thesis table and figure scripts
#' read it as is. There is no cohort function: `dplyr::bind_rows()` the rows of every fit.
#'
#' **Columns**, in order:
#' * **Factors** (12): `run_key` to `seed`, from `factors`; then `missingness_condition` and `fit_type`.
#' * **Truth**: `n_true_cases` (cases in the truth), `n_observed_cases` (the Analysis's Infections),
#'   `true_average_outdegree` (true edges among the Infections over their number), `true_unobserved_source_proportion`
#'   (the share of Infections whose true Parent set holds the Latent parent; see [pt_score_parent_sets()]).
#' * **Functionals**, each a posterior mean with a 95% interval (`_lower`, `_upper`; see [pt_interval()]), over the
#'   sampling Draws of every Run: `estimated_average_outdegree` ([pt_mean_out_degree()]),
#'   `estimated_unobserved_source_proportion` ([pt_mean_p_latent_parent()]); the composition rolled up from
#'   [pt_composition()]: `mean_observed_only_probability` (the `observed_k` classes), `mean_mixed_observed_unobserved_
#'   probability` (`observed_k_latent`), `mean_unobserved_only_probability` (`latent_only`) and
#'   `mean_no_parent_probability`, always 0 since a Parent set is never empty; `mean_multi_source_probability_augmented`
#'   and `_observed_only` ([pt_multi_source()]); `mean_parent_set_entropy` and `mean_effective_parent_sets`
#'   ([pt_entropy()]). The out-degree and unobserved-source intervals are named `average_outdegree_lower` and
#'   `unobserved_source_lower` (and `_upper`).
#' * **Edge recovery**: the `directed_` and `undirected_` columns of [pt_score_edges()] with direct truth.
#' * **Parent sets** at level 0.8, from [pt_score_parent_sets()]: `mean_credible_parent_set_size_80`,
#'   `median_credible_parent_set_size_80`, `credible_parent_set_coverage_80` (the share covered),
#'   `mean_top_parent_set_probability`.
#' * **Errors**: `average_outdegree_error` and `unobserved_source_error` (estimate minus truth), each with its
#'   `_abs_error` and `_sq_error`.
#'
#' **Fit types.** `"posterior"` scores a model fit and `"null_prior"` the null model's fit, both the same way.
#' `"dcifer"` scores a pairwise baseline instead: the edge-recovery columns come from [pt_score_pairs()] on `scores`
#' over the Analysis's Infections, and the functional, Parent-set and error columns are `NA`. Pass the Analysis of the
#' fit the baseline is compared with: it gives the Infections.
#'
#' **Missingness** labels the fit only. In a fit to part of the cases, true parents outside the Analysis count as the
#' Latent parent, which the truth columns and scores already account for.
#'
#' @inheritParams pt_score_edges
#' @param fit_type `"posterior"`, `"null_prior"` or `"dcifer"`.
#' @param missingness The missingness condition label, such as `"full"` or `"missing_50"`.
#' @param factors `NULL`, or a one-row data frame of factor columns bound in front of the other columns, such as the
#'   thesis factors parsed from a directory path by `tools/plasmotrackr/thesis-factors.R`. When it holds the 12
#'   thesis factor columns, the row has the 70-column layout; with `NULL` those columns are `NA`. Its
#'   `missingness_condition` and `fit_type` columns, if any, stand in for `missingness` and `fit_type` when those
#'   arguments are not given.
#' @param scores With `fit_type = "dcifer"` only: a data frame of `a`, `b`, `score`, as for [pt_score_pairs()].
#' @return A one-row tibble.
#' @seealso [pt_score_edges()], [pt_score_parent_sets()], [pt_interval()]; [pt_read_truth()].
#' @export
#' @examples
#' \dontrun{
#' truth <- pt_read_truth("sim_0_s7.json.gz")
#' full <- pt_read("full/sim_0_s7")
#' rows <- dplyr::bind_rows(
#'   pt_fit_summary(full, truth),
#'   pt_fit_summary(pt_read("null/sim_0_s7"), truth, fit_type = "null_prior"),
#'   pt_fit_summary(full, truth, fit_type = "dcifer", scores = read.csv("dcifer_ibd.csv"))
#' )
#' }
pt_fit_summary <- function(a, truth, fit_type = "posterior", missingness = "full", factors = NULL,
                           scores = NULL) {
  call <- rlang::current_env()
  factors <- check_fit_factors(factors, call)
  if (missing(fit_type) && "fit_type" %in% names(factors)) fit_type <- factors$fit_type
  if (missing(missingness) && "missingness_condition" %in% names(factors)) {
    missingness <- factors$missingness_condition
  }
  fit_type <- rlang::arg_match(fit_type, c("posterior", "null_prior", "dcifer"))
  if (!rlang::is_string(missingness)) cli::cli_abort("{.arg missingness} must be a single string.", call = call)
  if (fit_type == "dcifer" && is.null(scores)) {
    cli::cli_abort("{.code fit_type = \"dcifer\"} needs the pairwise {.arg scores}.", call = call)
  }
  if (fit_type != "dcifer" && !is.null(scores)) {
    cli::cli_abort("{.arg scores} is used only with {.code fit_type = \"dcifer\"}.", call = call)
  }

  ids <- truth_infection_ids(a, truth, call)
  true_sets <- true_parent_sets(a, truth, call)
  true_edges <- unique(truth_direct_edges(truth, ids))
  truth_columns <- tibble::tibble(
    n_true_cases = nrow(truth$infections),
    n_observed_cases = length(ids),
    true_average_outdegree = nrow(true_edges) / length(ids),
    true_unobserved_source_proportion = mean(true_sets$has_latent_parent)
  )

  if (fit_type == "dcifer") {
    edges <- pt_score_pairs(scores, truth, ids)
    functionals <- fit_summary_functional_columns(NULL)
    sets <- tibble::tibble(mean_credible_parent_set_size_80 = NA_real_, median_credible_parent_set_size_80 = NA_real_,
                           credible_parent_set_coverage_80 = NA_real_, mean_top_parent_set_probability = NA_real_)
  } else {
    runs <- runs_with_draws(a, call = call)
    edges <- score_edges(a, truth, "direct", FALSE, runs, call)
    functionals <- fit_summary_functional_columns(fit_summary_functionals(a, runs))
    credible <- score_parent_sets(a, truth, 0.8, runs, call)
    sets <- tibble::tibble(
      mean_credible_parent_set_size_80 = mean(credible$credible_size),
      median_credible_parent_set_size_80 = stats::median(as.numeric(credible$credible_size)),
      credible_parent_set_coverage_80 = mean(credible$covered),
      mean_top_parent_set_probability = mean(credible$top_probability)
    )
  }

  errors <- function(stem, estimate, truth) {
    error <- estimate - truth
    rlang::set_names(list(error, abs(error), error^2), paste0(stem, c("_error", "_abs_error", "_sq_error")))
  }
  dplyr::bind_cols(
    factors[setdiff(names(factors), c("missingness_condition", "fit_type"))],
    tibble::tibble(missingness_condition = missingness, fit_type = fit_type),
    truth_columns,
    functionals,
    edge_columns(edges),
    sets,
    tibble::as_tibble(c(
      errors("average_outdegree", functionals$estimated_average_outdegree, truth_columns$true_average_outdegree),
      errors("unobserved_source", functionals$estimated_unobserved_source_proportion,
             truth_columns$true_unobserved_source_proportion)
    ))
  )
}

# `factors` as a one-row tibble (the typed NA thesis factors when NULL); no name may clash with a computed column.
check_fit_factors <- function(factors, call) {
  if (is.null(factors)) return(fit_summary_factors)
  if (!is.data.frame(factors) || nrow(factors) != 1) {
    cli::cli_abort("{.arg factors} must be a one-row data frame.", call = call)
  }
  factors <- tibble::as_tibble(factors)
  for (column in intersect(names(factors), c("missingness_condition", "fit_type"))) {
    if (!rlang::is_string(factors[[column]])) {
      cli::cli_abort("{.arg factors}' {.field {column}} must be a string.", call = call)
    }
  }
  clash <- intersect(names(factors), fit_summary_computed_columns())
  if (length(clash)) {
    cli::cli_abort("{.arg factors} can't hold {.field {clash}}: {?it is a/they are} computed column{?s}.",
                   call = call)
  }
  factors
}

fit_summary_computed_columns <- function() {
  c("n_true_cases", "n_observed_cases", "true_average_outdegree", "true_unobserved_source_proportion",
    names(fit_summary_functional_columns(NULL)),
    paste0(rep(c("directed_", "undirected_"), each = 6),
           c("n_edge_universe", "n_true_edges", "n_distinct_scores", "average_precision", "pr_auc", "top_k_precision")),
    "mean_credible_parent_set_size_80", "median_credible_parent_set_size_80", "credible_parent_set_coverage_80",
    "mean_top_parent_set_probability",
    paste0(rep(c("average_outdegree", "unobserved_source"), each = 3), c("_error", "_abs_error", "_sq_error")))
}

# Every functional's posterior mean and 95% interval: a tibble of value (as in fit_summary_functional_names),
# estimate, lower, upper.
fit_summary_functionals <- function(a, runs) {
  # The composition rolled up by kind of class, per Draw. Every Draw has a row for every class of the universe.
  kind <- dbplyr::sql(paste(
    "CASE WHEN class = 'latent_only' THEN 'unobserved_only'",
    "WHEN suffix(class, '_latent') THEN 'mixed'",
    "WHEN prefix(class, 'observed_') THEN 'observed_only'",
    "ELSE 'no_parent' END"
  ))
  composition <- pt_composition(a, runs = runs) |>
    dplyr::mutate(kind = kind) |>
    dplyr::group_by(.data$run, .data$iteration) |>
    dplyr::summarise(
      observed_only = sum(ifelse(.data$kind == "observed_only", .data$share, 0), na.rm = TRUE),
      mixed = sum(ifelse(.data$kind == "mixed", .data$share, 0), na.rm = TRUE),
      unobserved_only = sum(ifelse(.data$kind == "unobserved_only", .data$share, 0), na.rm = TRUE),
      .groups = "drop"
    )
  multi <- pt_multi_source(a, runs = runs) |>
    dplyr::rename(multi_observed_only = "observed_only")
  intervals <- list(
    pt_interval(pt_mean_out_degree(a, runs = runs)) |> tibble::add_column(value = "mean_out_degree"),
    pt_interval(pt_mean_p_latent_parent(a, runs = runs)) |> tibble::add_column(value = "mean_p_latent_parent"),
    pt_interval(composition),
    pt_interval(multi),
    pt_interval(pt_entropy(a, runs = runs))
  )
  out <- dplyr::bind_rows(intervals)[c("value", "estimate", "lower", "upper")]
  # A Parent set is never empty: its weight is 0 at every Draw.
  dplyr::bind_rows(out, tibble::tibble(value = "no_parent", estimate = 0, lower = 0, upper = 0))
}

# The 30 functional columns in layout order, from fit_summary_functionals() rows (NA when NULL).
fit_summary_functional_columns <- function(values) {
  columns <- list()
  for (stem in names(fit_summary_functional_names)) {
    at <- if (is.null(values)) NA_integer_ else match(fit_summary_functional_names[[stem]], values$value)
    pick <- function(field) if (is.na(at)) NA_real_ else values[[field]][[at]]
    names <- switch(stem,
      estimated_average_outdegree = c(stem, "average_outdegree_lower", "average_outdegree_upper"),
      estimated_unobserved_source_proportion = c(stem, "unobserved_source_lower", "unobserved_source_upper"),
      paste0(stem, c("", "_lower", "_upper"))
    )
    columns[names] <- list(pick("estimate"), pick("lower"), pick("upper"))
  }
  tibble::as_tibble(columns)
}

# pt_score_edges()/pt_score_pairs() rows as one row of directed_ and undirected_ columns.
edge_columns <- function(edges) {
  row <- function(direction) {
    x <- edges[edges$direction == direction, setdiff(names(edges), "direction")]
    rlang::set_names(x, paste0(direction, "_", names(x)))
  }
  dplyr::bind_cols(row("directed"), row("undirected"))
}
