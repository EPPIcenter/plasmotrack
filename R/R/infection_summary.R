# The Infections table: one row per Infection with its key posterior summaries, each read off the function it
# summarises (pt_infection_draws, pt_class_draws, pt_parent_sets_pooled, pt_edge_marginals, pt_coi).

#' Infection summary
#'
#' One row per Infection with its key posterior summaries, pooled over the sampling Draws of every Run: the table of
#' Infections a paper reports. Each column summarises a function of the package, and agrees with it:
#'
#' * `observation_time`: from [pt_infections()] (`NA` without the input JSON).
#' * `p_latent_parent`, `expected_observed_parents`: means over Draws of [pt_infection_draws()]'s columns.
#'   `p_latent_parent` is the manuscript's q_i, the probability of an unobserved source.
#' * `share_<class>`, one column per Composition class: the mean over Draws of the class's share in
#'   [pt_class_draws()] (`share_latent_only`, `share_observed_1`, ...). The classes are those of any Parent set with
#'   weight in the Runs read, ordered by the number of observed parents, then without before with the Latent parent.
#' * `mean_entropy`, `mean_effective_size`: means over Draws of the per-Draw Parent-set entropy \eqn{-\sum w \ln w}
#'   and of its exponential (see [pt_infection_draws()]).
#' * `posterior_entropy`, `posterior_effective_size`: the entropy \eqn{-\sum \bar w \ln \bar w} of the pooled
#'   weights \eqn{\bar w} of [pt_parent_sets_pooled()], and its exponential.
#' * `map_members`, `map_class`, `map_weight`: the MAP Parent set, the rank-1 row of [pt_map_graph()]: its observed
#'   members (a list column), Composition class and pooled weight.
#' * `top_parent`, `top_parent_marginal`: the observed parent with the largest Edge marginal in
#'   [pt_edge_marginals()], ties going to the smaller id as a string (`"10"` before `"2"`). `NA` when the Infection
#'   has no observed parent with weight.
#' * `expected_out_degree`: the mean over Draws of [pt_out_degree_draws()], the sum of the Infection's Edge
#'   marginals as a parent.
#' * `coi_median`: the median of [pt_coi()] (`NA` for an Infection genotyped at no Locus).
#'
#' **Two entropies.** The per-Draw entropy measures how spread an Infection's Parent-set weights are given one Draw's
#' Ordering, genotypes and parameters. When each Draw allows a single set, it is 0 (effective size 1) even if the
#' Draws disagree on which set. The posterior entropy measures how spread the posterior over Parent sets is once the
#' Draws are pooled, which is the uncertainty a reader usually means. It is at least the mean per-Draw entropy (up to
#' the weight the floor drops).
#'
#' Means over Draws pool as everywhere in the package: sums over every Draw of every Run divided by the total Draw
#' count, so each Run weighs by its Draws. An Infection with no Parent set row at a Draw counts 0 there (entropy 0,
#' effective size 1). Weights are used as stored (see [pt_infection_draws()] for the floor). Runs with no sampling
#' Draws are skipped with a message.
#'
#' The summary has no R-hat or ESS. Those are per Infection in `pt_convergence()$infections`; join them on
#' `infection`:
#' `dplyr::left_join(pt_infection_summary(a), pt_convergence(a)$infections, by = "infection")`.
#'
#' @inheritParams pt_edge_marginals
#' @param by_run `FALSE` to pool every Run; `TRUE` for one row per Run per Infection, each from that Run's own
#'   Draws.
#' @return A tibble with `run` first when `by_run = TRUE`, then the columns above, one row per Infection (per Run),
#'   sorted by `run` (in [pt_runs()] order) then `infection` as strings.
#' @seealso [pt_infection_draws()], [pt_parent_sets_pooled()], [pt_edge_marginals()], [pt_coi()].
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_infection_summary(a)
#' pt_infection_summary(a, by_run = TRUE)
#' pt_close(a)
pt_infection_summary <- function(a, by_run = FALSE, runs = NULL) {
  check_flag(by_run)
  runs <- runs_with_draws(a, runs)
  draws <- draws_of(a, runs)
  ids <- sort(a$infections$infection, method = "radix")
  groups <- if (by_run) runs else ""
  out <- tibble::tibble(run = rep(groups, each = length(ids)), infection = rep(ids, length(groups)))
  keys <- paste(out$run, out$infection, sep = "\r")
  denominator <- if (by_run) unname(draws[out$run]) else rep(sum(draws), nrow(out))
  # Sum `values` of rows keyed by (run, infection) into the output rows, 0 where a row has none.
  total <- function(run, infection, values) {
    if (!by_run) run <- rep("", length(run))
    k <- factor(paste(run, infection, sep = "\r"), levels = keys)
    as.vector(rowsum(c(values, numeric(length(keys))), c(k, factor(keys, levels = keys)), reorder = TRUE))
  }

  out$observation_time <- a$infections$observation_time[match(out$infection, a$infections$infection)]

  values <- pt_infection_draws(a, runs, "sampling") |>
    dplyr::group_by(.data$run, .data$infection) |>
    dplyr::summarise(
      n = dplyr::n(),
      p_latent_parent = sum(.data$p_latent_parent, na.rm = TRUE),
      expected_observed_parents = sum(.data$expected_observed_parents, na.rm = TRUE),
      entropy = sum(.data$entropy, na.rm = TRUE),
      effective_size = sum(.data$effective_size, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::collect()
  mean_of <- function(column) total(values$run, values$infection, as.numeric(values[[column]])) / denominator
  out$p_latent_parent <- mean_of("p_latent_parent")
  out$expected_observed_parents <- mean_of("expected_observed_parents")

  shares <- weighted_sets(a, runs, "sampling") |>
    dplyr::group_by(.data$run, .data$infection, .data$class, .data$n_observed, .data$has_latent_parent) |>
    dplyr::summarise(share = sum(.data$prob, na.rm = TRUE), .groups = "drop") |>
    dplyr::collect()
  classes <- unique(shares[c("class", "n_observed", "has_latent_parent")])
  classes <- classes$class[order(classes$n_observed, classes$has_latent_parent)]
  for (class in classes) {
    mine <- shares$class == class
    out[[paste0("share_", class)]] <- total(shares$run[mine], shares$infection[mine], shares$share[mine]) / denominator
  }

  out$mean_entropy <- mean_of("entropy")
  # A Draw with no row counts effective size 1.
  out$mean_effective_size <- mean_of("effective_size") - mean_of("n") + 1

  sets <- pooled_sets(a, runs, by_run, 0.8)
  if (!by_run) sets$run <- rep("", nrow(sets))
  positive <- sets$weight > 0
  out$posterior_entropy <- -total(sets$run[positive], sets$infection[positive],
                                  sets$weight[positive] * log(sets$weight[positive]))
  out$posterior_effective_size <- exp(out$posterior_entropy)
  map <- sets[sets$rank == 1L, ]
  at <- match(keys, paste(map$run, map$infection, sep = "\r"))
  out$map_members <- lapply(at, function(i) if (is.na(i)) character() else map$members[[i]])
  out$map_class <- map$class[at]
  out$map_weight <- map$weight[at]

  edges <- pt_edge_marginals(a, by_run = by_run, runs = runs)
  edge_run <- if (by_run) edges$run else rep("", nrow(edges))
  top <- edges[edges$marginal > 0, ]
  top_run <- edge_run[edges$marginal > 0]
  top_key <- paste(top_run, top$child, sep = "\r")
  o <- order(top_key, -top$marginal, top$parent, method = "radix")
  top <- top[o, ]
  top_key <- top_key[o]
  first <- !duplicated(top_key)
  at <- match(keys, top_key[first])
  out$top_parent <- top$parent[first][at]
  out$top_parent_marginal <- top$marginal[first][at]
  out$expected_out_degree <- total(edge_run, edges$parent, edges$marginal)

  coi <- pt_coi(a, by_run = by_run, runs = runs)
  coi_run <- if (by_run) coi$run else rep("", nrow(coi))
  out$coi_median <- coi$median[match(keys, paste(coi_run, coi$infection, sep = "\r"))]

  if (!by_run) out$run <- NULL
  out
}
