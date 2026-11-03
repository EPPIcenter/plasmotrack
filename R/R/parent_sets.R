# Pooled Parent sets and what is read off them: the MAP graph, ψ_MAP, the size range (ψ_lo, ψ_hi) and co-inclusion.
#
# Internal building blocks:
#   pooled_sets(a, runs, by_run, level)  eager: the ranked sets of pt_parent_sets_pooled() (plus `run` when by_run).
#   set_totals(a, runs)   eager: run, infection, observed (members joined by "\x1f"; "" for none),
#                         has_latent_parent, total (the set's summed weight over the Run's sampling Draws).
#   rank_sets(...)        rank, cumulative mass and the Credible collection within each (run, infection).

#' Pooled Parent sets, MAP graph and Credible collections
#'
#' Each Infection's posterior over Parent sets, pooled over the sampling Draws of every Run, ranked, with its
#' Credible collection marked.
#'
#' A Parent set is identified by its observed members and whether it holds the Latent parent (the manuscript's
#' unobserved source), so the same set is matched across Runs whatever its `set_id` there. Its `weight` is the sum of
#' its stored weight over every sampling Draw of every Run, divided by the total Draw count, so each Run weighs by its
#' Draws; a set with no row at a Draw counts 0 there. Weights are used as stored: the Store drops sets below
#' `parent_set_floor` (default 1e-9) without renormalising, so an Infection's weights may sum to slightly less than 1.
#'
#' `rank` orders an Infection's sets by weight, largest first. Ties in weight go to the members list compared element
#' by element: the observed parents in id string order (`"10" < "2"`), then the Latent parent, which sorts after any
#' observed id; a list that is a prefix of another comes first. So `{10}` precedes `{10, Latent}`, which precedes
#' `{2}`, which precedes `{Latent}`.
#'
#' The Credible collection at `level` is the highest-ranked sets whose cumulative weight reaches `level` (up to
#' 1e-12), plus every set tied in weight with the last set needed. When an Infection's weights never reach `level`,
#' every one of its sets is credible.
#'
#' `pt_map_graph()` keeps the rank-1 rows: each Infection's MAP Parent set. Together they make the MAP graph.
#'
#' @inheritParams pt_edge_marginals
#' @param level Level of the Credible collection, in (0, 1].
#' @param by_run `FALSE` to pool over Runs; `TRUE` to rank each Run's own sets, one block of rows per Run (a set
#'   gets a row only in the Runs where it has weight).
#' @return A tibble with `run` first when `by_run = TRUE`, then `infection`, `members` (a list of the observed
#'   parent ids, in string order), `n_observed`, `has_latent_parent`, `class` (the Composition class, named as in
#'   [pt_class_draws()]), `weight`, `rank`, `cumulative_mass` (the weight of this set and every set ranked above it)
#'   and `credible`. Rows are sorted by `run` (in [pt_runs()] order), `infection` (string order), then `rank`.
#' @seealso [pt_map_out_degree()] and [pt_out_degree_range()] for ψ_MAP and the size range; [pt_class_draws()] for
#'   per-Draw class shares.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_parent_sets_pooled(a)
#' pt_parent_sets_pooled(a, level = 0.95, by_run = TRUE)
#' pt_map_graph(a)
#' pt_close(a)
pt_parent_sets_pooled <- function(a, level = 0.8, by_run = FALSE, runs = NULL) {
  check_credible_level(level)
  check_flag(by_run)
  pooled_sets(a, runs_with_draws(a, runs), by_run, level)
}

#' @rdname pt_parent_sets_pooled
#' @export
pt_map_graph <- function(a, by_run = FALSE, runs = NULL) {
  check_flag(by_run)
  sets <- pooled_sets(a, runs_with_draws(a, runs), by_run, 0.8)
  sets[sets$rank == 1L, ]
}

check_credible_level <- function(level, call = rlang::caller_env()) {
  if (!is.numeric(level) || length(level) != 1 || is.na(level) || level <= 0 || level > 1) {
    cli::cli_abort("{.arg level} must be a number in (0, 1].", call = call)
  }
}

pooled_sets <- function(a, runs, by_run, level) {
  totals <- set_totals(a, runs)
  draws <- stats::setNames(as.numeric(a$runs$draws), a$runs$run)[runs]

  # Sum each set's totals within its group: (run,) infection, observed members, has_latent_parent.
  group_run <- if (by_run) totals$run else rep("", nrow(totals))
  id <- paste(group_run, totals$infection, totals$observed, totals$has_latent_parent, sep = "\r")
  first <- !duplicated(id)
  sets <- tibble::tibble(
    run = group_run[first], infection = totals$infection[first], observed = totals$observed[first],
    has_latent_parent = totals$has_latent_parent[first],
    total = as.vector(rowsum(totals$total, factor(id, levels = id[first]), reorder = FALSE))
  )
  sets$weight <- sets$total / if (by_run) unname(draws[sets$run]) else sum(draws)
  members <- strsplit(sets$observed, "\x1f", fixed = TRUE)
  sets$n_observed <- as.numeric(lengths(members))

  sets <- rank_sets(sets, members, run_order = runs, level = level)
  out <- tibble::tibble(
    run = sets$run, infection = sets$infection, members = sets$members, n_observed = sets$n_observed,
    has_latent_parent = sets$has_latent_parent,
    class = composition_class(sets$n_observed, sets$has_latent_parent),
    weight = sets$weight, rank = sets$rank, cumulative_mass = sets$cumulative_mass, credible = sets$credible
  )
  if (!by_run) out$run <- NULL
  out
}

# Each set's summed stored weight over the Run's sampling Draws, with its members. DuckDB aggregates; one row per
# (run, infection, set_id) comes back.
set_totals <- function(a, runs) {
  totals <- pt_table(a, "parent_sets", runs = runs, phase = "sampling") |>
    dplyr::group_by(.data$run, .data$infection, .data$set_id) |>
    dplyr::summarise(total = sum(.data$prob, na.rm = TRUE), .groups = "drop")
  members <- pt_table(a, "parent_set_members", runs = runs) |>
    dplyr::inner_join(structure_infections(a), by = "infection") |>
    dplyr::group_by(.data$run, .data$infection, .data$set_id) |>
    dplyr::summarise(
      observed = dbplyr::sql("coalesce(string_agg(member, chr(31) ORDER BY member) FILTER (WHERE member <> latent_parent), '')"),
      has_latent_parent = sum(as.integer(.data$member == .data$latent_parent), na.rm = TRUE) > 0L,
      .groups = "drop"
    )
  rows <- totals |>
    dplyr::inner_join(members, by = c("run", "infection", "set_id")) |>
    dplyr::select("run", "infection", "observed", "has_latent_parent", "total") |>
    dplyr::collect()
  rows$total <- as.numeric(rows$total)
  # Canonical members string: ids in string (byte) order whatever the database collation.
  canonical <- unique(rows$observed)
  sorted <- vapply(strsplit(canonical, "\x1f", fixed = TRUE),
                   function(x) paste(sort(x, method = "radix"), collapse = "\x1f"), "")
  rows$observed <- sorted[match(rows$observed, canonical)]
  rows
}

# Rank sets within (run, infection): weight descending, then the members list element by element (observed ids in
# string order, the Latent parent after every id, a prefix first). Adds members (list), rank, cumulative_mass and
# credible, and sorts rows by run (in `run_order`), infection, rank.
rank_sets <- function(sets, members, run_order, level) {
  ids <- sort(unique(unlist(members)), method = "radix")
  latent_token <- length(ids) + 1L
  tokens <- lapply(seq_along(members), function(i) {
    c(match(members[[i]], ids), if (sets$has_latent_parent[[i]]) latent_token)
  })
  width <- max(1L, lengths(tokens))
  # Token matrix padded with 0: the end of a list sorts before any member.
  columns <- lapply(seq_len(width), function(k) vapply(tokens, function(t) if (k <= length(t)) t[[k]] else 0L, 0L))
  run_index <- match(sets$run, c("", run_order))
  infection_index <- match(sets$infection, sort(unique(sets$infection), method = "radix"))
  o <- do.call(order, c(list(run_index, infection_index, -sets$weight), columns, list(method = "radix")))
  sets <- sets[o, ]
  sets$members <- members[o]

  group <- paste(sets$run, sets$infection, sep = "\r")
  sets$rank <- stats::ave(seq_len(nrow(sets)), group, FUN = seq_along)
  sets$cumulative_mass <- stats::ave(sets$weight, group, FUN = cumsum)
  sets$credible <- stats::ave(seq_len(nrow(sets)), group, FUN = function(i) {
    cumulative <- sets$cumulative_mass[i]
    m <- which(cumulative >= level - 1e-12)
    if (length(m) == 0) return(rep(1L, length(i)))
    m <- m[[1]]
    as.integer(seq_along(i) <= m | sets$weight[i] == sets$weight[i][[m]])
  }) == 1L
  sets
}

composition_class <- function(n_observed, has_latent_parent) {
  ifelse(n_observed == 0,
         ifelse(has_latent_parent, "latent_only", "no_parent"),
         paste0("observed_", n_observed, ifelse(has_latent_parent, "_latent", "")))
}
