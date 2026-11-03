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

#' MAP out-degree and out-degree size range
#'
#' Summaries of how many children an Infection infected, read off the pooled Parent sets of
#' [pt_parent_sets_pooled()]. Both are means over every Infection as a source, so they are on the scale of
#' [pt_mean_out_degree()] (ψ_edge).
#'
#' * `pt_map_out_degree()` estimates the manuscript's ψ_MAP: the mean out-degree in the MAP graph, the number of
#'   observed parents over every Infection's MAP Parent set divided by the number of Infections. The Latent parent
#'   is never counted.
#' * `pt_out_degree_range()` gives the manuscript's size range (ψ_lo, ψ_hi): for each Infection j, the fewest and the
#'   most observed parents over the sets of its Credible collection at `level`, summed over every Infection and
#'   divided by the number of Infections. In the manuscript's notation, ψ_lo(S) = n_S⁻¹ Σ_j min over U in C_j of
#'   |U_obs ∩ S| with S every Infection, and ψ_hi the same with max. The MAP Parent set is always credible, so
#'   `lo <= pt_map_out_degree() <= hi`.
#'
#' Sets are pooled over the sampling Draws of every Run (or of each Run with `by_run = TRUE`); Runs with no sampling
#' Draws are skipped with a message.
#'
#' **Strata.** With `by`, S is one stratum's Infections, counted as sources: ψ_MAP(S) is the number of MAP-set
#' parents (of every recipient) that lie in S, over n_S, the number of Infections in S; ψ_lo(S) and ψ_hi(S) take, for
#' every recipient j, the min and max over its Credible collection of the observed parents in S, summed and divided by
#' n_S. A stratum is a value of the epi column `by`; Infections with an NA value, or missing from epi, form the NA
#' stratum, which is kept. ψ_MAP(S) weighted by n_S averages back to ψ_MAP, and `lo <= map_out_degree <= hi` in every
#' stratum. The size range does not average back: a sum of per-stratum minima is at most the overall minimum (and of
#' maxima at least the overall maximum).
#'
#' @inheritParams pt_parent_sets_pooled
#' @param by `NULL`, or the name of one epi column (a string; see [pt_read()]) to stratify the sources by.
#' @param by_run `FALSE` for one row pooled over Runs; `TRUE` for one row per Run, from that Run's own sets.
#' @return A tibble with `run` first when `by_run = TRUE`, then `stratum` with `by` (one row per stratum, sorted, NA
#'   last). `pt_map_out_degree()`: `map_out_degree` and `infections` (the number of Infections averaged over, n_S).
#'   `pt_out_degree_range()`: `lo`, `hi` and `level`.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_map_out_degree(a)
#' pt_out_degree_range(a)
#' pt_out_degree_range(a, level = 0.95, by_run = TRUE)
#' pt_close(a)
#'
#' epi <- data.frame(infection = c("1", "2"), district = c("north", "south"))
#' a <- pt_read(file.path(golden, "single-chain"), epi = epi, threads = 1)
#' pt_map_out_degree(a, by = "district")
#' pt_out_degree_range(a, by = "district")
#' pt_close(a)
pt_map_out_degree <- function(a, by = NULL, by_run = FALSE, runs = NULL) {
  check_by(a, by)
  check_flag(by_run)
  sets <- pooled_sets(a, runs_with_draws(a, runs), by_run, 0.8)
  map <- sets[sets$rank == 1L, ]
  strata <- source_strata(a, by)
  per_run(map, by_run, function(rows) {
    counts <- colSums(stratum_counts(rows$members, strata))
    stratum_rows(strata, by, map_out_degree = counts / strata$n, infections = strata$n)
  })
}

#' @rdname pt_map_out_degree
#' @export
pt_out_degree_range <- function(a, level = 0.8, by = NULL, by_run = FALSE, runs = NULL) {
  check_credible_level(level)
  check_by(a, by)
  check_flag(by_run)
  sets <- pooled_sets(a, runs_with_draws(a, runs), by_run, level)
  credible <- sets[sets$credible, ]
  strata <- source_strata(a, by)
  per_run(credible, by_run, function(rows) {
    # Per credible set, its observed parents in each stratum; min and max within each recipient.
    counts <- stratum_counts(rows$members, strata)
    child <- factor(rows$infection)
    lo <- vapply(seq_along(strata$n), function(k) sum(tapply(counts[, k], child, min)), 0)
    hi <- vapply(seq_along(strata$n), function(k) sum(tapply(counts[, k], child, max)), 0)
    stratum_rows(strata, by, lo = lo / strata$n, hi = hi / strata$n, level = level)
  })
}

# Apply `f` to the rows of each Run (by_run: one result row per Run, `run` first, in row order) or to all rows.
per_run <- function(rows, by_run, f) {
  if (!by_run) return(f(rows))
  runs <- unique(rows$run)
  parts <- lapply(runs, function(r) f(rows[rows$run == r, ]))
  tibble::tibble(run = rep(runs, vapply(parts, nrow, 1L)), dplyr::bind_rows(parts))
}

#' Co-inclusion of parents
#'
#' For each child, the posterior probability that two observed parents infected it together: the pooled weight of
#' the child's Parent sets that hold both. Only sets with at least two observed parents count, so the Latent parent
#' never appears. Weights are pooled as in [pt_parent_sets_pooled()].
#'
#' @inheritParams pt_parent_sets_pooled
#' @param infections Children to report, as Infection ids; `NULL` for every Infection.
#' @param by_run `FALSE` for probabilities pooled over Runs; `TRUE` for one row per Run per pair, where a pair seen in
#'   any Run has a row (0 if need be) in every Run.
#' @return A tibble of `child`, `parent_a`, `parent_b` (`parent_a < parent_b` in string order) and `probability`,
#'   after `run` when `by_run = TRUE`, sorted by `run` (in [pt_runs()] order), then `child`, `parent_a` and
#'   `parent_b` as strings. Pairs that never share a set have no row.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_co_inclusion(a)
#' pt_co_inclusion(a, infections = "3", by_run = TRUE)
#' pt_close(a)
pt_co_inclusion <- function(a, infections = NULL, by_run = FALSE, runs = NULL) {
  check_flag(by_run)
  if (!is.null(infections)) {
    unknown <- setdiff(infections, a$infections$infection)
    if (!is.character(infections) || length(unknown)) {
      cli::cli_abort("Unknown Infection{?s}: {.val {unknown}}.")
    }
  }
  runs <- runs_with_draws(a, runs)
  sets <- pooled_sets(a, runs, TRUE, 0.8)
  keep <- sets$n_observed >= 2 & (is.null(infections) | sets$infection %in% infections)
  sets <- sets[keep, ]
  pairs <- lapply(sets$members, function(m) utils::combn(m, 2))
  n <- vapply(pairs, ncol, 1L)
  rows <- tibble::tibble(
    run = rep(sets$run, n), child = rep(sets$infection, n),
    parent_a = as.character(unlist(lapply(pairs, function(p) p[1, ]))),
    parent_b = as.character(unlist(lapply(pairs, function(p) p[2, ]))),
    # A by-Run weight times the Run's Draws is the set's summed weight; sums are divided by the Draws of the Run, or
    # of every Run.
    total = rep(sets$weight * unname(draws_of(a, sets$run)), n)
  )
  if (!by_run) rows$run <- ""
  pair <- paste(rows$child, rows$parent_a, rows$parent_b, sep = "\r")
  keys <- unique(rows[c("child", "parent_a", "parent_b")])
  keys <- keys[order(keys$child, keys$parent_a, keys$parent_b, method = "radix"), ]
  key <- paste(keys$child, keys$parent_a, keys$parent_b, sep = "\r")
  columns <- if (by_run) runs else ""
  total <- matrix(0, nrow(keys), length(columns))
  for (j in seq_along(columns)) {
    mine <- rows$run == columns[[j]]
    total[, j] <- as.vector(rowsum(c(rows$total[mine], numeric(nrow(keys))),
                                   factor(c(pair[mine], key), levels = key), reorder = FALSE))
  }
  if (!by_run) {
    return(tibble::tibble(keys, probability = total[, 1] / sum(draws_of(a, runs))))
  }
  tibble::tibble(
    run = rep(runs, each = nrow(keys)), child = rep(keys$child, length(runs)),
    parent_a = rep(keys$parent_a, length(runs)), parent_b = rep(keys$parent_b, length(runs)),
    probability = as.vector(sweep(total, 2, draws_of(a, runs), "/"))
  )
}

draws_of <- function(a, runs) stats::setNames(as.numeric(a$runs$draws), a$runs$run)[runs]

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
