# Flows between strata: lazy, per Draw. Every flow is computed on integer stratum ids (`sid`, 1..k in
# stratum_levels() order, the NA stratum included), so the NA stratum joins and compares like any other; the values
# are joined back at the end.
#
#   flow_strata(a, by)            list(map = lazy infection, sid; levels = eager sid, stratum; k).
#   flow_edges(a, by, runs, phase)  lazy: pt_edge_draws() rows with the parent's (`psid`) and child's (`csid`) sid.
#   flow_grid(a, runs, phase, ...)  lazy: every Draw crossed with inline tables of keys.
#   flow_values(x, levels, column, sid)  lazy: `x` with the stratum of `sid` added as `column` (NA where no match).

flow_strata <- function(a, by) {
  strata <- infection_strata(a, by)
  levels <- stratum_levels(strata$stratum)
  map <- tibble::tibble(infection = strata$infection, sid = match(strata$stratum, levels))
  list(
    map = dbplyr::copy_inline(analysis_con(a), map),
    levels = tibble::tibble(sid = seq_along(levels), stratum = levels),
    k = length(levels)
  )
}

flow_edges <- function(a, strata, runs, phase) {
  observed_edges(a, runs, phase) |>
    dplyr::inner_join(dplyr::rename(strata$map, child = "infection", csid = "sid"), by = "child") |>
    dplyr::inner_join(dplyr::rename(strata$map, parent = "infection", psid = "sid"), by = "parent")
}

flow_grid <- function(a, runs, phase, ...) {
  grid <- draw_grid(a, runs, phase)
  for (keys in list(...)) grid <- dplyr::cross_join(grid, dbplyr::copy_inline(analysis_con(a), keys))
  grid
}

flow_values <- function(x, a, strata, column, sid) {
  levels <- dbplyr::copy_inline(analysis_con(a), strata$levels)
  dplyr::left_join(x, dplyr::rename(levels, !!sid := "sid", !!column := "stratum"), by = sid)
}

#' Flows between strata
#'
#' How transmission runs between the strata of one epi column, one value per Draw. A stratum is a value of the epi
#' column `by` (see [pt_read()]); Infections with an NA value, or missing from epi, form the NA stratum, which every
#' output keeps. Each function returns a lazy table with `run` and `iteration` first and one row per Draw per
#' stratum (or pair of strata), zeros included; [pt_interval()] turns it into estimates with intervals. Weights are
#' the stored Parent-set weights, as in [pt_edge_draws()]: an edge's mass at a Draw is its `prob` there.
#'
#' * `pt_attribution()`: `recipient`, `source`, `latent`, `share`. A recipient stratum's attributed parent mass at a
#'   Draw is its Infections' expected observed parents plus their Latent parent probabilities. `share` splits it by
#'   the source stratum of the observed parents, plus one row with `latent = TRUE` and `source = NA` for the Latent
#'   parent, normalised so that the shares of a recipient sum to 1 at each Draw. The NA stratum as a source has
#'   `latent = FALSE`.
#' * `pt_enrichment()`: `source`, `recipient`, `enrichment`: the cell's share of all observed edge mass at the Draw
#'   over its share of the candidate pairs. A candidate pair (parent, child) lies in the cell of the parent's and the
#'   child's strata. `enrichment` is NA where the cell has no candidate pair (and at a Draw with no observed edge);
#'   filter those rows before [pt_interval()]. Above 1, the cell carries more transmission than its candidate pairs
#'   alone would give it.
#' * `pt_net_flow()`: `stratum`, `net_flow`: the stratum's expected out-edges to other strata minus its in-edges
#'   from other strata. Edges within a stratum do not count. The net flows of a Draw sum to 0.
#'
#' @inheritParams pt_table
#' @param by The name of one epi column, a string.
#' @param candidates A data frame of `parent` and `child` Infection ids, the candidate pairs. The default,
#'   [pt_candidates()], is the input's `allowed_parents`. Pairs naming an id that is not an Infection of the
#'   Analysis are dropped, and repeated pairs count once.
#' @return A lazy table (`tbl_lazy`) with `run`, `iteration`, then the columns above. Stratum columns keep the epi
#'   column's type (a factor becomes character).
#' @seealso [pt_reassign_long_lags()]; [pt_mean_out_degree()] and the other functionals with `by`; [pt_interval()].
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' epi <- data.frame(infection = c("1", "2", "3"), district = c("north", "south", "south"))
#' a <- pt_read(store, epi = epi, threads = 1)
#' pt_attribution(a, by = "district") |> pt_interval(by = c(recipient, source, latent))
#' pt_enrichment(a, by = "district") |>
#'   dplyr::filter(!is.na(enrichment)) |>
#'   pt_interval(by = c(source, recipient))
#' pt_net_flow(a, by = "district") |> pt_interval(by = stratum)
#' pt_close(a)
pt_attribution <- function(a, by, runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  check_flow_by(a, by)
  strata <- flow_strata(a, by)
  observed <- flow_edges(a, strata, runs, phase) |>
    dplyr::group_by(.data$run, .data$iteration, rsid = .data$csid, ssid = .data$psid) |>
    dplyr::summarise(mass = sum(.data$prob, na.rm = TRUE), .groups = "drop")
  # The Latent parent is source id 0.
  latent <- pt_infection_draws(a, runs, phase) |>
    dplyr::inner_join(strata$map, by = "infection") |>
    dplyr::group_by(.data$run, .data$iteration, rsid = .data$sid) |>
    dplyr::summarise(mass = sum(.data$p_latent_parent, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(ssid = 0L)
  sources <- tibble::tibble(ssid = c(seq_len(strata$k), 0L), latent = c(rep(FALSE, strata$k), TRUE))
  flow_grid(a, runs, phase, tibble::tibble(rsid = seq_len(strata$k)), sources) |>
    dplyr::left_join(dplyr::union_all(observed, latent), by = c("run", "iteration", "rsid", "ssid")) |>
    dplyr::mutate(mass = dplyr::coalesce(.data$mass, 0)) |>
    dplyr::group_by(.data$run, .data$iteration, .data$rsid) |>
    dplyr::mutate(total = sum(.data$mass, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(share = dplyr::if_else(.data$total > 0, .data$mass / .data$total, NA)) |>
    flow_values(a, strata, "recipient", "rsid") |>
    flow_values(a, strata, "source", "ssid") |>
    dplyr::select("run", "iteration", "recipient", "source", "latent", "share")
}

#' @rdname pt_attribution
#' @export
pt_enrichment <- function(a, by, candidates = pt_candidates(a), runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  check_flow_by(a, by)
  strata <- flow_strata(a, by)
  cells <- candidate_cells(a, by, candidates, strata$k)
  mass <- flow_edges(a, strata, runs, phase) |>
    dplyr::group_by(.data$run, .data$iteration, ssid = .data$psid, rsid = .data$csid) |>
    dplyr::summarise(mass = sum(.data$prob, na.rm = TRUE), .groups = "drop")
  flow_grid(a, runs, phase, cells) |>
    dplyr::left_join(mass, by = c("run", "iteration", "ssid", "rsid")) |>
    dplyr::mutate(mass = dplyr::coalesce(.data$mass, 0)) |>
    dplyr::group_by(.data$run, .data$iteration) |>
    dplyr::mutate(total = sum(.data$mass, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(enrichment = dplyr::if_else(.data$candidate_share > 0 & .data$total > 0,
                                              .data$mass / .data$total / .data$candidate_share, NA)) |>
    flow_values(a, strata, "source", "ssid") |>
    flow_values(a, strata, "recipient", "rsid") |>
    dplyr::select("run", "iteration", "source", "recipient", "enrichment")
}

# Every (source, recipient) cell of stratum ids with its share of the candidate pairs (0 without any).
candidate_cells <- function(a, by, candidates, k, call = rlang::caller_env()) {
  if (!is.data.frame(candidates) || !all(c("parent", "child") %in% names(candidates))) {
    cli::cli_abort("{.arg candidates} must be a data frame with {.field parent} and {.field child} columns.",
                   call = call)
  }
  pairs <- unique(data.frame(parent = as.character(candidates$parent), child = as.character(candidates$child)))
  strata <- infection_strata(a, by)
  sid <- match(strata$stratum, stratum_levels(strata$stratum))
  ssid <- sid[match(pairs$parent, strata$infection)]
  rsid <- sid[match(pairs$child, strata$infection)]
  known <- !is.na(ssid) & !is.na(rsid)
  counts <- tabulate((ssid[known] - 1L) * k + rsid[known], k * k)
  total <- sum(counts)
  tibble::tibble(
    ssid = rep(seq_len(k), each = k),
    rsid = rep(seq_len(k), k),
    candidate_share = if (total > 0) counts / total else rep(0, k * k)
  )
}

check_flow_by <- function(a, by, call = rlang::caller_env()) {
  if (is.null(by)) cli::cli_abort("{.arg by} must name an epi column: flows are between strata.", call = call)
  check_by(a, by, call = call)
}
