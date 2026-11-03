#' Sampler acceptance rates
#'
#' Proposals accepted and made by each sampler type, or each sampler, summed from the Store's `acceptance` Block.
#' A Block has one row per sampler per Window and counts every Sampling iteration of the Window, thinned or not. In a
#' replica-exchange Run only the Hot chain is counted: at each Sampling iteration, the counts of whichever replica
#' held temperature 1.
#'
#' The right rate depends on the sampler, so there is no target band. As a guide:
#'
#' * `ConstrainedContinuousRandomWalk` (COI, transmission parameters, error rates, durations) and `SALTSampler`
#'   (Allele frequencies) adapt their proposal variance during burn-in toward a rate of 0.23. A sampling rate far
#'   from that, near 0 or 1, means adaptation did not settle.
#' * Discrete samplers (`RandomAllelesBitSetSampler*`, `StrainCountSampler`, `AlleleSwapSampler`,
#'   `ParentRelinkSampler`, the Ordering samplers) don't adapt. Low rates are common once the State settles; a rate
#'   of 0 over many proposals means that part of the State never moves.
#' * A rate of `NA` means nothing was proposed. `CloneOrderSwapSampler`, for example, proposes only when an Infection
#'   has a clone to swap with.
#'
#' Window numbers are per Run. With `per_window = TRUE` and `by_run = FALSE`, Windows with the same number in
#' different Runs are summed.
#'
#' @inheritParams pt_table
#' @param level `"type"` for one row per sampler type, `"sampler"` for one row per sampler id (`<type>/<target>`,
#'   as in the Store's `structure/samplers`).
#' @param phase `"sampling"`, `"burn_in"`, or both: the Windows summed.
#' @param by_run `FALSE` sums over Runs; `TRUE` gives one row per Run per key.
#' @param per_window `TRUE` gives one row per Window per key, with its `window` and `phase`.
#' @return A tibble: `run` (with `by_run`), `window` and `phase` (with `per_window`), `type`, `sampler` (with
#'   `level = "sampler"`), `accepted`, `proposed`, and `rate` = `accepted / proposed`, `NA` when `proposed` is 0.
#' @seealso [pt_swaps()] for replica-exchange swap rates; [pt_table()] with `"acceptance"` for the raw Block.
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_acceptance(a)
#' pt_acceptance(a, level = "sampler", per_window = TRUE)
#' pt_close(a)
pt_acceptance <- function(a, level = c("type", "sampler"), phase = "sampling", by_run = FALSE, per_window = FALSE) {
  level <- rlang::arg_match(level)
  phase <- check_phase(phase)
  check_flag(by_run)
  check_flag(per_window)
  keys <- c(block_keys(by_run, per_window), "type", if (level == "sampler") "sampler")
  rows <- pt_table(a, "acceptance", phase = phase) |>
    dplyr::mutate(type = dbplyr::sql("split_part(sampler, '/', 1)"))
  sum_block(a, rows, keys, c("accepted", "proposed"))
}

#' Replica-exchange swap rates
#'
#' Swaps accepted and attempted between each pair of neighbouring temperature slots, summed from the Store's
#' `swaps` Block. Temperature pair `k` is the swaps between slots `k` and `k + 1`; slot 0 is the Hot chain. A Run with
#' one chain has no swaps and adds no rows, so the result is empty when every Run has one chain.
#'
#' A pair whose rate is near 0 is a barrier: replicas rarely pass between those temperatures, so the hotter chains
#' do little for the Hot chain.
#'
#' Window numbers are per Run. With `per_window = TRUE` and `by_run = FALSE`, Windows with the same number in
#' different Runs are summed.
#'
#' @inheritParams pt_acceptance
#' @param by_run `TRUE` (the default) gives one row per Run per pair; `FALSE` sums over Runs.
#' @return A tibble: `run` (with `by_run`), `window` and `phase` (with `per_window`), `temperature_pair`,
#'   `accepted`, `attempted`, and `rate` = `accepted / attempted`, `NA` when `attempted` is 0.
#' @seealso [pt_acceptance()]; [pt_table()] with `"chains"` for each temperature slot's replica and temperature, and
#'   with `"events"` for ladder adaptation (`temperatures_adapted` payloads are JSON strings, read with
#'   [jsonlite::fromJSON()]).
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "replica-exchange", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_swaps(a)
#' pt_swaps(a, phase = c("burn_in", "sampling"), per_window = TRUE)
#' pt_close(a)
pt_swaps <- function(a, phase = "sampling", by_run = TRUE, per_window = FALSE) {
  phase <- check_phase(phase)
  check_flag(by_run)
  check_flag(per_window)
  keys <- c(block_keys(by_run, per_window), "temperature_pair")
  sum_block(a, pt_table(a, "swaps", phase = phase), keys, c("accepted", "attempted"))
}

# The leading keys every Block summary groups by.
block_keys <- function(by_run, per_window) {
  c(if (by_run) "run", if (per_window) c("window", "phase"))
}

# Sum a Block's count columns within `keys`, collect, and add `rate` = first count / second count (NA when the
# second is 0). Counts come back as doubles (DuckDB sums BIGINT to HUGEINT). Rows are ordered by the keys, Runs in
# the Analysis's order.
sum_block <- function(a, rows, keys, counts) {
  rows <- rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(counts), \(x) sum(x, na.rm = TRUE)), .groups = "drop") |>
    dplyr::collect()
  for (count in counts) rows[[count]] <- as.numeric(rows[[count]])
  rate <- rows[[counts[1]]] / rows[[counts[2]]]
  rate[rows[[counts[2]]] == 0] <- NA_real_
  rows$rate <- rate
  sort_keys <- as.list(rows[keys])
  if ("run" %in% keys) sort_keys$run <- match(rows$run, a$runs$run)
  rows <- rows[do.call(order, unname(sort_keys)), c(keys, counts, "rate")]
  tibble::as_tibble(rows)
}

check_flag <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!rlang::is_bool(x)) cli::cli_abort("{.arg {arg}} must be {.val {TRUE}} or {.val {FALSE}}.", call = call)
}
