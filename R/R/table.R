#' A raw Store table
#'
#' Any table of the Store, by its name in the Store, as a lazy dplyr table over the closed Windows of every Run,
#' unioned. Filter it before [dplyr::collect()]: the rows stay in the Parquet files until then.
#'
#' The tables and their columns are described in `doc/store-layout.md` in the plasmotrack repository
#' (<https://github.com/EPPIcenter/plasmotrack/blob/master/doc/store-layout.md>):
#'
#' * Traces, one row per written iteration: `run_scalars`, `infections`, `ordering`, `allele_frequencies`,
#'   `latent_genotypes`, `latent_parent_genotypes` (always empty), `parent_sets`, and `chains` (one row per
#'   temperature slot; the only table with tempered chains).
#' * The lookup table `parent_set_members`, read over every closed Window of both Phases whatever `phase` says,
#'   because a Parent set's members are written in the Window the set first appears in.
#' * Blocks, one row per Window: `acceptance` and `swaps`.
#' * Run events: `events`.
#'
#' Only closed Windows are read (the Store's reader rule): rows of a Window being written, or of one a killed Run
#' never finished, are never returned. A Run that lacks a table (a one-chain Run has no `swaps`) adds no rows.
#'
#' @inheritParams pt_runs
#' @param name A Store table name, listed above.
#' @param runs Run labels to read (see [pt_runs()]). `NULL` reads every Run.
#' @param phase `"sampling"`, `"burn_in"`, or both. Burn-in includes the starting State at iteration 0.
#' @return A lazy table (`tbl_lazy`) with `run` first, then the Store's columns. `run` and `iteration` are the Draw
#'   key of every Trace.
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_table(a, "run_scalars") |>
#'   dplyr::mutate(log_posterior = prior + likelihood) |>
#'   dplyr::select(run, iteration, log_posterior) |>
#'   dplyr::collect()
#' pt_close(a)
pt_table <- function(a, name, runs = NULL, phase = "sampling") {
  con <- analysis_con(a)
  if (!rlang::is_string(name) || !name %in% names(store_tables)) {
    cli::cli_abort(c(
      "No Store table {.val {name}}.",
      i = "Tables: {names(store_tables)}.",
      i = "See {.url https://github.com/EPPIcenter/plasmotrack/blob/master/doc/store-layout.md}."
    ))
  }
  phase <- check_phase(phase)
  selected <- select_runs(a, runs)

  parts <- character()
  for (i in which(selected)) {
    files <- window_files(a$runs[i, ], name, phase)
    if (length(files) == 0) next
    where <- if (store_tables[[name]]$phase) sprintf(" WHERE phase IN (%s)", paste(sql_string(con, phase),
                                                                                       collapse = ", ")) else ""
    parts <- c(parts, sprintf("SELECT %s AS run, * FROM read_parquet(%s)%s", sql_string(con, a$runs$run[[i]]),
                              sql_list(con, files), where))
  }
  sql <- if (length(parts)) paste(parts, collapse = "\nUNION ALL BY NAME\n") else empty_table_sql(con, name)
  dplyr::tbl(con, dbplyr::sql(sql))
}

# A query with no rows and the table's columns, for a table no selected Run has files for.
empty_table_sql <- function(con, name) {
  columns <- store_tables[[name]]$columns
  types <- column_types[columns]
  types[is.na(types)] <- "VARCHAR"
  sprintf("SELECT %s WHERE false", paste(
    c("CAST(NULL AS VARCHAR) AS run", sprintf("CAST(NULL AS %s) AS %s", types, DBI::dbQuoteIdentifier(con, columns))),
    collapse = ", "
  ))
}

check_phase <- function(phase, call = rlang::caller_env()) {
  phases <- c("burn_in", "sampling")
  if (!is.character(phase) || length(phase) == 0 || !all(phase %in% phases)) {
    cli::cli_abort("{.arg phase} must be {.or {.val {phases}}}, or both.", call = call)
  }
  unique(phase)
}

# Logical index of the Runs `runs` names; every Run when NULL.
select_runs <- function(a, runs, call = rlang::caller_env()) {
  if (is.null(runs)) return(rep(TRUE, nrow(a$runs)))
  unknown <- setdiff(runs, a$runs$run)
  if (!is.character(runs) || length(unknown)) {
    cli::cli_abort(c("Unknown Run{?s}: {.val {unknown}}.", i = "Runs: {.val {a$runs$run}}."), call = call)
  }
  a$runs$run %in% runs
}

# The Runs a pooled summary or diagnostic uses: those `runs` names (every Run when NULL) that have sampling Draws.
# Runs with 0 sampling Draws stay in the Run table but are skipped here, with a message; none left is an error.
# Returns Run labels in Run-table order.
runs_with_draws <- function(a, runs = NULL, call = rlang::caller_env()) {
  selected <- select_runs(a, runs, call = call)
  empty <- selected & a$runs$draws == 0
  if (any(empty)) {
    cli::cli_inform("Skipping Run{?s} {.val {a$runs$run[empty]}}: no sampling Draws.")
  }
  if (!any(selected & !empty)) {
    cli::cli_abort("No sampling Draws in the selected Runs.", call = call)
  }
  a$runs$run[selected & !empty]
}
