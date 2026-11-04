# Latent genotypes and COI: per-Draw lazy tables and pooled summaries, from the Store's `latent_genotypes` Trace.
#
# A Latent genotype is an Allele bitstring: character k (0-based, from the left) is Allele k, "1" means present.
# Only Loci the Infection was genotyped at have rows. Pooled summaries count bitstrings in DuckDB first and expand
# the few distinct ones in R: joining every row to its Alleles first took 25 GB on a real 457-Infection Store.

# The Alleles present in `alleles`, as a DuckDB list of 0-based indices.
present_alleles_sql <- "list_filter(range(length(alleles)), lambda k: substr(alleles, CAST(k AS INTEGER) + 1, 1) = '1')"

# The number of Alleles present in `alleles`.
allele_count_sql <- "length(alleles) - length(replace(alleles, '1', ''))"

genotype_rows <- function(a, runs = NULL, phase = "sampling") {
  pt_table(a, "latent_genotypes", runs = runs, phase = phase) |>
    dplyr::select("run", "iteration", "infection", "locus", "alleles")
}

#' Per-Draw Latent genotypes and COI
#'
#' Lazy tables with one row per Draw (`run`, `iteration`) per key, from the Store's `latent_genotypes` Trace. Filter
#' them before [dplyr::collect()].
#'
#' A Latent genotype is stored as an Allele bitstring per Infection and Locus, such as `00101000`: character k
#' (0-based, from the left) is Allele k, numbered as in the `allele` column of `pt_table(a, "allele_frequencies")`,
#' and `1` means present. Only Loci an Infection was genotyped at have Latent genotypes; other Loci have no rows.
#'
#' * `pt_alleles()`: `infection`, `locus`, `allele`, one row per Allele present. The bitstrings are expanded in
#'   DuckDB.
#' * `pt_coi_draws()`: `infection`, `coi`, the complexity of infection: the largest number of Alleles present at any
#'   Locus the Infection was genotyped at. An Infection genotyped at no Locus has no rows.
#'
#' @inheritParams pt_table
#' @return A lazy table (`tbl_lazy`) with `run` and `iteration` first, then the columns above.
#' @seealso [pt_coi()], [pt_genotype_presence()] and [pt_genotype_mode()] for summaries pooled over Draws.
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' pt_alleles(a) |> dplyr::filter(infection == "3") |> dplyr::collect()
#' pt_coi_draws(a) |> dplyr::collect()
#' pt_close(a)
pt_alleles <- function(a, runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  genotype_rows(a, runs, phase) |>
    dplyr::mutate(allele = dbplyr::sql(paste0("UNNEST(", present_alleles_sql, ")"))) |>
    dplyr::select("run", "iteration", "infection", "locus", "allele")
}

#' @rdname pt_alleles
#' @export
pt_coi_draws <- function(a, runs = NULL, phase = "sampling") {
  phase <- check_phase(phase)
  select_runs(a, runs)
  genotype_rows(a, runs, phase) |>
    dplyr::mutate(n_present = dbplyr::sql(allele_count_sql)) |>
    dplyr::group_by(.data$run, .data$iteration, .data$infection) |>
    dplyr::summarise(coi = max(.data$n_present, na.rm = TRUE), .groups = "drop") |>
    dplyr::select("run", "iteration", "infection", "coi")
}

# Count rows of a lazy per-Draw table by `keys` in DuckDB (with `run` first when by_run), collect, and attach the
# Draw total each count is a fraction of: the Run's Draws, or every selected Run's. With `count`, the name of a
# logical column, `n` counts only the rows where it is TRUE (0 for a group with none). Rows sorted by Run (Run-table
# order) then keys as strings.
tally_draws <- function(a, rows, keys, by_run, runs, count = NULL) {
  groups <- c(if (by_run) "run", keys)
  n <- if (is.null(count)) rlang::quo(dplyr::n()) else {
    rlang::quo(sum(dplyr::if_else(.data[[!!count]], 1, 0), na.rm = TRUE))
  }
  counts <- rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
    dplyr::summarise(n = !!n, .groups = "drop") |>
    dplyr::collect()
  counts$n <- as.numeric(counts$n)
  counts$draws <- if (by_run) unname(draws_of(a, counts$run)) else rep(sum(draws_of(a, runs)), nrow(counts))
  order_by <- lapply(groups, function(g) if (g == "run") match(counts$run, a$runs$run) else counts[[g]])
  counts[do.call(order, c(unname(order_by), method = "radix")), ]
}

# The type-7 quantile at `p` of a sample given as distinct values and their counts.
quantile7_counts <- function(values, counts, p) {
  o <- order(values)
  values <- values[o]
  ends <- cumsum(counts[o])
  h <- (sum(counts) - 1) * p + 1
  at <- function(i) values[[which(ends >= i)[[1]]]]
  lo <- at(floor(h))
  hi <- at(ceiling(h))
  if (hi == lo) return(lo)
  # stats::quantile()'s own form, so the result matches it to the last bit.
  w <- h - floor(h)
  (1 - w) * lo + w * hi
}

#' Pooled COI
#'
#' Each Infection's complexity of infection (COI, see [pt_coi_draws()]) over the sampling Draws: its median and its
#' distribution.
#'
#' Draws pool as everywhere in the package: the fraction of a COI value is the number of Draws with that value,
#' summed over Runs, divided by the total Draw count, so each Run weighs by its Draws. The median is the type-7
#' quantile (R's default, see [stats::quantile()]) of the pooled per-Draw values, so it can fall halfway between two
#' values. Runs with no sampling Draws are skipped with a message. An Infection genotyped at no Locus has no row.
#'
#' @inheritParams pt_edge_marginals
#' @param by_run `FALSE` to pool every Run; `TRUE` for one row per Run per Infection.
#' @return A tibble of `infection`, `median` and `distribution`, after `run` when `by_run = TRUE`, sorted by `run`
#'   (in [pt_runs()] order), then `infection` as strings. `distribution` is a list column of tibbles of `coi` and
#'   `fraction`, one row per COI value seen, in increasing order.
#' @seealso [pt_coi_draws()] for the per-Draw values.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_coi(a)
#' pt_coi(a)$distribution[[1]]
#' pt_close(a)
pt_coi <- function(a, by_run = FALSE, runs = NULL) {
  check_flag(by_run)
  runs <- runs_with_draws(a, runs)
  counts <- tally_draws(a, pt_coi_draws(a, runs), c("infection", "coi"), by_run, runs)
  keys <- c(if (by_run) "run", "infection")
  group <- do.call(row_key, unname(as.list(counts[keys])))
  first <- !duplicated(group)
  rows <- split(seq_len(nrow(counts)), factor(group, levels = unique(group)))
  out <- counts[first, keys]
  out$median <- vapply(rows, function(i) quantile7_counts(counts$coi[i], counts$n[i], 0.5), numeric(1),
                       USE.NAMES = FALSE)
  out$distribution <- lapply(rows, function(i) {
    tibble::tibble(coi = counts$coi[i], fraction = counts$n[i] / counts$draws[i])
  }) |> unname()
  tibble::as_tibble(out)
}

# Draws of each Latent genotype bitstring per (infection, locus), with `run` first when by_run, counted in DuckDB.
bitstring_counts <- function(a, by_run, runs) {
  tally_draws(a, genotype_rows(a, runs, "sampling"), c("infection", "locus", "alleles"), by_run, runs)
}

#' Pooled Latent genotypes
#'
#' Each Infection's Latent genotype at each Locus it was genotyped at, summarised over the sampling Draws. Loci an
#' Infection was not genotyped at have no rows. Bitstrings and Allele numbers are as in [pt_alleles()].
#'
#' * `pt_genotype_presence()`: `infection`, `locus`, `allele`, `presence`, the fraction of Draws at which the
#'   Allele is in the Infection's Latent genotype. Every Allele of the Locus has a row, 0 when it never appears.
#' * `pt_genotype_mode()`: `infection`, `locus`, `alleles`, `fraction`, `agreement`. `alleles` is the posterior
#'   mode: the bitstring held at the most Draws, ties going to the smaller bitstring as a string (`"0010"` before
#'   `"0100"`). `fraction` is the share of Draws holding it. `agreement` is the share of the Draws with a Latent
#'   genotype at the Locus that hold it; the two differ only if some Draw lacks the Infection's row at the Locus.
#'
#' Draws pool as everywhere in the package: counts sum over Runs and divide by the total Draw count, so each Run
#' weighs by its Draws. Runs with no sampling Draws are skipped with a message. Bitstrings are counted in DuckDB
#' before they are expanded to Alleles, so the summaries stay cheap on large Stores.
#'
#' @inheritParams pt_coi
#' @return A tibble with the columns above, after `run` when `by_run = TRUE`, sorted by `run` (in [pt_runs()]
#'   order), then `infection` and `locus` as strings, then `allele`.
#' @seealso [pt_alleles()] for the per-Draw Alleles.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_genotype_presence(a) |> dplyr::filter(presence > 0)
#' pt_genotype_mode(a)
#' pt_genotype_mode(a, by_run = TRUE)
#' pt_close(a)
pt_genotype_presence <- function(a, by_run = FALSE, runs = NULL) {
  check_flag(by_run)
  runs <- runs_with_draws(a, runs)
  counts <- bitstring_counts(a, by_run, runs)
  keys <- c(if (by_run) "run", "infection", "locus")
  # One row per counted bitstring and Allele of its Locus; presence sums the Draws of the bitstrings holding it.
  n_alleles <- a$loci$n_alleles[match(counts$locus, a$loci$locus)]
  at <- rep(seq_len(nrow(counts)), n_alleles)
  allele <- sequence(n_alleles) - 1
  held <- substr(counts$alleles[at], allele + 1, allele + 1) == "1"
  group <- do.call(row_key, c(unname(as.list(counts[at, keys])), list(allele)))
  total <- rowsum(counts$n[at] * held, factor(group, levels = unique(group)), reorder = FALSE)
  first <- at[!duplicated(group)]
  out <- counts[first, keys]
  out$allele <- allele[!duplicated(group)]
  out$presence <- as.vector(total) / counts$draws[first]
  tibble::as_tibble(out)
}

#' @rdname pt_genotype_presence
#' @export
pt_genotype_mode <- function(a, by_run = FALSE, runs = NULL) {
  check_flag(by_run)
  runs <- runs_with_draws(a, runs)
  counts <- bitstring_counts(a, by_run, runs)
  keys <- c(if (by_run) "run", "infection", "locus")
  group <- do.call(row_key, unname(as.list(counts[keys])))
  group <- match(group, unique(group))
  seen <- as.vector(rowsum(counts$n, group, reorder = FALSE))[group]
  # Within each group, most Draws first, then the smaller bitstring: the first row is the mode.
  ord <- order(group, -counts$n, counts$alleles, method = "radix")
  best <- ord[!duplicated(group[ord])]
  mode <- counts[best, ]
  mode$seen <- seen[best]
  out <- mode[keys]
  out$alleles <- mode$alleles
  out$fraction <- mode$n / mode$draws
  out$agreement <- mode$n / mode$seen
  tibble::as_tibble(out)
}
