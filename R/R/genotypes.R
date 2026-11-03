# Latent genotypes and COI: per-Draw lazy tables and pooled summaries, from the Store's `latent_genotypes` Trace.
#
# A Latent genotype is an Allele bitstring: character k (0-based, from the left) is Allele k, "1" means present.
# Only Loci the Infection was genotyped at have rows. Pooled summaries count bitstrings in DuckDB first and expand
# the few distinct ones in R: joining every row to its Alleles first cost 25 GB on Eswatini.

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
