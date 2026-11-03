# Scale check (manual): every exported function on a real Store under a 4 GB DuckDB memory cap, with timings and
# peak memory. Not part of R CMD check or CI (tests/manual is Rbuildignored); skipped unless
# PLASMOTRACKR_SCALE_STORE is set.
#
# Environment:
#   PLASMOTRACKR_SCALE_STORE    the Store: one parent directory whose subdirectories are Runs, or Run directories
#                               separated by ":" (required).
#   PLASMOTRACKR_SCALE_INPUT    the input JSON, passed to pt_read(input =) (optional: the Store's copy or the
#                               manifest's path otherwise).
#   PLASMOTRACKR_SCALE_TRUTH    a simulator truth JSON (optional). Without it the scoring functions run against a
#                               stand-in truth built from the Store itself: the MAP graph's edges, the posterior-mode
#                               genotypes as one strain per Infection, and the input's observation times.
#   PLASMOTRACKR_SCALE_THREADS  DuckDB threads (default 2).
#   PLASMOTRACKR_SCALE_OUT      a CSV file for the timings table (optional).
#
# Run from the repository root, with the package installed:
#   PLASMOTRACKR_SCALE_STORE=/path/to/store Rscript -e \
#     'testthat::test_file("R/tests/manual/test-scale.R", package = "plasmotrackr", load_package = "installed")'
#
# Each function is called once, in order, on one Analysis opened with memory_limit = "4GB". A lazy result is
# materialised in DuckDB (`compute()` into a temporary table, then dropped), so its query runs in full under the cap
# and `rows` is its row count; nothing that size is pulled into R. Eager results report their row count (edges for an
# igraph, variables for a draws_array, quantities for a pt_convergence).
#
# Memory: `peak_rss_mb` is the R process's peak resident set during the call (Linux only: the peak is reset through
# /proc/self/clear_refs before each call, and DuckDB's buffers count since they live in the process); `r_heap_mb` is
# R's own peak heap from gc(). Each call must finish without error and peak below 4 GiB of resident memory.

scale_store <- Sys.getenv("PLASMOTRACKR_SCALE_STORE")
cap <- "4GB"
cap_mb <- 4096

scale_paths <- function(store) {
  paths <- strsplit(store, ":", fixed = TRUE)[[1]]
  if (length(paths) == 1) paths else stats::setNames(paths, basename(paths))
}

reset_peak <- function() {
  if (file.exists("/proc/self/clear_refs")) try(writeLines("5", "/proc/self/clear_refs"), silent = TRUE)
  invisible(gc(reset = TRUE, full = TRUE))
}

peak_rss_mb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  hwm <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  if (!length(hwm)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", hwm)) / 1024
}

r_heap_mb <- function() sum(gc()[, 6])

# Runs a lazy result's query in full inside DuckDB and returns its row count.
materialise <- function(x) {
  con <- dbplyr::remote_con(x)
  dplyr::compute(x, name = "pt_scale_result", temporary = TRUE, overwrite = TRUE)
  on.exit(DBI::dbRemoveTable(con, "pt_scale_result"), add = TRUE)
  DBI::dbGetQuery(con, "SELECT count(*) AS n FROM pt_scale_result")$n
}

size_of <- function(x) {
  if (inherits(x, "tbl_lazy")) return(materialise(x))
  if (inherits(x, "pt_convergence")) return(nrow(x$quantities))
  if (inherits(x, "pt_analysis")) return(sum(pt_runs(x)$draws))
  if (inherits(x, "pt_truth")) return(nrow(x$edges))
  if (inherits(x, "draws")) return(posterior::nvariables(x))
  if (inherits(x, "igraph")) return(igraph::ecount(x))
  if (is.data.frame(x)) return(nrow(x))
  NA_real_
}

# A truth JSON for an Analysis with no simulation behind it: the MAP graph's observed edges, each Infection's
# posterior-mode genotype as its one strain, the input's observation times. Scoring then runs in full on real sizes.
write_stand_in_truth <- function(a, path) {
  infections <- pt_infections(a)
  loci <- pt_loci(a)
  map <- pt_map_graph(a)
  edges <- do.call(rbind, Map(function(child, members) {
    if (length(members)) data.frame(src = members, dest = child) else NULL
  }, map$infection, map$members))
  mode <- pt_genotype_mode(a)
  bits <- function(s) as.integer(strsplit(s, "")[[1]])
  nodes <- lapply(seq_len(nrow(infections)), function(i) {
    id <- infections$infection[[i]]
    rows <- mode[mode$infection == id, ]
    strain <- lapply(seq_len(nrow(rows)), function(j) list(locus = rows$locus[[j]], genotype = I(bits(rows$alleles[[j]]))))
    time <- infections$observation_time[[i]]
    list(id = id, latent_genotype = if (length(strain)) list(strain) else list(), infection_duration = 10,
         infection_time = if (is.na(time)) 0 else time - 10, observation_time = if (is.na(time)) 10 else time,
         symptomatic = TRUE)
  })
  json <- list(
    network = if (is.null(edges)) list() else lapply(seq_len(nrow(edges)), function(i) as.list(edges[i, ])),
    nodes = nodes,
    loci = lapply(seq_len(nrow(loci)), function(i) {
      list(locus = loci$locus[[i]], allele_freqs = I(rep(1 / loci$n_alleles[[i]], loci$n_alleles[[i]])),
           num_alleles = loci$n_alleles[[i]])
    })
  )
  jsonlite::write_json(json, path, auto_unbox = TRUE, digits = NA)
  path
}

test_that("every exported function completes on a real Store under a 4 GB cap", {
  skip_if(!nzchar(scale_store), "PLASMOTRACKR_SCALE_STORE is not set")
  skip_on_ci()

  paths <- scale_paths(scale_store)
  input <- Sys.getenv("PLASMOTRACKR_SCALE_INPUT")
  input <- if (nzchar(input)) input else NULL
  threads <- as.integer(Sys.getenv("PLASMOTRACKR_SCALE_THREADS", "2"))
  temp <- withr::local_tempdir("pt-scale-")

  results <- list()
  failures <- character()
  run <- function(name, expr) {
    reset_peak()
    start <- proc.time()[["elapsed"]]
    value <- tryCatch(force(expr), error = function(e) e)
    rows <- if (inherits(value, "error")) NA_real_ else tryCatch(size_of(value), error = function(e) e)
    seconds <- proc.time()[["elapsed"]] - start
    if (inherits(rows, "error")) {
      value <- rows
      rows <- NA_real_
    }
    peak <- peak_rss_mb()
    error <- if (inherits(value, "error")) conditionMessage(value) else NA_character_
    results[[length(results) + 1]] <<- data.frame(
      function_name = name, seconds = round(seconds, 2), peak_rss_mb = round(peak), r_heap_mb = round(r_heap_mb()),
      rows = rows, error = error
    )
    cat(sprintf("%-40s %9.2f s %8s MB %12s rows%s\n", name, seconds, format(round(peak)), format(rows),
                if (is.na(error)) "" else paste0("  ERROR: ", error)))
    if (!is.na(error)) failures[[length(failures) + 1]] <<- paste0(name, ": ", error)
    if (!is.na(peak) && peak > cap_mb) failures[[length(failures) + 1]] <<- sprintf("%s: peak %.0f MB", name, peak)
    invisible(if (inherits(value, "error")) NULL else value)
  }
  open <- function(epi = NULL) {
    suppressWarnings(pt_read(paths, epi = epi, input = input, memory_limit = cap, temp_directory = temp,
                             threads = threads))
  }

  first <- run("pt_read", open())
  if (is.null(first)) fail(paste(failures, collapse = "\n"))
  ids <- pt_infections(first)$infection
  run("pt_close", pt_close(first))

  # Stand-in strata: Infections dealt round-robin into three groups.
  epi <- data.frame(infection = ids, region = c("a", "b", "c")[(seq_along(ids) - 1) %% 3 + 1])
  a <- open(epi)
  withr::defer(pt_close(a))
  expect_true(grepl("GiB|GB", a$engine$settings$memory_limit))

  run("pt_runs", pt_runs(a))
  run("pt_infections", pt_infections(a))
  run("pt_loci", pt_loci(a))
  run("pt_candidates", pt_candidates(a))
  tables <- c("run_scalars", "infections", "ordering", "allele_frequencies", "latent_genotypes",
              "latent_parent_genotypes", "parent_sets", "parent_set_members", "chains", "acceptance", "swaps",
              "events")
  for (table in tables) run(sprintf("pt_table(%s)", table), pt_table(a, table))

  # Blocks
  run("pt_acceptance", pt_acceptance(a))
  run("pt_swaps", pt_swaps(a))

  # Per-Draw Parent sets and Edge marginals
  run("pt_edge_draws", pt_edge_draws(a))
  run("pt_infection_draws", pt_infection_draws(a))
  run("pt_class_draws", pt_class_draws(a))
  run("pt_out_degree_draws", pt_out_degree_draws(a))
  run("pt_edge_marginals", pt_edge_marginals(a))
  run("pt_edge_marginals(by_run)", pt_edge_marginals(a, by_run = TRUE))

  # Pooled Parent sets
  run("pt_parent_sets_pooled", pt_parent_sets_pooled(a))
  run("pt_map_graph", pt_map_graph(a))
  run("pt_map_out_degree", pt_map_out_degree(a))
  run("pt_map_out_degree(by)", pt_map_out_degree(a, by = "region"))
  run("pt_out_degree_range", pt_out_degree_range(a))
  run("pt_out_degree_range(by)", pt_out_degree_range(a, by = "region"))
  run("pt_co_inclusion", pt_co_inclusion(a))

  # Functionals
  run("pt_mean_out_degree", pt_mean_out_degree(a))
  run("pt_mean_out_degree(by)", pt_mean_out_degree(a, by = "region"))
  run("pt_mean_p_latent_parent", pt_mean_p_latent_parent(a))
  run("pt_composition", pt_composition(a))
  run("pt_multi_source", pt_multi_source(a))
  run("pt_entropy", pt_entropy(a))
  functional <- dplyr::collect(pt_mean_out_degree(a))
  run("pt_interval", pt_interval(functional))

  # Out-degree posteriors
  run("pt_out_degree_dist", pt_out_degree_dist(a))
  run("pt_out_degree_dist(pmf)", pt_out_degree_dist(a, pmf = TRUE))
  run("pt_realized_out_degree", pt_realized_out_degree(a))
  run("pt_realized_out_degree(by)", pt_realized_out_degree(a, by = "region"))
  set.seed(1)
  run("pt_sample_networks(100)", pt_sample_networks(a, 100))

  # Flows
  run("pt_attribution", pt_attribution(a, "region"))
  run("pt_enrichment", pt_enrichment(a, "region"))
  run("pt_net_flow", pt_net_flow(a, "region"))
  run("pt_reassign_long_lags", pt_reassign_long_lags(a, max_lag = 60))

  # Genotypes and COI
  run("pt_alleles", pt_alleles(a))
  run("pt_coi_draws", pt_coi_draws(a))
  run("pt_coi", pt_coi(a))
  run("pt_genotype_presence", pt_genotype_presence(a))
  run("pt_genotype_mode", pt_genotype_mode(a))

  # Ordering
  run("pt_order_pair_draws", pt_order_pair_draws(a))
  run("pt_order_pairs", pt_order_pairs(a))
  run("pt_infection_times", pt_infection_times(a))

  # Summaries, convergence, export
  run("pt_infection_summary", pt_infection_summary(a))
  run("pt_convergence", pt_convergence(a))
  run("pt_as_draws", pt_as_draws(a))
  if (requireNamespace("igraph", quietly = TRUE)) {
    run("pt_as_igraph", pt_as_igraph(a))
    run("pt_as_igraph(map)", pt_as_igraph(a, graph = "map"))
  }
  if (requireNamespace("tidygraph", quietly = TRUE)) run("pt_as_tbl_graph", pt_as_tbl_graph(a))

  # Truth and scoring
  truth_path <- Sys.getenv("PLASMOTRACKR_SCALE_TRUTH")
  if (!nzchar(truth_path)) truth_path <- write_stand_in_truth(a, file.path(temp, "truth.json"))
  truth <- run("pt_read_truth", pt_read_truth(truth_path))
  run("pt_score_edges", pt_score_edges(a, truth))
  run("pt_score_edges(transitive, curve)", pt_score_edges(a, truth, edges = "transitive", curve = TRUE))
  marginals <- pt_edge_marginals(a)
  scores <- data.frame(a = marginals$parent, b = marginals$child, score = marginals$marginal)
  run("pt_score_pairs", pt_score_pairs(scores, truth, ids))
  run("pt_score_parent_sets", pt_score_parent_sets(a, truth))
  run("pt_score_map", pt_score_map(a, truth))
  run("pt_score_genotypes", pt_score_genotypes(a, truth))
  run("pt_score_coi", pt_score_coi(a, truth))
  run("pt_fit_summary", pt_fit_summary(a, truth))

  timings <- do.call(rbind, results)
  out <- Sys.getenv("PLASMOTRACKR_SCALE_OUT")
  if (nzchar(out)) utils::write.csv(timings, out, row.names = FALSE)
  cat(sprintf("\nTotal: %.1f s, highest peak %s MB\n", sum(timings$seconds),
              format(max(timings$peak_rss_mb, na.rm = FALSE))))

  exported <- getNamespaceExports("plasmotrackr")
  called <- unique(sub("\\(.*", "", timings$function_name))
  expect_setequal(called, exported)
  expect(length(failures) == 0, paste(c("Scale check failures:", failures), collapse = "\n"))
})
