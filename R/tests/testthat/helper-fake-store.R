# fake_store(): a small, valid Run store written from R data frames, for edge cases the golden Stores lack.
#
# Every table the package requires is written, with defaults for anything not given, so pt_read() accepts the Store.
# Parquet is written through DuckDB COPY with the engine's column types (int64 counts, double values, string ids).
#
# Iterations follow the engine's grid: iteration 0 (the starting State), burn-in at thin, 2 thin, ..., burn_in, and
# sampling at burn_in + thin, ..., burn_in + samples. Windows hold `window` iterations counted from each Phase's
# start; burn-in Window 0 also holds iteration 0, so a Run always has at least one burn-in Window. Every Window gets a
# `_closed/` marker unless listed in `unclosed`.
#
# Arguments (all but `path` optional):
#   path         Run directory to create. Default: a fresh directory under a temp dir removed when `.local_envir`
#                exits. Put several Runs under one parent to read them with pt_read(parent).
#   infections   Infection ids, in any order. Default: every id named in `parent_sets`.
#   parent_sets  Data frame of iteration, infection, parents, latent, prob: one row per Parent set with weight at
#                a written iteration. `parents` is the observed members, either a list column of character vectors
#                or one string with ids separated by "," ("" for none). `latent` (default FALSE) is whether the set
#                includes the Infection's Latent parent. Set ids are assigned per Infection in order of first
#                appearance (by iteration, then row), and the members are written in the Window the set first
#                appears in, as the engine does. Rows are written as given: drop rows yourself to mimic the
#                floor. Every written iteration with no row here (iteration 0, burn-in, or all of them when
#                `parent_sets` is NULL) gets the default: each Infection has its Latent parent alone, prob 1.
#   symptomatic  Logical, recycled over `infections`. Default TRUE.
#   loci         Named integer vector: Locus id -> Allele count. Default c(L1 = 2).
#   burn_in, samples, thin, window
#                The Run's -b, -s, -t and --window. `samples` defaults to the last iteration of `parent_sets`
#                minus `burn_in` (or 4 with no `parent_sets`). `window` defaults to `samples` (one sampling Window).
#   tables       Named list of data frames that replace the default rows of other Traces (run_scalars, infections,
#                ordering, allele_frequencies, latent_genotypes, chains). Give the Store's columns except phase and
#                window, which are filled in from `iteration`.
#   manifest     Named list of manifest fields to override (such as null_model = TRUE, input_sha256 = list(...)).
#   unclosed     Window numbers written without a marker (the reader rule must skip them).
#
# Returns the Run directory, invisibly. Example:
#   sets <- data.frame(iteration = c(1, 1, 2), infection = c("2", "2", "2"), parents = c("10", "", "10"),
#                      latent = c(FALSE, TRUE, FALSE), prob = c(0.7, 0.3, 1))
#   path <- fake_store(parent_sets = sets, infections = c("10", "2"))
fake_store <- function(path = NULL, infections = NULL, parent_sets = NULL, symptomatic = TRUE, loci = c(L1 = 2L),
                       burn_in = 0, samples = NULL, thin = 1, window = NULL, tables = list(), manifest = list(),
                       unclosed = integer(), .local_envir = parent.frame()) {
  if (is.null(path)) path <- file.path(withr::local_tempdir(.local_envir = .local_envir), "fake")
  if (!is.null(parent_sets)) parent_sets <- fake_parent_set_rows(parent_sets)
  if (is.null(infections)) {
    stopifnot(!is.null(parent_sets))
    infections <- unique(c(parent_sets$infection, unlist(parent_sets$parents)))
  }
  infections <- as.character(infections)
  if (is.null(samples)) samples <- if (is.null(parent_sets)) 4 else max(parent_sets$iteration) - burn_in
  if (is.null(window)) window <- samples
  stopifnot(burn_in %% thin == 0, samples %% thin == 0, window %% thin == 0, samples > 0)

  iterations <- c(0, seq_len(burn_in / thin) * thin, burn_in + seq_len(samples / thin) * thin)
  grid <- fake_window_grid(iterations, burn_in, window)
  # Iterations without any parent_sets row get the default: every Infection has its Latent parent alone.
  bare <- setdiff(iterations, parent_sets$iteration)
  if (length(bare)) {
    parent_sets <- rbind(parent_sets, fake_parent_set_rows(data.frame(
      iteration = rep(bare, each = length(infections)), infection = rep(infections, length(bare)),
      parents = "", latent = TRUE, prob = 1
    )))
  }
  stopifnot(all(parent_sets$iteration %in% iterations))

  dir.create(path, recursive = TRUE)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  write <- function(rows, table, file) {
    dir.create(file.path(path, table), showWarnings = FALSE)
    fake_write_parquet(con, rows, file.path(path, table, file))
  }
  # A windowed table: phase and window from the iteration grid, one file per Window with rows.
  write_trace <- function(rows, table) {
    rows <- as.data.frame(rows)
    at <- match(rows$iteration, grid$iteration)
    stopifnot(!anyNA(at))
    rows <- cbind(rows["iteration"], phase = grid$phase[at], window = grid$window[at],
                  rows[setdiff(names(rows), "iteration")])
    for (w in unique(rows$window)) {
      write(rows[rows$window == w, , drop = FALSE], table, sprintf("window-%06d.parquet", w))
    }
  }
  each_iteration <- function(...) {
    cols <- list(...)
    n <- length(cols[[1]])
    cbind(iteration = rep(iterations, each = n), as.data.frame(lapply(cols, rep, times = length(iterations))))
  }

  locus_rows <- data.frame(locus = rep(names(loci), loci), allele = unlist(lapply(loci, function(k) seq_len(k) - 1)))
  defaults <- list(
    run_scalars = data.frame(iteration = iterations, mean_coi = 1.5, mean_strains_tx = 1.2,
                             parent_set_size_prob = 0.5, prior = -1, likelihood = -2, observation = -0.5,
                             transmission = -1.5, source = 0),
    infections = each_iteration(infection = infections, expected_false_positives = 0.01,
                                expected_false_negatives = 0.01, infection_duration = 10),
    ordering = each_iteration(infection = infections, rank = seq_along(infections) - 1),
    allele_frequencies = each_iteration(locus = locus_rows$locus, allele = locus_rows$allele,
                                        frequency = 1 / loci[locus_rows$locus]),
    latent_genotypes = each_iteration(
      infection = rep(infections, each = length(loci)), locus = rep(names(loci), length(infections)),
      alleles = rep(vapply(loci, function(k) paste0("1", strrep("0", k - 1)), character(1)), length(infections))
    )
  )
  defaults[names(tables)] <- tables
  for (table in names(defaults)) write_trace(defaults[[table]], table)

  # Parent sets: ids per Infection by first appearance; members in the Window of first appearance.
  parent_sets <- parent_sets[order(parent_sets$iteration, seq_len(nrow(parent_sets))), ]
  key <- vapply(seq_len(nrow(parent_sets)), function(i) {
    paste(c(parent_sets$infection[[i]], sort(parent_sets$parents[[i]], method = "radix"),
            if (parent_sets$latent[[i]]) "\001latent"), collapse = "\002")
  }, character(1))
  first <- !duplicated(key)
  set_id <- integer(nrow(parent_sets))
  for (inf in unique(parent_sets$infection)) {
    own <- parent_sets$infection == inf
    set_id[own] <- match(key[own], unique(key[own])) - 1L
  }
  write_trace(data.frame(iteration = parent_sets$iteration, infection = parent_sets$infection, set_id = set_id,
                         prob = parent_sets$prob), "parent_sets")
  members <- do.call(rbind, lapply(which(first), function(i) {
    m <- c(parent_sets$parents[[i]], if (parent_sets$latent[[i]]) paste0(parent_sets$infection[[i]], "_copy"))
    data.frame(infection = rep(parent_sets$infection[[i]], length(m)), set_id = rep(set_id[[i]], length(m)),
               member = m, window = rep(grid$window[match(parent_sets$iteration[[i]], grid$iteration)], length(m)))
  }))
  for (w in unique(members$window)) {
    write(members[members$window == w, c("infection", "set_id", "member")], "parent_set_members",
          sprintf("window-%06d.parquet", w))
  }

  # Blocks and events, one Window at a time.
  windows <- unique(grid[c("window", "phase")])
  acceptance <- data.frame(window = windows$window, phase = windows$phase, sampler = "FakeSampler/mean_coi",
                           accepted = 1, proposed = 2)
  for (w in windows$window) write(acceptance[acceptance$window == w, ], "acceptance", sprintf("window-%06d.parquet", w))
  last_of <- tapply(grid$iteration, grid$window, max)
  events <- rbind(
    data.frame(iteration = 0, phase = "burn_in", window = 0, type = "run_started",
               payload = sprintf('{"burn_in":%d,"samples":%d}', as.integer(burn_in), as.integer(samples))),
    data.frame(iteration = as.numeric(last_of), phase = windows$phase, window = windows$window, type = "window_closed",
               payload = '{"sampling_iterations":1,"seconds":1,"iterations_per_second":1}'),
    data.frame(iteration = max(iterations), phase = "sampling", window = max(windows$window), type = "run_ended",
               payload = '{"reason":"completed","seconds":1}')
  )
  events$wall_time <- "2026-11-03T00:00:00.000Z"
  events <- events[c("iteration", "phase", "window", "wall_time", "type", "payload")]
  for (w in windows$window) write(events[events$window == w, ], "events", sprintf("window-%06d.parquet", w))

  # structure/, markers, manifest.
  write(data.frame(infection = infections, latent_parent = paste0(infections, "_copy"),
                   symptomatic = rep_len(as.logical(symptomatic), length(infections))), "structure",
        "infections.parquet")
  write(locus_rows, "structure", "loci.parquet")
  write(data.frame(sampler = "FakeSampler/mean_coi", type = "FakeSampler", target = "mean_coi"), "structure",
        "samplers.parquet")
  dir.create(file.path(path, "_closed"))
  for (w in setdiff(windows$window, unclosed)) {
    its <- grid$iteration[grid$window == w]
    first_it <- if (w == 0) 0 else max(grid$iteration[grid$window < w]) + 1
    jsonlite::write_json(list(first_iteration = first_it, last_iteration = max(its), phase = windows$phase[windows$window == w],
                              window = w), file.path(path, "_closed", sprintf("window-%06d", w)), auto_unbox = TRUE)
  }
  fields <- list(
    burn_in = burn_in, chains = 1, cli_arguments = list("<cli>"), engine_version = "fake", git_commit = "fake",
    gradient = 0, input_sha256 = list(asymptomatic_idp = "fake", input = "fake", symptomatic_idp = "fake"),
    null_model = FALSE, parent_set_floor = 1e-9, run_id = "fake", samples = samples, seed = 1,
    started_at = "2026-11-03T00:00:00.000Z", thin = thin, window = window
  )
  fields[names(manifest)] <- manifest
  jsonlite::write_json(fields, file.path(path, "manifest.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE)
  invisible(path)
}

# parent_sets rows with `parents` as a list of character vectors and `latent` logical.
fake_parent_set_rows <- function(rows) {
  rows <- as.data.frame(rows)
  if (is.null(rows$latent)) rows$latent <- FALSE
  if (!is.list(rows$parents)) {
    rows$parents <- lapply(strsplit(as.character(rows$parents), ",", fixed = TRUE), function(p) p[p != ""])
  }
  rows$parents <- lapply(rows$parents, as.character)
  rows$infection <- as.character(rows$infection)
  rows$latent <- as.logical(rows$latent)
  rows
}

# Phase and Window of each iteration of the grid.
fake_window_grid <- function(iterations, burn_in, window) {
  phase <- ifelse(iterations <= burn_in, "burn_in", "sampling")
  burn_windows <- max(1, ceiling(burn_in / window))
  w <- ifelse(phase == "burn_in", pmax(ceiling(iterations / window), 1) - 1,
              burn_windows + ceiling((iterations - burn_in) / window) - 1)
  data.frame(iteration = iterations, phase = phase, window = w)
}

# Write a data frame to Parquet through DuckDB COPY, casting columns to the Store's types.
fake_write_parquet <- function(con, rows, file) {
  rows <- as.data.frame(rows)
  name <- paste0("fake_", as.integer(stats::runif(1, 0, 1e9)))
  duckdb::duckdb_register(con, name, rows)
  on.exit(duckdb::duckdb_unregister(con, name))
  int64 <- c("iteration", "window", "set_id", "allele", "rank", "temperature_slot", "replica", "accepted",
             "proposed", "temperature_pair", "attempted")
  cols <- vapply(names(rows), function(column) {
    id <- as.character(DBI::dbQuoteIdentifier(con, column))
    type <- if (column %in% int64) "BIGINT" else if (is.logical(rows[[column]])) "BOOLEAN"
            else if (is.numeric(rows[[column]])) "DOUBLE" else "VARCHAR"
    sprintf("CAST(%s AS %s) AS %s", id, type, id)
  }, character(1))
  DBI::dbExecute(con, sprintf("COPY (SELECT %s FROM %s) TO %s (FORMAT parquet)", paste(cols, collapse = ", "),
                              name, as.character(DBI::dbQuoteString(con, file))))
}
