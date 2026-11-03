# The test oracle: an independent computation of every quantity from the raw Parquet tables.
#
# Raw tables are loaded with a plain DuckDB `SELECT * FROM read_parquet([...files...])` into data frames (no dbplyr,
# no package code: nanoparquet 0.5.2 misreads parent_sets$set_id). Everything after that is base R, written as
# plain loops over the spec's formulas, so it shares nothing with the package's SQL but the Parquet decoder.
#
#   o <- oracle_read(paths)        same `paths` as pt_read(): named Run directories, or unnamed (label = basename).
#   oracle_table(o, run, table, phase)   data frame of one Run's raw rows over its closed Windows (reader rule
#                                  reimplemented from the _closed/ markers).
#   oracle_draws(o, run, phase)    the Run's Draws: sorted distinct iterations of run_scalars in `phase`.
#   oracle_parent_sets(o, run, phase)  parent_sets rows with each set's observed members (list), n_observed and
#                                  has_latent_parent attached; the building block of every Parent-set quantity.
#   oracle_<quantity>(o, ...)      one function per package function, returning a data frame with the package's
#                                  columns. Add one per ticket, next to the ones here.
#   expect_oracle_equal(actual, expected, keys)  collects `actual`, sorts both by `keys` and compares every numeric
#                                  column element by element at relative 1e-12, other columns exactly.

oracle_read <- function(paths) {
  labels <- names(paths)
  paths <- unname(paths)
  if (is.null(labels)) labels <- basename(paths)
  labels[labels == ""] <- basename(paths)[labels == ""]
  runs <- lapply(seq_along(paths), function(i) {
    markers <- list.files(file.path(paths[[i]], "_closed"), pattern = "^window-\\d+$", full.names = TRUE)
    windows <- do.call(rbind, lapply(markers, function(f) {
      m <- jsonlite::read_json(f)
      data.frame(window = m$window, phase = m$phase)
    }))
    list(label = labels[[i]], path = paths[[i]], windows = windows,
         infections = oracle_parquet(file.path(paths[[i]], "structure", "infections.parquet")))
  })
  names(runs) <- labels
  list(runs = runs, cache = new.env())
}

oracle_parquet <- function(files) {
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  quoted <- paste(vapply(files, function(f) as.character(DBI::dbQuoteString(con, f)), ""), collapse = ", ")
  DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet([%s])", quoted))
}

# One Run's rows of a windowed table over its closed Windows; parent_set_members over every closed Window, rows of
# other tables restricted to `phase`. NULL when the Run has no file for the table.
oracle_table <- function(o, run, table, phase = "sampling") {
  key <- paste(run, table, paste(phase, collapse = "+"))
  if (!is.null(o$cache[[key]])) return(o$cache[[key]])
  r <- o$runs[[run]]
  windows <- if (table == "parent_set_members") r$windows$window else r$windows$window[r$windows$phase %in% phase]
  files <- file.path(r$path, table, sprintf("window-%06d.parquet", windows))
  files <- files[file.exists(files)]
  rows <- if (length(files)) oracle_parquet(files) else NULL
  if (!is.null(rows) && "phase" %in% names(rows)) rows <- rows[rows$phase %in% phase, , drop = FALSE]
  o$cache[[key]] <- rows
  rows
}

oracle_runs <- function(o, runs = NULL) if (is.null(runs)) names(o$runs) else runs

oracle_draws <- function(o, run, phase = "sampling") {
  rows <- oracle_table(o, run, "run_scalars", phase)
  if (is.null(rows)) return(numeric())
  sort(unique(rows$iteration))
}

oracle_infection_ids <- function(o) o$runs[[1]]$infections$infection

# parent_sets rows of one Run in `phase`, each with its set's observed members, n_observed and has_latent_parent.
oracle_parent_sets <- function(o, run, phase = "sampling") {
  sets <- oracle_table(o, run, "parent_sets", phase)
  members <- oracle_table(o, run, "parent_set_members")
  infections <- o$runs[[run]]$infections
  latent_of <- setNames(infections$latent_parent, infections$infection)
  if (is.null(sets) || nrow(sets) == 0) return(NULL)
  observed <- list()
  latent <- logical()
  for (i in seq_len(nrow(members))) {
    k <- paste(members$infection[[i]], members$set_id[[i]], sep = "\r")
    if (is.null(observed[[k]])) observed[[k]] <- character()
    if (is.na(latent[k])) latent[k] <- FALSE
    if (members$member[[i]] == latent_of[[members$infection[[i]]]]) {
      latent[k] <- TRUE
    } else {
      observed[[k]] <- c(observed[[k]], members$member[[i]])
    }
  }
  keys <- paste(sets$infection, sets$set_id, sep = "\r")
  stopifnot(all(keys %in% names(latent)))
  sets$observed <- lapply(keys, function(k) observed[[k]])
  sets$n_observed <- vapply(sets$observed, length, numeric(1))
  sets$has_latent_parent <- unname(latent[keys])
  sets
}

oracle_class_name <- function(n_observed, has_latent_parent) {
  if (n_observed == 0) return(if (has_latent_parent) "latent_only" else "no_parent")
  paste0("observed_", n_observed, if (has_latent_parent) "_latent" else "")
}

oracle_edge_draws <- function(o, runs = NULL, phase = "sampling") {
  out <- list()
  for (run in oracle_runs(o, runs)) {
    sets <- oracle_parent_sets(o, run, phase)
    if (is.null(sets)) next
    total <- list()
    for (i in seq_len(nrow(sets))) {
      for (p in sets$observed[[i]]) {
        k <- paste(sets$iteration[[i]], sets$infection[[i]], p, sep = "\r")
        total[[k]] <- (if (is.null(total[[k]])) 0 else total[[k]]) + sets$prob[[i]]
      }
    }
    if (length(total) == 0) next
    parts <- do.call(rbind, strsplit(names(total), "\r", fixed = TRUE))
    out[[run]] <- data.frame(run = run, iteration = as.numeric(parts[, 1]), child = parts[, 2], parent = parts[, 3],
                             prob = unlist(total, use.names = FALSE))
  }
  oracle_bind(out, data.frame(run = character(), iteration = numeric(), child = character(), parent = character(),
                              prob = numeric()))
}

oracle_infection_draws <- function(o, runs = NULL, phase = "sampling") {
  out <- list()
  for (run in oracle_runs(o, runs)) {
    sets <- oracle_parent_sets(o, run, phase)
    if (is.null(sets)) next
    groups <- split(seq_len(nrow(sets)), paste(sets$iteration, sets$infection, sep = "\r"))
    rows <- lapply(groups, function(idx) {
      q <- 0; e <- 0; h <- 0; aug <- 0; obs <- 0
      for (i in idx) {
        w <- sets$prob[[i]]
        n <- sets$n_observed[[i]]
        l <- sets$has_latent_parent[[i]]
        if (l) q <- q + w
        e <- e + w * n
        if (w > 0) h <- h - w * log(w)
        if (n + l > 1) aug <- aug + w
        if (n > 1) obs <- obs + w
      }
      data.frame(run = run, iteration = sets$iteration[[idx[1]]], infection = sets$infection[[idx[1]]],
                 p_latent_parent = q, expected_observed_parents = e, entropy = h, effective_size = exp(h),
                 multi_source_augmented = aug, multi_source_observed = obs)
    })
    out[[run]] <- do.call(rbind, rows)
  }
  oracle_bind(out, NULL)
}

# Every class of any set with weight in the Runs read, for every Infection at every Draw.
oracle_class_draws <- function(o, runs = NULL, phase = "sampling") {
  runs <- oracle_runs(o, runs)
  all_sets <- lapply(setNames(runs, runs), oracle_parent_sets, o = o, phase = phase)
  classes <- data.frame(class = character(), n_observed = numeric(), has_latent_parent = logical())
  for (sets in all_sets) {
    for (i in seq_len(NROW(sets))) {
      name <- oracle_class_name(sets$n_observed[[i]], sets$has_latent_parent[[i]])
      if (!name %in% classes$class) {
        classes <- rbind(classes, data.frame(class = name, n_observed = sets$n_observed[[i]],
                                             has_latent_parent = sets$has_latent_parent[[i]]))
      }
    }
  }
  out <- list()
  for (run in runs) {
    sets <- all_sets[[run]]
    share <- list()
    for (i in seq_len(NROW(sets))) {
      k <- paste(sets$iteration[[i]], sets$infection[[i]],
                 oracle_class_name(sets$n_observed[[i]], sets$has_latent_parent[[i]]), sep = "\r")
      share[[k]] <- (if (is.null(share[[k]])) 0 else share[[k]]) + sets$prob[[i]]
    }
    rows <- list()
    for (it in oracle_draws(o, run, phase)) {
      for (inf in oracle_infection_ids(o)) {
        for (c in seq_len(nrow(classes))) {
          k <- paste(it, inf, classes$class[[c]], sep = "\r")
          rows[[length(rows) + 1]] <- data.frame(
            run = run, iteration = it, infection = inf, class = classes$class[[c]],
            n_observed = classes$n_observed[[c]], has_latent_parent = classes$has_latent_parent[[c]],
            share = if (is.null(share[[k]])) 0 else share[[k]]
          )
        }
      }
    }
    out[[run]] <- do.call(rbind, rows)
  }
  oracle_bind(out, NULL)
}

oracle_out_degree_draws <- function(o, runs = NULL, phase = "sampling") {
  edges <- oracle_edge_draws(o, runs, phase)
  out <- list()
  for (run in oracle_runs(o, runs)) {
    for (it in oracle_draws(o, run, phase)) {
      for (inf in oracle_infection_ids(o)) {
        mine <- edges$run == run & edges$iteration == it & edges$parent == inf
        out[[length(out) + 1]] <- data.frame(run = run, iteration = it, infection = inf,
                                             out_degree = sum(edges$prob[mine]))
      }
    }
  }
  oracle_bind(out, NULL)
}

# Pooled Edge marginals: sum of an edge's per-Draw weight over every sampling Draw, divided by the Draws (pooled:
# of every Run with Draws; by_run: of the Run). Lag from `observation_time`, a named vector (NA by default).
oracle_edge_marginals <- function(o, floor = 0, by_run = FALSE, runs = NULL, observation_time = NULL) {
  runs <- oracle_runs(o, runs)
  draws <- vapply(runs, function(r) length(oracle_draws(o, r)), numeric(1))
  runs <- runs[draws > 0]
  draws <- draws[draws > 0]
  edges <- oracle_edge_draws(o, runs)
  keys <- unique(paste(edges$child, edges$parent, sep = "\r"))
  total <- matrix(0, length(keys), length(runs), dimnames = list(keys, runs))
  for (i in seq_len(nrow(edges))) {
    k <- paste(edges$child[[i]], edges$parent[[i]], sep = "\r")
    total[k, edges$run[[i]]] <- total[k, edges$run[[i]]] + edges$prob[[i]]
  }
  lag_of <- function(child, parent) {
    if (is.null(observation_time)) return(NA_real_)
    unname(observation_time[child] - observation_time[parent])
  }
  out <- list()
  for (k in keys) {
    parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
    per_run <- total[k, ] / draws
    if (by_run) {
      if (max(per_run) < floor) next
      for (j in seq_along(runs)) {
        out[[length(out) + 1]] <- data.frame(run = runs[[j]], child = parts[[1]], parent = parts[[2]],
                                             marginal = per_run[[j]], lag = lag_of(parts[[1]], parts[[2]]))
      }
    } else {
      marginal <- sum(total[k, ]) / sum(draws)
      if (marginal < floor) next
      out[[length(out) + 1]] <- data.frame(child = parts[[1]], parent = parts[[2]], marginal = marginal,
                                           lag = lag_of(parts[[1]], parts[[2]]))
    }
  }
  oracle_bind(out, NULL)
}

# Pooled Parent sets (pt_parent_sets_pooled): sets keyed by sorted observed members and has_latent_parent, weight =
# summed stored weight / Draws, ranked by a pairwise comparison (weight, then the members list with the Latent parent
# after every id and a prefix first), Credible collection by a forward scan.
oracle_parent_sets_pooled <- function(o, level = 0.8, by_run = FALSE, runs = NULL) {
  runs <- oracle_runs(o, runs)
  draws <- vapply(runs, function(r) length(oracle_draws(o, r)), numeric(1))
  runs <- runs[draws > 0]
  draws <- draws[draws > 0]
  sets <- list()
  for (j in seq_along(runs)) {
    rows <- oracle_parent_sets(o, runs[[j]])
    group <- if (by_run) runs[[j]] else ""
    for (i in seq_len(NROW(rows))) {
      observed <- sort(rows$observed[[i]], method = "radix")
      k <- paste(group, rows$infection[[i]], paste(observed, collapse = ","), rows$has_latent_parent[[i]], sep = "|")
      if (is.null(sets[[k]])) {
        sets[[k]] <- list(run = group, infection = rows$infection[[i]], members = observed,
                          has_latent_parent = rows$has_latent_parent[[i]], total = 0,
                          draws = if (by_run) draws[[j]] else sum(draws))
      }
      sets[[k]]$total <- sets[[k]]$total + rows$prob[[i]]
    }
  }
  list_before <- function(x, y) {
    x <- c(x$members, if (x$has_latent_parent) NA)
    y <- c(y$members, if (y$has_latent_parent) NA)
    for (k in seq_len(min(length(x), length(y)))) {
      if (is.na(x[k]) && is.na(y[k])) next
      if (is.na(x[k])) return(FALSE)
      if (is.na(y[k])) return(TRUE)
      if (x[k] != y[k]) return(sort(c(x[k], y[k]), method = "radix")[1] == x[k])
    }
    length(x) < length(y)
  }
  for (k in names(sets)) sets[[k]]$weight <- sets[[k]]$total / sets[[k]]$draws
  out <- list()
  for (g in unique(vapply(sets, function(s) paste(s$run, s$infection, sep = "|"), ""))) {
    mine <- Filter(function(s) paste(s$run, s$infection, sep = "|") == g, sets)
    rank <- vapply(mine, function(s) {
      1 + sum(vapply(mine, function(t) t$weight > s$weight || (t$weight == s$weight && list_before(t, s)), TRUE))
    }, 1)
    mine <- mine[order(rank)]
    cumulative <- 0
    reached <- NA
    rows <- list()
    for (r in seq_along(mine)) {
      cumulative <- cumulative + mine[[r]]$weight
      if (is.na(reached) && cumulative >= level - 1e-12) reached <- r
      rows[[r]] <- data.frame(run = mine[[r]]$run, infection = mine[[r]]$infection, members = I(list(mine[[r]]$members)),
                              n_observed = length(mine[[r]]$members), has_latent_parent = mine[[r]]$has_latent_parent,
                              class = oracle_class_name(length(mine[[r]]$members), mine[[r]]$has_latent_parent),
                              weight = mine[[r]]$weight, rank = r, cumulative_mass = cumulative)
    }
    rows <- do.call(rbind, rows)
    rows$credible <- if (is.na(reached)) TRUE else rows$rank <= reached | rows$weight == rows$weight[reached]
    out[[g]] <- rows
  }
  out <- oracle_bind(out, NULL)
  out$members <- unclass(out$members)
  if (!by_run) out$run <- NULL
  out
}

oracle_bind <- function(parts, empty) {
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) return(empty)
  rows <- do.call(rbind, unname(parts))
  rownames(rows) <- NULL
  rows
}

# Compare a package result (lazy or eager) with the oracle's: same columns in the same order, same rows once both
# are sorted by `keys` (string order), numbers equal element by element at relative `tolerance`.
expect_oracle_equal <- function(actual, expected, keys, tolerance = 1e-12) {
  actual <- as.data.frame(dplyr::collect(actual))
  expected <- as.data.frame(expected)
  expect_equal(names(actual), names(expected))
  expect_equal(nrow(actual), nrow(expected))
  sort_rows <- function(x) {
    x <- x[do.call(order, c(unname(as.list(x[keys])), method = "radix")), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  actual <- sort_rows(actual)
  expected <- sort_rows(expected)
  for (column in names(expected)) {
    a <- actual[[column]]
    e <- expected[[column]]
    if (is.numeric(e)) {
      expect_true(is.numeric(a), label = column)
      ok <- (is.na(a) & is.na(e)) | (!is.na(a) & !is.na(e) & (a == e | abs(a - e) <= tolerance * abs(e)))
      expect_true(all(ok), label = sprintf("%s within relative %g of the oracle (worst row %s)", column, tolerance,
                                           paste(which(!ok)[1])))
    } else {
      expect_identical(as.vector(a), as.vector(e), label = column)
    }
  }
  invisible(actual)
}
