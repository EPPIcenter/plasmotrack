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

# pt_score_edges() with direct truth: `true_edges` is a data frame of parent, child among the Analysis's Infections,
# written by hand (the oracle does not read truth files). Scores from oracle_edge_marginals(); every metric by loops.
oracle_score_edges <- function(o, true_edges) {
  ids <- oracle_infection_ids(o)
  marginals <- oracle_edge_marginals(o)
  marginal <- function(parent, child) {
    hit <- marginals$parent == parent & marginals$child == child
    if (any(hit)) marginals$marginal[hit] else 0
  }
  is_true <- function(parent, child) any(true_edges$parent == parent & true_edges$child == child)
  normalise <- function(s) round(min(1, max(0, s)), 10)
  pairs <- list(directed = list(), undirected = list())
  for (i in ids) {
    for (j in ids) {
      if (i == j) next
      pairs$directed[[length(pairs$directed) + 1]] <- c(score = normalise(marginal(i, j)), true = is_true(i, j))
      if (i < j) {
        pairs$undirected[[length(pairs$undirected) + 1]] <- c(
          score = normalise(marginal(i, j) + marginal(j, i)), true = is_true(i, j) || is_true(j, i)
        )
      }
    }
  }
  rows <- lapply(names(pairs), function(direction) {
    score <- vapply(pairs[[direction]], function(p) p[["score"]], numeric(1))
    true <- vapply(pairs[[direction]], function(p) p[["true"]] == 1, logical(1))
    thresholds <- sort(unique(score), decreasing = TRUE)
    n_true <- sum(true)
    ap <- 0
    auc <- 0
    last_recall <- 0
    last_precision <- 1
    for (t in thresholds) {
      tp <- 0
      reached <- 0
      for (k in seq_along(score)) {
        if (score[[k]] >= t) {
          reached <- reached + 1
          if (true[[k]]) tp <- tp + 1
        }
      }
      recall <- tp / n_true
      precision <- tp / reached
      ap <- ap + (recall - last_recall) * precision
      auc <- auc + (recall - last_recall) * (precision + last_precision) / 2
      last_recall <- recall
      last_precision <- precision
    }
    kth <- sort(score, decreasing = TRUE)[[min(n_true, length(score))]]
    top <- sum(true & score >= kth) / sum(score >= kth)
    if (n_true == 0) ap <- auc <- top <- NA_real_
    data.frame(direction = direction, n_edge_universe = length(score), n_true_edges = n_true,
               n_distinct_scores = length(thresholds), average_precision = ap, pr_auc = auc, top_k_precision = top)
  })
  do.call(rbind, rows)
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

# ψ_MAP: observed members of every Infection's MAP set, over the number of Infections.
oracle_map_out_degree <- function(o, by_run = FALSE, runs = NULL) {
  sets <- oracle_parent_sets_pooled(o, by_run = by_run, runs = runs)
  if (!by_run) sets$run <- ""
  n <- length(oracle_infection_ids(o))
  out <- list()
  for (run in unique(sets$run)) {
    count <- 0
    for (i in which(sets$run == run & sets$rank == 1)) count <- count + length(sets$members[[i]])
    out[[run]] <- data.frame(run = run, map_out_degree = count / n, infections = n)
  }
  out <- oracle_bind(out, NULL)
  if (!by_run) out$run <- NULL
  out
}

# ψ_lo and ψ_hi: per Infection the fewest and most observed members over its Credible collection, summed over
# Infections (0 for one with no set), over the number of Infections.
oracle_out_degree_range <- function(o, level = 0.8, by_run = FALSE, runs = NULL) {
  sets <- oracle_parent_sets_pooled(o, level = level, by_run = by_run, runs = runs)
  if (!by_run) sets$run <- ""
  n <- length(oracle_infection_ids(o))
  out <- list()
  for (run in unique(sets$run)) {
    lo <- 0
    hi <- 0
    for (child in oracle_infection_ids(o)) {
      sizes <- numeric()
      for (i in which(sets$run == run & sets$infection == child & sets$credible)) {
        sizes <- c(sizes, length(sets$members[[i]]))
      }
      if (length(sizes)) {
        lo <- lo + min(sizes)
        hi <- hi + max(sizes)
      }
    }
    out[[run]] <- data.frame(run = run, lo = lo / n, hi = hi / n, level = level)
  }
  out <- oracle_bind(out, NULL)
  if (!by_run) out$run <- NULL
  out
}

# Co-inclusion: per child, the pooled weight of its sets holding both parents of each pair (a < b in string order).
# With by_run, a pair seen in any Run has a row in every Run (0 where it is absent).
oracle_co_inclusion <- function(o, infections = NULL, by_run = FALSE, runs = NULL) {
  sets <- oracle_parent_sets_pooled(o, by_run = by_run, runs = runs)
  if (!by_run) sets$run <- "pooled"  # not "": a list element named "" is never matched
  total <- list()
  for (i in seq_len(nrow(sets))) {
    m <- sets$members[[i]]
    if (length(m) < 2 || (!is.null(infections) && !sets$infection[[i]] %in% infections)) next
    for (x in seq_along(m)) for (y in seq_along(m)) {
      if (x >= y) next
      k <- paste(sets$infection[[i]], m[[x]], m[[y]], sep = "\r")
      if (is.null(total[[k]])) total[[k]] <- list()
      total[[k]][[sets$run[[i]]]] <- (if (is.null(total[[k]][[sets$run[[i]]]])) 0 else total[[k]][[sets$run[[i]]]]) +
        sets$weight[[i]]
    }
  }
  out <- list()
  for (k in names(total)) {
    parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
    for (run in unique(sets$run)) {
      p <- if (is.null(total[[k]][[run]])) 0 else total[[k]][[run]]
      out[[length(out) + 1]] <- data.frame(run = run, child = parts[[1]], parent_a = parts[[2]], parent_b = parts[[3]],
                                           probability = p)
    }
  }
  out <- oracle_bind(out, data.frame(run = character(), child = character(), parent_a = character(),
                                     parent_b = character(), probability = numeric()))
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

# ---- Latent genotypes and COI (#121) ----
# Bitstrings are split character by character in R: character k (1-based here) is Allele k - 1.

oracle_bits <- function(alleles) strsplit(alleles, "", fixed = TRUE)[[1]] == "1"

oracle_alleles <- function(o, runs = NULL, phase = "sampling") {
  out <- list()
  for (run in oracle_runs(o, runs)) {
    rows <- oracle_table(o, run, "latent_genotypes", phase)
    for (i in seq_len(NROW(rows))) {
      for (k in which(oracle_bits(rows$alleles[[i]]))) {
        out[[length(out) + 1]] <- data.frame(run = run, iteration = rows$iteration[[i]],
                                             infection = rows$infection[[i]], locus = rows$locus[[i]],
                                             allele = k - 1)
      }
    }
  }
  oracle_bind(out, NULL)
}

oracle_coi_draws <- function(o, runs = NULL, phase = "sampling") {
  out <- list()
  for (run in oracle_runs(o, runs)) {
    rows <- oracle_table(o, run, "latent_genotypes", phase)
    coi <- list()
    for (i in seq_len(NROW(rows))) {
      k <- paste(rows$iteration[[i]], rows$infection[[i]], sep = "\r")
      coi[[k]] <- max(if (is.null(coi[[k]])) 0 else coi[[k]], sum(oracle_bits(rows$alleles[[i]])))
    }
    for (k in names(coi)) {
      parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
      out[[length(out) + 1]] <- data.frame(run = run, iteration = as.numeric(parts[[1]]), infection = parts[[2]],
                                           coi = coi[[k]])
    }
  }
  oracle_bind(out, NULL)
}

# The sampling Runs a pooled summary uses (those with Draws), and the groups it pools: every Run together, or each.
oracle_pool_groups <- function(o, runs, by_run) {
  runs <- oracle_runs(o, runs)
  runs <- runs[vapply(runs, function(r) length(oracle_draws(o, r)) > 0, logical(1))]
  if (by_run) setNames(as.list(runs), runs) else list(pooled = runs)
}

# COI median (type 7 over every Draw's COI) and distribution (long: infection, coi, fraction of the Draws).
oracle_coi <- function(o, by_run = FALSE, runs = NULL) {
  medians <- list()
  dist <- list()
  for (group in oracle_pool_groups(o, runs, by_run)) {
    draws <- sum(vapply(group, function(r) length(oracle_draws(o, r)), numeric(1)))
    per_draw <- oracle_coi_draws(o, group)
    for (inf in sort(unique(per_draw$infection), method = "radix")) {
      values <- per_draw$coi[per_draw$infection == inf]
      lead <- if (by_run) data.frame(run = group) else data.frame(row.names = 1)
      medians[[length(medians) + 1]] <- cbind(lead, infection = inf,
                                              median = unname(stats::quantile(values, 0.5, type = 7)))
      for (v in sort(unique(values))) {
        dist[[length(dist) + 1]] <- cbind(lead, infection = inf, coi = v, fraction = sum(values == v) / draws)
      }
    }
  }
  list(median = oracle_bind(medians, NULL), distribution = oracle_bind(dist, NULL))
}

# Per (infection, locus) seen at a Draw, the count of Draws of each bitstring, over the Runs of one group.
oracle_bitstring_counts <- function(o, group) {
  counts <- list()
  for (run in group) {
    rows <- oracle_table(o, run, "latent_genotypes", "sampling")
    for (i in seq_len(NROW(rows))) {
      k <- paste(rows$infection[[i]], rows$locus[[i]], sep = "\r")
      tally <- counts[[k]]
      if (is.null(tally)) tally <- numeric()
      b <- rows$alleles[[i]]
      tally[b] <- (if (is.na(tally[b])) 0 else tally[b]) + 1
      counts[[k]] <- tally
    }
  }
  counts
}

oracle_genotype_presence <- function(o, by_run = FALSE, runs = NULL) {
  loci <- oracle_parquet(file.path(o$runs[[1]]$path, "structure", "loci.parquet"))
  out <- list()
  for (group in oracle_pool_groups(o, runs, by_run)) {
    draws <- sum(vapply(group, function(r) length(oracle_draws(o, r)), numeric(1)))
    counts <- oracle_bitstring_counts(o, group)
    for (k in names(counts)) {
      parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
      n_alleles <- sum(loci$locus == parts[[2]])
      present <- numeric(n_alleles)
      tally <- counts[[k]]
      for (b in names(tally)) present <- present + oracle_bits(b) * tally[[b]]
      for (allele in seq_len(n_alleles)) {
        out[[length(out) + 1]] <- cbind(if (by_run) data.frame(run = group) else data.frame(row.names = 1),
                                        infection = parts[[1]], locus = parts[[2]], allele = allele - 1,
                                        presence = present[[allele]] / draws)
      }
    }
  }
  oracle_bind(out, NULL)
}

oracle_genotype_mode <- function(o, by_run = FALSE, runs = NULL) {
  out <- list()
  for (group in oracle_pool_groups(o, runs, by_run)) {
    draws <- sum(vapply(group, function(r) length(oracle_draws(o, r)), numeric(1)))
    counts <- oracle_bitstring_counts(o, group)
    for (k in names(counts)) {
      parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
      tally <- counts[[k]]
      best <- names(tally)[tally == max(tally)]
      mode <- sort(best, method = "radix")[[1]]
      out[[length(out) + 1]] <- cbind(if (by_run) data.frame(run = group) else data.frame(row.names = 1),
                                      infection = parts[[1]], locus = parts[[2]], alleles = mode,
                                      fraction = tally[[mode]] / draws, agreement = tally[[mode]] / sum(tally))
    }
  }
  oracle_bind(out, NULL)
}

# Functionals (#119): at every Draw, the mean over every Infection of a per-Infection value from an oracle per-Draw
# table; an Infection with no row at the Draw takes `fill`. `columns` maps output names to the table's columns.
oracle_infection_mean <- function(o, rows, columns, fill, runs = NULL, phase = "sampling") {
  out <- list()
  ids <- oracle_infection_ids(o)
  for (run in oracle_runs(o, runs)) {
    for (it in oracle_draws(o, run, phase)) {
      row <- data.frame(run = run, iteration = it)
      for (name in names(columns)) {
        total <- 0
        for (inf in ids) {
          hit <- which(rows$run == run & rows$iteration == it & rows$infection == inf)
          stopifnot(length(hit) <= 1)
          total <- total + if (length(hit)) rows[[columns[[name]]]][[hit]] else fill[[name]]
        }
        row[[name]] <- total / length(ids)
      }
      out[[length(out) + 1]] <- row
    }
  }
  oracle_bind(out, NULL)
}

oracle_mean_out_degree <- function(o, runs = NULL, phase = "sampling") {
  oracle_infection_mean(o, oracle_out_degree_draws(o, runs, phase), c(mean_out_degree = "out_degree"),
                        c(mean_out_degree = 0), runs, phase)
}

oracle_mean_p_latent_parent <- function(o, runs = NULL, phase = "sampling") {
  oracle_infection_mean(o, oracle_infection_draws(o, runs, phase), c(mean_p_latent_parent = "p_latent_parent"),
                        c(mean_p_latent_parent = 0), runs, phase)
}

oracle_multi_source <- function(o, runs = NULL, phase = "sampling") {
  oracle_infection_mean(o, oracle_infection_draws(o, runs, phase),
                        c(augmented = "multi_source_augmented", observed_only = "multi_source_observed"),
                        c(augmented = 0, observed_only = 0), runs, phase)
}

oracle_entropy <- function(o, runs = NULL, phase = "sampling") {
  oracle_infection_mean(o, oracle_infection_draws(o, runs, phase),
                        c(mean_entropy = "entropy", mean_effective_size = "effective_size"),
                        c(mean_entropy = 0, mean_effective_size = 1), runs, phase)
}

# Long: one row per Draw per class of oracle_class_draws() (already a full zero-filled grid).
oracle_composition <- function(o, runs = NULL, phase = "sampling") {
  shares <- oracle_class_draws(o, runs, phase)
  n <- length(oracle_infection_ids(o))
  out <- list()
  for (k in unique(paste(shares$run, shares$iteration, shares$class, sep = "\r"))) {
    parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
    mine <- shares$run == parts[[1]] & shares$iteration == as.numeric(parts[[2]]) & shares$class == parts[[3]]
    total <- 0
    for (s in shares$share[mine]) total <- total + s
    out[[length(out) + 1]] <- data.frame(run = parts[[1]], iteration = as.numeric(parts[[2]]), class = parts[[3]],
                                         share = total / n)
  }
  oracle_bind(out, data.frame(run = character(), iteration = numeric(), class = character(), share = numeric()))
}

# Type-7 quantile written out: h = 1 + (n - 1) p; x[floor h] moved (h - floor h) of the way to x[ceiling h].
oracle_quantile7 <- function(x, p) {
  x <- sort(x)
  h <- 1 + (length(x) - 1) * p
  lo <- floor(h)
  hi <- ceiling(h)
  if (hi == lo || x[[hi]] == x[[lo]]) return(x[[lo]])
  (1 - (h - lo)) * x[[lo]] + (h - lo) * x[[hi]]
}

# pt_interval() by hand: `x` a collected per-Draw data frame, `by` and `value` character vectors (value: one
# column, or several with a `value` column naming each). Rows grouped by the `by` keys; pooled mean over the
# group's rows, equal-tailed type-7 interval.
oracle_interval <- function(x, by = character(), value, level = 0.95) {
  x <- as.data.frame(x)
  groups <- if (length(by)) split(seq_len(nrow(x)), lapply(x[by], as.character), drop = TRUE)
            else list(seq_len(nrow(x)))
  out <- list()
  for (idx in groups) {
    for (v in value) {
      y <- x[[v]][idx]
      total <- 0
      for (yi in y) total <- total + yi
      row <- x[idx[1], by, drop = FALSE]
      if (length(value) > 1) row$value <- v
      row$estimate <- total / length(y)
      row$lower <- oracle_quantile7(y, (1 - level) / 2)
      row$upper <- oracle_quantile7(y, (1 + level) / 2)
      row$level <- level
      row$draws <- length(y)
      out[[length(out) + 1]] <- row
    }
  }
  oracle_bind(out, NULL)
}

# ---- Infection summary (#123) ----
# One row per Infection (per Run with by_run), each column from the oracle of the function it summarises, summed
# over the group's Draws in loops. `observation_time`: named vector (NA by default).
oracle_infection_summary <- function(o, by_run = FALSE, runs = NULL, observation_time = NULL) {
  groups <- oracle_pool_groups(o, runs, by_run)
  all_runs <- unlist(groups, use.names = FALSE)
  ids <- sort(oracle_infection_ids(o), method = "radix")
  shares <- oracle_class_draws(o, all_runs)
  classes <- unique(shares[c("class", "n_observed", "has_latent_parent")])
  classes <- classes$class[order(classes$n_observed, classes$has_latent_parent)]
  values <- oracle_infection_draws(o, all_runs)
  degrees <- oracle_out_degree_draws(o, all_runs)
  sets <- oracle_parent_sets_pooled(o, by_run = by_run, runs = all_runs)
  if (!by_run) sets$run <- "pooled"
  edges <- oracle_edge_marginals(o, by_run = by_run, runs = all_runs)
  if (!by_run) edges$run <- "pooled"
  coi <- oracle_coi(o, by_run = by_run, runs = all_runs)$median
  if (!by_run) coi$run <- "pooled"
  out <- list()
  for (g in names(groups)) {
    group <- groups[[g]]
    draws <- 0
    for (r in group) draws <- draws + length(oracle_draws(o, r))
    for (inf in ids) {
      row <- if (by_run) data.frame(run = g, infection = inf) else data.frame(infection = inf)
      row$observation_time <- if (is.null(observation_time)) NA_real_ else unname(observation_time[inf])
      q <- 0; e <- 0; h <- 0; es <- 0
      for (r in group) for (it in oracle_draws(o, r)) {
        hit <- which(values$run == r & values$iteration == it & values$infection == inf)
        if (length(hit)) {
          q <- q + values$p_latent_parent[[hit]]
          e <- e + values$expected_observed_parents[[hit]]
          h <- h + values$entropy[[hit]]
          es <- es + values$effective_size[[hit]]
        } else {
          es <- es + 1
        }
      }
      row$p_latent_parent <- q / draws
      row$expected_observed_parents <- e / draws
      for (class in classes) {
        s <- 0
        for (i in which(shares$run %in% group & shares$infection == inf & shares$class == class)) {
          s <- s + shares$share[[i]]
        }
        row[[paste0("share_", class)]] <- s / draws
      }
      row$mean_entropy <- h / draws
      row$mean_effective_size <- es / draws
      ph <- 0
      mine <- which(sets$run == g & sets$infection == inf)
      for (i in mine) if (sets$weight[[i]] > 0) ph <- ph - sets$weight[[i]] * log(sets$weight[[i]])
      row$posterior_entropy <- ph
      row$posterior_effective_size <- exp(ph)
      map <- mine[sets$rank[mine] == 1]
      row$map_members <- I(list(if (length(map)) sets$members[[map]] else character()))
      row$map_class <- if (length(map)) sets$class[[map]] else NA_character_
      row$map_weight <- if (length(map)) sets$weight[[map]] else NA_real_
      top <- NA_character_
      best <- 0
      for (i in which(edges$run == g & edges$child == inf)) {
        m <- edges$marginal[[i]]
        p <- edges$parent[[i]]
        if (m > best || (m == best && m > 0 && sort(c(p, top), method = "radix")[[1]] == p)) {
          top <- p
          best <- m
        }
      }
      row$top_parent <- top
      row$top_parent_marginal <- if (is.na(top)) NA_real_ else best
      d <- 0
      for (i in which(degrees$run %in% group & degrees$infection == inf)) d <- d + degrees$out_degree[[i]]
      row$expected_out_degree <- d / draws
      hit <- which(coi$run == g & coi$infection == inf)
      row$coi_median <- if (length(hit)) coi$median[[hit]] else NA_real_
      out[[length(out) + 1]] <- row
    }
  }
  out <- oracle_bind(out, NULL)
  out$map_members <- unclass(out$map_members)
  out
}

# ---- Infection times and the Ordering (#122) ----

# Every pair i < j (string order) at every Draw, from each Draw's ranks: before = rank_i < rank_j.
oracle_order_pair_draws <- function(o, runs = NULL, phase = "sampling") {
  out <- list()
  for (run in oracle_runs(o, runs)) {
    rows <- oracle_table(o, run, "ordering", phase)
    for (it in sort(unique(rows$iteration))) {
      at <- rows[rows$iteration == it, ]
      rank <- setNames(at$rank, at$infection)
      ids <- sort(at$infection, method = "radix")
      for (x in seq_along(ids)) for (y in seq_along(ids)) {
        if (x >= y) next
        out[[length(out) + 1]] <- data.frame(run = run, iteration = it, i = ids[[x]], j = ids[[y]],
                                             before = rank[[ids[[x]]]] < rank[[ids[[y]]]])
      }
    }
  }
  oracle_bind(out, NULL)
}

# Share of the group's Draws with i before j; `infections` keeps pairs touching any of them.
oracle_order_pairs <- function(o, infections = NULL, by_run = FALSE, runs = NULL) {
  out <- list()
  for (group in oracle_pool_groups(o, runs, by_run)) {
    draws <- sum(vapply(group, function(r) length(oracle_draws(o, r)), numeric(1)))
    pairs <- oracle_order_pair_draws(o, group)
    if (!is.null(infections)) pairs <- pairs[pairs$i %in% infections | pairs$j %in% infections, ]
    keys <- unique(paste(pairs$i, pairs$j, sep = "\r"))
    for (k in keys) {
      ij <- strsplit(k, "\r", fixed = TRUE)[[1]]
      n <- 0
      for (r in which(pairs$i == ij[[1]] & pairs$j == ij[[2]])) n <- n + pairs$before[[r]]
      lead <- if (by_run) data.frame(run = group) else data.frame(row.names = 1)
      out[[length(out) + 1]] <- cbind(lead, i = ij[[1]], j = ij[[2]], p_before = n / draws, draws = draws)
    }
  }
  oracle_bind(out, NULL)
}

# Observation time (named vector id -> time, NA when missing) minus oracle_quantile7() of every Draw's duration.
oracle_infection_times <- function(o, observation_time, level = 0.95, by_run = FALSE, runs = NULL) {
  out <- list()
  for (group in oracle_pool_groups(o, runs, by_run)) {
    rows <- oracle_bind(lapply(group, function(r) oracle_table(o, r, "infections")), NULL)
    for (inf in sort(unique(rows$infection), method = "radix")) {
      d <- rows$infection_duration[rows$infection == inf]
      time <- if (inf %in% names(observation_time)) unname(observation_time[[inf]]) else NA_real_
      lead <- if (by_run) data.frame(run = group) else data.frame(row.names = 1)
      out[[length(out) + 1]] <- cbind(lead, infection = inf,
                                      lower = time - oracle_quantile7(d, (1 + level) / 2),
                                      median = time - oracle_quantile7(d, 0.5),
                                      upper = time - oracle_quantile7(d, (1 - level) / 2))
    }
  }
  oracle_bind(out, NULL)
}

# ---- Fit summary (#131) ----
# pt_fit_summary() for a posterior fit, every column from the oracles above. Truth is written by hand:
# `true_edges` (parent, child among the Infections), `true_sets` (a named list per Infection of list(members,
# latent)), `n_true_cases`. Factor columns are left out: compare from n_true_cases on.
oracle_fit_summary <- function(o, true_edges, true_sets, n_true_cases) {
  ids <- oracle_infection_ids(o)
  row <- list(n_true_cases = n_true_cases, n_observed_cases = length(ids),
              true_average_outdegree = nrow(true_edges) / length(ids))
  latent <- 0
  for (id in ids) if (true_sets[[id]]$latent) latent <- latent + 1
  row$true_unobserved_source_proportion <- latent / length(ids)

  interval <- function(x, value) oracle_interval(x, value = value)
  put <- function(names, x) {
    row[[names[[1]]]] <<- x$estimate
    row[[names[[2]]]] <<- x$lower
    row[[names[[3]]]] <<- x$upper
  }
  put(c("estimated_average_outdegree", "average_outdegree_lower", "average_outdegree_upper"),
      interval(oracle_mean_out_degree(o), "mean_out_degree"))
  put(c("estimated_unobserved_source_proportion", "unobserved_source_lower", "unobserved_source_upper"),
      interval(oracle_mean_p_latent_parent(o), "mean_p_latent_parent"))
  shares <- oracle_composition(o)
  draws <- unique(shares[c("run", "iteration")])
  for (kind in c("observed_only", "mixed", "unobserved_only")) {
    draws[[kind]] <- 0
  }
  for (k in seq_len(nrow(shares))) {
    class <- shares$class[[k]]
    kind <- if (class == "latent_only") "unobserved_only"
            else if (grepl("^observed_[0-9]+_latent$", class)) "mixed"
            else if (grepl("^observed_[0-9]+$", class)) "observed_only"
            else stop("unexpected class ", class)
    at <- which(draws$run == shares$run[[k]] & draws$iteration == shares$iteration[[k]])
    draws[[kind]][at] <- draws[[kind]][at] + shares$share[[k]]
  }
  stem <- function(s) paste0(s, c("", "_lower", "_upper"))
  put(stem("mean_observed_only_probability"), interval(draws, "observed_only"))
  put(stem("mean_mixed_observed_unobserved_probability"), interval(draws, "mixed"))
  put(stem("mean_unobserved_only_probability"), interval(draws, "unobserved_only"))
  put(stem("mean_no_parent_probability"), list(estimate = 0, lower = 0, upper = 0))
  multi <- oracle_multi_source(o)
  put(stem("mean_multi_source_probability_augmented"), interval(multi, "augmented"))
  put(stem("mean_multi_source_probability_observed_only"), interval(multi, "observed_only"))
  entropy <- oracle_entropy(o)
  put(stem("mean_parent_set_entropy"), interval(entropy, "mean_entropy"))
  put(stem("mean_effective_parent_sets"), interval(entropy, "mean_effective_size"))

  edges <- oracle_score_edges(o, true_edges)
  for (direction in c("directed", "undirected")) {
    for (column in setdiff(names(edges), "direction")) {
      row[[paste0(direction, "_", column)]] <- edges[[column]][edges$direction == direction]
    }
  }

  sets <- oracle_parent_sets_pooled(o, level = 0.8)
  size <- covered <- top <- numeric()
  for (id in ids) {
    mine <- sets[sets$infection == id, ]
    credible <- mine[mine$credible, ]
    size[[id]] <- nrow(credible)
    hit <- FALSE
    for (r in seq_len(nrow(credible))) {
      if (identical(sort(credible$members[[r]]), sort(true_sets[[id]]$members)) &&
          credible$has_latent_parent[[r]] == true_sets[[id]]$latent) hit <- TRUE
    }
    covered[[id]] <- hit
    top[[id]] <- if (nrow(mine)) mine$weight[mine$rank == 1] else 0
  }
  sorted <- sort(size)
  n <- length(sorted)
  row$mean_credible_parent_set_size_80 <- sum(size) / n
  row$median_credible_parent_set_size_80 <- if (n %% 2) sorted[[(n + 1) / 2]] else (sorted[[n / 2]] + sorted[[n / 2 + 1]]) / 2
  row$credible_parent_set_coverage_80 <- sum(covered) / n
  row$mean_top_parent_set_probability <- sum(top) / n

  for (pair in list(c("average_outdegree", "estimated_average_outdegree", "true_average_outdegree"),
                    c("unobserved_source", "estimated_unobserved_source_proportion",
                      "true_unobserved_source_proportion"))) {
    error <- row[[pair[[2]]]] - row[[pair[[3]]]]
    row[[paste0(pair[[1]], "_error")]] <- error
    row[[paste0(pair[[1]], "_abs_error")]] <- abs(error)
    row[[paste0(pair[[1]], "_sq_error")]] <- error * error
  }
  as.data.frame(row)
}
