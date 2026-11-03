# Oracle for out-degree posteriors (#128), apart from helper-oracle.R so parallel tickets don't collide.
#
# At every Draw it enumerates every combination of one Parent set per child, weighted by the product of the
# weights, and counts each quantity in each combination. A child's weights are used as stored; the weight they
# leave (what the floor dropped) is one more option with no member. A child with no row at a Draw has that option
# alone. The posterior is the mean of the per-Draw distributions over every Draw of every Run with Draws.
#
#   oracle_count_draws(o, count)   per Draw, a list of pmfs (one per quantity `count` returns: a named vector of
#                                  counts given the observed members of every child's chosen set).
#   oracle_count_summary(pmfs, level)  mixture of a list of per-Draw pmfs: pmf, estimate, lower, upper, var_within,
#                                  var_between.
#   oracle_out_degree_dist(o, infections, level, pmf), oracle_realized_out_degree(o, strata, level)
#   expect_pmf_equal(actual, expected, keys)   pmfs equal within 1e-12 absolute, counts missing on one side = 0.

oracle_count_draws <- function(o, count) {
  ids <- oracle_infection_ids(o)
  out <- list()
  for (run in names(o$runs)) {
    sets <- oracle_parent_sets(o, run)
    for (it in oracle_draws(o, run)) {
      options <- lapply(ids, function(child) {
        rows <- if (is.null(sets)) integer() else which(sets$iteration == it & sets$infection == child)
        members <- lapply(rows, function(i) sets$observed[[i]])
        weights <- vapply(rows, function(i) sets$prob[[i]], 0)
        left <- 1 - sum(weights)
        if (left > 0) {
          members <- c(members, list(character()))
          weights <- c(weights, left)
        }
        list(members = members, weights = weights)
      })
      combos <- expand.grid(lapply(options, function(x) seq_along(x$weights)))
      pmfs <- list()
      for (r in seq_len(nrow(combos))) {
        weight <- 1
        chosen <- list()
        for (j in seq_along(ids)) {
          k <- combos[r, j]
          weight <- weight * options[[j]]$weights[[k]]
          chosen[[j]] <- options[[j]]$members[[k]]
        }
        counts <- count(chosen)
        for (q in names(counts)) {
          p <- pmfs[[q]]
          if (is.null(p)) p <- numeric()
          n <- counts[[q]] + 1
          if (length(p) < n) p <- c(p, numeric(n - length(p)))
          p[[n]] <- p[[n]] + weight
          pmfs[[q]] <- p
        }
      }
      out[[length(out) + 1]] <- pmfs
    }
  }
  out
}

oracle_count_summary <- function(pmfs, level) {
  width <- max(vapply(pmfs, length, 0))
  pmfs <- lapply(pmfs, function(p) c(p, numeric(width - length(p))))
  k <- seq_len(width) - 1
  means <- vapply(pmfs, function(p) sum(k * p), 0)
  variances <- vapply(seq_along(pmfs), function(d) sum((k - means[[d]])^2 * pmfs[[d]]), 0)
  mixed <- numeric(width)
  for (p in pmfs) mixed <- mixed + p
  mixed <- mixed / length(pmfs)
  estimate <- mean(means)
  quantile_at <- function(p) {
    cdf <- 0
    for (j in seq_len(width)) {
      cdf <- cdf + mixed[[j]]
      if (cdf >= p - 1e-12) return(j - 1)
    }
    width - 1
  }
  list(pmf = mixed, estimate = estimate, lower = quantile_at((1 - level) / 2), upper = quantile_at((1 + level) / 2),
       var_within = mean(variances), var_between = mean((means - estimate)^2))
}

oracle_out_degree_dist <- function(o, infections = NULL, level = 0.95, pmf = FALSE) {
  ids <- sort(if (is.null(infections)) oracle_infection_ids(o) else infections, method = "radix")
  draws <- oracle_count_draws(o, function(chosen) {
    counts <- setNames(numeric(length(ids)), ids)
    for (members in chosen) for (m in members) if (m %in% ids) counts[[m]] <- counts[[m]] + 1
    counts
  })
  rows <- lapply(ids, function(i) {
    s <- oracle_count_summary(lapply(draws, function(d) d[[i]]), level)
    if (pmf) return(data.frame(infection = i, k = seq_along(s$pmf) - 1L, probability = s$pmf))
    data.frame(infection = i, estimate = s$estimate, lower = s$lower, upper = s$upper, level = level,
               var_within = s$var_within, var_between = s$var_between, p_zero = s$pmf[[1]], draws = length(draws))
  })
  do.call(rbind, rows)
}

# `strata`: the stratum of every Infection (oracle_strata()), or NULL for every Infection in one stratum.
oracle_realized_out_degree <- function(o, strata = NULL, level = 0.95) {
  ids <- oracle_infection_ids(o)
  key <- if (is.null(strata)) rep("all", length(ids)) else ifelse(is.na(strata), "\001NA", as.character(strata))
  levels <- unique(key)
  draws <- oracle_count_draws(o, function(chosen) {
    counts <- setNames(numeric(length(levels)), levels)
    for (members in chosen) for (m in members) counts[[key[[match(m, ids)]]]] <- counts[[key[[match(m, ids)]]]] + 1
    counts
  })
  rows <- lapply(levels, function(s) {
    n <- sum(key == s)
    x <- oracle_count_summary(lapply(draws, function(d) d[[s]]), level)
    row <- data.frame(n = n, estimate = x$estimate / n, lower = x$lower / n, upper = x$upper / n, level = level,
                      var_within = x$var_within / n^2, var_between = x$var_between / n^2, draws = length(draws))
    if (is.null(strata)) row else cbind(data.frame(stratum = strata[match(s, key)]), row)
  })
  do.call(rbind, rows)
}

expect_pmf_equal <- function(actual, expected, tolerance = 1e-12) {
  actual <- as.data.frame(actual)
  expect_equal(names(actual), names(expected))
  key <- function(x) paste(x$infection, x$k, sep = "\r")
  all <- union(key(actual), key(expected))
  a <- setNames(numeric(length(all)), all)
  e <- a
  a[key(actual)] <- actual$probability
  e[key(expected)] <- expected$probability
  expect_true(all(abs(a - e) <= tolerance), label = sprintf("pmf within %g of the oracle (worst %s)", tolerance,
                                                            names(which.max(abs(a - e)))))
  expect_setequal(unique(actual$infection), unique(expected$infection))
}
