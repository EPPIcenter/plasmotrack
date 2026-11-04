# Out-degree posteriors: the out-degree's full distribution, mixed over Draws, and sampled networks.
#
# Given a Draw, children's Parent sets are independent, so a count that sums over children (|Pa(j) ∩ S| for a set S
# of sources) has an exact distribution given the Draw: the convolution of the children's. The posterior is the
# mixture over every Draw of every Run, each Draw weighing equally.
#
#   draw_index(a, runs)           eager: run, iteration, draw (1..D) over the sampling Draws, Runs in `runs` order.
#   group_count_rows(a, runs, groups)  eager: run, iteration, infection, grp, c, weight: per child and Draw, the
#                                 weight of its sets with c >= 1 observed members in group `grp` (one DuckDB query,
#                                 whatever the number of groups).
#   child_count_pmfs(rows, draws, groups)  per group, a list of D x (max c + 1) matrices, one per child: the child's
#                                 count distribution at every Draw, count 0 holding the weight the others leave.
#   count_posterior(pmfs, draws, level)  the mixed pmf, mean, var_within, var_between, lower, upper.

out_degree_tail_sds <- 12
out_degree_tolerance <- 1e-12

#' Out-degree posteriors
#'
#' The posterior distribution of out-degree, the number of observed children an Infection infected. An interval from
#' [pt_interval()] over [pt_mean_out_degree()] or [pt_out_degree_draws()] covers only the variation of the expected
#' out-degree between Draws. Given a Draw the out-degree still varies over the Parent sets the weights allow, and that
#' variation can be the larger part. These functions give the full distribution.
#'
#' Given a Draw, children's Parent sets are independent, so the distribution given the Draw is exact:
#'
#' * `pt_out_degree_dist()`: an Infection's out-degree is Poisson-binomial over its children, each child holding it
#'   as a parent with probability the edge's probability at the Draw (its per-Draw Edge marginal, see
#'   [pt_edge_draws()]). An Infection that is no Infection's parent has out-degree 0.
#' * `pt_realized_out_degree()`: the realized mean out-degree of a stratum's Infections, counted as sources as
#'   [pt_mean_out_degree()] counts them. Given a Draw, the stratum's total out-degree is the sum over children of the
#'   number of observed parents in the stratum, independent across children. Values are on the mean scale (the total
#'   divided by `n`, the stratum's size). Without `by`, the stratum is every Infection. Each Run is read once, whatever
#'   the number of strata.
#'
#' **Mixture.** The posterior is the mean of the distributions given each Draw, over every sampling Draw of every
#' Run, so each Run weighs by its Draws, as all pooling does. Runs with no sampling Draws are skipped with a message.
#'
#' **Weights** are used as stored. The Store drops Parent sets below `parent_set_floor` without renormalising, and a
#' child's count is 0 with whatever weight its other counts leave. So `estimate` equals the pooled mean of
#' [pt_mean_out_degree()] (or of [pt_out_degree_draws()] for one Infection) exactly, and each distribution moves by at
#' most the dropped weight.
#'
#' **Intervals** are equal-tailed quantiles of the discrete mixture: the smallest count whose cumulative probability
#' is at least \eqn{p - 10^{-12}}, for \eqn{p = (1 - level) / 2} and \eqn{(1 + level) / 2}. These are not type-7
#' sample quantiles, and an interval need not contain the estimate: a stratum with a mean of 0.003 children can have
#' the interval \[0, 0\].
#'
#' **Variance.** `var_within` is the mean over Draws of the variance given the Draw; `var_between` is the variance
#' over Draws (divisor: the number of Draws) of the mean given the Draw, which is all that [pt_interval()] sees.
#' Their sum is the posterior variance.
#'
#' **Support.** The convolution is capped at 12 conditional standard deviations plus 12 above the largest per-Draw
#' mean, with the mass beyond kept in the last count, so a pmf sums to 1. Means and variances come from the children's
#' distributions, not the capped support.
#'
#' The Latent parent is never counted: out-degree is to observed children.
#'
#' @param a An Analysis from [pt_read()].
#' @param infections `NULL` for every Infection, or a character vector of Infection ids.
#' @param level Interval probability, in (0, 1].
#' @param pmf `FALSE` for one row per Infection; `TRUE` for the mixed distribution, one row per Infection per count.
#' @param by `NULL`, or the name of one epi column (a string; see [pt_read()]) whose values are the strata of sources.
#'   Infections with an NA value, or missing from epi, form the NA stratum, which is kept.
#' @return A tibble.
#'   * `pt_out_degree_dist()`: `infection`, `estimate` (the posterior mean), `lower`, `upper`, `level`,
#'     `var_within`, `var_between`, `p_zero` (the probability of no observed child) and `draws`, one row per
#'     Infection in string order. With `pmf = TRUE`: `infection`, `k` and `probability`, `k` from 0 to the capped
#'     support's last count.
#'   * `pt_realized_out_degree()`: `stratum` (with `by` only, sorted, NA last), `n`, `estimate`, `lower`, `upper`,
#'     `level`, `var_within`, `var_between` and `draws`, all on the mean scale.
#' @seealso [pt_mean_out_degree()] and [pt_interval()] for the between-Draw interval; [pt_sample_networks()] for
#'   functionals that are not sums over children.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' pt_out_degree_dist(a)
#' pt_out_degree_dist(a, infections = "1", pmf = TRUE)
#' pt_realized_out_degree(a)
#' pt_close(a)
#'
#' epi <- data.frame(infection = c("1", "2"), district = c("north", "south"))
#' a <- pt_read(file.path(golden, "single-chain"), epi = epi, threads = 1)
#' pt_realized_out_degree(a, by = "district")
#' pt_close(a)
pt_out_degree_dist <- function(a, infections = NULL, level = 0.95, pmf = FALSE) {
  infections <- check_infection_subset(a, infections)
  check_credible_level(level)
  check_flag(pmf)
  runs <- runs_with_draws(a)
  ids <- sort(if (is.null(infections)) a$infections$infection else infections, method = "radix")
  index <- draw_index(a, runs)
  rows <- group_count_rows(a, runs, tibble::tibble(member = ids, grp = seq_along(ids)), index)
  pmfs <- child_count_pmfs(rows, nrow(index), length(ids))
  posteriors <- lapply(pmfs, count_posterior, draws = nrow(index), level = level)
  if (pmf) {
    support <- lapply(posteriors, `[[`, "pmf")
    return(tibble::tibble(
      infection = rep(ids, lengths(support)),
      k = unlist(lapply(support, function(p) seq_along(p) - 1L)),
      probability = unlist(support)
    ))
  }
  field <- function(name) vapply(posteriors, `[[`, 0, name)
  tibble::tibble(
    infection = ids, estimate = field("mean"), lower = field("lower"), upper = field("upper"), level = level,
    var_within = field("var_within"), var_between = field("var_between"),
    p_zero = vapply(posteriors, function(p) p$pmf[[1]], 0), draws = nrow(index)
  )
}

#' @rdname pt_out_degree_dist
#' @export
pt_realized_out_degree <- function(a, by = NULL, level = 0.95) {
  check_by(a, by)
  check_credible_level(level)
  runs <- runs_with_draws(a)
  strata <- source_strata(a, by)
  index <- draw_index(a, runs)
  groups <- tibble::tibble(member = strata$ids, grp = strata$index)
  pmfs <- child_count_pmfs(group_count_rows(a, runs, groups, index), nrow(index), length(strata$n))
  posteriors <- lapply(pmfs, count_posterior, draws = nrow(index), level = level)
  field <- function(name) vapply(posteriors, `[[`, 0, name)
  n <- strata$n
  stratum_rows(
    strata, by, n = n, estimate = field("mean") / n, lower = field("lower") / n, upper = field("upper") / n,
    level = level, var_within = field("var_within") / n^2, var_between = field("var_between") / n^2,
    draws = nrow(index)
  )
}

# Every sampling Draw of `runs`, numbered 1..D in Run then iteration order.
draw_index <- function(a, runs) {
  grid <- dplyr::collect(draw_grid(a, runs, "sampling"))
  grid <- grid[order(match(grid$run, runs), grid$iteration), c("run", "iteration")]
  grid$draw <- seq_len(nrow(grid))
  grid
}

# Per child, Draw, group and count c >= 1: the weight of the child's sets with c observed members in the group.
# `groups` maps member ids to group numbers (an id in at most one group). One query over every Run.
group_count_rows <- function(a, runs, groups, index) {
  con <- analysis_con(a)
  counted <- pt_table(a, "parent_set_members", runs = runs) |>
    dplyr::inner_join(structure_infections(a), by = "infection") |>
    dplyr::filter(.data$member != .data$latent_parent) |>
    dplyr::inner_join(dbplyr::copy_inline(con, groups), by = "member") |>
    dplyr::group_by(.data$run, .data$infection, .data$set_id, .data$grp) |>
    dplyr::summarise(c = dplyr::n(), .groups = "drop")
  rows <- pt_table(a, "parent_sets", runs = runs, phase = "sampling") |>
    dplyr::select("run", "iteration", "infection", "set_id", "prob") |>
    dplyr::inner_join(counted, by = c("run", "infection", "set_id")) |>
    dplyr::group_by(.data$run, .data$iteration, .data$infection, .data$grp, .data$c) |>
    dplyr::summarise(weight = sum(.data$prob, na.rm = TRUE), .groups = "drop") |>
    dplyr::collect()
  rows$draw <- index$draw[match(row_key(rows$run, rows$iteration),
                                row_key(index$run, index$iteration))]
  rows <- rows[!is.na(rows$draw), ]
  rows$grp <- as.integer(rows$grp)
  rows$c <- as.integer(rows$c)
  rows$weight <- as.numeric(rows$weight)
  rows
}

# Per group, one D x (max c + 1) matrix per child with a counted set at some Draw. Count 0 takes the weight the
# others leave (sets without a member in the group, and weight the floor dropped).
child_count_pmfs <- function(rows, draws, n_groups) {
  out <- rep(list(list()), n_groups)
  # Rows are unique per (group, child, Draw, count), so each lands in its own cell.
  keys <- split(seq_len(nrow(rows)), row_key(rows$grp, rows$infection))
  for (key in keys) {
    g <- rows$grp[[key[[1]]]]
    p <- matrix(0, draws, max(rows$c[key]) + 1L)
    p[cbind(rows$draw[key], rows$c[key] + 1L)] <- rows$weight[key]
    p[, -1] <- pmin(p[, -1], 1)
    p[, 1] <- pmax(1 - rowSums(p[, -1, drop = FALSE]), 0)
    out[[g]] <- c(out[[g]], list(p))
  }
  out
}

# The posterior of a sum of independent counts, each given by a D x support matrix of per-Draw pmfs, mixed over the
# D equally weighted Draws. A sum of no counts is 0.
count_posterior <- function(pmfs, draws, level) {
  mean <- numeric(draws)
  variance <- numeric(draws)
  for (p in pmfs) {
    s <- seq_len(ncol(p)) - 1
    m <- as.vector(p %*% s)
    mean <- mean + m
    variance <- variance + as.vector(p %*% (s * s)) - m * m
  }
  largest <- sum(vapply(pmfs, ncol, 0L) - 1L)
  bound <- max(mean + out_degree_tail_sds * sqrt(pmax(variance, 0))) + out_degree_tail_sds
  cap <- as.integer(min(largest, ceiling(bound)))
  mixed <- colMeans(convolve_counts(pmfs, draws, cap))
  estimate <- mean(mean)
  bounds <- count_quantiles(mixed, level)
  list(pmf = mixed, mean = estimate, var_within = mean(variance), var_between = mean((mean - estimate)^2),
       lower = bounds[[1]], upper = bounds[[2]])
}

# Per Draw (row), the distribution of the sum: a D x (cap + 1) matrix whose last column holds the mass at cap and
# above. Counts are non-negative, so capping as it goes is exact below the cap.
convolve_counts <- function(pmfs, draws, cap) {
  total <- matrix(0, draws, cap + 1L)
  total[, 1] <- 1
  for (p in pmfs) {
    out <- total * p[, 1]
    for (s in seq_len(ncol(p) - 1L)) {
      w <- p[, s + 1L]
      if (s <= cap) {
        out[, (s + 1L):(cap + 1L)] <- out[, (s + 1L):(cap + 1L)] + total[, 1:(cap + 1L - s), drop = FALSE] * w
        out[, cap + 1L] <- out[, cap + 1L] + rowSums(total[, (cap + 2L - s):(cap + 1L), drop = FALSE]) * w
      } else {
        out[, cap + 1L] <- out[, cap + 1L] + rowSums(total) * w
      }
    }
    total <- out
  }
  total
}

# Equal-tailed quantiles of a pmf over 0, 1, ...: the smallest counts whose cumulative probability reaches
# (1 - level) / 2 and (1 + level) / 2, within the tolerance.
count_quantiles <- function(pmf, level) {
  cdf <- cumsum(pmf)
  tail <- (1 - level) / 2
  at <- function(p) {
    k <- which(cdf >= p - out_degree_tolerance)
    if (length(k)) k[[1]] - 1 else length(pmf) - 1
  }
  c(at(tail), at(1 - tail))
}

#' Sample transmission networks
#'
#' Draws whole transmission networks from the posterior, for functionals that are not sums over children (such as
#' an offspring distribution fitted per network). Each network picks one Draw uniformly from every sampling Draw of
#' every Run, then one Parent set per child from the child's weights at that Draw. The stored weights are
#' normalised to draw from them: the Store drops sets below `parent_set_floor`, so they can sum to slightly less
#' than 1. Runs with no sampling Draws are skipped with a message.
#'
#' It uses R's random number generator, so [set.seed()] reproduces it. Over many networks, the share holding an
#' edge estimates its Edge marginal ([pt_edge_marginals()]).
#'
#' @inheritParams pt_out_degree_dist
#' @param n Number of networks, a whole number.
#' @return A tibble with one row per network per child with a Parent set at the Draw: `sample` (1 to `n`), `run` and
#'   `iteration` (the Draw), `infection` (the child), `members` (a list of its observed parent ids, in string order)
#'   and `has_latent_parent`, as in [pt_parent_sets_pooled()]. Rows are sorted by `sample`, then `infection` (string
#'   order). The Latent parent's `<id>_copy` id never appears.
#' @seealso [pt_out_degree_dist()] for exact out-degree distributions.
#' @export
#' @examples
#' store <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' a <- pt_read(store, threads = 1)
#' set.seed(1)
#' networks <- pt_sample_networks(a, 100)
#' # Out-degree of each Infection in each network.
#' table(unlist(networks$members))
#' pt_close(a)
pt_sample_networks <- function(a, n) {
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 0 || n != round(n)) {
    cli::cli_abort("{.arg n} must be a whole number of networks, 0 or more.")
  }
  runs <- runs_with_draws(a)
  index <- draw_index(a, runs)
  picked <- index[sample.int(nrow(index), n, replace = TRUE), ]
  sets <- picked_sets(a, picked)

  draw_key <- row_key(sets$run, sets$iteration)
  child_key <- row_key(draw_key, sets$infection)
  by_child <- split(seq_len(nrow(sets)), factor(child_key, levels = unique(child_key)))
  children_of <- split(names(by_child), factor(sub("\r[^\r]*$", "", names(by_child)), levels = unique(draw_key)))

  # One uniform per (network, child), networks in order, children in string order.
  rows <- lapply(seq_len(n), function(s) {
    children <- children_of[[row_key(picked$run[[s]], picked$iteration[[s]])]]
    u <- stats::runif(length(children))
    vapply(seq_along(children), function(k) {
      i <- by_child[[children[[k]]]]
      cumulative <- cumsum(sets$prob[i])
      chosen <- which(cumulative > u[[k]] * cumulative[[length(i)]])
      i[[if (length(chosen)) chosen[[1]] else length(i)]]
    }, 0L)
  })
  chosen <- unlist(rows)
  observed <- sets$observed[chosen]
  tibble::tibble(
    sample = rep(seq_len(n), lengths(rows)), run = sets$run[chosen], iteration = sets$iteration[chosen],
    infection = sets$infection[chosen],
    members = lapply(observed, function(x) if (x == "") character() else strsplit(x, "\x1f", fixed = TRUE)[[1]]),
    has_latent_parent = sets$has_latent_parent[chosen]
  )
}

# The Parent sets with weight at the picked Draws: run, iteration, infection, set_id, prob, observed (canonical),
# has_latent_parent, sorted by run, iteration, infection (string order) and set_id.
picked_sets <- function(a, picked) {
  con <- analysis_con(a)
  draws <- unique(picked[c("run", "iteration")])
  if (nrow(draws) == 0) {
    return(tibble::tibble(run = character(), iteration = numeric(), infection = character(), set_id = numeric(),
                          prob = numeric(), observed = character(), has_latent_parent = logical()))
  }
  runs <- unique(draws$run)
  rows <- pt_table(a, "parent_sets", runs = runs, phase = "sampling") |>
    dplyr::select("run", "iteration", "infection", "set_id", "prob") |>
    dplyr::inner_join(dbplyr::copy_inline(con, draws), by = c("run", "iteration")) |>
    dplyr::inner_join(set_members(a, runs), by = c("run", "infection", "set_id")) |>
    dplyr::collect()
  rows <- rows[order(match(rows$run, runs), rows$iteration, rows$infection, rows$set_id, method = "radix"), ]
  rows$prob <- as.numeric(rows$prob)
  rows$observed <- canonical_observed(rows$observed)
  rows
}
