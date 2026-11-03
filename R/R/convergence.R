# Convergence across Runs: the Draw set, per-quantity diagnostics (#124) and the tree diagnostics (#125).
#
#   convergence_draw_set(a, runs, iterations)  the Draw set: each Run's first n sampling Draws (n the shortest Run's
#                                count, within `iterations`). list(runs, n, iterations = named list, counts, used).
#   in_draw_set(x, ds)           lazy: the rows of a per-Draw table whose Draw is in the Draw set.
#   edge_run_means(a, ds)        eager: run, child, parent, mean: every edge's Run mean over the Draw set,
#                                zero-filled (a Run without the edge gets 0). The edge floor and #125's spread use it.
#   convergence_groups           quantity groups, in report order.
#   group_array(a, ds, group, edge_floor)  n x Runs x variables array of one group, variables in name order.
#   diagnose_array(x, group, rhat, ess)    the quantities rows of one group.

convergence_groups <- c("log_posterior", "global", "expected_false_positives", "expected_false_negatives",
                        "infection_duration", "allele_frequency", "edge", "p_latent_parent", "functional")

#' Convergence across Runs
#'
#' Checks whether the Runs of an Analysis agree, quantity by quantity. Every diagnostic uses the same Draws: each
#' Run's first *n* sampling Draws, where *n* is the shortest Run's count (after `iterations`, when given). Pooled
#' estimates elsewhere in the package still use every Draw.
#'
#' **Quantities.** Each row of `quantities` is one quantity with a value at every Draw, named as its variable in
#' [pt_as_draws()]:
#'
#' * `log_posterior`: `log_posterior`, the log prior plus the log likelihood.
#' * `global`: `mean_coi`, `mean_strains_tx`, `parent_set_size_prob`.
#' * `expected_false_positives`, `expected_false_negatives`, `infection_duration`: one per Infection, such as
#'   `infection_duration[12]`. The per-Allele error rate is the expected count over the Infection's Alleles, a
#'   constant factor, so its diagnostics are the count's.
#' * `allele_frequency`: one per Locus and Allele (0-based), such as `allele_frequency[AS1,0]`.
#' * `edge`: the per-Draw Edge weight, zero-filled, such as `edge[12,7]` for child 12 and parent 7. Only edges whose
#'   mean in some Run is at least `edge_floor` are diagnosed.
#' * `p_latent_parent`: one per Infection, the Latent parent probability (the manuscript's q).
#' * `functional`: `mean_out_degree` (ψ_edge), `mean_p_latent_parent`, `multi_source_augmented`,
#'   `multi_source_observed`, `mean_entropy`, `mean_effective_size`, and one `composition[<class>]` per composition
#'   class (see [pt_mean_out_degree()]).
#'
#' **Status.** Each quantity is classified before `posterior` is called:
#'
#' * `constant`: one value over every Draw of every Run, such as an edge that never appears. `rhat`, `ess_bulk` and
#'   `ess_tail` are `NA`, and the quantity is never flagged. `NA` R-hat means only this.
#' * `stuck`: every Run is constant, but the Runs differ. `rhat` is `Inf`, and the quantity is always flagged.
#' * `mixing`: everything else. `rhat` is [posterior::rhat()], or [posterior::rhat_basic()] when that is `NA`:
#'   when the Draws are exactly half one value and half another, `rhat()`'s folded R-hat is undefined.
#'
#' `ess_bulk` and `ess_tail` are [posterior::ess_bulk()] and [posterior::ess_tail()]. `ess_tail` is `NA` for most 0/1
#' quantities and for edges that are 0 in more than 5% of the Draws: the 5% or 95% quantile is then the smallest or
#' largest value, and the tail indicator is constant.
#'
#' A quantity is `flagged` when it is stuck, or mixing with `rhat` above `rhat` or `ess_bulk` below `ess`. The thresholds
#' follow Vehtari et al. (2021): R-hat 1.01, and 100 effective Draws per Run.
#'
#' @inheritParams pt_table
#' @param runs Run labels to check (see [pt_runs()]). `NULL` checks every Run. Runs with no sampling Draws are
#'   skipped with a message.
#' @param iterations `c(first, last)`: use only sampling Draws with iteration in this range (inclusive), such as the
#'   second half of a Run. `NULL` uses every sampling Draw.
#' @param edge_floor Diagnose an edge only when its mean in some Run is at least this.
#' @param rhat Flag quantities with R-hat above this.
#' @param ess Flag quantities with `ess_bulk` below this. `NULL` is 100 times the number of Runs.
#' @param spread Flag Infections and order pairs whose Spread (max minus min of the Run means) is above this.
#' @return A `pt_convergence`: a list of tibbles.
#'   * `runs`: `run`, `draws` (the Draws used, the same in every Run), `first_iteration` and `last_iteration` of
#'     those Draws, `log_posterior` (the Run's mean log posterior over them) and `gap` (the best Run's mean minus
#'     this Run's, in log units).
#'   * `quantities`: `group`, `name`, `status`, `rhat`, `ess_bulk`, `ess_tail` and `flagged`, in the group order
#'     above, names sorted by their keys within a group.
#' @seealso [pt_as_draws()] for the Draws themselves, to use any `posterior` diagnostic.
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Bürkner, P.-C. (2021). Rank-normalization,
#'   folding, and localization: an improved R-hat for assessing convergence of MCMC. *Bayesian Analysis* 16(2).
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' input <- system.file("extdata", "example_input.json", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c("single-chain", "replica-exchange")), input = input, threads = 1)
#' conv <- pt_convergence(a)
#' conv
#' conv$quantities[conv$quantities$flagged, ]
#' pt_close(a)
pt_convergence <- function(a, runs = NULL, iterations = NULL, edge_floor = 0.01, rhat = 1.01, ess = NULL,
                           spread = 0.2) {
  check_fraction(edge_floor, "edge_floor")
  check_positive(rhat, "rhat")
  if (!is.null(ess)) check_positive(ess, "ess")
  check_fraction(spread, "spread")
  ds <- convergence_draw_set(a, runs, iterations)
  if (is.null(ess)) ess <- 100 * length(ds$runs)

  log_posterior <- group_array(a, ds, "log_posterior")
  quantities <- lapply(convergence_groups, function(group) {
    x <- if (group == "log_posterior") log_posterior else group_array(a, ds, group, edge_floor)
    diagnose_array(x, group, rhat, ess)
  })

  means <- colMeans(matrix(log_posterior, ds$n, length(ds$runs)))
  out <- list(
    runs = tibble::tibble(
      run = ds$runs, draws = ds$n, first_iteration = ds$used$first_iteration,
      last_iteration = ds$used$last_iteration, log_posterior = means, gap = max(means) - means
    ),
    quantities = dplyr::bind_rows(quantities)
  )
  structure(out, class = "pt_convergence", settings = list(
    iterations = iterations, edge_floor = edge_floor, rhat = rhat, ess = ess, spread = spread,
    counts = ds$counts
  ))
}

#' Convergence Draws as a posterior draws array
#'
#' The quantities of [pt_convergence()] as a [posterior::draws_array()]: iterations by Runs by variables, for any
#' `posterior` or `bayesplot` function. Runs are cut to the shortest Run's Draw count, with a message, since
#' `posterior` needs Runs of equal length. Chains are the Runs in [pt_runs()] order (the `runs` table of
#' [pt_convergence()]).
#'
#' R* is not in [pt_convergence()] because it needs `caret`. With `caret` and `randomForest` installed,
#' `posterior::rstar(pt_as_draws(a, groups = "functional"))` gives it.
#'
#' @inheritParams pt_convergence
#' @param groups Quantity groups to include (see [pt_convergence()]). `NULL` includes every group.
#' @return A `draws_array`. Variables are named as `name` in `pt_convergence()$quantities`.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' input <- system.file("extdata", "example_input.json", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c("single-chain", "replica-exchange")), input = input, threads = 1)
#' d <- pt_as_draws(a, groups = c("log_posterior", "functional"))
#' posterior::summarise_draws(d)
#' pt_close(a)
pt_as_draws <- function(a, groups = NULL, runs = NULL, iterations = NULL, edge_floor = 0.01) {
  if (is.null(groups)) groups <- convergence_groups
  if (!is.character(groups) || length(groups) == 0 || !all(groups %in% convergence_groups)) {
    cli::cli_abort("{.arg groups} must be one or more of {.val {convergence_groups}}.")
  }
  check_fraction(edge_floor, "edge_floor")
  ds <- convergence_draw_set(a, runs, iterations)
  if (length(unique(ds$counts)) > 1) {
    cli::cli_inform(c(
      "Truncating every Run to the shortest Run's {ds$n} Draw{?s}.",
      i = "Draws per Run: {paste(names(ds$counts), ds$counts, sep = ' ', collapse = ', ')}."
    ))
  }
  arrays <- lapply(convergence_groups[convergence_groups %in% groups], group_array, a = a, ds = ds,
                   edge_floor = edge_floor)
  variables <- unlist(lapply(arrays, function(x) dimnames(x)[[3]]))
  out <- array(NA_real_, c(ds$n, length(ds$runs), length(variables)),
               dimnames = list(iteration = NULL, chain = NULL, variable = variables))
  at <- 0
  for (x in arrays) {
    k <- dim(x)[[3]]
    if (k) out[, , at + seq_len(k)] <- x
    at <- at + k
  }
  posterior::as_draws_array(out)
}

# ---- The Draw set ----

convergence_draw_set <- function(a, runs = NULL, iterations = NULL, call = rlang::caller_env()) {
  runs <- runs_with_draws(a, runs, call = call)
  thin <- stats::setNames(vapply(a$runs$manifest, function(m) as.numeric(m$thin), numeric(1)), a$runs$run)[runs]
  if (length(unique(thin)) > 1) {
    cli::cli_abort(c(
      "Runs have different thinning strides; their Draws are not comparable lag by lag.",
      rlang::set_names(cli_escape(sprintf("Run %s: thin %s", runs, format(thin, trim = TRUE))), "*"),
      i = "Check Runs of one stride together with {.arg runs}."
    ), call = call)
  }
  grid <- draw_grid(a, runs)
  if (!is.null(iterations)) {
    if (!is.numeric(iterations) || length(iterations) != 2 || anyNA(iterations) || iterations[[1]] > iterations[[2]]) {
      cli::cli_abort("{.arg iterations} must be {.code c(first, last)} with {.code first <= last}.", call = call)
    }
    first <- iterations[[1]]
    last <- iterations[[2]]
    grid <- dplyr::filter(grid, .data$iteration >= !!first, .data$iteration <= !!last)
  }
  grid <- dplyr::collect(grid)
  its <- lapply(stats::setNames(runs, runs), function(r) sort(as.numeric(grid$iteration[grid$run == r])))
  empty <- lengths(its) == 0
  if (any(empty)) {
    if (all(empty)) {
      cli::cli_abort("No sampling Draws with iteration in {iterations[[1]]}-{iterations[[2]]}.", call = call)
    }
    cli::cli_inform("Skipping Run{?s} {.val {runs[empty]}}: no sampling Draws with iteration in
                     {iterations[[1]]}-{iterations[[2]]}.")
    its <- its[!empty]
  }
  counts <- lengths(its)
  n <- min(counts)
  its <- lapply(its, function(x) x[seq_len(n)])
  list(
    runs = names(its), n = n, iterations = its, counts = counts,
    used = tibble::tibble(run = names(its), first_iteration = vapply(its, function(x) x[[1]], numeric(1)),
                          last_iteration = vapply(its, function(x) x[[n]], numeric(1)))
  )
}

# The rows of a lazy per-Draw table (run, iteration, ...) whose Draw is in the Draw set. A Run's Draws are its
# sorted iterations, so its first n are those up to its n-th.
in_draw_set <- function(x, ds) {
  bounds <- dbplyr::copy_inline(dbplyr::remote_con(x), data.frame(
    run = ds$used$run, draw_set_first = ds$used$first_iteration, draw_set_last = ds$used$last_iteration
  ))
  x |>
    dplyr::inner_join(bounds, by = "run") |>
    dplyr::filter(.data$iteration >= .data$draw_set_first, .data$iteration <= .data$draw_set_last) |>
    dplyr::select(-"draw_set_first", -"draw_set_last")
}

edge_run_means <- function(a, ds) {
  totals <- in_draw_set(observed_edges(a, ds$runs), ds) |>
    dplyr::group_by(.data$run, .data$child, .data$parent) |>
    dplyr::summarise(total = sum(.data$prob, na.rm = TRUE), .groups = "drop") |>
    dplyr::collect()
  edges <- unique(totals[c("child", "parent")])
  edges <- edges[order(edges$child, edges$parent, method = "radix"), ]
  k <- nrow(edges)
  means <- tibble::tibble(run = rep(ds$runs, each = k), child = rep(edges$child, length(ds$runs)),
                          parent = rep(edges$parent, length(ds$runs)), mean = 0)
  at <- match(paste(totals$run, totals$child, totals$parent, sep = "\r"),
              paste(means$run, means$child, means$parent, sep = "\r"))
  means$mean[at] <- totals$total / ds$n
  means
}

# ---- Quantity arrays ----

# One group's values: an n x Runs x variables array, Runs in Draw-set order, variables in name order.
group_array <- function(a, ds, group, edge_floor = 0.01) {
  g <- group_values(a, ds, group, edge_floor)
  x <- array(g$fill, c(ds$n, length(ds$runs), length(g$names)), dimnames = list(NULL, ds$runs, g$names))
  rows <- g$rows
  if (nrow(rows)) {
    draw <- unlist(lapply(ds$iterations, seq_along))
    names(draw) <- paste(rep(ds$runs, each = ds$n), unlist(ds$iterations), sep = "\r")
    x[cbind(draw[paste(rows$run, as.numeric(rows$iteration), sep = "\r")], match(rows$run, ds$runs),
            match(rows$name, g$names))] <- rows$value
  }
  x
}

# Long rows (run, iteration, name, value) of one group over the Draw set, the group's variable names in order, and
# the value of a Draw with no row (0 for edges, NA otherwise: every other Trace has a row at every Draw).
group_values <- function(a, ds, group, edge_floor) {
  runs <- ds$runs
  sorted <- function(x) x[order(x, method = "radix")]
  ids <- sorted(a$infections$infection)
  long <- function(x, columns, names = columns) {
    x <- dplyr::collect(in_draw_set(x, ds))
    dplyr::bind_rows(lapply(seq_along(columns), function(k) {
      tibble::tibble(run = x$run, iteration = x$iteration, name = names[[k]], value = as.numeric(x[[columns[[k]]]]))
    }))
  }
  keyed <- function(x, column, key, prefix) {
    x <- dplyr::collect(in_draw_set(x, ds))
    tibble::tibble(run = x$run, iteration = x$iteration, name = sprintf("%s[%s]", prefix, x[[key]]),
                   value = as.numeric(x[[column]]))
  }
  out <- function(rows, names, fill = NA_real_) list(rows = rows, names = names, fill = fill)

  switch(group,
    log_posterior = out(
      long(dplyr::mutate(pt_table(a, "run_scalars", runs = runs), log_posterior = .data$prior + .data$likelihood),
           "log_posterior"),
      "log_posterior"
    ),
    global = {
      columns <- c("mean_coi", "mean_strains_tx", "parent_set_size_prob")
      out(long(pt_table(a, "run_scalars", runs = runs), columns), sorted(columns))
    },
    expected_false_positives = ,
    expected_false_negatives = ,
    infection_duration = out(
      keyed(pt_table(a, "infections", runs = runs), group, "infection", group),
      sprintf("%s[%s]", group, ids)
    ),
    allele_frequency = {
      x <- dplyr::collect(in_draw_set(pt_table(a, "allele_frequencies", runs = runs), ds))
      loci <- a$loci[order(a$loci$locus, method = "radix"), ]
      names <- sprintf("allele_frequency[%s,%d]", rep(loci$locus, loci$n_alleles),
                       unlist(lapply(loci$n_alleles, function(k) seq_len(k) - 1L)))
      out(tibble::tibble(run = x$run, iteration = x$iteration,
                         name = sprintf("allele_frequency[%s,%d]", x$locus, as.integer(x$allele)),
                         value = as.numeric(x$frequency)), names)
    },
    edge = {
      means <- edge_run_means(a, ds)
      top <- tapply(means$mean, paste(means$child, means$parent, sep = "\r"), max)
      kept <- unique(means[c("child", "parent")])
      kept <- kept[top[paste(kept$child, kept$parent, sep = "\r")] >= edge_floor, ]
      names <- sprintf("edge[%s,%s]", kept$child, kept$parent)
      rows <- tibble::tibble(run = character(), iteration = numeric(), name = character(), value = numeric())
      if (nrow(kept)) {
        edges <- observed_edges(a, runs) |>
          dplyr::semi_join(dbplyr::copy_inline(analysis_con(a), as.data.frame(kept)), by = c("child", "parent"))
        x <- dplyr::collect(in_draw_set(edges, ds))
        rows <- tibble::tibble(run = x$run, iteration = x$iteration, name = sprintf("edge[%s,%s]", x$child, x$parent),
                               value = as.numeric(x$prob))
      }
      out(rows, names, fill = 0)
    },
    p_latent_parent = out(
      keyed(infection_values(a, runs, "sampling"), "p_latent_parent", "infection", "p_latent_parent"),
      sprintf("p_latent_parent[%s]", ids)
    ),
    functional = {
      composition <- dplyr::collect(in_draw_set(pt_composition(a, runs), ds))
      rows <- dplyr::bind_rows(
        long(pt_mean_out_degree(a, runs), "mean_out_degree"),
        long(pt_mean_p_latent_parent(a, runs), "mean_p_latent_parent"),
        long(pt_multi_source(a, runs), c("augmented", "observed_only"),
             c("multi_source_augmented", "multi_source_observed")),
        long(pt_entropy(a, runs), c("mean_entropy", "mean_effective_size")),
        tibble::tibble(run = composition$run, iteration = composition$iteration,
                       name = sprintf("composition[%s]", composition$class), value = composition$share)
      )
      out(rows, sorted(unique(c(
        "mean_out_degree", "mean_p_latent_parent", "multi_source_augmented", "multi_source_observed",
        "mean_entropy", "mean_effective_size", sprintf("composition[%s]", unique(composition$class))
      ))))
    }
  )
}

# ---- Diagnostics ----

diagnose_array <- function(x, group, rhat, ess) {
  n <- dim(x)[[1]]
  chains <- dim(x)[[2]]
  rows <- lapply(seq_len(dim(x)[[3]]), function(k) diagnose_quantity(matrix(x[, , k], n, chains), rhat, ess))
  tibble::tibble(
    group = rep(group, length(rows)), name = if (is.null(dimnames(x)[[3]])) character() else dimnames(x)[[3]],
    status = vapply(rows, `[[`, character(1), "status"),
    rhat = vapply(rows, `[[`, numeric(1), "rhat"),
    ess_bulk = vapply(rows, `[[`, numeric(1), "ess_bulk"),
    ess_tail = vapply(rows, `[[`, numeric(1), "ess_tail"),
    flagged = vapply(rows, `[[`, logical(1), "flagged")
  )
}

# Status first: constant (one value everywhere), stuck (each Run constant, the Runs not), else mixing. A matrix
# with a missing or non-finite value is mixing, and posterior gives it NA.
diagnose_quantity <- function(m, rhat, ess) {
  finite <- all(is.finite(m))
  status <- if (finite && all(m == m[[1]])) "constant"
            else if (finite && all(apply(m, 2, function(col) all(col == col[[1]])))) "stuck"
            else "mixing"
  if (status == "constant") {
    return(list(status = status, rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_, flagged = FALSE))
  }
  r <- if (status == "stuck") Inf else posterior::rhat(m)
  if (is.na(r)) r <- posterior::rhat_basic(m)
  bulk <- posterior::ess_bulk(m)
  flagged <- status == "stuck" || isTRUE(r > rhat) || isTRUE(bulk < ess)
  list(status = status, rhat = r, ess_bulk = bulk, ess_tail = posterior::ess_tail(m), flagged = flagged)
}

# ---- Argument checks ----

check_fraction <- function(x, arg, call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < 0 || x > 1) {
    cli::cli_abort("{.arg {arg}} must be a number between 0 and 1.", call = call)
  }
}

check_positive <- function(x, arg, call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
    cli::cli_abort("{.arg {arg}} must be a positive number.", call = call)
  }
}
