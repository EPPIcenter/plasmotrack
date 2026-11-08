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
#   tree_infections(), tree_pairs()        the infections and pairs tables (Spread, top parents, classes).

convergence_groups <- c("log_posterior", "global", "expected_false_positives", "expected_false_negatives",
                        "infection_duration", "allele_frequency", "edge", "order_pair", "p_latent_parent", "functional")

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
#' * `global`: `mean_coi`, `mean_strains_tx`, `parent_set_size_prob`, and the Run scalars a Leak prior or the
#'   hierarchical Error model adds when the Runs have them: `transmission_leak`, `efp_mean`, `efp_shape`, `efn_mean`,
#'   `efn_shape`.
#' * `expected_false_positives`, `expected_false_negatives`, `infection_duration`: one per Infection, such as
#'   `infection_duration[12]`. The per-Allele error rate is the expected count over the Infection's Alleles, a
#'   constant factor, so its diagnostics are the count's.
#' * `allele_frequency`: one per Locus and Allele (0-based), such as `allele_frequency[AS1,0]`.
#' * `edge`: the edge's probability at each Draw (the per-Draw Edge marginal), zero-filled, such as `edge[12,7]` for child 12 and parent 7. Only edges whose
#'   mean in some Run is at least `edge_floor` are diagnosed.
#' * `order_pair`: the 0/1 indicator that `i` was acquired before `j` (see [pt_order_pair_draws()]), such as
#'   `order_pair[10,2]` with `i` the smaller id as a string. Only pairs joined by a diagnosed edge, in either
#'   direction, are diagnosed: the full set is quadratic in the Infections.
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
#' **Tree diagnostics.** A quantity's *Run mean* is its mean over one Run's Draws, and its *Spread* is the largest Run
#' mean minus the smallest. An Infection's parent keys are every Infection that is its parent at some Draw of some Run,
#' plus its Latent parent; an edge's Run mean is its Edge marginal over the Run's Draws, 0 in a Run that never has it,
#' and the Latent parent's is the Run mean of `p_latent_parent`. Edges below `edge_floor` count here. The Infection's
#' spread is the largest spread over its keys, and its *top parent* in a Run is the key with the largest Run mean.
#' Ties go to the smaller id in string order (`"10"` before `"2"`), and an observed parent wins a tie against the
#' Latent parent. An Infection is flagged when its spread is above `spread`, and then gets the first disagreement
#' class that matches, with *t* = `spread`:
#'
#' 1. `same_top`: every Run has the same top parent.
#' 2. `root_flip`: some Run's top is the Latent parent, and in that Run the Infection infects another Run's top
#'    parent with Edge marginal above *t*.
#' 3. `latent_vs_observed`: one Run's top is the Latent parent, another's is an observed parent, and rule 2 doesn't
#'    apply.
#' 4. `direction_flip`: in some Run, the Infection infects another Run's top parent with Edge marginal above *t*.
#' 5. `linked_parents`: the top parents differ, and in some Run one of them infects another with Edge marginal above
#'    *t*.
#' 6. `different_parent`: everything else.
#'
#' The Log-posterior gap (the spread of the log posterior) is reported but not flagged on its own: the log
#' posterior's R-hat flags it.
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
#'   * `infections`: one row per Infection, in id (string) order. `infection`, `spread`, `class` (the disagreement
#'     class above; `NA` when not flagged), `top_parent` (a list column: per row, a character vector of top parents
#'     named by Run, `"latent"` for the Latent parent; an Infection with id `"latent"` is an error, since the two
#'     would be indistinguishable), `worst_rhat` (the largest `rhat` of the quantities the
#'     Infection owns: its `expected_false_positives`, `expected_false_negatives`, `infection_duration`,
#'     `p_latent_parent` and the edges into it; `NA` when all are constant) and `flagged`.
#'   * `pairs`: one row per diagnosed order pair (the `order_pair` quantities). `i`, `j`, `p_before` (a list column:
#'     P(i acquired before j) named by Run), `spread` and `flagged` (spread above `spread`).
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
#' conv$infections
#' pt_close(a)
pt_convergence <- function(a, runs = NULL, iterations = NULL, edge_floor = 0.01, rhat = 1.01, ess = NULL,
                           spread = 0.2) {
  check_fraction(edge_floor, "edge_floor")
  check_positive(rhat, "rhat")
  if (!is.null(ess)) check_positive(ess, "ess")
  check_fraction(spread, "spread")
  check_no_latent_id(a)
  ds <- convergence_draw_set(a, runs, iterations)
  if (is.null(ess)) ess <- 100 * length(ds$runs)

  arrays <- lapply(stats::setNames(convergence_groups, convergence_groups), group_array, a = a, ds = ds,
                   edge_floor = edge_floor)
  quantities <- dplyr::bind_rows(lapply(convergence_groups, function(group) {
    diagnose_array(arrays[[group]], group, rhat, ess)
  }))

  means <- colMeans(matrix(arrays$log_posterior, ds$n, length(ds$runs)))
  out <- list(
    runs = tibble::tibble(
      run = ds$runs, draws = ds$n, first_iteration = ds$used$first_iteration,
      last_iteration = ds$used$last_iteration, log_posterior = means, gap = max(means) - means
    ),
    quantities = quantities,
    infections = tree_infections(a, ds, run_means(arrays$p_latent_parent), quantities, spread),
    pairs = tree_pairs(run_means(arrays$order_pair), spread)
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
      cli::cli_abort("No sampling Draws with iteration in {paste(iterations, collapse = '-')}.", call = call)
    }
    skipped <- runs[empty]
    range <- paste(iterations, collapse = "-")
    cli::cli_inform(c("Skipping {cli::qty(skipped)}Run{?s} {.val {skipped}}:",
                       " " = "no sampling Draws with iteration in {range}."))
    its <- its[!empty]
  }
  counts <- lengths(its)
  n <- min(counts)
  its <- lapply(its, function(x) x[seq_len(n)])
  list(
    runs = names(its), n = n, iterations = its, counts = counts,
    used = tibble::tibble(run = names(its), first_iteration = unname(vapply(its, function(x) x[[1]], numeric(1))),
                          last_iteration = unname(vapply(its, function(x) x[[n]], numeric(1))))
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
  at <- match(row_key(totals$run, totals$child, totals$parent),
              row_key(means$run, means$child, means$parent))
  means$mean[at] <- totals$total / ds$n
  means
}

# The edges diagnosed: those whose Run mean over the Draw set is at least edge_floor in some Run, sorted.
kept_edges <- function(a, ds, edge_floor) {
  means <- edge_run_means(a, ds)
  top <- tapply(means$mean, row_key(means$child, means$parent), max)
  kept <- unique(means[c("child", "parent")])
  as.data.frame(kept[top[row_key(kept$child, kept$parent)] >= edge_floor, ])
}

# ---- Quantity arrays ----

# One group's values: an n x Runs x variables array, Runs in Draw-set order, variables in name order.
group_array <- function(a, ds, group, edge_floor = 0.01) {
  g <- group_values(a, ds, group, edge_floor)
  x <- array(g$fill, c(ds$n, length(ds$runs), length(g$names)), dimnames = list(NULL, ds$runs, g$names))
  rows <- g$rows
  if (nrow(rows)) {
    draw <- unlist(lapply(ds$iterations, seq_along))
    names(draw) <- row_key(rep(ds$runs, each = ds$n), unlist(ds$iterations))
    x[cbind(draw[row_key(rows$run, as.numeric(rows$iteration))], match(rows$run, ds$runs),
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
      scalars <- pt_table(a, "run_scalars", runs = runs)
      columns <- c("mean_coi", "mean_strains_tx", "parent_set_size_prob",
                   intersect(model_choice_scalars, colnames(scalars)))
      out(long(scalars, columns), sorted(columns))
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
      kept <- kept_edges(a, ds, edge_floor)
      names <- sprintf("edge[%s,%s]", kept$child, kept$parent)
      rows <- tibble::tibble(run = character(), iteration = numeric(), name = character(), value = numeric())
      if (nrow(kept)) {
        edges <- observed_edges(a, runs) |>
          dplyr::semi_join(dbplyr::copy_inline(analysis_con(a), kept), by = c("child", "parent"))
        x <- dplyr::collect(in_draw_set(edges, ds))
        rows <- tibble::tibble(run = x$run, iteration = x$iteration, name = sprintf("edge[%s,%s]", x$child, x$parent),
                               value = as.numeric(x$prob))
      }
      out(rows, names, fill = 0)
    },
    order_pair = {
      kept <- kept_edges(a, ds, edge_floor)
      # i < j in string (byte) order, as in pt_order_pair_draws().
      swap <- vapply(seq_len(nrow(kept)), function(k) {
        order(c(kept$child[[k]], kept$parent[[k]]), method = "radix")[[1]] == 2L
      }, logical(1))
      pairs <- unique(data.frame(i = ifelse(swap, kept$parent, kept$child), j = ifelse(swap, kept$child, kept$parent)))
      pairs <- pairs[order(pairs$i, pairs$j, method = "radix"), , drop = FALSE]
      rows <- tibble::tibble(run = character(), iteration = numeric(), name = character(), value = numeric())
      if (nrow(pairs)) {
        x <- order_pairs_lazy(a, runs, "sampling", infections = unique(c(pairs$i, pairs$j))) |>
          dplyr::semi_join(dbplyr::copy_inline(analysis_con(a), pairs), by = c("i", "j")) |>
          in_draw_set(ds = ds) |>
          dplyr::collect()
        rows <- tibble::tibble(run = x$run, iteration = x$iteration, name = sprintf("order_pair[%s,%s]", x$i, x$j),
                               value = as.numeric(x$before))
      }
      out(rows, sprintf("order_pair[%s,%s]", pairs$i, pairs$j))
    },
    p_latent_parent = out(
      keyed(infection_values(a, runs, "sampling"), "p_latent_parent", "infection", "p_latent_parent"),
      sprintf("p_latent_parent[%s]", ids)
    ),
    functional = {
      composition <- dplyr::collect(in_draw_set(pt_composition(a, runs = runs), ds))
      rows <- dplyr::bind_rows(
        long(pt_mean_out_degree(a, runs = runs), "mean_out_degree"),
        long(pt_mean_p_latent_parent(a, runs = runs), "mean_p_latent_parent"),
        long(pt_multi_source(a, runs = runs), c("augmented", "observed_only"),
             c("multi_source_augmented", "multi_source_observed")),
        long(pt_entropy(a, runs = runs), c("mean_entropy", "mean_effective_size")),
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

# ---- Tree diagnostics ----

# The key of the Latent parent among an Infection's parent keys (top_parent values).
latent_key <- "latent"

# `top_parent` names the Latent parent `latent_key`, so an Infection with that id would read as the Latent parent.
check_no_latent_id <- function(a, call = rlang::caller_env()) {
  if (latent_key %in% a$infections$infection) {
    cli::cli_abort(c(
      "Infection id {.val {latent_key}} collides with the Latent parent's key in {.field top_parent}.",
      i = "Rename that Infection in the input so the two can be told apart."
    ), call = call)
  }
}

# Run means of an n x Runs x variables array: a Runs x variables matrix.
run_means <- function(x) {
  matrix(colMeans(x), dim(x)[[2]], dim(x)[[3]], dimnames = dimnames(x)[2:3])
}

# One row per Infection (string order): spread over its parent keys (zero-filled Edge marginals over the Draw set
# plus the Latent parent), each Run's top parent, the disagreement class of flagged Infections, the worst R-hat of
# the quantities it owns. `latent_means`: Runs x p_latent_parent[id] matrix.
tree_infections <- function(a, ds, latent_means, quantities, spread) {
  runs <- ds$runs
  ids <- a$infections$infection[order(a$infections$infection, method = "radix")]
  edges <- edge_run_means(a, ds)
  edge_mean <- stats::setNames(edges$mean, row_key(edges$run, edges$child, edges$parent))
  edge <- function(run, child, parent) {
    if (length(child) == 0 || length(parent) == 0) return(numeric())
    v <- edge_mean[row_key(run, child, parent)]
    ifelse(is.na(v), 0, unname(v))
  }
  by_child <- split(edges, edges$child)
  rows <- lapply(ids, function(id) {
    mine <- by_child[[id]]
    parents <- if (is.null(mine)) character() else unique(mine$parent)
    parents <- parents[order(parents, method = "radix")]
    # keys x Runs: observed parents in id order, then the Latent parent.
    m <- rbind(
      matrix(vapply(runs, function(r) edge(r, id, parents), numeric(length(parents))), length(parents),
             length(runs)),
      unname(latent_means[, sprintf("p_latent_parent[%s]", id)])
    )
    keys <- c(parents, latent_key)
    s <- max(apply(m, 1, function(v) max(v) - min(v)))
    # Ties (to rounding) go to the first key: the smaller id, observed before the Latent parent.
    tops <- vapply(seq_along(runs), function(k) keys[[which(m[, k] >= max(m[, k]) - 1e-9)[[1]]]], character(1))
    list(spread = s, top = stats::setNames(tops, runs),
         class = if (s > spread) disagreement_class(id, tops, runs, edge, spread) else NA_character_)
  })
  owner <- quantity_owner(quantities)
  worst <- vapply(ids, function(id) {
    r <- quantities$rhat[owner %in% id & !is.na(quantities$rhat)]
    if (length(r)) max(r) else NA_real_
  }, numeric(1))
  spreads <- vapply(rows, `[[`, numeric(1), "spread")
  tibble::tibble(
    infection = ids, spread = spreads, class = vapply(rows, `[[`, character(1), "class"),
    top_parent = lapply(rows, `[[`, "top"), worst_rhat = unname(worst), flagged = spreads > !!spread
  )
}

# The first class, in `disagreement_classes` order, whose rule matches (the rules are documented on pt_convergence()).
# `edge(run, child, parent)` is a Run mean; `tops` each Run's top parent, in `runs` order.
disagreement_class <- function(child, tops, runs, edge, t) {
  if (length(unique(tops)) == 1) return("same_top")
  observed <- setdiff(unique(tops), latent_key)
  if (latent_key %in% tops && length(observed)) {
    for (k in which(tops == latent_key)) {
      if (any(edge(runs[[k]], observed, child) > t)) return("root_flip")
    }
    return("latent_vs_observed")
  }
  for (k in seq_along(runs)) {
    others <- setdiff(observed, tops[[k]])
    if (any(edge(runs[[k]], others, child) > t)) return("direction_flip")
  }
  pairs <- expand.grid(x = observed, y = observed, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$x != pairs$y, ]
  for (run in runs) {
    if (any(edge(run, pairs$x, pairs$y) > t)) return("linked_parents")
  }
  "different_parent"
}

# The Infection that owns each quantity (NA for none): the id of a per-Infection quantity, the child of an edge.
quantity_owner <- function(quantities) {
  name <- quantities$name
  owner <- rep(NA_character_, length(name))
  single <- quantities$group %in% c("expected_false_positives", "expected_false_negatives", "infection_duration",
                                    "p_latent_parent")
  owner[single] <- sub("^[^[]*\\[(.*)\\]$", "\\1", name[single])
  edge <- quantities$group == "edge"
  owner[edge] <- sub("^edge\\[(.*),[^,]*\\]$", "\\1", name[edge])
  owner
}

# One row per diagnosed order pair: each Run's P(i before j), its spread, flagged. `means`: Runs x order_pair[i,j].
tree_pairs <- function(means, spread) {
  names <- colnames(means)
  if (is.null(names)) names <- character()
  ends <- strsplit(sub("^order_pair\\[(.*)\\]$", "\\1", names), ",", fixed = TRUE)
  spreads <- if (length(names)) apply(means, 2, function(v) max(v) - min(v)) else numeric()
  tibble::tibble(
    i = vapply(ends, `[[`, character(1), 1), j = vapply(ends, `[[`, character(1), 2),
    p_before = lapply(seq_along(names), function(k) stats::setNames(means[, k], rownames(means))),
    spread = unname(spreads), flagged = unname(spreads) > !!spread
  )
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

# ---- Print ----

#' @export
print.pt_convergence <- function(x, ...) {
  for (section in convergence_print_sections) section(x)
  invisible(x)
}

# The report, one function per section, in order.
convergence_print_sections <- list(
  function(x) {
    s <- attr(x, "settings")
    runs <- x$runs
    n <- runs$draws[[1]]
    cat(sprintf("<pt_convergence> %d Run%s, %s Draw%s each%s\n", nrow(runs), if (nrow(runs) == 1) "" else "s",
                format_count(n), if (n == 1) "" else "s",
                if (is.null(s$iterations)) "" else sprintf(", iterations %s", paste(format_count(s$iterations),
                                                                                    collapse = "-"))))
    if (length(unique(s$counts)) > 1) {
      cat(sprintf("Runs cut to the shortest: %s Draws\n",
                  paste(names(s$counts), format_count(s$counts), sep = " ", collapse = ", ")))
    }
    cat_columns(data.frame(
      run = runs$run, draws = format_count(runs$draws),
      iterations = paste0(format_count(runs$first_iteration), "-", format_count(runs$last_iteration)),
      log_posterior = formatC(runs$log_posterior, format = "f", digits = 2),
      gap = formatC(runs$gap, format = "f", digits = 2)
    ), left = c("run", "iterations"))
    cat(sprintf("Log-posterior gap: %s (reported, not flagged)\n", formatC(max(runs$gap), format = "f", digits = 2)))
  },
  function(x) {
    s <- attr(x, "settings")
    q <- x$quantities
    cat(sprintf("Quantities: %s, %s flagged (stuck, R-hat > %s or ess_bulk < %s)\n", format_count(nrow(q)),
                format_count(sum(q$flagged)), format(s$rhat), format(s$ess)))
    groups <- convergence_groups
    by_group <- lapply(groups, function(g) q[q$group == g, ])
    worst <- function(v, f) {
      v <- v[!is.na(v)]
      if (length(v)) formatC(f(v), format = "g", digits = 3) else "-"
    }
    cat_columns(data.frame(
      group = groups,
      quantities = format_count(vapply(by_group, nrow, integer(1))),
      flagged = format_count(vapply(by_group, function(g) sum(g$flagged), integer(1))),
      stuck = format_count(vapply(by_group, function(g) sum(g$status == "stuck"), integer(1))),
      constant = format_count(vapply(by_group, function(g) sum(g$status == "constant"), integer(1))),
      max_rhat = vapply(by_group, function(g) worst(g$rhat, max), character(1)),
      min_ess_bulk = vapply(by_group, function(g) worst(g$ess_bulk, min), character(1))
    ), left = "group")
  },
  function(x) {
    s <- attr(x, "settings")
    inf <- x$infections
    cat(sprintf("Infections: %s, %s flagged (spread > %s)\n", format_count(nrow(inf)), format_count(sum(inf$flagged)),
                format(s$spread)))
    if (any(inf$flagged)) {
      cat_columns(data.frame(
        class = disagreement_classes,
        infections = format_count(vapply(disagreement_classes, function(k) sum(inf$class %in% k), integer(1)))
      ), left = "class")
    }
    worst <- inf[order(-inf$spread, inf$infection, method = "radix"), ][seq_len(min(10, nrow(inf))), ]
    cat(sprintf("%s by spread:\n", if (nrow(inf) > 10) "10 worst Infections" else "Infections"))
    cat_columns(data.frame(
      infection = worst$infection,
      spread = formatC(worst$spread, format = "f", digits = 3),
      class = ifelse(is.na(worst$class), "-", worst$class),
      worst_rhat = ifelse(is.na(worst$worst_rhat), "-", formatC(worst$worst_rhat, format = "g", digits = 3)),
      top_parent = vapply(worst$top_parent, function(t) paste(names(t), t, sep = ": ", collapse = ", "),
                          character(1))
    ), left = c("infection", "class", "top_parent"))
    pairs <- x$pairs
    cat(sprintf("Order pairs: %s, %s flagged (spread > %s)\n", format_count(nrow(pairs)),
                format_count(sum(pairs$flagged)), format(s$spread)))
  }
)

disagreement_classes <- c("same_top", "root_flip", "latent_vs_observed", "direction_flip", "linked_parents",
                          "different_parent")

# Print a data frame of strings as aligned columns: `left` columns left-aligned, the rest right-aligned.
cat_columns <- function(table, left = character()) {
  lines <- vapply(names(table), function(column) {
    values <- c(column, table[[column]])
    formatC(values, width = if (column %in% left) -max(nchar(values)) else max(nchar(values)))
  }, character(nrow(table) + 1))
  lines <- matrix(lines, ncol = ncol(table))
  cat(paste0("  ", apply(lines, 1, paste, collapse = "  "), "\n"), sep = "")
}
