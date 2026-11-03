# ---- Convergence tree diagnostics (#125) ----
# Oracle for pt_convergence()$infections and $pairs, in loops over the raw Parent-set and Ordering rows.
# `draws` from oracle_convergence_draws(), `m` from oracle_convergence_matrices(), `q` from
# oracle_convergence_quantities(m). Gotcha: inside tibble(), `spread` is the column just made, not the argument.

# Run means over the Draw set: list run -> child -> named vector parent -> mean, with "latent" for the Latent parent.
# Every edge seen in any Run is present (0) in every Run.
oracle_tree_means <- function(o, draws) {
  runs <- names(draws)
  n <- length(draws[[1]])
  edges <- oracle_edge_draws(o, runs)
  infection <- oracle_infection_draws(o, runs)
  ids <- oracle_infection_ids(o)
  out <- list()
  for (run in runs) {
    out[[run]] <- list()
    for (id in ids) out[[run]][[id]] <- c(latent = 0)
  }
  for (k in seq_len(nrow(edges))) {
    for (run in runs) {
      if (is.na(out[[run]][[edges$child[[k]]]][edges$parent[[k]]])) out[[run]][[edges$child[[k]]]][edges$parent[[k]]] <- 0
    }
  }
  for (k in seq_len(nrow(edges))) {
    run <- edges$run[[k]]
    if (!edges$iteration[[k]] %in% draws[[run]]) next
    out[[run]][[edges$child[[k]]]][edges$parent[[k]]] <- out[[run]][[edges$child[[k]]]][edges$parent[[k]]] +
      edges$prob[[k]] / n
  }
  for (k in seq_len(nrow(infection))) {
    run <- infection$run[[k]]
    if (!infection$iteration[[k]] %in% draws[[run]]) next
    id <- infection$infection[[k]]
    out[[run]][[id]]["latent"] <- out[[run]][[id]]["latent"] + infection$p_latent_parent[[k]] / n
  }
  out
}

oracle_tree_top <- function(v) {
  keys <- names(v)
  observed <- keys[keys != "latent"]
  best <- max(v)
  cands <- sort(observed[abs(v[observed] - best) <= 1e-9], method = "radix")
  if (length(cands)) cands[[1]] else "latent"
}

oracle_tree_class <- function(child, tops, means, t) {
  get <- function(run, c, p) {
    v <- means[[run]][[c]][p]
    if (is.null(v) || is.na(v)) 0 else unname(v)
  }
  runs <- names(tops)
  if (length(unique(tops)) == 1) return("same_top")
  observed <- unique(tops[tops != "latent"])
  if (any(tops == "latent") && length(observed)) {
    for (run in runs[tops == "latent"]) for (p in observed) if (get(run, p, child) > t) return("root_flip")
    return("latent_vs_observed")
  }
  for (run in runs) for (p in observed) if (p != tops[[run]] && get(run, p, child) > t) return("direction_flip")
  for (run in runs) for (x in observed) for (y in observed) if (x != y && get(run, x, y) > t) return("linked_parents")
  "different_parent"
}

oracle_convergence_infections <- function(o, draws, q, spread = 0.2) {
  means <- oracle_tree_means(o, draws)
  runs <- names(draws)
  ids <- sort(oracle_infection_ids(o), method = "radix")
  t <- spread
  rows <- list()
  for (id in ids) {
    keys <- names(means[[runs[[1]]]][[id]])
    s <- 0
    for (key in keys) {
      v <- vapply(runs, function(r) unname(means[[r]][[id]][key]), numeric(1))
      s <- max(s, max(v) - min(v))
    }
    tops <- vapply(runs, function(r) oracle_tree_top(means[[r]][[id]]), character(1))
    own <- c(sprintf("%s[%s]", c("expected_false_positives", "expected_false_negatives", "infection_duration",
                                  "p_latent_parent"), id))
    mine <- q$name %in% own | startsWith(q$name, sprintf("edge[%s,", id))
    r <- q$rhat[mine & !is.na(q$rhat)]
    rows[[id]] <- tibble::tibble(
      infection = id, spread = s, class = if (s > t) oracle_tree_class(id, tops, means, t) else NA_character_,
      top_parent = list(tops), worst_rhat = if (length(r)) max(r) else NA_real_, flagged = s > t
    )
  }
  dplyr::bind_rows(rows)
}

oracle_convergence_pairs <- function(m, spread = 0.2) {
  t <- spread
  rows <- list()
  for (name in m$groups$order_pair) {
    ends <- strsplit(sub("^order_pair\\[(.*)\\]$", "\\1", name), ",", fixed = TRUE)[[1]]
    p <- colMeans(m$matrices[[name]])
    rows[[name]] <- tibble::tibble(i = ends[[1]], j = ends[[2]], p_before = list(p), spread = max(p) - min(p),
                                   flagged = max(p) - min(p) > t)
  }
  dplyr::bind_rows(rows)
}
