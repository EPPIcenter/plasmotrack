# Oracle for flows between strata (#127), apart from helper-oracle.R so parallel tickets don't collide. Base R loops
# over the spec (doc/r-package.md, "Strata and flows"), per Draw. `strata` is oracle_strata()'s vector (one value per
# oracle_infection_ids(), NA kept); its levels are the sorted distinct values, NA last.
#
#   oracle_attribution(o, strata, runs, phase)
#   oracle_enrichment(o, strata, candidates, runs, phase)   candidates: data frame of parent, child.
#   oracle_net_flow(o, strata, runs, phase)
#   oracle_reassign_long_lags(o, observation_time, max_lag, strata = NULL, runs, phase)
#                                   observation_time: named vector (as golden_observation_times()).

oracle_flow_levels <- function(strata) sort(unique(strata), na.last = TRUE, method = "radix")

# Stratum of an Infection id, as a level index.
oracle_flow_index <- function(o, strata) {
  levels <- oracle_flow_levels(strata)
  index <- vapply(strata, function(s) which(vapply(levels, identical, logical(1), s)), integer(1))
  setNames(index, oracle_infection_ids(o))
}

oracle_attribution <- function(o, strata, runs = NULL, phase = "sampling") {
  levels <- oracle_flow_levels(strata)
  of <- oracle_flow_index(o, strata)
  edges <- oracle_edge_draws(o, runs, phase)
  infection <- oracle_infection_draws(o, runs, phase)
  out <- list()
  for (run in oracle_runs(o, runs)) {
    for (it in oracle_draws(o, run, phase)) {
      for (r in seq_along(levels)) {
        mass <- numeric(length(levels) + 1)
        for (i in which(edges$run == run & edges$iteration == it)) {
          if (of[[edges$child[[i]]]] == r) {
            s <- of[[edges$parent[[i]]]]
            mass[[s]] <- mass[[s]] + edges$prob[[i]]
          }
        }
        for (i in which(infection$run == run & infection$iteration == it)) {
          if (of[[infection$infection[[i]]]] == r) {
            mass[[length(mass)]] <- mass[[length(mass)]] + infection$p_latent_parent[[i]]
          }
        }
        for (s in seq_along(mass)) {
          row <- data.frame(run = run, iteration = it)
          row$recipient <- levels[[r]]
          row$source <- if (s <= length(levels)) levels[[s]] else levels[NA_integer_]
          row$latent <- s > length(levels)
          row$share <- if (sum(mass) > 0) mass[[s]] / sum(mass) else NA_real_
          out[[length(out) + 1]] <- row
        }
      }
    }
  }
  oracle_bind(out, NULL)
}

oracle_enrichment <- function(o, strata, candidates, runs = NULL, phase = "sampling") {
  levels <- oracle_flow_levels(strata)
  of <- oracle_flow_index(o, strata)
  k <- length(levels)
  count <- matrix(0, k, k)
  candidates <- unique(data.frame(parent = as.character(candidates$parent), child = as.character(candidates$child)))
  for (i in seq_len(nrow(candidates))) {
    p <- candidates$parent[[i]]
    c <- candidates$child[[i]]
    if (p %in% names(of) && c %in% names(of)) count[of[[p]], of[[c]]] <- count[of[[p]], of[[c]]] + 1
  }
  edges <- oracle_edge_draws(o, runs, phase)
  out <- list()
  for (run in oracle_runs(o, runs)) {
    for (it in oracle_draws(o, run, phase)) {
      mass <- matrix(0, k, k)
      for (i in which(edges$run == run & edges$iteration == it)) {
        s <- of[[edges$parent[[i]]]]
        r <- of[[edges$child[[i]]]]
        mass[s, r] <- mass[s, r] + edges$prob[[i]]
      }
      for (s in seq_len(k)) {
        for (r in seq_len(k)) {
          row <- data.frame(run = run, iteration = it)
          row$source <- levels[[s]]
          row$recipient <- levels[[r]]
          row$enrichment <- if (count[s, r] > 0 && sum(mass) > 0) {
            (mass[s, r] / sum(mass)) / (count[s, r] / sum(count))
          } else {
            NA_real_
          }
          out[[length(out) + 1]] <- row
        }
      }
    }
  }
  oracle_bind(out, NULL)
}

oracle_net_flow <- function(o, strata, runs = NULL, phase = "sampling") {
  levels <- oracle_flow_levels(strata)
  of <- oracle_flow_index(o, strata)
  edges <- oracle_edge_draws(o, runs, phase)
  out <- list()
  for (run in oracle_runs(o, runs)) {
    for (it in oracle_draws(o, run, phase)) {
      for (s in seq_along(levels)) {
        out_edges <- 0
        in_edges <- 0
        for (i in which(edges$run == run & edges$iteration == it)) {
          p <- of[[edges$parent[[i]]]]
          c <- of[[edges$child[[i]]]]
          if (p == c) next
          if (p == s) out_edges <- out_edges + edges$prob[[i]]
          if (c == s) in_edges <- in_edges + edges$prob[[i]]
        }
        row <- data.frame(run = run, iteration = it)
        row$stratum <- levels[[s]]
        row$net_flow <- out_edges - in_edges
        out[[length(out) + 1]] <- row
      }
    }
  }
  oracle_bind(out, NULL)
}

oracle_reassign_long_lags <- function(o, observation_time, max_lag, strata = NULL, runs = NULL, phase = "sampling") {
  ids <- oracle_infection_ids(o)
  edges <- oracle_edge_draws(o, runs, phase)
  infection <- oracle_infection_draws(o, runs, phase)
  lag <- unname(observation_time[edges$child] - observation_time[edges$parent])
  moved <- !is.na(lag) & lag > max_lag
  groups <- if (is.null(strata)) list(ids) else lapply(oracle_flow_levels(strata), function(s) {
    ids[vapply(strata, identical, logical(1), s)]
  })
  out <- list()
  for (run in oracle_runs(o, runs)) {
    for (it in oracle_draws(o, run, phase)) {
      out_degree <- setNames(numeric(length(ids)), ids)
      q <- setNames(numeric(length(ids)), ids)
      for (i in which(infection$run == run & infection$iteration == it)) {
        q[[infection$infection[[i]]]] <- infection$p_latent_parent[[i]]
      }
      for (i in which(edges$run == run & edges$iteration == it)) {
        if (moved[[i]]) {
          q[[edges$child[[i]]]] <- q[[edges$child[[i]]]] + edges$prob[[i]]
        } else {
          out_degree[[edges$parent[[i]]]] <- out_degree[[edges$parent[[i]]]] + edges$prob[[i]]
        }
      }
      q <- pmin(q, 1)
      for (g in seq_along(groups)) {
        row <- data.frame(run = run, iteration = it)
        if (!is.null(strata)) row$stratum <- oracle_flow_levels(strata)[[g]]
        row$mean_out_degree <- sum(out_degree[groups[[g]]]) / length(groups[[g]])
        row$mean_p_latent_parent <- sum(q[groups[[g]]]) / length(groups[[g]])
        out[[length(out) + 1]] <- row
      }
    }
  }
  oracle_bind(out, NULL)
}
