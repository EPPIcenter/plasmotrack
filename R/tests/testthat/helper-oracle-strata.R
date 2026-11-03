# Oracle for strata (#126), apart from helper-oracle.R so parallel tickets don't collide. Base R loops over the spec:
# a stratum S is the Infections whose epi value (column `by`) is S; Infections missing from epi or with NA form the
# NA stratum. n_S = |S|.
#
#   oracle_strata(o, epi, by)   stratum of every Infection (a vector in oracle_infection_ids() order, NA kept).
#   oracle_stratified_mean(o, rows, columns, fill, strata, runs, phase)   per Draw and stratum, the mean over the
#                               stratum's Infections (sources or recipients: `rows$infection` is whichever counts).
#   oracle_stratified_composition(), oracle_stratified_map_out_degree(), oracle_stratified_out_degree_range().

oracle_strata <- function(o, epi, by) {
  value <- epi[[by]]
  if (is.factor(value)) value <- as.character(value)
  ids <- oracle_infection_ids(o)
  out <- value[0][rep(NA_integer_, length(ids))]
  for (k in seq_along(ids)) {
    hit <- which(as.character(epi$infection) == ids[[k]])
    if (length(hit)) out[[k]] <- value[[hit]]
  }
  out
}

oracle_stratum_members <- function(o, strata, s) oracle_infection_ids(o)[strata %in% s]

oracle_stratified_mean <- function(o, rows, columns, fill, strata, runs = NULL, phase = "sampling") {
  out <- list()
  for (run in oracle_runs(o, runs)) {
    for (it in oracle_draws(o, run, phase)) {
      for (s in unique(strata)) {
        ids <- oracle_stratum_members(o, strata, s)
        row <- data.frame(run = run, iteration = it)
        row$stratum <- s
        for (name in names(columns)) {
          total <- 0
          for (inf in ids) {
            hit <- which(rows$run == run & rows$iteration == it & rows$infection == inf)
            total <- total + if (length(hit)) rows[[columns[[name]]]][[hit]] else fill[[name]]
          }
          row[[name]] <- total / length(ids)
        }
        out[[length(out) + 1]] <- row
      }
    }
  }
  oracle_bind(out, NULL)
}

oracle_stratified_composition <- function(o, strata, runs = NULL, phase = "sampling") {
  shares <- oracle_class_draws(o, runs, phase)
  out <- list()
  for (run in oracle_runs(o, runs)) {
    for (it in oracle_draws(o, run, phase)) {
      for (s in unique(strata)) {
        ids <- oracle_stratum_members(o, strata, s)
        for (class in unique(shares$class)) {
          total <- 0
          for (i in which(shares$run == run & shares$iteration == it & shares$class == class &
                            shares$infection %in% ids)) {
            total <- total + shares$share[[i]]
          }
          row <- data.frame(run = run, iteration = it)
          row$stratum <- s
          row$class <- class
          row$share <- total / length(ids)
          out[[length(out) + 1]] <- row
        }
      }
    }
  }
  oracle_bind(out, NULL)
}

# ψ_MAP(S): MAP members (over every recipient) that lie in S, over n_S.
oracle_stratified_map_out_degree <- function(o, strata, by_run = FALSE, runs = NULL) {
  sets <- oracle_parent_sets_pooled(o, by_run = by_run, runs = runs)
  if (!by_run) sets$run <- "pooled"
  out <- list()
  for (run in unique(sets$run)) {
    for (s in unique(strata)) {
      ids <- oracle_stratum_members(o, strata, s)
      count <- 0
      for (i in which(sets$run == run & sets$rank == 1)) count <- count + sum(sets$members[[i]] %in% ids)
      row <- data.frame(run = run)
      row$stratum <- s
      row$map_out_degree <- count / length(ids)
      row$infections <- length(ids)
      out[[length(out) + 1]] <- row
    }
  }
  out <- oracle_bind(out, NULL)
  if (!by_run) out$run <- NULL
  out
}

# ψ_lo(S), ψ_hi(S): per recipient j the min and max over its Credible collection of |U_obs ∩ S|, summed, over n_S.
oracle_stratified_out_degree_range <- function(o, strata, level = 0.8, by_run = FALSE, runs = NULL) {
  sets <- oracle_parent_sets_pooled(o, level = level, by_run = by_run, runs = runs)
  if (!by_run) sets$run <- "pooled"
  out <- list()
  for (run in unique(sets$run)) {
    for (s in unique(strata)) {
      ids <- oracle_stratum_members(o, strata, s)
      lo <- 0
      hi <- 0
      for (child in oracle_infection_ids(o)) {
        sizes <- numeric()
        for (i in which(sets$run == run & sets$infection == child & sets$credible)) {
          sizes <- c(sizes, sum(sets$members[[i]] %in% ids))
        }
        if (length(sizes)) {
          lo <- lo + min(sizes)
          hi <- hi + max(sizes)
        }
      }
      row <- data.frame(run = run)
      row$stratum <- s
      row$lo <- lo / length(ids)
      row$hi <- hi / length(ids)
      row$level <- level
      out[[length(out) + 1]] <- row
    }
  }
  out <- oracle_bind(out, NULL)
  if (!by_run) out$run <- NULL
  out
}
