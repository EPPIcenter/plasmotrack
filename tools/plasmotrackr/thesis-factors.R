#!/usr/bin/env Rscript
# Thesis factors from a thesis simulation directory path: the 14 factor columns of the fit summary, for
# plasmotrackr::pt_fit_summary(factors =). Not part of the package.
#
# Source it and call thesis_factors(paths), or run it: Rscript thesis-factors.R <path>... prints a CSV.
#
# A path names, in order (anything before or after is ignored, such as the Run directory, named chainN in the rerun
# layout):
#   <model>      full | partial | null (the reruns) or output_multichain | partially_obs_multichain (the thesis)
#   <scenario>   [simple_]total_nodes_<n>_founder_prop_<p>_super_infected_prop_<p>_r0_<r>_mean_moi_<m>_asymptomatic_prop_<p>
#   <network>    [simple_]tx_network_<panel>_div, or the bare panel (low | med | high)
#   <simulation> sim_<index>_s<seed>, optionally with a file suffix such as _to_run.json.gz or .tar.gz
# The model gives missingness_condition and fit_type: full -> full, posterior; partial -> missing_50, posterior;
# null -> full, null_prior. A simple_ prefix on the scenario or network makes the simulation model "simplified".
# For a dcifer row pass fit_type = "dcifer" to pt_fit_summary(): explicit arguments win over these columns.

thesis_factors <- function(paths) {
  rows <- lapply(paths, thesis_factors_one)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

thesis_factors_one <- function(path) {
  parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1]]
  scenario_pattern <- paste0(
    "^(simple_)?total_nodes_([0-9.]+)_founder_prop_([0-9.]+)_super_infected_prop_([0-9.]+)",
    "_r0_([0-9.]+)_mean_moi_([0-9.]+)_asymptomatic_prop_([0-9.]+)$"
  )
  at <- grep(scenario_pattern, parts)
  if (length(at) != 1 || at < 2 || at + 2 > length(parts)) {
    stop(sprintf("Not a thesis simulation path (<model>/<scenario>/<network>/sim_<i>_s<seed>): %s", path),
         call. = FALSE)
  }
  scenario <- regmatches(parts[[at]], regexec(scenario_pattern, parts[[at]]))[[1]]
  network <- regmatches(parts[[at + 1]],
                        regexec("^(simple_)?(?:tx_network_)?(low|med|high)(?:_div)?$", parts[[at + 1]], perl = TRUE))[[1]]
  simulation <- regmatches(parts[[at + 2]], regexec("^sim_([0-9]+)_s([0-9]+)(?:[_.].*)?$", parts[[at + 2]],
                                                     perl = TRUE))[[1]]
  if (length(network) == 0 || length(simulation) == 0) {
    stop(sprintf("Not a thesis simulation path (<model>/<scenario>/<network>/sim_<i>_s<seed>): %s", path),
         call. = FALSE)
  }
  model <- parts[[at - 1]]
  condition <- switch(model,
    full = , output_multichain = c("full", "posterior"),
    partial = , partially_obs_multichain = c("missing_50", "posterior"),
    null = c("full", "null_prior"),
    stop(sprintf("Unknown model directory '%s' (full, partial, null, output_multichain or partially_obs_multichain): %s",
                 model, path), call. = FALSE)
  )

  number <- function(i) as.numeric(scenario[[i]])
  out <- data.frame(
    run_key = NA_character_,
    simulation_model = if (nzchar(scenario[[2]]) || nzchar(network[[2]])) "simplified" else "biological",
    total_nodes = number(3),
    founder_probability = number(4),
    superinfection_probability = number(5),
    secondary_case_rate = number(6),
    founder_coi_rate = number(7),
    asymptomatic_probability = number(8),
    panel_name = network[[3]],
    panel_alleles_per_locus = c(low = 5, med = 10, high = 20)[[network[[3]]]],
    simulation_index = as.integer(simulation[[2]]),
    seed = as.integer(simulation[[3]]),
    missingness_condition = condition[[1]],
    fit_type = condition[[2]],
    stringsAsFactors = FALSE
  )
  # The thesis run key: numbers as R prints them (0.05, 1, 1.5).
  out$run_key <- paste(out$simulation_model, paste0("fp", out$founder_probability),
                       paste0("sp", out$superinfection_probability), paste0("r0", out$secondary_case_rate),
                       paste0("coi", out$founder_coi_rate), out$panel_name, paste0("sim", out$simulation_index),
                       paste0("s", out$seed), sep = "__")
  out
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    stop("Usage: Rscript thesis-factors.R <path>...", call. = FALSE)
  }
  utils::write.csv(thesis_factors(args), stdout(), row.names = FALSE)
}
