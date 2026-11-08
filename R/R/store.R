# The Store format, as plasmotrackr needs it: tables, required columns, and the reader rule.
# doc/store-layout.md in the plasmotrack repository is the reference.

# Windowed Store tables in the layout doc's order. `required` tables must have files in a Run with closed Windows;
# `columns` are the columns the package needs (a Store may have more). `phase` is FALSE for the one table read over
# every closed Window whatever the Phase asked for.
store_tables <- list(
  run_scalars = list(required = TRUE, phase = TRUE, columns = c(
    "iteration", "phase", "window", "mean_coi", "mean_strains_tx", "parent_set_size_prob", "prior", "likelihood",
    "observation", "transmission"
  )),
  infections = list(required = TRUE, phase = TRUE, columns = c(
    "iteration", "phase", "window", "infection", "expected_false_positives", "expected_false_negatives",
    "infection_duration"
  )),
  ordering = list(required = TRUE, phase = TRUE, columns = c("iteration", "phase", "window", "infection", "rank")),
  allele_frequencies = list(required = TRUE, phase = TRUE, columns = c(
    "iteration", "phase", "window", "locus", "allele", "frequency"
  )),
  latent_genotypes = list(required = TRUE, phase = TRUE, columns = c(
    "iteration", "phase", "window", "infection", "locus", "alleles"
  )),
  latent_parent_genotypes = list(required = FALSE, phase = TRUE, columns = c(
    "iteration", "phase", "window", "latent_parent", "locus", "alleles"
  )),
  parent_sets = list(required = TRUE, phase = TRUE, columns = c(
    "iteration", "phase", "window", "infection", "set_id", "prob"
  )),
  parent_set_members = list(required = TRUE, phase = FALSE, columns = c("infection", "set_id", "member")),
  chains = list(required = FALSE, phase = TRUE, columns = c(
    "iteration", "phase", "window", "temperature_slot", "replica", "temperature", "prior", "likelihood"
  )),
  acceptance = list(required = TRUE, phase = TRUE, columns = c(
    "window", "phase", "sampler", "accepted", "proposed"
  )),
  swaps = list(required = FALSE, phase = TRUE, columns = c(
    "window", "phase", "temperature_pair", "accepted", "attempted"
  )),
  events = list(required = TRUE, phase = TRUE, columns = c(
    "iteration", "phase", "window", "wall_time", "type", "payload"
  ))
)

# DuckDB types of the columns above, for the empty table of a name no Run has rows for.
column_types <- c(
  iteration = "BIGINT", window = "BIGINT", set_id = "BIGINT", allele = "BIGINT", rank = "BIGINT",
  temperature_slot = "BIGINT", replica = "BIGINT", accepted = "BIGINT", proposed = "BIGINT",
  temperature_pair = "BIGINT", attempted = "BIGINT",
  mean_coi = "DOUBLE", mean_strains_tx = "DOUBLE", parent_set_size_prob = "DOUBLE", prior = "DOUBLE",
  likelihood = "DOUBLE", observation = "DOUBLE", transmission = "DOUBLE", expected_false_positives = "DOUBLE",
  expected_false_negatives = "DOUBLE", infection_duration = "DOUBLE", frequency = "DOUBLE", prob = "DOUBLE",
  temperature = "DOUBLE"
)

# Run scalars only some Runs have: the Transmission leak with a Leak prior, and the rates' population means and shapes
# under the hierarchical Error model. Runs pooled in one Analysis have the same ones (check_comparable()).
model_choice_scalars <- c("transmission_leak", "efp_mean", "efp_shape", "efn_mean", "efn_shape")

structure_tables <- list(
  infections = c("infection", "latent_parent", "symptomatic"),
  loci = c("locus", "allele"),
  samplers = c("sampler", "type", "target")
)

# The closed Windows of a Run, from its `_closed/` markers: a tibble of window, phase, first_iteration and
# last_iteration, ordered by window. The reader rule starts here: a Window without a marker does not exist.
closed_windows <- function(path) {
  markers <- list.files(file.path(path, "_closed"), pattern = "^window-\\d+$", full.names = TRUE)
  rows <- lapply(markers, function(f) {
    m <- jsonlite::read_json(f)
    tibble::tibble(
      window = as.numeric(m$window), phase = as.character(m$phase),
      first_iteration = as.numeric(m$first_iteration), last_iteration = as.numeric(m$last_iteration)
    )
  })
  windows <- dplyr::bind_rows(
    tibble::tibble(window = numeric(), phase = character(), first_iteration = numeric(), last_iteration = numeric()),
    rows
  )
  windows[order(windows$window), ]
}

# Files of one table in a Run's closed Windows of the given Phases: built from the markers, never a glob, since a
# glob also reads `.in-flight.parquet` and the files of Windows without a marker.
window_files <- function(run, table, phase) {
  windows <- run$windows[[1]]
  if (store_tables[[table]]$phase) windows <- windows[windows$phase %in% phase, ]
  files <- file.path(run$path, table, sprintf("window-%06.0f.parquet", windows$window))
  files[file.exists(files)]
}

read_manifest <- function(path) {
  jsonlite::read_json(file.path(path, "manifest.json"))
}
