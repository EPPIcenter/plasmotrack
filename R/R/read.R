#' Open an Analysis
#'
#' An Analysis is the Runs of one input. `pt_read()` opens their Stores, checks that they can be pooled, and returns
#' a `pt_analysis` that every other `pt_` function takes.
#'
#' `pt_read()` reads only the Window markers, the manifests, `structure/` and the Draw count of each Run; it reads no
#' Traces. Tables are read later, lazily, from the Windows that are closed, so a Run that is still running, or was
#' killed, reads up to its last closed Window.
#'
#' Each Analysis has its own DuckDB connection, closed by [pt_close()] or when the Analysis is garbage-collected. The
#' memory cap bounds each Analysis, not the R session: two open Analyses may use twice the cap. Set the engine once
#' with R options, and override them for one Analysis with the arguments:
#'
#' * `plasmotrackr.memory_limit`: DuckDB's memory cap, such as `"8GB"`. Default: DuckDB's own (80% of RAM).
#' * `plasmotrackr.temp_directory`: where DuckDB spills when it reaches the cap. Default: [tempdir()].
#' * `plasmotrackr.threads`: DuckDB threads. Default: DuckDB's own (every core).
#'
#' @param paths Run directories (a Store each, holding `manifest.json`), or one parent directory whose
#'   subdirectories are Runs. Names, when given, become the Run labels; otherwise each Run is labelled by its
#'   directory name.
#' @param epi Optional data frame of epidemiological data with an `infection` column of Infection ids, one row per
#'   Infection.
#' @param memory_limit,temp_directory,threads DuckDB settings for this Analysis. `NULL` uses the R option, then the
#'   default.
#' @return A `pt_analysis`: an S3 list holding the DuckDB connection, the Run table, the Infections, the Loci and
#'   `epi`.
#' @seealso [pt_runs()], [pt_table()], [pt_close()]
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", package = "plasmotrackr")
#' a <- pt_read(file.path(golden, c(one = "single-chain", four = "replica-exchange")), threads = 1)
#' a
#' pt_runs(a)
#' pt_close(a)
pt_read <- function(paths, epi = NULL, memory_limit = NULL, temp_directory = NULL, threads = NULL) {
  runs <- resolve_runs(paths)
  runs$manifest <- lapply(runs$path, read_manifest)
  check_comparable(runs)
  runs$windows <- lapply(runs$path, closed_windows)
  epi <- check_epi(epi)

  engine <- open_engine(
    memory_limit = engine_option(memory_limit, "plasmotrackr.memory_limit"),
    temp_directory = engine_option(temp_directory, "plasmotrackr.temp_directory", tempdir()),
    threads = engine_option(threads, "plasmotrackr.threads")
  )
  ok <- FALSE
  on.exit(if (!ok) close_engine(engine))
  con <- engine$con

  check_store_tables(con, runs)
  structure <- read_structure(con, runs$path[[1]])
  runs <- count_draws(con, runs)
  runs <- summarise_events(con, runs)
  ok <- TRUE
  warn_runs(runs)

  structure(
    list(engine = engine, runs = runs, infections = structure$infections, loci = structure$loci, epi = epi),
    class = "pt_analysis"
  )
}

# Run directories and labels from pt_read()'s `paths`.
resolve_runs <- function(paths, call = rlang::caller_env()) {
  if (!is.character(paths) || length(paths) == 0 || anyNA(paths)) {
    cli::cli_abort("{.arg paths} must be Run directories, or one directory of Runs.", call = call)
  }
  labels <- names(paths)
  paths <- path.expand(unname(paths))
  is_run <- file.exists(file.path(paths, "manifest.json"))
  if (length(paths) == 1 && is.null(labels) && dir.exists(paths) && !is_run) {
    subdirs <- sort(list.dirs(paths, recursive = FALSE))
    subdirs <- subdirs[file.exists(file.path(subdirs, "manifest.json"))]
    if (length(subdirs) == 0) {
      cli::cli_abort("{.path {paths}} is not a Run directory and holds no Run directories.", call = call)
    }
    paths <- subdirs
    is_run <- rep(TRUE, length(paths))
  }
  if (!all(is_run)) {
    cli::cli_abort(c(
      "Not a Run directory (no {.file manifest.json}):",
      rlang::set_names(cli_escape(paths[!is_run]), "x")
    ), call = call)
  }
  paths <- normalizePath(paths, winslash = "/")
  if (is.null(labels)) labels <- rep("", length(paths))
  labels[is.na(labels) | labels == ""] <- basename(paths)[is.na(labels) | labels == ""]
  if (anyDuplicated(labels)) {
    cli::cli_abort(c(
      "Run labels must be unique: {.val {unique(labels[duplicated(labels)])}} repeats.",
      i = "Name the paths to label the Runs, as in {.code pt_read(c(a = \"run1\", b = \"run2\"))}."
    ), call = call)
  }
  tibble::tibble(run = labels, label = labels, path = paths)
}

# Runs of one Analysis must share their input and model.
check_comparable <- function(runs, call = rlang::caller_env()) {
  differs <- function(field, value, describe) {
    values <- vapply(runs$manifest, function(m) jsonlite::toJSON(value(m), auto_unbox = TRUE), character(1))
    if (length(unique(values)) > 1) {
      groups <- split(runs$run, factor(values, levels = unique(values)))
      first <- match(unique(values), values)
      lines <- paste0(describe(runs$manifest[first]), ": ", vapply(groups, paste, character(1), collapse = ", "))
      cli::cli_abort(c(
        "Runs differ in {.field {field}}; they cannot be pooled in one Analysis.",
        rlang::set_names(cli_escape(lines), "*")
      ), call = call)
    }
  }
  differs("input_sha256", function(m) m$input_sha256[sort(names(m$input_sha256))],
          function(ms) vapply(ms, function(m) paste0("input ", substr(m$input_sha256$input, 1, 12)), character(1)))
  differs("null_model", function(m) isTRUE(m$null_model),
          function(ms) vapply(ms, function(m) paste0("null_model ", tolower(isTRUE(m$null_model))), character(1)))
}

check_epi <- function(epi, call = rlang::caller_env()) {
  if (is.null(epi)) return(NULL)
  if (!is.data.frame(epi) || !"infection" %in% names(epi)) {
    cli::cli_abort("{.arg epi} must be a data frame with an {.field infection} column.", call = call)
  }
  epi <- tibble::as_tibble(epi)
  epi$infection <- as.character(epi$infection)
  if (anyNA(epi$infection) || anyDuplicated(epi$infection)) {
    cli::cli_abort("{.arg epi} must have one row per Infection: {.field infection} has NA or repeated ids.",
                   call = call)
  }
  epi
}

# Every required table and column, in every Run. Reads Parquet metadata only.
check_store_tables <- function(con, runs, call = rlang::caller_env()) {
  missing <- character()
  for (i in seq_len(nrow(runs))) {
    run <- runs[i, ]
    for (name in names(structure_tables)) {
      file <- file.path(run$path, "structure", paste0(name, ".parquet"))
      if (!file.exists(file)) {
        missing <- c(missing, sprintf("Run %s: table structure/%s", run$run, name))
        next
      }
      have <- query(con, sprintf("SELECT name FROM parquet_schema(%s)", sql_string(con, file)))$name
      lost <- setdiff(structure_tables[[name]], have)
      if (length(lost)) {
        missing <- c(missing, sprintf("Run %s: structure/%s column%s %s", run$run, name,
                                      if (length(lost) > 1) "s" else "", paste(lost, collapse = ", ")))
      }
    }
    if (nrow(run$windows[[1]]) == 0) next
    for (name in names(store_tables)) {
      files <- window_files(run, name, c("burn_in", "sampling"))
      if (length(files) == 0) {
        if (store_tables[[name]]$required) missing <- c(missing, sprintf("Run %s: table %s", run$run, name))
        next
      }
      have <- query(con, sprintf("SELECT name FROM parquet_schema(%s)", sql_string(con, files[[1]])))$name
      lost <- setdiff(store_tables[[name]]$columns, have)
      if (length(lost)) {
        missing <- c(missing, sprintf("Run %s: %s column%s %s", run$run, name,
                                      if (length(lost) > 1) "s" else "", paste(lost, collapse = ", ")))
      }
    }
  }
  if (length(missing)) {
    cli::cli_abort(c(
      "Stores are missing tables or columns plasmotrackr needs:",
      rlang::set_names(cli_escape(missing), "x")
    ), call = call)
  }
}

read_structure <- function(con, path) {
  read <- function(name) {
    query(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(con, file.path(path, "structure",
                                                                                   paste0(name, ".parquet")))))
  }
  infections <- read("infections")
  infections <- tibble::tibble(
    infection = as.character(infections$infection),
    observation_time = NA_real_,
    symptomatic = as.logical(infections$symptomatic),
    latent_parent = as.character(infections$latent_parent)
  )
  infections <- infections[order(infections$infection, method = "radix"), ]
  loci <- read("loci")
  loci <- tibble::tibble(locus = as.character(loci$locus), allele = loci$allele)
  loci <- dplyr::count(loci, .data$locus, name = "n_alleles")
  loci$n_alleles <- as.integer(loci$n_alleles)
  loci <- loci[order(loci$locus, method = "radix"), ]
  list(infections = infections, loci = loci)
}

# D_r and the first and last sampling iteration of each Run, from run_scalars in closed sampling Windows.
count_draws <- function(con, runs) {
  counts <- lapply(seq_len(nrow(runs)), function(i) {
    files <- window_files(runs[i, ], "run_scalars", "sampling")
    if (length(files) == 0) return(tibble::tibble(draws = 0, first_iteration = NA_real_, last_iteration = NA_real_))
    query(con, sprintf(paste(
      "SELECT count(DISTINCT iteration) AS draws, min(iteration) AS first_iteration,",
      "max(iteration) AS last_iteration FROM read_parquet(%s) WHERE phase = 'sampling'"
    ), sql_list(con, files)))
  })
  counts <- dplyr::bind_rows(counts)
  runs$draws <- as.integer(counts$draws)
  runs$first_iteration <- as.numeric(counts$first_iteration)
  runs$last_iteration <- as.numeric(counts$last_iteration)
  runs
}

# Status, resumes, sampling time and non-finite values of each Run, from its events in closed Windows of both Phases.
summarise_events <- function(con, runs) {
  summaries <- lapply(seq_len(nrow(runs)), function(i) {
    run <- runs[i, ]
    files <- window_files(run, "events", c("burn_in", "sampling"))
    events <- if (length(files)) {
      query(con, sprintf(paste(
        "SELECT \"window\", type, payload FROM read_parquet(%s, file_row_number = true)",
        "ORDER BY \"window\", file_row_number"
      ), sql_list(con, files)))
    } else {
      tibble::tibble(window = numeric(), type = character(), payload = character())
    }
    payload <- function(rows) lapply(events$payload[rows], jsonlite::parse_json)

    # The last start, resume or end decides: a Run whose events stop after a start or resume was killed or is
    # still running. window_closed rows follow run_ended in its Window, so they are not looked at.
    lifecycle <- which(events$type %in% c("run_started", "run_resumed", "run_ended"))
    status <- "incomplete"
    if (length(lifecycle) && events$type[[max(lifecycle)]] == "run_ended") {
      reason <- payload(max(lifecycle))[[1]]$reason
      if (identical(reason, "completed")) status <- "finished"
      if (identical(reason, "interrupted")) status <- "interrupted"
    }

    sampling <- run$windows[[1]]$window[run$windows[[1]]$phase == "sampling"]
    closed <- payload(which(events$type == "window_closed" & events$window %in% sampling))
    seconds <- sum(vapply(closed, function(p) as.numeric(p$seconds), numeric(1)))
    iterations <- sum(vapply(closed, function(p) as.numeric(p$sampling_iterations), numeric(1)))

    tibble::tibble(
      status = status,
      resumes = sum(events$type == "run_resumed"),
      sampling_seconds = seconds,
      iterations_per_second = if (seconds > 0) iterations / seconds else NA_real_,
      non_finite = sum(events$type == "non_finite_value")
    )
  })
  summaries <- dplyr::bind_rows(summaries)
  runs$status <- summaries$status
  runs$resumes <- as.integer(summaries$resumes)
  runs$started_from <- vapply(runs$manifest, function(m) {
    id <- m$started_from$run_id
    if (is.null(id)) NA_character_ else as.character(id)
  }, character(1))
  runs$sampling_seconds <- summaries$sampling_seconds
  runs$iterations_per_second <- summaries$iterations_per_second
  runs$non_finite <- as.integer(summaries$non_finite)
  runs
}

# pt_read()'s warnings: Runs that are not finished, Draw counts that differ, and non-finite values.
warn_runs <- function(runs, call = rlang::caller_env()) {
  unfinished <- runs$status != "finished"
  if (any(unfinished)) {
    lines <- sprintf("%s: %s", runs$run[unfinished], runs$status[unfinished])
    cli::cli_warn(c(
      "{sum(unfinished)} Run{?s} {?is/are} not finished; {?it is/they are} read up to the last closed Window.",
      rlang::set_names(cli_escape(lines), "*")
    ), call = call)
  }
  if (length(unique(runs$draws)) > 1) {
    lines <- sprintf("%s: %s", runs$run, format_count(runs$draws))
    cli::cli_warn(c(
      "Draw counts differ between Runs; pooled values weigh each Run by its Draws.",
      rlang::set_names(cli_escape(lines), "*")
    ), call = call)
  }
  bad <- runs$non_finite > 0
  if (any(bad)) {
    lines <- sprintf("%s: %d non_finite_value event%s", runs$run[bad], runs$non_finite[bad],
                     ifelse(runs$non_finite[bad] == 1, "", "s"))
    cli::cli_warn(c(
      "{sum(bad)} Run{?s} hit a non-finite prior or likelihood; see {.code pt_table(a, \"events\")}.",
      rlang::set_names(cli_escape(lines), "*")
    ), call = call)
  }
}

#' The Runs of an Analysis
#'
#' @param a A `pt_analysis` from [pt_read()].
#' @return A tibble with one row per Run:
#'   * `run`: the Run's label, the value of the `run` column in every table.
#'   * `label`: the Run's label, from the names given to [pt_read()] or the directory name.
#'   * `path`: the Run directory.
#'   * `draws`: \eqn{D_r}, the number of sampling Draws in closed Windows. Every per-Run mean divides by it.
#'   * `first_iteration`, `last_iteration`: the first and last sampling Draw's iteration, `NA` with no Draws.
#'   * `windows`: closed Windows, both Phases.
#'   * `status`: from the Run's events. `"finished"` when the last start, resume or end of the Run is a `run_ended`
#'     that `completed`; `"interrupted"` when it is a `run_ended` after an interrupt (the Run can be resumed);
#'     otherwise `"incomplete"`: the Run is still running, or was killed, and is read up to its last closed Window.
#'   * `resumes`: how many times the Run was resumed (`--resume`), the count of `run_resumed` events. A resumed Run
#'     is one Run with continuous Draws.
#'   * `started_from`: for a Run started with `--start-from`, the `run_id` of the Run whose checkpoint it started
#'     from; otherwise `NA`. The two Runs stay separate: `pt_read()` never follows this link.
#'   * `thin`, `burn_in`: the Run's `-t` and `-b`.
#'   * `sampling_seconds`: wall time spent on sampling Windows, summed from their `window_closed` events.
#'   * `iterations_per_second`: sampling iterations (thinned or not) per second over those Windows; `NA` when no
#'     sampling time was recorded.
#'   * `non_finite`: the count of `non_finite_value` events, raised when a chain's prior or likelihood became NaN
#'     or infinite. Anything above 0 deserves a look.
#' @export
pt_runs <- function(a) {
  runs <- a$runs
  tibble::tibble(
    run = runs$run, label = runs$label, path = runs$path, draws = runs$draws,
    first_iteration = runs$first_iteration, last_iteration = runs$last_iteration,
    windows = vapply(runs$windows, nrow, integer(1)),
    status = runs$status, resumes = runs$resumes, started_from = runs$started_from,
    thin = vapply(runs$manifest, function(m) as.integer(m$thin), integer(1)),
    burn_in = vapply(runs$manifest, function(m) as.integer(m$burn_in), integer(1)),
    sampling_seconds = runs$sampling_seconds, iterations_per_second = runs$iterations_per_second,
    non_finite = runs$non_finite
  )
}

#' The Infections and Loci of an Analysis
#'
#' The input's Infections and Loci, from the Store's `structure/` tables. They are the same in every Run of an
#' Analysis.
#'
#' The Store does not record observation times, so `observation_time` is `NA` for now.
#'
#' @inheritParams pt_runs
#' @return `pt_infections()`: a tibble of `infection`, `observation_time`, `symptomatic` (which infection-duration
#'   prior the Infection uses) and `latent_parent` (the id of the Infection's Latent parent). `pt_loci()`: a tibble
#'   of `locus` and `n_alleles`. Both are sorted by id as strings.
#' @export
pt_infections <- function(a) a$infections

#' @rdname pt_infections
#' @export
pt_loci <- function(a) a$loci

#' Close an Analysis
#'
#' Closes the Analysis's DuckDB connection. Lazy tables made from it stop working. An Analysis that is not closed is
#' closed when it is garbage-collected.
#'
#' @inheritParams pt_runs
#' @return `a`, invisibly.
#' @export
pt_close <- function(a) {
  close_engine(a$engine)
  invisible(a)
}

#' @export
print.pt_analysis <- function(x, ...) {
  runs <- pt_runs(x)
  cat(sprintf("<pt_analysis> %d Run%s, %d Infections, %d Loci\n", nrow(runs), if (nrow(runs) == 1) "" else "s",
              nrow(x$infections), nrow(x$loci)))
  iterations <- ifelse(is.na(runs$first_iteration), "-",
                       paste0(format_count(runs$first_iteration), "-", format_count(runs$last_iteration)))
  table <- data.frame(
    run = runs$run, draws = format_count(runs$draws), iterations = iterations,
    windows = format(runs$windows), thin = format(runs$thin)
  )
  lines <- vapply(names(table), function(column) {
    values <- c(column, table[[column]])
    if (column %in% c("run", "iterations")) formatC(values, width = -max(nchar(values)))
    else formatC(values, width = max(nchar(values)))
  }, character(nrow(table) + 1))
  cat(paste0("  ", apply(matrix(lines, ncol = ncol(table)), 1, paste, collapse = "  "), "\n"), sep = "")
  if (is.null(x$epi)) {
    cat("Epi: none\n")
  } else {
    matched <- sum(x$infections$infection %in% x$epi$infection)
    columns <- setdiff(names(x$epi), "infection")
    cat(sprintf("Epi: %d of %d Infections matched; columns %s\n", matched, nrow(x$infections),
                if (length(columns)) paste(columns, collapse = ", ") else "none"))
  }
  if (is.null(x$engine$con)) {
    cat("Engine: closed\n")
  } else {
    s <- x$engine$settings
    cat(sprintf("Engine: DuckDB, memory_limit %s, threads %s, temp_directory %s\n", s$memory_limit, s$threads,
                s$temp_directory))
  }
  invisible(x)
}

format_count <- function(x) formatC(x, format = "d", big.mark = ",")
