# The input JSON: the Infections' observation times and allowed parents. New Stores embed a copy; old Stores don't.

# Parse an input JSON file, gzipped or not, into its Infections' observation times and candidate pairs.
read_input <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con))
  nodes <- jsonlite::parse_json(paste(readLines(con, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))$nodes
  id <- vapply(nodes, function(n) as.character(n$id), character(1))
  time <- vapply(nodes, function(n) if (is.null(n$observation_time)) NA_real_ else as.numeric(n$observation_time),
                 numeric(1))
  parents <- lapply(nodes, function(n) as.character(unlist(n$allowed_parents)))
  candidates <- tibble::tibble(parent = as.character(unlist(parents)), child = rep(id, lengths(parents)))
  candidates <- unique(candidates[order(candidates$child, candidates$parent, method = "radix"), ])
  list(times = tibble::tibble(infection = id, observation_time = time), candidates = candidates)
}

#' Candidate pairs of an Analysis
#'
#' The (parent, child) pairs the model may join: each Infection's `allowed_parents` in the input JSON. They come from
#' the input [pt_read()] found: the Store's embedded copy in new Stores, or the input file for old ones. The Latent parent is never a candidate here. Pass
#' them to `pt_enrichment()`, which uses them by default.
#'
#' @inheritParams pt_runs
#' @return A tibble of `parent` and `child` Infection ids, one row per candidate pair, sorted by child then parent as
#'   strings. Errors when [pt_read()] did not find the input.
#' @export
#' @examples
#' golden <- system.file("extdata", "golden", "single-chain", package = "plasmotrackr")
#' input <- system.file("extdata", "example_input.json", package = "plasmotrackr")
#' a <- pt_read(golden, input = input, threads = 1)
#' pt_candidates(a)
#' pt_close(a)
pt_candidates <- function(a) {
  if (is.null(a$input)) {
    cli::cli_abort(c(
      "This Analysis has no input JSON, which holds the candidate pairs.",
      i = "Give it with {.code pt_read(..., input = )}."
    ))
  }
  a$input$candidates
}

# The input of an Analysis: list(path, times, candidates), or NULL when it is not found.
#
# With `input` given, it is checked against every Run's manifest and a mismatch errors. Without it, pt_read() tries
# the first Run's own copy (input.json or input.json.gz in the Store root, engine >= #140), then the input path in
# its manifest's cli_arguments, and uses the first that exists and whose hash matches every Run's. When none does,
# it says once what it tried and how to supply the input, and observation times stay NA.
find_input <- function(input, runs, call = rlang::caller_env()) {
  if (!is.null(input)) return(use_input(input, runs, call = call))
  run <- runs$run[[1]]
  copies <- file.path(runs$path[[1]], c("input.json", "input.json.gz"))
  tried <- character()
  for (copy in copies[file.exists(copies)]) {
    if (all(input_matches(cli::hash_file_sha256(copy), runs))) return(input_from(copy))
    tried <- c(tried, sprintf("Run %s's copy %s: its SHA-256 differs from the manifest's", run, basename(copy)))
  }
  if (!length(tried)) tried <- sprintf("Run %s holds no copy of its input", run)

  args <- as.character(unlist(runs$manifest[[1]]$cli_arguments))
  named <- cli_input_path(args)
  if (is.na(named)) {
    tried <- c(tried, sprintf("Run %s's cli_arguments name no input", run))
  } else {
    found <- input_locations(named, cli_output_path(args), runs$path[[1]])
    found <- found[file.exists(found) & !dir.exists(found)]
    if (!length(found)) {
      tried <- c(tried, sprintf("%s, named in Run %s's cli_arguments: not found", named, run))
    } else if (all(input_matches(cli::hash_file_sha256(found[[1]]), runs))) {
      return(input_from(found[[1]]))
    } else {
      tried <- c(tried, sprintf("%s, named in Run %s's cli_arguments: its SHA-256 differs from the manifest's",
                                found[[1]], run))
    }
  }
  cli::cli_inform(c(
    "No input JSON found; observation times are NA and {.fn pt_candidates} errors.",
    rlang::set_names(cli_escape(tried), "x"),
    i = "Give the input with {.code pt_read(..., input = )}."
  ))
  NULL
}

input_from <- function(path) c(list(path = normalizePath(path, winslash = "/")), read_input(path))

# The input given to pt_read(): checked against every Run's manifest, then parsed. A mismatch errors.
use_input <- function(input, runs, call = rlang::caller_env()) {
  if (!is.character(input) || length(input) != 1 || is.na(input)) {
    cli::cli_abort("{.arg input} must be the path of the input JSON the Runs were given.", call = call)
  }
  input <- path.expand(input)
  if (!file.exists(input) || dir.exists(input)) {
    cli::cli_abort("Input file {.path {cli_escape(input)}} does not exist.", call = call)
  }
  hash <- cli::hash_file_sha256(input)
  wrong <- !input_matches(hash, runs)
  if (any(wrong)) {
    cli::cli_abort(c(
      "{.arg input} is not the input of {?Run/Runs} {.val {runs$run[wrong]}}: its SHA-256 differs.",
      x = "{.path {cli_escape(input)}}: {substr(hash, 1, 12)}",
      i = "Run manifests: {substr(runs$manifest[[which(wrong)[[1]]]]$input_sha256$input, 1, 12)}"
    ), call = call)
  }
  input_from(input)
}

# Whether `hash` is each Run's input_sha256.input.
input_matches <- function(hash, runs) {
  vapply(runs$manifest, function(m) identical(m$input_sha256$input, hash), logical(1))
}

# Options of the engine's command line that take no value; every other option takes one.
cli_flags <- c("--help", "--version", "--null-model")

# The input path in an engine command line (argv, program first): `-i`/`--input`, else the first positional argument.
# NA when there is none.
cli_input_path <- function(args) {
  value <- cli_option(args, c("-i", "--input"))
  if (!is.na(value)) return(value)
  cli_positionals(args)[1]
}

# The output directory in an engine command line: `-o`/`--output-dir`, the second positional, or `--resume`.
cli_output_path <- function(args) {
  value <- cli_option(args, c("-o", "--output-dir"))
  if (!is.na(value)) return(value)
  value <- cli_positionals(args)[2]
  if (!is.na(value)) return(value)
  cli_option(args, "--resume")
}

# The value of an option, given as `-x value`, `-xvalue`, `--name value` or `--name=value`. NA when absent.
cli_option <- function(args, names) {
  args <- args[-1]
  for (i in seq_along(args)) {
    a <- args[[i]]
    if (a %in% names) return(if (i < length(args)) args[[i + 1]] else NA_character_)
    for (name in names) {
      long <- startsWith(name, "--")
      prefix <- if (long) paste0(name, "=") else name
      if (startsWith(a, prefix) && nchar(a) > nchar(prefix) && (long || !startsWith(a, "--"))) {
        return(substring(a, nchar(prefix) + 1))
      }
    }
  }
  NA_character_
}

# Positional arguments: neither an option nor an option's value.
cli_positionals <- function(args) {
  args <- args[-1]
  out <- character()
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--") return(c(out, args[-seq_len(i)]))
    if (startsWith(a, "-") && nchar(a) > 1) {
      takes_next <- if (startsWith(a, "--")) !grepl("=", a, fixed = TRUE) && !a %in% cli_flags else nchar(a) == 2
      i <- i + if (takes_next) 2 else 1
    } else {
      out <- c(out, a)
      i <- i + 1
    }
  }
  out
}

# Where a recorded input path may be now. A relative path is relative to the directory the engine was started in:
# the Run directory minus a relative output path, when the Run directory ends with it; then R's working directory.
input_locations <- function(input, output, run_path) {
  input <- path.expand(input)
  if (is_absolute(input)) return(input)
  places <- character()
  if (!is.na(output) && !is_absolute(output)) {
    parts <- function(p) {
      p <- strsplit(gsub("\\\\", "/", p), "/", fixed = TRUE)[[1]]
      p[p != "" & p != "."]
    }
    out <- parts(output)
    run <- parts(run_path)
    n <- length(out)
    if (n > 0 && !".." %in% out && length(run) > n && identical(utils::tail(run, n), out)) {
      root <- if (startsWith(run_path, "/")) "/" else ""
      places <- paste0(root, paste(utils::head(run, -n), collapse = "/"))
    }
  }
  c(file.path(places, input), input)
}

is_absolute <- function(path) grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
