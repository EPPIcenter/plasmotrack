# pt_example(): download the example simulation once into the user's R cache directory and return its paths.
#
#   example_url          the release asset's URL. The asset's file name (minus .tar.gz) names its cache directory,
#                        so a new asset version gets a fresh cache entry.
#   example_files        the archive's documented layout, checked before anything is cached.

example_url <- paste0(
  "https://github.com/EPPIcenter/plasmotrack/releases/download/",
  "plasmotrackr-example-v1/plasmotrackr-example-v1.tar.gz"
)

example_root <- "plasmotrackr-example"

#' Download the example simulation
#'
#' Downloads a simulated outbreak analysed with plasmotrack, about 23 MB, the first time you call it, and returns
#' the paths of its files. Later calls reuse the copy in your R user cache directory,
#' `tools::R_user_dir("plasmotrackr", "cache")`, and download nothing. The articles on the package website use
#' this example.
#'
#' The example is one simulation of 214 Infections with a mean COI of 2 and medium genetic diversity. It holds:
#'
#' * four Runs of the full model (`full`), which open together as one Analysis.
#' * one Run of the null model (`null`), on the same input. Open it as its own Analysis: [pt_read()] errors when
#'   full-model and null-model Runs are mixed.
#' * the simulator's truth (`truth`), for [pt_read_truth()].
#' * the input JSON the Runs were started from (`input`), for `pt_read(input = )`.
#' * a synthetic epi table (`epi`), for `pt_read(epi = )`: `district` (character) and `travel` (logical). These
#'   values are random draws, not part of the simulation, so the strata show the methods, not a real pattern.
#'
#' The Runs were written without checkpoints, so they cannot be resumed.
#'
#' The archive unpacks into a directory `plasmotrackr-example/` that holds `full/` (one subdirectory per Run),
#' `null/` (one subdirectory), `truth.json` or `truth.json.gz`, `input.json` or `input.json.gz`, and `epi.csv`
#' (columns `infection`, `district`, `travel`).
#'
#' @param refresh `TRUE` to download the example again and replace the cached copy.
#' @param url Where to download the archive from. The default is the release asset on GitHub; give another URL,
#'   such as a `file://` URL, to use a copy you host yourself.
#' @return A list:
#'   * `full`: the directory of the four full-model Runs, for [pt_read()].
#'   * `null`: the directory of the null-model Run, for [pt_read()].
#'   * `truth`: the path of the truth JSON.
#'   * `input`: the path of the input JSON.
#'   * `epi`: a tibble with columns `infection` (character), `district` and `travel`.
#' @export
#' @examples
#' \dontrun{
#' ex <- pt_example()
#' a <- pt_read(ex$full, epi = ex$epi, input = ex$input)
#' truth <- pt_read_truth(ex$truth)
#' }
pt_example <- function(refresh = FALSE, url = example_url) {
  call <- rlang::current_env()
  if (!rlang::is_string(url) || !nzchar(url)) cli::cli_abort("{.arg url} must be a single URL.", call = call)
  check_flag(refresh)

  name <- sub("\\.tar\\.gz$|\\.tgz$", "", basename(url))
  cache <- tools::R_user_dir("plasmotrackr", "cache")
  target <- file.path(cache, name)
  if (refresh || !dir.exists(target)) {
    download_example(url, cache, target, call)
  }
  example_paths(file.path(target, example_root), call)
}

# Download and unpack the archive into a staging directory beside `target`, check its layout, then move it into
# place. Nothing is left in the cache when any step fails.
download_example <- function(url, cache, target, call) {
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile("download-", tmpdir = cache)
  dir.create(staging)
  on.exit({
    unlink(staging, recursive = TRUE)
    if (length(list.files(cache, all.files = TRUE, no.. = TRUE)) == 0) unlink(cache, recursive = TRUE)
  })

  cli::cli_inform(c(
    "Downloading the plasmotrackr example from {.url {cli_escape(url)}}.",
    i = "Caching it in {.path {cli_escape(target)}}."
  ))
  archive <- file.path(staging, "example.tar.gz")
  old <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old), add = TRUE, after = FALSE)
  status <- tryCatch(
    utils::download.file(url, archive, mode = "wb", quiet = TRUE),
    error = function(e) e
  )
  if (inherits(status, "error") || !identical(status, 0L) || !file.exists(archive)) {
    reason <- if (inherits(status, "error")) conditionMessage(status) else "the download did not complete"
    cli::cli_abort(
      c("Could not download the example from {.url {cli_escape(url)}}.", x = "{cli_escape(reason)}"),
      call = call
    )
  }
  unpacked <- file.path(staging, "unpacked")
  if (!identical(utils::untar(archive, exdir = unpacked), 0L)) {
    cli::cli_abort("Could not unpack the example archive from {.url {cli_escape(url)}}.", call = call)
  }
  example_paths(file.path(unpacked, example_root), call)

  unlink(target, recursive = TRUE)
  if (!file.rename(unpacked, target)) {
    cli::cli_abort("Could not move the example into {.path {cli_escape(target)}}.", call = call)
  }
}

# The example's paths, or an error naming every part of the layout that is missing.
example_paths <- function(root, call) {
  find_one <- function(stem) {
    path <- file.path(root, paste0(stem, c(".json", ".json.gz")))
    path[file.exists(path)][1]
  }
  paths <- list(
    full = file.path(root, "full"),
    null = file.path(root, "null"),
    truth = find_one("truth"),
    input = find_one("input"),
    epi = file.path(root, "epi.csv")
  )
  expected <- c(full = "full/", null = "null/", truth = "truth.json", input = "input.json", epi = "epi.csv")
  present <- vapply(paths, function(p) !is.na(p) && file.exists(p), logical(1))
  if (!all(present)) {
    cli::cli_abort(
      c(
        "The example archive does not have the documented layout.",
        x = "Missing from {.path {example_root}/}: {.file {expected[!present]}}.",
        i = "Run {.code pt_example(refresh = TRUE)} to download it again."
      ),
      call = call
    )
  }
  epi <- utils::read.csv(paths$epi, colClasses = c(infection = "character"), stringsAsFactors = FALSE)
  paths$epi <- tibble::as_tibble(epi)
  paths
}
