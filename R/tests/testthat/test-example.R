# pt_example(): the downloaded example, served here from a local file URL. No test touches the network; every test
# points R's user cache at a temporary directory.

# Build an example archive with the documented layout from the golden Stores. `drop` leaves files out, `epi` replaces
# the epi table. Returns the archive's file URL.
local_example_asset <- function(drop = character(), epi = NULL, .local_envir = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = .local_envir)
  root <- file.path(dir, "build", "plasmotrackr-example")
  dir.create(file.path(root, "full"), recursive = TRUE)
  dir.create(file.path(root, "null"))
  file.copy(golden_store("single-chain"), file.path(root, "full"), recursive = TRUE)
  file.copy(golden_store("replica-exchange"), file.path(root, "full"), recursive = TRUE)
  file.copy(golden_store("single-chain"), file.path(root, "null"), recursive = TRUE)
  file.copy(golden_input(), file.path(root, "input.json"))
  write_truth(list(c("1", "2")), list(truth_node("1"), truth_node("2")), path = file.path(root, "truth.json"))
  if (is.null(epi)) {
    epi <- data.frame(infection = c("1", "2", "3"), district = c("North", "North", "South"), travel = c(FALSE, TRUE, FALSE))
  }
  utils::write.csv(epi, file.path(root, "epi.csv"), row.names = FALSE)
  unlink(file.path(root, drop), recursive = TRUE)
  archive <- file.path(dir, "plasmotrackr-example-test.tar.gz")
  withr::with_dir(file.path(dir, "build"), utils::tar(archive, "plasmotrackr-example", compression = "gzip"))
  file_url(archive)
}

file_url <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  paste0("file://", if (!startsWith(path, "/")) "/", path)
}

local_cache <- function(.local_envir = parent.frame()) {
  cache <- withr::local_tempdir(.local_envir = .local_envir)
  withr::local_envvar(R_USER_CACHE_DIR = cache, .local_envir = .local_envir)
  cache
}

test_that("pt_example() downloads the asset once, then reuses the cached copy", {
  cache <- local_cache()
  url <- local_example_asset()

  expect_message(ex <- pt_example(url = url), "Downloading")
  expect_named(ex, c("full", "null", "truth", "input", "epi"))
  expect_true(startsWith(normalizePath(ex$full), normalizePath(cache)))
  expect_setequal(list.files(ex$full), c("single-chain", "replica-exchange"))
  expect_identical(list.files(ex$null), "single-chain")
  expect_true(file.exists(ex$truth))
  expect_true(file.exists(ex$input))
  expect_identical(ex$epi$infection, c("1", "2", "3"))
  expect_identical(ex$epi$district, c("North", "North", "South"))
  expect_identical(ex$epi$travel, c(FALSE, TRUE, FALSE))

  a <- pt_read(ex$full, epi = ex$epi, input = ex$input)
  expect_identical(pt_runs(a)$run, c("replica-exchange", "single-chain"))
  expect_s3_class(pt_read_truth(ex$truth), "pt_truth")

  unlink(sub("^file://", "", url))
  expect_no_message(again <- pt_example(url = url))
  expect_identical(again, ex)
})

test_that("refresh = TRUE downloads the asset again", {
  local_cache()
  url <- local_example_asset()
  suppressMessages(pt_example(url = url))

  changed <- data.frame(infection = c("1", "2", "3"), district = "East", travel = TRUE)
  file.copy(sub("^file://", "", local_example_asset(epi = changed)), sub("^file://", "", url), overwrite = TRUE)
  expect_identical(pt_example(url = url)$epi$district, c("North", "North", "South"))
  expect_message(ex <- pt_example(url = url, refresh = TRUE), "Downloading")
  expect_identical(ex$epi$district, rep("East", 3))
})

test_that("a failed download errors and caches nothing", {
  cache <- local_cache()
  missing <- file_url(file.path(withr::local_tempdir(), "plasmotrackr-example-none.tar.gz"))
  expect_error(suppressWarnings(suppressMessages(pt_example(url = missing))), "download")
  expect_length(list.files(cache, recursive = TRUE), 0)
})

test_that("an archive without the documented layout errors and caches nothing", {
  cache <- local_cache()
  url <- local_example_asset(drop = c("truth.json", "null"))
  expect_error(suppressMessages(pt_example(url = url)), "truth")
  expect_length(list.files(cache, recursive = TRUE), 0)
})

test_that("pt_example() checks its arguments", {
  local_cache()
  expect_error(pt_example(url = 1), "url")
  expect_error(pt_example(url = local_example_asset(), refresh = NA), "refresh")
})
