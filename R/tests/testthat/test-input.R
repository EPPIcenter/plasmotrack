# Observation times and candidate pairs come from the input JSON: the Store's own copy (engine #140 on), the input given,
# or the one named in cli_arguments.
# Expected values are the literals in examples/example_input.json.

test_that("pt_read(input =) gives the input's observation times", {
  a <- pt_read(golden_store(), input = golden_input())
  on.exit(pt_close(a))
  infections <- pt_infections(a)
  expect_identical(infections$infection, c("1", "2", "3"))
  expect_identical(infections$observation_time, c(234.0737, 306.2297, 350.5))
})

test_that("a gzipped input with the same bytes inside is read when its hash matches the manifests", {
  # The engine hashes the file as given, so a gzipped input matches only Runs that were given the gzipped file.
  dir <- withr::local_tempdir()
  gz <- file.path(dir, "input.json.gz")
  con <- gzfile(gz, "wb")
  writeBin(readBin(golden_input(), "raw", file.size(golden_input())), con)
  close(con)
  run <- local_golden_copy()
  manifest <- jsonlite::read_json(file.path(run, "manifest.json"))
  hashes <- manifest$input_sha256
  hashes$input <- cli::hash_file_sha256(gz)
  edit_manifest(run, input_sha256 = hashes)
  a <- pt_read(run, input = gz)
  on.exit(pt_close(a))
  expect_identical(pt_infections(a)$observation_time, c(234.0737, 306.2297, 350.5))
})

test_that("an input whose SHA-256 differs from the manifests' errors", {
  dir <- withr::local_tempdir()
  wrong <- file.path(dir, "input.json")
  writeBin(c(readBin(golden_input(), "raw", file.size(golden_input())), charToRaw("\n")), wrong)
  expect_snapshot(pt_read(golden_store(), input = wrong), error = TRUE,
                  transform = function(x) gsub(dir, "<dir>", x, fixed = TRUE))
  expect_error(pt_read(golden_store(), input = file.path(dir, "absent.json")), "does not exist")
})

# Point a Store's manifest at `input` the way the engine records it: argv with the input as the first positional.
set_cli_input <- function(run, input, args = c("<cli>", input, "<out>", "-s", "50")) {
  edit_manifest(run, cli_arguments = as.list(args))
}

test_that("without input =, a golden Store gives its own copy's observation times", {
  expect_no_message(a <- pt_read(golden_store()))
  on.exit(pt_close(a))
  expect_identical(pt_infections(a)$observation_time, c(234.0737, 306.2297, 350.5))
})

test_that("without input =, pt_read() uses the input named in the first Run's cli_arguments", {
  run <- local_golden_copy_without_input()
  set_cli_input(run, golden_input())
  expect_no_message(a <- pt_read(run))
  on.exit(pt_close(a))
  expect_identical(pt_infections(a)$observation_time, c(234.0737, 306.2297, 350.5))
})

test_that("cli_arguments may name the input with -i or --input", {
  for (args in list(c("<cli>", "-i", golden_input(), "-o", "<out>"),
                    c("<cli>", "--seed", "1", paste0("--input=", golden_input())),
                    c("<cli>", "-s", "50", "--input", golden_input(), "<out>"))) {
    run <- local_golden_copy_without_input()
    edit_manifest(run, cli_arguments = as.list(args))
    a <- pt_read(run)
    expect_identical(pt_infections(a)$observation_time, c(234.0737, 306.2297, 350.5), info = args[[3]])
    pt_close(a)
  }
})

test_that("a missing or mismatched input in cli_arguments gives NA times and one message", {
  # The golden manifests record "<input>", a path that does not exist; the copy is removed.
  expect_snapshot(a <- pt_read(local_golden_copy_without_input()))
  expect_identical(pt_infections(a)$observation_time, rep(NA_real_, 3))
  pt_close(a)

  run <- local_golden_copy_without_input()
  wrong <- file.path(dirname(run), "input.json")
  writeBin(c(readBin(golden_input(), "raw", file.size(golden_input())), charToRaw("\n")), wrong)
  set_cli_input(run, wrong)
  expect_message(b <- pt_read(run), "SHA-256 differs", fixed = TRUE)
  expect_identical(pt_infections(b)$observation_time, rep(NA_real_, 3))
  pt_close(b)
})

test_that("a relative input in cli_arguments is found from the directory the engine was started in", {
  launch <- withr::local_tempdir()
  dir.create(file.path(launch, "runs"))
  dir.create(file.path(launch, "data"))
  file.copy(golden_store(), file.path(launch, "runs"), recursive = TRUE)
  file.copy(golden_input(), file.path(launch, "data", "input.json"))
  run <- file.path(launch, "runs", "single-chain")
  unlink(file.path(run, "input.json"))
  set_cli_input(run, NULL, c("<cli>", "data/input.json", "./runs/single-chain/", "-b", "20"))
  expect_no_message(a <- pt_read(run))
  on.exit(pt_close(a))
  expect_identical(pt_infections(a)$observation_time, c(234.0737, 306.2297, 350.5))
})

test_that("pt_candidates() returns the input's allowed parents as parent, child pairs", {
  a <- pt_read(golden_store(), input = golden_input())
  on.exit(pt_close(a))
  expect_identical(pt_candidates(a), tibble::tibble(parent = c("1", "1", "2"), child = c("2", "3", "3")))
})

test_that("pt_candidates() errors when pt_read() found no input", {
  a <- suppressMessages(pt_read(local_golden_copy_without_input()))
  on.exit(pt_close(a))
  expect_snapshot(pt_candidates(a), error = TRUE)
})

test_that("without input =, pt_read() uses the Store's own copy of the input before cli_arguments", {
  run <- local_golden_copy()
  # cli_arguments name a different file: the Store's copy comes first.
  wrong <- file.path(dirname(run), "other.json")
  writeBin(c(readBin(golden_input(), "raw", file.size(golden_input())), charToRaw("\n")), wrong)
  set_cli_input(run, wrong)
  expect_no_message(a <- pt_read(run))
  on.exit(pt_close(a))
  expect_identical(pt_infections(a)$observation_time, c(234.0737, 306.2297, 350.5))
})

test_that("the Store's copy may be gzipped, and a copy with the wrong hash gives NA times and one message", {
  run <- local_golden_copy_without_input()
  con <- gzfile(file.path(run, "input.json.gz"), "wb")
  writeBin(readBin(golden_input(), "raw", file.size(golden_input())), con)
  close(con)
  # The engine hashes the input it was given, not its copy, so a gzipped copy matches only a gzipped input.
  hashes <- jsonlite::read_json(file.path(run, "manifest.json"))$input_sha256
  hashes$input <- cli::hash_file_sha256(file.path(run, "input.json.gz"))
  edit_manifest(run, input_sha256 = hashes)
  expect_no_message(a <- pt_read(run))
  expect_identical(pt_infections(a)$observation_time, c(234.0737, 306.2297, 350.5))
  pt_close(a)

  other <- local_golden_copy()
  writeLines(c(readLines(golden_input()), ""), file.path(other, "input.json"))
  expect_snapshot(b <- pt_read(other))
  expect_identical(pt_infections(b)$observation_time, rep(NA_real_, 3))
  pt_close(b)
})

test_that("print() says where the input came from, or that it was not found", {
  a <- pt_read(c(one = golden_store()), input = golden_input(), memory_limit = "1GiB")
  on.exit(pt_close(a))
  b <- suppressMessages(pt_read(c(one = local_golden_copy_without_input()), memory_limit = "1GiB"))
  on.exit(pt_close(b), add = TRUE)
  expect_snapshot({
    print(a)
    print(b)
  }, transform = function(x) {
    x <- sub("temp_directory .*$", "temp_directory <temp>", x)
    gsub(normalizePath(golden_input(), winslash = "/"), "<golden input>", x, fixed = TRUE)
  })
})

test_that("fake_store() Stores carry their own input, with the observation times and candidates given", {
  path <- fake_store(infections = c("10", "2", "3"), observation_time = c("2" = 5, "10" = 1.5),
                     allowed_parents = list("2" = c("3", "10"), "3" = "10"))
  expect_no_message(a <- pt_read(c(r = path)))
  on.exit(pt_close(a))
  expect_identical(pt_infections(a)$observation_time, c(1.5, 5, NA))
  expect_identical(pt_candidates(a), tibble::tibble(parent = c("10", "3", "10"), child = c("2", "2", "3")))
})
