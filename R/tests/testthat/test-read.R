test_that("pt_read() opens a golden Store as a pt_analysis labelled by its directory", {
  a <- pt_read(golden_store("single-chain"))
  on.exit(pt_close(a))
  expect_s3_class(a, "pt_analysis")
  expect_equal(pt_runs(a)$run, "single-chain")
  expect_equal(pt_runs(a)$label, "single-chain")
})

test_that("names given to pt_read() become Run labels", {
  a <- pt_read(c(one = golden_store("single-chain"), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$run, c("one", "four"))
  expect_equal(pt_runs(a)$label, c("one", "four"))
  expect_equal(pt_runs(a)$path, normalizePath(c(golden_store("single-chain"), golden_store("replica-exchange")),
                                              winslash = "/"))
})

test_that("pt_read() on a parent directory finds its Run subdirectories", {
  a <- pt_read(dirname(golden_store("single-chain")))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$run, c("replica-exchange", "single-chain"))
})

test_that("pt_read() errors on a directory that is not a Run and holds none", {
  expect_error(pt_read(withr::local_tempdir()), "not a Run directory")
  expect_error(pt_read(c(golden_store(), withr::local_tempdir())), "Not a Run directory")
})

test_that("pt_runs() counts Draws, Windows and iterations from the closed Windows", {
  a <- pt_read(golden_store("single-chain"))
  on.exit(pt_close(a))
  runs <- pt_runs(a)
  # The golden single-chain Run: -b 20 -s 50 -t 5 --window 25; Windows 0 (burn-in, 0-20), 1 (21-45), 2 (46-70).
  expect_equal(runs$draws, 10L)
  expect_equal(runs$first_iteration, 25)
  expect_equal(runs$last_iteration, 70)
  expect_equal(runs$windows, 3L)
  expect_equal(runs$thin, 5L)
  expect_equal(runs$burn_in, 20L)
  expect_named(runs, c("run", "label", "path", "draws", "first_iteration", "last_iteration", "windows", "thin",
                       "burn_in"))
})

test_that("pt_infections() and pt_loci() come from structure/", {
  a <- pt_read(golden_store("single-chain"))
  on.exit(pt_close(a))
  expect_equal(pt_infections(a), tibble::tibble(
    infection = c("1", "2", "3"), observation_time = NA_real_, symptomatic = c(TRUE, TRUE, FALSE),
    latent_parent = c("1_copy", "2_copy", "3_copy")
  ))
  expect_equal(pt_loci(a), tibble::tibble(locus = c("AS1", "AS11"), n_alleles = c(8L, 12L)))
})

test_that("pt_read() errors when Runs come from different inputs", {
  other <- local_golden_copy("replica-exchange")
  manifest <- jsonlite::read_json(file.path(other, "manifest.json"))
  edit_manifest(other, input_sha256 = utils::modifyList(manifest$input_sha256, list(input = strrep("0", 64))))
  expect_error(pt_read(c(a = golden_store(), b = other)), "input_sha256")
})

test_that("pt_read() errors when Runs mix the null and full models", {
  other <- local_golden_copy("replica-exchange")
  edit_manifest(other, null_model = TRUE)
  expect_error(pt_read(c(a = golden_store(), b = other)), "null_model")
})

test_that("pt_read() errors on a missing required table, naming it", {
  store <- local_golden_copy()
  unlink(file.path(store, "parent_sets"), recursive = TRUE)
  expect_error(pt_read(store), "table parent_sets")
})

test_that("pt_read() errors on a missing required column, naming it", {
  store <- local_golden_copy()
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  for (f in list.files(file.path(store, "ordering"), full.names = TRUE)) {
    DBI::dbExecute(con, sprintf("COPY (SELECT * EXCLUDE (rank) FROM read_parquet('%s')) TO '%s.tmp' (FORMAT parquet)",
                                f, f))
    file.rename(paste0(f, ".tmp"), f)
  }
  expect_error(pt_read(store), "ordering column rank")
})

test_that("pt_read() accepts Runs that differ only in chains, thin or engine version", {
  other <- local_golden_copy("replica-exchange")
  edit_manifest(other, thin = 10L, engine_version = "9.9.9")
  a <- pt_read(c(a = golden_store(), b = other))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$thin, c(5L, 10L))
})
