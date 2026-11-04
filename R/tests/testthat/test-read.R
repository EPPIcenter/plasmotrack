test_that("pt_read() opens a golden Store as a pt_analysis labelled by its directory", {
  a <- read_golden(golden_store("single-chain"))
  on.exit(pt_close(a))
  expect_s3_class(a, "pt_analysis")
  expect_equal(pt_runs(a)$run, "single-chain")
  expect_equal(pt_runs(a)$label, "single-chain")
})

test_that("names given to pt_read() become Run labels", {
  a <- read_golden(c(one = golden_store("single-chain"), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$run, c("one", "four"))
  expect_equal(pt_runs(a)$label, c("one", "four"))
  expect_equal(pt_runs(a)$path, normalizePath(c(golden_store("single-chain"), golden_store("replica-exchange")),
                                              winslash = "/"))
})

test_that("pt_read() on a parent directory finds its Run subdirectories", {
  a <- read_golden(dirname(golden_store("single-chain")))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$run, c("replica-exchange", "single-chain"))
})

test_that("pt_read() errors on a directory that is not a Run and holds none", {
  expect_error(read_golden(withr::local_tempdir()), "not a Run directory")
  expect_error(read_golden(c(golden_store(), withr::local_tempdir())), "Not a Run directory")
})

test_that("pt_runs() counts Draws, Windows and iterations from the closed Windows", {
  a <- read_golden(golden_store("single-chain"))
  on.exit(pt_close(a))
  runs <- pt_runs(a)
  # The golden single-chain Run: -b 20 -s 50 -t 5 --window 25; Windows 0 (burn-in, 0-20), 1 (21-45), 2 (46-70).
  expect_equal(runs$draws, 10L)
  expect_equal(runs$first_iteration, 25)
  expect_equal(runs$last_iteration, 70)
  expect_equal(runs$windows, 3L)
  expect_equal(runs$thin, 5L)
  expect_equal(runs$burn_in, 20L)
  expect_named(runs, c("run", "label", "path", "draws", "first_iteration", "last_iteration", "windows", "status",
                       "resumes", "started_from", "thin", "burn_in", "sampling_seconds", "iterations_per_second",
                       "non_finite"))
})

test_that("pt_runs() reports status, lineage, timing and health from the golden Runs' events", {
  a <- read_golden(c(one = golden_store("single-chain"), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  runs <- pt_runs(a)
  expect_equal(runs$status, c("finished", "finished"))
  expect_equal(runs$resumes, c(0L, 0L))
  expect_equal(runs$started_from, c(NA_character_, NA_character_))
  expect_equal(runs$non_finite, c(0L, 0L))
  # window_closed seconds of sampling Windows 1 and 2: 0.034737 + 0.024007, and 0.15717 + 0.103418; 25 + 25 iterations.
  expect_equal(runs$sampling_seconds, c(0.058744, 0.260588), tolerance = 1e-12)
  expect_equal(runs$iterations_per_second, c(50 / 0.058744, 50 / 0.260588), tolerance = 1e-12)
})

test_that("pt_read() opens the golden Stores without a warning", {
  expect_no_warning(a <- read_golden(c(one = golden_store("single-chain"), four = golden_store("replica-exchange"))))
  pt_close(a)
})

test_that("pt_infections() and pt_loci() come from structure/, observation times from the input", {
  a <- read_golden(golden_store("single-chain"))
  on.exit(pt_close(a))
  expect_equal(pt_infections(a), tibble::tibble(
    infection = c("1", "2", "3"), observation_time = c(234.0737, 306.2297, 350.5), symptomatic = c(TRUE, TRUE, FALSE),
    latent_parent = c("1_copy", "2_copy", "3_copy")
  ))
  expect_equal(pt_loci(a), tibble::tibble(locus = c("AS1", "AS11"), n_alleles = c(8L, 12L)))
})

test_that("pt_read() errors when Runs come from different inputs", {
  other <- local_golden_copy("replica-exchange")
  manifest <- jsonlite::read_json(file.path(other, "manifest.json"))
  edit_manifest(other, input_sha256 = utils::modifyList(manifest$input_sha256, list(input = strrep("0", 64))))
  expect_error(read_golden(c(a = golden_store(), b = other)), "input_sha256")
})

test_that("pt_read() errors when Runs mix the null and full models", {
  other <- local_golden_copy("replica-exchange")
  edit_manifest(other, null_model = TRUE)
  expect_error(read_golden(c(a = golden_store(), b = other)), "null_model")
})

test_that("pt_read() errors on a missing required table, naming it", {
  store <- local_golden_copy()
  unlink(file.path(store, "parent_sets"), recursive = TRUE)
  expect_error(read_golden(store), "table parent_sets")
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
  expect_error(read_golden(store), "ordering column rank")
})

test_that("pt_read() accepts Runs that differ only in chains, thin or engine version", {
  other <- local_golden_copy("replica-exchange")
  edit_manifest(other, thin = 10L, engine_version = "9.9.9")
  a <- read_golden(c(a = golden_store(), b = other))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$thin, c(5L, 10L))
})

engine_settings <- function(a) {
  con <- dbplyr::remote_con(pt_table(a, "run_scalars"))
  DBI::dbGetQuery(con, paste(
    "SELECT current_setting('memory_limit') AS memory_limit, current_setting('threads') AS threads,",
    "current_setting('temp_directory') AS temp_directory"
  ))
}

test_that("the connection uses the plasmotrackr engine options", {
  spill <- withr::local_tempdir()
  withr::local_options(plasmotrackr.memory_limit = "512MiB", plasmotrackr.threads = 2,
                       plasmotrackr.temp_directory = spill)
  a <- read_golden(golden_store())
  on.exit(pt_close(a))
  settings <- engine_settings(a)
  expect_equal(settings$memory_limit, "512.0 MiB")
  expect_equal(settings$threads, 2)
  expect_equal(normalizePath(settings$temp_directory), normalizePath(spill))
})

test_that("pt_read() arguments override the engine options for one Analysis", {
  spill <- withr::local_tempdir()
  withr::local_options(plasmotrackr.memory_limit = "512MiB", plasmotrackr.threads = 2)
  a <- read_golden(golden_store(), memory_limit = "256MiB", threads = 1, temp_directory = spill)
  on.exit(pt_close(a))
  settings <- engine_settings(a)
  expect_equal(settings$memory_limit, "256.0 MiB")
  expect_equal(settings$threads, 1)
  expect_equal(normalizePath(settings$temp_directory), normalizePath(spill))
})

test_that("temp_directory defaults to tempdir()", {
  withr::local_options(plasmotrackr.temp_directory = NULL)
  a <- read_golden(golden_store())
  on.exit(pt_close(a))
  expect_equal(normalizePath(engine_settings(a)$temp_directory), normalizePath(tempdir()))
})

test_that("pt_close() closes the connection", {
  a <- read_golden(golden_store())
  con <- dbplyr::remote_con(pt_table(a, "run_scalars"))
  pt_close(a)
  expect_false(DBI::dbIsValid(con))
  expect_error(pt_table(a, "run_scalars"), "closed")
})

test_that("pt_read() checks epi", {
  expect_error(read_golden(golden_store(), epi = data.frame(id = "1")), "infection")
  expect_error(read_golden(golden_store(), epi = data.frame(infection = c("1", "1"))), "one row per Infection")
})

# Scrub the temp directory and the installed golden input's path, which differ between machines.
scrub_temp <- function(x) {
  x <- sub("temp_directory .*$", "temp_directory <temp>", x)
  gsub(normalizePath(golden_input(), winslash = "/"), "<golden input>", x, fixed = TRUE)
}

test_that("print() shows one line per Run, the epi match and the engine", {
  spill <- withr::local_tempdir()
  epi <- data.frame(infection = c(1, 2, 7), district = c("a", "b", "c"), travel = c(TRUE, FALSE, NA))
  a <- read_golden(c(one = golden_store(), four = golden_store("replica-exchange")), epi = epi, memory_limit = "1GiB",
               threads = 1, temp_directory = spill)
  on.exit(pt_close(a))
  expect_snapshot(print(a), transform = scrub_temp)
  b <- read_golden(golden_store(), memory_limit = "1GiB", threads = 1, temp_directory = spill)
  pt_close(b)
  expect_snapshot(print(b), transform = scrub_temp)
})

test_that("an Analysis is closed when garbage-collected", {
  a <- read_golden(golden_store())
  con <- dbplyr::remote_con(pt_table(a, "run_scalars"))
  rm(a)
  gc()
  expect_false(DBI::dbIsValid(con))
})

test_that("a Run killed before its last Window is incomplete, and pt_read() warns about it and its Draws", {
  killed <- drop_windows(local_golden_copy(), 2)
  expect_warning(
    expect_warning(a <- read_golden(c(one = golden_store(), killed = killed)), "not finished"),
    "Draw counts differ"
  )
  on.exit(pt_close(a))
  runs <- pt_runs(a)
  expect_equal(runs$status, c("finished", "incomplete"))
  # Window 1 alone: iterations 25-45 thinned by 5, and its window_closed event.
  expect_equal(runs$draws, c(10L, 5L))
  expect_equal(runs$sampling_seconds[[2]], 0.034737, tolerance = 1e-12)
  expect_equal(runs$iterations_per_second[[2]], 25 / 0.034737, tolerance = 1e-12)
})

test_that("a Run that ended on an interrupt is interrupted, and pt_read() warns", {
  store <- local_golden_copy()
  events <- read_events(store)
  ended <- events$type == "run_ended"
  events$payload[ended] <- '{"reason":"interrupted","seconds":0.085819}'
  events <- rbind(events[seq_len(which(ended) - 1), ],
                  event(2, "interrupted", list(signal = "SIGINT", second = FALSE), iteration = 70),
                  events[which(ended):nrow(events), ])
  write_events(store, events)
  expect_warning(a <- read_golden(store), "not finished")
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$status, "interrupted")
})

test_that("a Run whose events stop at its start, or end with an unknown reason, is incomplete", {
  started <- local_golden_copy()
  write_events(started, event(0, "run_started", list(burn_in = 20, samples = 50), phase = "burn_in"))
  odd <- local_golden_copy("replica-exchange")
  events <- read_events(odd)
  events$payload[events$type == "run_ended"] <- '{"reason":"crashed","seconds":1}'
  write_events(odd, events)
  expect_warning(a <- read_golden(c(started = started, odd = odd)), "2 Runs are not finished")
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$status, c("incomplete", "incomplete"))
  expect_equal(pt_runs(a)$sampling_seconds[[1]], 0)
  expect_equal(pt_runs(a)$iterations_per_second[[1]], NA_real_)
})

test_that("non_finite counts non_finite_value events, and pt_read() warns", {
  store <- local_golden_copy()
  events <- read_events(store)
  bad <- function(term) event(1, "non_finite_value", list(temperature_slot = 0, term = term, value = "nan"),
                              iteration = 30)
  events <- rbind(events[events$window == 0, ], bad("prior"), bad("likelihood"), events[events$window > 0, ])
  write_events(store, events)
  expect_warning(a <- read_golden(c(one = golden_store(), bad = store)), "non-finite")
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$non_finite, c(0L, 2L))
  expect_equal(pt_runs(a)$status, c("finished", "finished"))
})

test_that("resumes counts run_resumed events and a resumed Run that completed is finished", {
  store <- local_golden_copy()
  events <- read_events(store)
  # Interrupted in Window 0, resumed; interrupted in Window 1, resumed; completed in Window 2.
  stop <- function(window, iteration) {
    rbind(event(window, "run_ended", list(reason = "interrupted", seconds = 1), iteration = iteration),
          event(window + 1, "run_resumed", list(iteration = iteration, window = window, checkpoint = "c",
                                                samples = 50), iteration = iteration))
  }
  events <- rbind(events[events$window == 0, ], stop(0, 20), events[events$window == 1, ], stop(1, 45),
                  events[events$window == 2, ])
  write_events(store, events)
  expect_no_warning(a <- read_golden(store))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$resumes, 2L)
  expect_equal(pt_runs(a)$status, "finished")
})

test_that("a resumed Run killed after the resume is incomplete", {
  store <- local_golden_copy()
  events <- read_events(store)
  events <- rbind(events[events$window <= 1, ],
                  event(1, "run_ended", list(reason = "interrupted", seconds = 1), iteration = 45),
                  event(2, "run_resumed", list(iteration = 45, window = 1, checkpoint = "c", samples = 50),
                        iteration = 45),
                  events[events$window == 2 & events$type == "window_closed", ])
  write_events(store, events)
  expect_warning(a <- read_golden(store), "not finished")
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$resumes, 1L)
  expect_equal(pt_runs(a)$status, "incomplete")
})

test_that("started_from carries the source run_id of a --start-from Run", {
  store <- local_golden_copy("replica-exchange")
  edit_manifest(store, started_from = list(run_id = "0123456789abcdef0123456789abcdef", engine_version = "1.0.0",
                                           temperatures = list(1, 0.8, 0.7, 0.6)))
  a <- read_golden(c(one = golden_store(), from = store))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$started_from, c(NA, "0123456789abcdef0123456789abcdef"))
})

test_that("a Run with 0 sampling Draws stays in the Run table with draws = 0", {
  burn_in_only <- drop_windows(local_golden_copy(), c(1, 2))
  expect_warning(expect_warning(a <- read_golden(c(one = golden_store(), early = burn_in_only)), "not finished"),
                 "Draw counts differ")
  on.exit(pt_close(a))
  runs <- pt_runs(a)
  expect_equal(runs$run, c("one", "early"))
  expect_equal(runs$draws, c(10L, 0L))
  expect_equal(runs$windows, c(3L, 1L))
  expect_equal(runs$first_iteration[[2]], NA_real_)
  expect_equal(runs$sampling_seconds[[2]], 0)
  expect_equal(runs$iterations_per_second[[2]], NA_real_)
})

# runs_with_draws() is internal: the helper pooled summaries and diagnostics call to skip Runs with no Draws.
test_that("runs_with_draws() skips Runs with 0 sampling Draws, with a message", {
  burn_in_only <- drop_windows(local_golden_copy(), c(1, 2))
  suppressWarnings(a <- read_golden(c(one = golden_store(), early = burn_in_only)))
  on.exit(pt_close(a))
  expect_message(expect_equal(runs_with_draws(a), "one"), "Skipping Run .*early.* no sampling Draws")
  expect_no_message(expect_equal(runs_with_draws(a, runs = "one"), "one"))
  expect_error(suppressMessages(runs_with_draws(a, runs = "early")), "No sampling Draws")
  expect_error(runs_with_draws(a, runs = "nope"), "Unknown Run")
})

test_that("print() shows each Run's status", {
  killed <- drop_windows(local_golden_copy(), 2)
  a <- suppressWarnings(read_golden(c(one = golden_store(), killed = killed), threads = 1))
  pt_close(a)
  expect_snapshot(print(a), transform = scrub_temp)
})
