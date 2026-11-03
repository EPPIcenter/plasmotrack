# Rows of one table in the given Windows of a golden Store, read with nanoparquet (no DuckDB). nanoparquet 0.5.2
# misreads DELTA_BINARY_PACKED columns past their first 32 values, which is how the engine writes
# parent_sets$set_id: don't use this for parent_sets.
golden_rows <- function(store, table, windows) {
  files <- file.path(store, table, sprintf("window-%06d.parquet", windows))
  rows <- do.call(rbind, lapply(files[file.exists(files)], function(f) as.data.frame(nanoparquet::read_parquet(f))))
  rows[do.call(order, unname(as.list(rows))), , drop = FALSE]
}

collect_sorted <- function(tbl) {
  rows <- as.data.frame(dplyr::collect(tbl))
  rownames(rows) <- NULL
  rows <- rows[do.call(order, unname(as.list(rows))), , drop = FALSE]
  rownames(rows) <- NULL
  rows
}

expect_same_rows <- function(tbl, expected) {
  rownames(expected) <- NULL
  expect_equal(collect_sorted(tbl), expected, ignore_attr = TRUE)
}

test_that("pt_table() is a lazy table of sampling rows, run first, then the Store's columns", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  tbl <- pt_table(a, "run_scalars")
  expect_s3_class(tbl, "tbl_lazy")
  expect_equal(colnames(tbl)[1:4], c("run", "iteration", "phase", "window"))
  expect_same_rows(tbl, cbind(run = "one", golden_rows(golden_store(), "run_scalars", 1:2)))
  expect_equal(nrow(dplyr::collect(tbl)), 10)
})

test_that("pt_table() unions every Run", {
  a <- read_golden(c(one = golden_store(), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  expected <- rbind(
    cbind(run = "four", golden_rows(golden_store("replica-exchange"), "infections", 1:2)),
    cbind(run = "one", golden_rows(golden_store(), "infections", 1:2))
  )
  expect_same_rows(pt_table(a, "infections"), expected)
})

test_that("phase selects burn-in rows, starting State included, or both Phases", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  burn_in <- dplyr::collect(pt_table(a, "ordering", phase = "burn_in"))
  expect_equal(sort(unique(burn_in$iteration)), c(0, 5, 10, 15, 20))
  expect_equal(unique(burn_in$phase), "burn_in")
  expect_same_rows(pt_table(a, "ordering", phase = c("burn_in", "sampling")),
                   cbind(run = "one", golden_rows(golden_store(), "ordering", 0:2)))
})

test_that("parent_set_members is read over every closed Window whatever the phase", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  # The golden Store writes every Parent set in burn-in Window 0.
  expect_same_rows(pt_table(a, "parent_set_members"),
                   cbind(run = "one", golden_rows(golden_store(), "parent_set_members", 0)))
})

test_that("runs narrows the Runs read", {
  a <- read_golden(c(one = golden_store(), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  expect_equal(unique(dplyr::pull(pt_table(a, "run_scalars", runs = "four"), run)), "four")
  expect_error(pt_table(a, "run_scalars", runs = "five"), "five")
})

test_that("a table no Run has rows for is empty, with the Store's columns", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  swaps <- dplyr::collect(pt_table(a, "swaps"))
  expect_equal(nrow(swaps), 0)
  expect_equal(names(swaps), c("run", "window", "phase", "temperature_pair", "accepted", "attempted"))
})

test_that("a table only some Runs have reads those Runs", {
  a <- read_golden(c(one = golden_store(), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  expect_same_rows(pt_table(a, "swaps"), cbind(run = "four", golden_rows(golden_store("replica-exchange"), "swaps", 1:2)))
})

test_that("an unknown table name errors and lists the tables", {
  a <- read_golden(golden_store())
  on.exit(pt_close(a))
  expect_snapshot(pt_table(a, "parent_set"), error = TRUE)
})

test_that("an unknown phase errors", {
  a <- read_golden(golden_store())
  on.exit(pt_close(a))
  expect_error(pt_table(a, "run_scalars", phase = "sample"), "phase")
})
