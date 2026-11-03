all_rows <- function(a, table) {
  rows <- as.data.frame(dplyr::collect(pt_table(a, table, phase = c("burn_in", "sampling"))))
  rows[do.call(order, unname(as.list(rows))), , drop = FALSE]
}

test_that("an unmarked Window and .in-flight.parquet files are never read", {
  killed <- add_unclosed_window(local_golden_copy("replica-exchange"))
  # The fixture really adds rows: a glob over each table directory would read them.
  expect_true(file.exists(file.path(killed, "parent_sets", "window-000003.parquet")))
  expect_true(file.exists(file.path(killed, "parent_set_members", ".in-flight.parquet")))

  a <- read_golden(c(r = killed))
  on.exit(pt_close(a), add = TRUE)
  golden <- read_golden(c(r = golden_store("replica-exchange")))
  on.exit(pt_close(golden), add = TRUE)

  for (table in c("run_scalars", "infections", "ordering", "allele_frequencies", "latent_genotypes", "parent_sets",
                  "parent_set_members", "chains", "acceptance", "swaps", "events")) {
    expect_equal(all_rows(a, table), all_rows(golden, table), info = table)
  }
  expect_equal(pt_runs(a)[c("draws", "last_iteration", "windows")],
               pt_runs(golden)[c("draws", "last_iteration", "windows")])
})
