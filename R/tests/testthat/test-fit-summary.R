# The fit summary: one row per fit in the thesis layout (fixtures/fit-summary-columns.csv: the committed names, order
# and types), values against the oracle on the golden Stores with a hand-written truth.

fit_summary_columns <- function() {
  utils::read.csv(test_path("fixtures", "fit-summary-columns.csv"), stringsAsFactors = FALSE)
}

expect_fit_summary_layout <- function(row) {
  columns <- fit_summary_columns()
  expect_s3_class(row, "tbl_df")
  expect_equal(nrow(row), 1L)
  expect_identical(names(row), columns$name)
  expect_identical(unname(vapply(row, typeof, "")), columns$type)
}

# Golden Infections 1, 2, 3. Truth: 1 -> 3, 3 -> 2, and 4 (a case outside the Analysis) -> 2. True Parent sets: 1 a
# root ({Latent}), 2 {3, Latent} (4 is outside), 3 {1}.
golden_fit_truth <- function(.local_envir = parent.frame()) {
  pt_read_truth(write_truth(edges = list(c(1, 3), c(3, 2), c(4, 2)), nodes = lapply(1:4, truth_node),
                            .local_envir = .local_envir))
}

golden_fit_paths <- function() c(one = golden_store(), four = golden_store("replica-exchange"))

test_that("pt_fit_summary() returns one row in the committed 70-column layout", {
  a <- read_golden(golden_fit_paths())
  on.exit(pt_close(a))
  row <- pt_fit_summary(a, golden_fit_truth())
  expect_fit_summary_layout(row)
  expect_identical(row$missingness_condition, "full")
  expect_identical(row$fit_type, "posterior")
  expect_true(all(is.na(row[, 1:12])))
})

test_that("pt_fit_summary() counts cases and scores the true network by hand", {
  a <- read_golden(golden_fit_paths())
  on.exit(pt_close(a))
  row <- pt_fit_summary(a, golden_fit_truth())
  expect_identical(row$n_true_cases, 4L)
  expect_identical(row$n_observed_cases, 3L)
  # Edges 1 -> 3 and 3 -> 2 among the 3 Infections; 1 and 2 need the Latent parent.
  expect_equal(row$true_average_outdegree, 2 / 3)
  expect_equal(row$true_unobserved_source_proportion, 2 / 3)
  expect_identical(row$mean_no_parent_probability, 0)
  expect_identical(row$mean_no_parent_probability_lower, 0)
  expect_identical(row$mean_no_parent_probability_upper, 0)
})
