# Convergence report (#124): the Draw set, per-quantity status, R-hat and ESS, pt_as_draws().

golden_pair <- function() golden_store(c("single-chain", "replica-exchange"))

test_that("pt_convergence() quantities equal posterior on the oracle's matrices (golden Runs)", {
  paths <- golden_pair()
  a <- read_golden(paths)
  on.exit(pt_close(a))
  conv <- pt_convergence(a)
  expect_s3_class(conv, "pt_convergence")
  o <- oracle_read(paths)
  expected <- oracle_convergence_quantities(oracle_convergence_matrices(o, oracle_convergence_draws(o)))
  expect_s3_class(conv$quantities, "tbl_df")
  expect_identical(conv$quantities$name, expected$name)
  expect_oracle_equal(conv$quantities, expected, keys = c("group", "name"))
})
