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

test_that("pt_fit_summary() matches the oracle on both golden Stores", {
  paths <- golden_fit_paths()
  a <- read_golden(paths)
  on.exit(pt_close(a))
  expected <- oracle_fit_summary(
    oracle_read(paths),
    true_edges = data.frame(parent = c("1", "3"), child = c("3", "2")),
    true_sets = list(`1` = list(members = character(), latent = TRUE),
                     `2` = list(members = "3", latent = TRUE),
                     `3` = list(members = "1", latent = FALSE)),
    n_true_cases = 4
  )
  for (fit_type in c("posterior", "null_prior")) {
    row <- pt_fit_summary(a, golden_fit_truth(), fit_type = fit_type, missingness = "missing_50")
    expect_identical(row$fit_type, fit_type)
    expect_identical(row$missingness_condition, "missing_50")
    expect_oracle_equal(row[names(expected)], expected, character())
  }
})

test_that("fit_type = \"dcifer\" scores pairwise scores and leaves the functionals NA", {
  a <- read_golden(golden_fit_paths())
  on.exit(pt_close(a))
  truth <- golden_fit_truth()
  ibd <- data.frame(a = c("1", "2", "3", "1"), b = c("3", "3", "4", "2"), score = c(0.9, 0.4, 1, 0.1))
  row <- pt_fit_summary(a, truth, fit_type = "dcifer", missingness = "missing_50", scores = ibd)
  expect_fit_summary_layout(row)
  expect_identical(row$fit_type, "dcifer")
  pairs <- pt_score_pairs(ibd, truth, c("1", "2", "3"))
  for (direction in c("directed", "undirected")) {
    for (column in setdiff(names(pairs), "direction")) {
      expect_identical(row[[paste0(direction, "_", column)]], pairs[[column]][pairs$direction == direction])
    }
  }
  # Pairs {1, 3} and {2, 3} are the true ones and score highest: perfect recovery.
  expect_identical(row$undirected_average_precision, 1)
  expect_identical(row$n_observed_cases, 3L)
  expect_equal(row$true_average_outdegree, 2 / 3)
  na_columns <- names(row)[match("estimated_average_outdegree", names(row)):match("mean_effective_parent_sets_upper",
                                                                                   names(row))]
  na_columns <- c(na_columns, names(row)[61:70])
  expect_true(all(is.na(unlist(row[na_columns]))))
  expect_false(anyNA(unlist(row[c(15:18, 49:60)])))
})

test_that("pt_fit_summary() checks fit_type and scores", {
  a <- read_golden(golden_fit_paths())
  on.exit(pt_close(a))
  truth <- golden_fit_truth()
  expect_error(pt_fit_summary(a, truth, fit_type = "dcfier"), "fit_type")
  expect_error(pt_fit_summary(a, truth, fit_type = "dcifer"), "scores")
  expect_error(pt_fit_summary(a, truth, scores = data.frame(a = "1", b = "2", score = 1)), "dcifer")
})

# A made-up path in the thesis rerun layout: <model>/<scenario>/<network>/<simulation>/<Run>.
thesis_path <- function(model = "null", simple = "simple_") {
  file.path("/data/thesis-reruns", model,
            "total_nodes_200_founder_prop_0.2_super_infected_prop_0.1_r0_1.5_mean_moi_3_asymptomatic_prop_0.5",
            paste0(simple, "tx_network_med_div"), "sim_0_s7", "chain2")
}

thesis_factors_script <- function() {
  path <- test_path("..", "..", "..", "tools", "plasmotrackr", "thesis-factors.R")
  skip_if_not(file.exists(path), "tools/ is not in the package")
  env <- new.env()
  sys.source(path, envir = env)
  env$thesis_factors
}

test_that("the thesis factor script parses a directory path into the 14 factor columns", {
  thesis_factors <- thesis_factors_script()
  columns <- fit_summary_columns()[1:14, ]
  f <- thesis_factors(thesis_path())
  expect_s3_class(f, "data.frame")
  expect_identical(names(f), columns$name)
  expect_identical(unname(vapply(f, typeof, "")), columns$type)
  expect_identical(f$run_key, "simplified__fp0.2__sp0.1__r01.5__coi3__med__sim0__s7")
  expect_identical(f$simulation_model, "simplified")
  expect_identical(c(f$total_nodes, f$founder_probability, f$superinfection_probability, f$secondary_case_rate,
                     f$founder_coi_rate, f$asymptomatic_probability, f$panel_alleles_per_locus),
                   c(200, 0.2, 0.1, 1.5, 3, 0.5, 10))
  expect_identical(c(f$panel_name, f$missingness_condition, f$fit_type), c("med", "full", "null_prior"))
  expect_identical(c(f$simulation_index, f$seed), c(0L, 7L))

  full <- thesis_factors(dirname(thesis_path("full", "")))
  expect_identical(c(full$simulation_model, full$missingness_condition, full$fit_type),
                   c("biological", "full", "posterior"))
  expect_identical(full$run_key, "biological__fp0.2__sp0.1__r01.5__coi3__med__sim0__s7")
  expect_identical(thesis_factors(thesis_path("partial"))$missingness_condition, "missing_50")
  # The thesis's own layout: tarballs under output_multichain, simple_ on the scenario, the panel bare.
  old <- thesis_factors(paste0("remote/partially_obs_multichain/simple_total_nodes_200_founder_prop_0.05_",
                               "super_infected_prop_0_r0_1_mean_moi_1_asymptomatic_prop_0.5/high/sim_0_s3/x.tar.gz"))
  expect_identical(old$run_key, "simplified__fp0.05__sp0__r01__coi1__high__sim0__s3")
  expect_identical(c(old$missingness_condition, old$fit_type, old$panel_name), c("missing_50", "posterior", "high"))
  expect_identical(nrow(thesis_factors(c(thesis_path(), thesis_path("full")))), 2L)
  expect_error(thesis_factors("runs/sim_0_s7"), "thesis")
  expect_error(thesis_factors(sub("null", "elsewhere", thesis_path())), "elsewhere")
})

test_that("factors are bound in front; their fit_type and missingness stand in for missing arguments", {
  thesis_factors <- thesis_factors_script()
  a <- read_golden(golden_fit_paths())
  on.exit(pt_close(a))
  truth <- golden_fit_truth()
  factors <- thesis_factors(thesis_path())
  row <- pt_fit_summary(a, truth, factors = factors)
  expect_fit_summary_layout(row)
  expect_identical(as.data.frame(row[1:14]), factors)
  # Explicit arguments win: the dcifer row of a fit labelled by its path.
  dcifer <- pt_fit_summary(a, truth, fit_type = "dcifer", missingness = "missing_50", factors = factors,
                           scores = data.frame(a = "1", b = "3", score = 1))
  expect_identical(c(dcifer$fit_type, dcifer$missingness_condition), c("dcifer", "missing_50"))
  expect_identical(dcifer$run_key, factors$run_key)
})

test_that("factors of another shape are bound in front as given", {
  a <- read_golden(golden_fit_paths())
  on.exit(pt_close(a))
  truth <- golden_fit_truth()
  row <- pt_fit_summary(a, truth, factors = data.frame(scenario = "a", replicate = 2L))
  expect_identical(names(row)[1:4], c("scenario", "replicate", "missingness_condition", "fit_type"))
  expect_identical(ncol(row), 60L)
  expect_error(pt_fit_summary(a, truth, factors = data.frame(scenario = c("a", "b"))), "one-row")
  expect_error(pt_fit_summary(a, truth, factors = data.frame(n_true_cases = 1)), "n_true_cases")
})
