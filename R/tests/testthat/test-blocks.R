test_that("pt_acceptance() sums each sampler type's sampling counts over Windows", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  acc <- pt_acceptance(a)
  expect_s3_class(acc, "tbl_df")
  expect_equal(names(acc), c("type", "accepted", "proposed", "rate"))
  # Hand sums of the golden Store's SALTSampler rows (AS1, AS11) in sampling Windows 1 and 2.
  salt <- acc[acc$type == "SALTSampler", ]
  expect_equal(salt$accepted, 6 + 7 + 10 + 12)
  expect_equal(salt$proposed, 18 + 18 + 27 + 27)
  expect_equal(salt$rate, 35 / 90)
  expect_setequal(acc$type, c(
    "AlleleSwapSampler", "CloneOrderSwapSampler", "ConstrainedContinuousRandomWalk", "GuidedGroupSampler", "GuidedPairSampler",
    "InformedOrderJumpSampler", "PairLocusGibbsSampler", "ParentRelinkSampler", "RandomAllelesBitSetSampler3",
    "RandomAllelesBitSetSampler4", "SALTSampler", "StrainCountSampler"
  ))
})

test_that("rate is NA when nothing was proposed", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  acc <- pt_acceptance(a, level = "sampler")
  # The golden Run's CloneOrderSwapSampler never proposes: one Infection order, no clone pairs.
  clone <- acc[acc$type == "CloneOrderSwapSampler", ]
  expect_equal(nrow(clone), 3)
  expect_equal(clone$proposed, c(0, 0, 0))
  expect_equal(clone$rate, rep(NA_real_, 3))
  expect_true(is.na(pt_acceptance(a)$rate[pt_acceptance(a)$type == "CloneOrderSwapSampler"]))
})

test_that("level = \"sampler\" gives one row per sampler id, after its type", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  acc <- pt_acceptance(a, level = "sampler")
  expect_equal(names(acc), c("type", "sampler", "accepted", "proposed", "rate"))
  expect_equal(nrow(acc), 39)
  as11 <- acc[acc$sampler == "SALTSampler/allele_frequencies/AS11", ]
  expect_equal(as11$type, "SALTSampler")
  expect_equal(c(as11$accepted, as11$proposed), c(10 + 12, 27 + 27))
  expect_equal(as11$rate, 22 / 54)
  coi <- acc[acc$sampler == "ConstrainedContinuousRandomWalk/mean_coi", ]
  expect_equal(c(coi$accepted, coi$proposed, coi$rate), c(40, 50, 40 / 50))
})

test_that("counts are summed over Runs, or kept per Run with by_run, Runs in Analysis order", {
  a <- read_golden(c(one = golden_store(), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  # Hand sums of SALTSampler sampling rows: one 35 / 90, four 10 + 10 + 10 + 9 = 39 / 90.
  pooled <- pt_acceptance(a)
  salt <- pooled[pooled$type == "SALTSampler", ]
  expect_equal(c(salt$accepted, salt$proposed, salt$rate), c(74, 180, 74 / 180))

  per_run <- pt_acceptance(a, level = "sampler", by_run = TRUE)
  expect_equal(names(per_run), c("run", "type", "sampler", "accepted", "proposed", "rate"))
  expect_equal(unique(per_run$run), c("one", "four"))
  coi <- per_run[per_run$sampler == "ConstrainedContinuousRandomWalk/mean_coi", ]
  expect_equal(coi$run, c("one", "four"))
  expect_equal(coi$accepted, c(22 + 18, 11 + 13))
  expect_equal(coi$proposed, c(50, 50))
  expect_equal(coi$rate, c(40 / 50, 24 / 50))
})

test_that("phase picks the Windows summed, and per_window keeps one row per Window", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  # Burn-in Window 0 of the golden Store: SALTSampler AS1 34 / 46, AS11 43 / 69.
  burn_in <- pt_acceptance(a, phase = "burn_in")
  salt <- burn_in[burn_in$type == "SALTSampler", ]
  expect_equal(c(salt$accepted, salt$proposed), c(77, 115))

  windows <- pt_acceptance(a, phase = c("burn_in", "sampling"), per_window = TRUE)
  expect_equal(names(windows), c("window", "phase", "type", "accepted", "proposed", "rate"))
  salt <- windows[windows$type == "SALTSampler", ]
  expect_equal(salt$window, c(0, 1, 2))
  expect_equal(salt$phase, c("burn_in", "sampling", "sampling"))
  expect_equal(salt$accepted, c(34 + 43, 6 + 10, 7 + 12))
  expect_equal(salt$proposed, c(46 + 69, 18 + 27, 18 + 27))
  expect_equal(salt$rate, c(77 / 115, 16 / 45, 19 / 45))
})

test_that("per_window with by_run keeps run first, then the Window", {
  a <- read_golden(c(one = golden_store(), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  windows <- pt_acceptance(a, level = "sampler", by_run = TRUE, per_window = TRUE)
  expect_equal(names(windows), c("run", "window", "phase", "type", "sampler", "accepted", "proposed", "rate"))
  coi <- windows[windows$sampler == "ConstrainedContinuousRandomWalk/mean_coi", ]
  expect_equal(coi$run, c("one", "one", "four", "four"))
  expect_equal(coi$window, c(1, 2, 1, 2))
  expect_equal(coi$accepted, c(22, 18, 11, 13))
  expect_equal(coi$rate, c(22, 18, 11, 13) / 25)
  # Pooled over Runs, Windows with the same number are summed.
  pooled <- pt_acceptance(a, level = "sampler", per_window = TRUE)
  coi <- pooled[pooled$sampler == "ConstrainedContinuousRandomWalk/mean_coi", ]
  expect_equal(c(coi$accepted, coi$proposed), c(33, 31, 50, 50))
})

test_that("pt_swaps() sums each temperature pair's sampling swaps per Run", {
  a <- read_golden(c(four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  swaps <- pt_swaps(a)
  expect_s3_class(swaps, "tbl_df")
  expect_equal(names(swaps), c("run", "temperature_pair", "accepted", "attempted", "rate"))
  # Hand sums of the golden Store's sampling Windows 1 and 2.
  expect_equal(swaps$run, rep("four", 3))
  expect_equal(swaps$temperature_pair, c(0, 1, 2))
  expect_equal(swaps$accepted, c(12 + 10, 9 + 11, 8 + 7))
  expect_equal(swaps$attempted, c(13 + 12, 12 + 13, 13 + 12))
  expect_equal(swaps$rate, c(22, 20, 15) / 25)
})

test_that("pt_swaps() is empty for one-chain Runs", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  swaps <- pt_swaps(a)
  expect_equal(nrow(swaps), 0)
  expect_equal(names(swaps), c("run", "temperature_pair", "accepted", "attempted", "rate"))
  expect_type(swaps$rate, "double")
  expect_equal(nrow(pt_swaps(a, by_run = FALSE, per_window = TRUE)), 0)
})

test_that("pt_swaps() pools over Runs, skipping one-chain Runs, and keeps Windows with per_window", {
  a <- read_golden(c(one = golden_store(), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  expect_equal(unique(pt_swaps(a)$run), "four")
  pooled <- pt_swaps(a, by_run = FALSE)
  expect_equal(names(pooled), c("temperature_pair", "accepted", "attempted", "rate"))
  expect_equal(pooled$accepted, c(22, 20, 15))

  windows <- pt_swaps(a, phase = c("burn_in", "sampling"), per_window = TRUE)
  expect_equal(names(windows), c("run", "window", "phase", "temperature_pair", "accepted", "attempted", "rate"))
  pair0 <- windows[windows$temperature_pair == 0, ]
  expect_equal(pair0$window, c(0, 1, 2))
  expect_equal(pair0$phase, c("burn_in", "sampling", "sampling"))
  expect_equal(pair0$accepted, c(7, 12, 10))
  expect_equal(pair0$rate, c(7 / 10, 12 / 13, 10 / 12))
})

test_that("by_run, per_window and level are checked", {
  a <- read_golden(c(one = golden_store()))
  on.exit(pt_close(a))
  expect_error(pt_acceptance(a, by_run = NA), "by_run")
  expect_error(pt_swaps(a, per_window = "yes"), "per_window")
  expect_error(pt_acceptance(a, level = "kind"), "level")
  expect_error(pt_acceptance(a, phase = "sample"), "phase")
})
