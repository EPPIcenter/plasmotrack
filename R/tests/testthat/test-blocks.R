test_that("pt_acceptance() sums each sampler type's sampling counts over Windows", {
  a <- pt_read(c(one = golden_store()))
  on.exit(pt_close(a))
  acc <- pt_acceptance(a)
  expect_s3_class(acc, "tbl_df")
  expect_equal(names(acc), c("type", "accepted", "proposed", "rate"))
  # Hand sums of the golden Store's SALTSampler rows (AS1, AS11) in sampling Windows 1 and 2.
  salt <- acc[acc$type == "SALTSampler", ]
  expect_equal(salt$accepted, 50 + 44 + 75 + 73)
  expect_equal(salt$proposed, 72 + 64 + 108 + 96)
  expect_equal(salt$rate, 242 / 340)
  expect_setequal(acc$type, c(
    "AlleleSwapSampler", "CloneOrderSwapSampler", "ConstrainedContinuousRandomWalk", "InformedOrderJumpSampler",
    "ParentRelinkSampler", "RandomAllelesBitSetSampler3", "RandomAllelesBitSetSampler4", "SALTSampler",
    "StrainCountSampler"
  ))
})

test_that("rate is NA when nothing was proposed", {
  a <- pt_read(c(one = golden_store()))
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
  a <- pt_read(c(one = golden_store()))
  on.exit(pt_close(a))
  acc <- pt_acceptance(a, level = "sampler")
  expect_equal(names(acc), c("type", "sampler", "accepted", "proposed", "rate"))
  expect_equal(nrow(acc), 36)
  as11 <- acc[acc$sampler == "SALTSampler/allele_frequencies/AS11", ]
  expect_equal(as11$type, "SALTSampler")
  expect_equal(c(as11$accepted, as11$proposed), c(75 + 73, 108 + 96))
  expect_equal(as11$rate, 148 / 204)
  coi <- acc[acc$sampler == "ConstrainedContinuousRandomWalk/mean_coi", ]
  expect_equal(c(coi$accepted, coi$proposed, coi$rate), c(24, 50, 24 / 50))
})

test_that("counts are summed over Runs, or kept per Run with by_run, Runs in Analysis order", {
  a <- pt_read(c(one = golden_store(), four = golden_store("replica-exchange")))
  on.exit(pt_close(a))
  # Hand sums of SALTSampler sampling rows: one 242 / 340, four 52 + 43 + 81 + 68 = 244 / 340.
  pooled <- pt_acceptance(a)
  salt <- pooled[pooled$type == "SALTSampler", ]
  expect_equal(c(salt$accepted, salt$proposed, salt$rate), c(486, 680, 486 / 680))

  per_run <- pt_acceptance(a, level = "sampler", by_run = TRUE)
  expect_equal(names(per_run), c("run", "type", "sampler", "accepted", "proposed", "rate"))
  expect_equal(unique(per_run$run), c("one", "four"))
  coi <- per_run[per_run$sampler == "ConstrainedContinuousRandomWalk/mean_coi", ]
  expect_equal(coi$run, c("one", "four"))
  expect_equal(coi$accepted, c(12 + 12, 13 + 7))
  expect_equal(coi$proposed, c(50, 50))
  expect_equal(coi$rate, c(24 / 50, 20 / 50))
})
