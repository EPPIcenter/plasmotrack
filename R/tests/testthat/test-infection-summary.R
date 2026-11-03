# pt_infection_summary(): one row per Infection, each column checked by hand on a small fake Store, against the
# function it summarises, and against the oracle on the golden Stores.

# Infections "10", "2", "3", "4" over 2 sampling Draws:
#   "2":  Draw 1 {10} 0.5, {3} 0.5;  Draw 2 {10} 1
#   "3":  Draw 1 {10} 0.6, {L} 0.4;  Draw 2 {L} 1
#   "4":  {3, 10} 1 at both Draws (top parent tie between "10" and "3")
#   "10": {L} 1 at both Draws
# COI at L1: "10" 1 and 1, "2" 2 and 3, "3" 2 and 2; "4" genotyped nowhere.
local_summary_store <- function(.local_envir = parent.frame()) {
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  genotypes <- data.frame(
    iteration = rep(0:2, each = 3), infection = rep(c("10", "2", "3"), 3), locus = "L1",
    alleles = c("100", "100", "100", "100", "110", "011", "100", "111", "011")
  )
  path <- fake_store(
    infections = c("10", "2", "3", "4"), loci = c(L1 = 3L),
    parent_sets = rbind(
      sets(1, "2", "10", FALSE, 0.5), sets(1, "2", "3", FALSE, 0.5), sets(2, "2", "10", FALSE, 1),
      sets(1, "3", "10", FALSE, 0.6), sets(1, "3", "", TRUE, 0.4), sets(2, "3", "", TRUE, 1),
      sets(1:2, "4", "3,10", FALSE, 1),
      sets(1:2, "10", "", TRUE, 1)
    ),
    tables = list(latent_genotypes = genotypes),
    observation_time = c("10" = 1, "2" = 5, "3" = 3, "4" = 7),
    .local_envir = .local_envir
  )
  a <- pt_read(c(r = path))
  withr::defer(pt_close(a), envir = .local_envir)
  a
}

summary_columns <- function(classes) {
  c("infection", "observation_time", "p_latent_parent", "expected_observed_parents", paste0("share_", classes),
    "mean_entropy", "mean_effective_size", "posterior_entropy", "posterior_effective_size",
    "map_members", "map_class", "map_weight", "top_parent", "top_parent_marginal", "expected_out_degree",
    "coi_median")
}

test_that("pt_infection_summary() gives one row per Infection with hand-computed values", {
  a <- local_summary_store()
  s <- pt_infection_summary(a)
  expect_s3_class(s, "tbl_df")
  expect_equal(names(s), summary_columns(c("latent_only", "observed_1", "observed_2")))
  expect_equal(s$infection, c("10", "2", "3", "4"))
  expect_equal(s$observation_time, c(1, 5, 3, 7))
  expect_equal(s$p_latent_parent, c(1, 0, 0.7, 0))
  expect_equal(s$expected_observed_parents, c(0, 1, 0.3, 2))
  expect_equal(s$share_latent_only, c(1, 0, 0.7, 0))
  expect_equal(s$share_observed_1, c(0, 1, 0.3, 0))
  expect_equal(s$share_observed_2, c(0, 0, 0, 1))

  h <- function(w) -sum(w * log(w))
  expect_equal(s$mean_entropy, c(0, h(c(0.5, 0.5)) / 2, h(c(0.6, 0.4)) / 2, 0))
  expect_equal(s$mean_effective_size, c(1, (2 + 1) / 2, (exp(h(c(0.6, 0.4))) + 1) / 2, 1))
  expect_equal(s$posterior_entropy, c(0, h(c(0.75, 0.25)), h(c(0.7, 0.3)), 0))
  expect_equal(s$posterior_effective_size, exp(s$posterior_entropy))

  expect_equal(s$map_members, list(character(), "10", character(), c("10", "3")))
  expect_equal(s$map_class, c("latent_only", "observed_1", "latent_only", "observed_2"))
  expect_equal(s$map_weight, c(1, 0.75, 0.7, 1))

  # "4" has "10" and "3" at marginal 1 each: the tie goes to the smaller id as a string, "10".
  expect_equal(s$top_parent, c(NA, "10", "10", "10"))
  expect_equal(s$top_parent_marginal, c(NA, 0.75, 0.3, 1))
  expect_equal(s$expected_out_degree, c(0.75 + 0.3 + 1, 0, 0.25 + 1, 0))
  expect_equal(s$coi_median, c(1, 2.5, 2, NA))
})
