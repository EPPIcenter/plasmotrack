# Pooled Parent sets, the MAP graph, ψ_MAP, the size range and co-inclusion: hand-computed on fake Stores, oracle on
# fake and golden Stores.

pset <- function(iteration, infection, parents, latent, prob) {
  data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
}

# Two Runs over Infections "10", "2", "3". Run a (2 Draws) and Run b (1 Draw) meet Infection 2's sets in different
# orders, so their set ids differ: sets must be matched by members, not set_id.
local_two_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  ids <- c("10", "2", "3")
  fake_store(file.path(parent, "a"), infections = ids, parent_sets = rbind(
    pset(1, "2", "10", FALSE, 0.5), pset(1, "2", "3", TRUE, 0.25), pset(1, "2", "", TRUE, 0.25),
    pset(2, "2", "10", FALSE, 0.5), pset(2, "2", "", TRUE, 0.4), pset(2, "2", "3,10", FALSE, 0.1),
    pset(1:2, "10", "", TRUE, 1), pset(1:2, "3", "", TRUE, 1)
  ))
  fake_store(file.path(parent, "b"), infections = ids, parent_sets = rbind(
    pset(1, "2", "", TRUE, 0.25), pset(1, "2", "10", FALSE, 0.75),
    pset(1, "10", "", TRUE, 1), pset(1, "3", "", TRUE, 1)
  ))
  file.path(parent, c("a", "b"))
}

read_two_runs <- function(paths) {
  expect_warning(a <- pt_read(paths), "Draw counts differ")
  a
}

# One Run, one Draw over Infections "10", "2", "3", "4", built for ties:
#   2: {10} .6, then {3}, {4}, {Latent} tied at .1 (the third reaches 0.8, the tied fourth is admitted), then {10,L}
#      and {3,4} tied at .05 and left out.
#   4: {3} and {10} tied at .4 (MAP {10}: "10" < "3"), {Latent} .2.
#   3: {2} and {Latent} tied at .35 (MAP {2}: the Latent parent sorts last); the weights never reach 0.8.
#  10: {Latent} alone.
local_ties <- function(.local_envir = parent.frame()) {
  fake_store(infections = c("10", "2", "3", "4"), .local_envir = .local_envir, parent_sets = rbind(
    pset(1, "2", c("10", "3", "4", "", "4,3", "10"), c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE),
         c(0.6, 0.1, 0.1, 0.1, 0.05, 0.05)),
    pset(1, "4", c("3", "10", ""), c(FALSE, FALSE, TRUE), c(0.4, 0.4, 0.2)),
    pset(1, "3", c("", "2"), c(TRUE, FALSE), c(0.35, 0.35)),
    pset(1, "10", "", TRUE, 1)
  ))
}

test_that("pt_parent_sets_pooled() pools sets matched by members across Runs", {
  a <- read_two_runs(local_two_runs())
  on.exit(pt_close(a))
  p <- pt_parent_sets_pooled(a)
  expect_s3_class(p, "tbl_df")
  expect_equal(names(p), c("infection", "members", "n_observed", "has_latent_parent", "class", "weight", "rank",
                           "cumulative_mass", "credible"))
  two <- p[p$infection == "2", ]
  expect_equal(two$members, list("10", character(), "3", c("10", "3")))
  expect_equal(two$has_latent_parent, c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(two$n_observed, c(1, 0, 1, 2))
  expect_equal(two$class, c("observed_1", "latent_only", "observed_1_latent", "observed_2"))
  expect_equal(two$weight, c(1.75, 0.9, 0.25, 0.1) / 3, tolerance = 1e-12)
  expect_equal(two$rank, 1:4)
  expect_equal(two$cumulative_mass, cumsum(c(1.75, 0.9, 0.25, 0.1)) / 3, tolerance = 1e-12)
  expect_equal(two$credible, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(p$infection, c("10", rep("2", 4), "3"))
  expect_equal(p$class[p$infection != "2"], c("latent_only", "latent_only"))
  expect_equal(p$weight[p$infection != "2"], c(1, 1))
})

test_that("pt_parent_sets_pooled(by_run = TRUE) ranks each Run's own sets", {
  a <- read_two_runs(local_two_runs())
  on.exit(pt_close(a))
  p <- pt_parent_sets_pooled(a, by_run = TRUE)
  expect_equal(names(p)[1:2], c("run", "infection"))
  two <- p[p$infection == "2", ]
  expect_equal(two$run, c("a", "a", "a", "a", "b", "b"))
  expect_equal(two$members, list("10", character(), "3", c("10", "3"), "10", character()))
  expect_equal(two$weight, c(0.5, 0.325, 0.125, 0.05, 0.75, 0.25), tolerance = 1e-12)
  expect_equal(two$rank, c(1:4, 1:2))
  expect_equal(two$credible, c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
  expect_equal(unique(p$run), c("a", "b"))
  b <- pt_parent_sets_pooled(a, runs = "b")
  expect_equal(b$weight[b$infection == "2"], c(0.75, 0.25))
})

test_that("weight ties go to the members list, Latent parent last", {
  a <- pt_read(c(r = local_ties()))
  on.exit(pt_close(a))
  p <- pt_parent_sets_pooled(a)
  two <- p[p$infection == "2", ]
  expect_equal(two$members, list("10", "3", "4", character(), "10", c("3", "4")))
  expect_equal(two$has_latent_parent, c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(two$rank, 1:6)
  four <- p[p$infection == "4", ]
  expect_equal(four$members, list("10", "3", character()))
  three <- p[p$infection == "3", ]
  expect_equal(three$members, list("2", character()))
  expect_equal(three$rank, 1:2)
  map <- pt_map_graph(a)
  expect_equal(map$infection, c("10", "2", "3", "4"))
  expect_equal(map$members, list(character(), "10", "2", "10"))
  expect_equal(map$rank, rep(1L, 4))
  expect_equal(names(map), names(p))
})

test_that("the Credible collection admits sets tied at the cutoff and holds every set when the level is not reached", {
  a <- pt_read(c(r = local_ties()))
  on.exit(pt_close(a))
  p <- pt_parent_sets_pooled(a)
  # 0.6 + 0.1 + 0.1 reaches 0.8 at rank 3; rank 4 ties it at 0.1.
  expect_equal(p$credible[p$infection == "2"], c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(p$credible[p$infection == "4"], c(TRUE, TRUE, FALSE))
  # Infection 3's weights sum to 0.7.
  expect_equal(p$credible[p$infection == "3"], c(TRUE, TRUE))
  expect_equal(pt_parent_sets_pooled(a, level = 0.5)$credible[p$infection == "2"], c(TRUE, FALSE, FALSE, FALSE,
                                                                                     FALSE, FALSE))
  expect_true(all(pt_parent_sets_pooled(a, level = 1)$credible))
  expect_error(pt_parent_sets_pooled(a, level = 0), "level")
  expect_error(pt_parent_sets_pooled(a, level = 1.5), "level")
})

golden_paths <- function() c(one = golden_store(), four = golden_store("replica-exchange"))

test_that("pt_parent_sets_pooled() and pt_map_graph() match the oracle on fake and golden Stores", {
  paths <- local_two_runs()
  a <- read_two_runs(paths)
  on.exit(pt_close(a))
  o <- oracle_read(paths)
  keys <- c("infection", "rank")
  expect_oracle_equal(pt_parent_sets_pooled(a), oracle_parent_sets_pooled(o), keys)
  expect_oracle_equal(pt_parent_sets_pooled(a, by_run = TRUE, level = 0.5),
                      oracle_parent_sets_pooled(o, by_run = TRUE, level = 0.5), c("run", keys))

  ties <- local_ties()
  t <- pt_read(c(r = ties))
  on.exit(pt_close(t), add = TRUE)
  expect_oracle_equal(pt_parent_sets_pooled(t), oracle_parent_sets_pooled(oracle_read(c(r = ties))), keys)

  g <- pt_read(golden_paths())
  on.exit(pt_close(g), add = TRUE)
  og <- oracle_read(golden_paths())
  pooled <- oracle_parent_sets_pooled(og)
  expect_oracle_equal(pt_parent_sets_pooled(g), pooled, keys)
  expect_oracle_equal(pt_parent_sets_pooled(g, by_run = TRUE, level = 0.95),
                      oracle_parent_sets_pooled(og, by_run = TRUE, level = 0.95), c("run", keys))
  expect_oracle_equal(pt_parent_sets_pooled(g, runs = "four"), oracle_parent_sets_pooled(og, runs = "four"), keys)
  expect_oracle_equal(pt_map_graph(g), pooled[pooled$rank == 1, ], keys)
})
