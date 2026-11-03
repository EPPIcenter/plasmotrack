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

# Two fake Runs of unequal Draw count: "r" as in local_summary_store(), "s" with 3 Draws, a class "r" lacks
# (observed_1_latent), a set dropped by the floor, and no COI rows for "3".
local_two_summary_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  ids <- c("10", "2", "3", "4")
  times <- c("10" = 1, "2" = 5, "3" = 3, "4" = 7)
  fake_store(file.path(parent, "r"), infections = ids, loci = c(L1 = 3L), observation_time = times,
             parent_sets = rbind(
               sets(1, "2", "10", FALSE, 0.5), sets(1, "2", "3", FALSE, 0.5), sets(2, "2", "10", FALSE, 1),
               sets(1, "3", "10", FALSE, 0.6), sets(1, "3", "", TRUE, 0.4), sets(2, "3", "", TRUE, 1),
               sets(1:2, "4", "3,10", FALSE, 1), sets(1:2, "10", "", TRUE, 1)
             ),
             tables = list(latent_genotypes = data.frame(
               iteration = rep(0:2, each = 3), infection = rep(c("10", "2", "3"), 3), locus = "L1",
               alleles = c("100", "100", "100", "100", "110", "011", "100", "111", "011")
             )))
  fake_store(file.path(parent, "s"), infections = ids, loci = c(L1 = 3L), observation_time = times,
             parent_sets = rbind(
               sets(1:3, "2", "3", TRUE, 1),
               sets(1, "3", "10", FALSE, 0.7), sets(1, "3", "", TRUE, 0.3 - 1e-10), sets(2:3, "3", "10", FALSE, 1),
               sets(1:3, "4", "2", FALSE, 0.5), sets(1:3, "4", "3", FALSE, 0.5), sets(1:3, "10", "", TRUE, 1)
             ),
             tables = list(latent_genotypes = data.frame(
               iteration = rep(0:3, each = 2), infection = rep(c("10", "2"), 4), locus = "L1",
               alleles = rep(c("110", "111"), 4)
             )))
  paths <- file.path(parent, c("r", "s"))
  names(paths) <- c("r", "s")
  paths
}

test_that("pt_infection_summary() matches the oracle on fake and golden Stores, pooled and by Run", {
  paths <- local_two_summary_runs()
  expect_warning(a <- pt_read(paths), "Draw counts differ")
  on.exit(pt_close(a))
  o <- oracle_read(paths)
  times <- c("10" = 1, "2" = 5, "3" = 3, "4" = 7)
  pooled <- pt_infection_summary(a)
  expect_equal(names(pooled), summary_columns(c("latent_only", "observed_1", "observed_1_latent", "observed_2")))
  expect_oracle_equal(pooled, oracle_infection_summary(o, observation_time = times), "infection")
  expect_oracle_equal(pt_infection_summary(a, by_run = TRUE),
                      oracle_infection_summary(o, by_run = TRUE, observation_time = times), c("run", "infection"))
  expect_oracle_equal(pt_infection_summary(a, runs = "s"),
                      oracle_infection_summary(o, runs = "s", observation_time = times), "infection")
  # "4" in "s": "2" and "3" tie at 0.5; "2" is the smaller id.
  by_run <- pt_infection_summary(a, by_run = TRUE)
  expect_equal(by_run$run, rep(c("r", "s"), each = 4))
  expect_equal(by_run$top_parent[by_run$run == "s" & by_run$infection == "4"], "2")

  golden <- c(one = golden_store(), four = golden_store("replica-exchange"))
  g <- read_golden(golden)
  on.exit(pt_close(g), add = TRUE)
  og <- oracle_read(golden)
  gt <- golden_observation_times()
  expect_oracle_equal(pt_infection_summary(g), oracle_infection_summary(og, observation_time = gt), "infection")
  expect_oracle_equal(pt_infection_summary(g, by_run = TRUE),
                      oracle_infection_summary(og, by_run = TRUE, observation_time = gt), c("run", "infection"))
})

test_that("pt_infection_summary() agrees with the functions it summarises", {
  g <- read_golden(c(one = golden_store(), four = golden_store("replica-exchange")))
  on.exit(pt_close(g))
  for (by_run in c(FALSE, TRUE)) {
    s <- pt_infection_summary(g, by_run = by_run)
    keys <- c(if (by_run) "run", "infection")
    key <- function(x) do.call(paste, c(unname(as.list(x[keys])), sep = "\r"))

    map <- pt_map_graph(g, by_run = by_run)
    at <- match(key(s), key(map))
    expect_identical(s$map_members, map$members[at])
    expect_identical(s$map_class, map$class[at])
    expect_identical(s$map_weight, map$weight[at])

    coi <- pt_coi(g, by_run = by_run)
    expect_identical(s$coi_median, coi$median[match(key(s), key(coi))])

    edges <- pt_edge_marginals(g, by_run = by_run)
    edges$infection <- edges$child
    for (i in seq_len(nrow(s))) {
      mine <- edges[key(edges) == key(s[i, ]) & edges$marginal > 0, ]
      if (nrow(mine) == 0) {
        expect_true(is.na(s$top_parent[[i]]))
        next
      }
      expect_equal(s$top_parent_marginal[[i]], max(mine$marginal))
      expect_equal(s$top_parent[[i]], min(mine$parent[mine$marginal == max(mine$marginal)]))
    }

    draws <- dplyr::collect(pt_out_degree_draws(g))
    infection_draws <- dplyr::collect(pt_infection_draws(g))
    shares <- dplyr::collect(pt_class_draws(g))
    if (!by_run) draws$run <- infection_draws$run <- shares$run <- ""
    for (i in seq_len(nrow(s))) {
      run <- if (by_run) s$run[[i]] else ""
      inf <- s$infection[[i]]
      expect_equal(s$expected_out_degree[[i]], mean(draws$out_degree[draws$run == run & draws$infection == inf]))
      mine <- infection_draws$run == run & infection_draws$infection == inf
      expect_equal(s$p_latent_parent[[i]], mean(infection_draws$p_latent_parent[mine]))
      expect_equal(s$mean_entropy[[i]], mean(infection_draws$entropy[mine]))
      expect_equal(s$mean_effective_size[[i]], mean(infection_draws$effective_size[mine]))
      for (class in unique(shares$class)) {
        expect_equal(s[[paste0("share_", class)]][[i]],
                     mean(shares$share[shares$run == run & shares$infection == inf & shares$class == class]))
      }
    }
  }
})

test_that("pt_infection_summary() has no R-hat or ESS columns and checks its arguments", {
  a <- local_summary_store()
  s <- pt_infection_summary(a)
  expect_false(any(grepl("rhat|ess", names(s), ignore.case = TRUE)))
  expect_error(pt_infection_summary(a, by_run = NA), "by_run")
  expect_error(pt_infection_summary(a, runs = "nope"), "nope")
})
