# Per-Draw Edges and pooled Edge marginals, on fake Stores (hand-computed) and golden Stores (oracle).

test_that("fake_store() writes a Store pt_read() accepts", {
  sets <- data.frame(iteration = c(1, 1, 2), infection = "2", parents = c("10", "", "10"),
                     latent = c(FALSE, TRUE, FALSE), prob = c(0.7, 0.3, 1))
  path <- fake_store(parent_sets = sets, infections = c("10", "2"), burn_in = 0, thin = 1)
  a <- pt_read(c(r = path))
  on.exit(pt_close(a))
  expect_equal(pt_runs(a)$draws, 2)
  expect_equal(pt_infections(a)$infection, c("10", "2"))
  expect_equal(nrow(dplyr::collect(pt_table(a, "parent_set_members"))), 2)
})

test_that("pt_edge_draws() gives observed parents only, run and iteration first", {
  sets <- data.frame(
    iteration = c(1, 1, 1, 2, 2),
    infection = c("2", "2", "2", "2", "3"),
    parents = c("10", "10,3", "", "3", ""),
    latent = c(FALSE, TRUE, TRUE, FALSE, TRUE),
    prob = c(0.5, 0.2, 0.3, 1, 1)
  )
  path <- fake_store(parent_sets = sets, infections = c("10", "2", "3"))
  a <- pt_read(c(r = path))
  on.exit(pt_close(a))
  draws <- pt_edge_draws(a)
  expect_s3_class(draws, "tbl_lazy")
  expect_equal(colnames(draws), c("run", "iteration", "child", "parent", "prob"))
  got <- dplyr::arrange(dplyr::collect(draws), iteration, child, parent)
  expect_equal(as.data.frame(got), data.frame(
    run = "r", iteration = c(1, 1, 2), child = "2", parent = c("10", "3", "3"), prob = c(0.7, 0.2, 1)
  ))
})

golden_paths <- function() c(one = golden_store(), four = golden_store("replica-exchange"))

test_that("per-Draw tables match the oracle on both golden Stores", {
  a <- pt_read(golden_paths())
  on.exit(pt_close(a))
  o <- oracle_read(golden_paths())
  expect_oracle_equal(pt_edge_draws(a), oracle_edge_draws(o), c("run", "iteration", "child", "parent"))
  expect_oracle_equal(pt_infection_draws(a), oracle_infection_draws(o), c("run", "iteration", "infection"))
  expect_oracle_equal(pt_class_draws(a), oracle_class_draws(o), c("run", "iteration", "infection", "class"))
  expect_oracle_equal(pt_out_degree_draws(a), oracle_out_degree_draws(o), c("run", "iteration", "infection"))
})

test_that("per-Draw tables take runs and phase", {
  a <- pt_read(golden_paths())
  on.exit(pt_close(a))
  o <- oracle_read(golden_paths())
  both <- c("burn_in", "sampling")
  expect_oracle_equal(pt_edge_draws(a, runs = "four", phase = both), oracle_edge_draws(o, "four", both),
                      c("run", "iteration", "child", "parent"))
  expect_oracle_equal(pt_infection_draws(a, phase = "burn_in"), oracle_infection_draws(o, phase = "burn_in"),
                      c("run", "iteration", "infection"))
  expect_oracle_equal(pt_class_draws(a, runs = "one", phase = both), oracle_class_draws(o, "one", both),
                      c("run", "iteration", "infection", "class"))
  burn_in <- dplyr::collect(pt_out_degree_draws(a, phase = "burn_in"))
  expect_true(0 %in% burn_in$iteration)
  expect_oracle_equal(burn_in, oracle_out_degree_draws(o, phase = "burn_in"), c("run", "iteration", "infection"))
  expect_error(pt_edge_draws(a, runs = "five"), "five")
  expect_error(pt_class_draws(a, phase = "warmup"), "phase")
})
