# Flows between strata (#127): attribution, enrichment, net flow, long-lag reassignment.

# One Run, Draws 1 and 2, Infections 1-4. Strata: x = {1, 2}, y = {3}, NA = {4} (missing from epi).
# Times 1 = 0, 2 = 2, 3 = 10, 4 = NA, so lags: 2<-1 is 2, 3<-1 is 10, 3<-2 is 8, 4<-3 is NA.
#   Draw 1: 1 {L} 1; 2 {1} 0.6, {L} 0.4; 3 {1, L} 0.7, {2} 0.3; 4 {3} 1.
#   Draw 2: 1 {L} 1; 2 {1} 1;            3 {L} 0.5, {1, 2} 0.5; 4 {L} 1.
local_flow_store <- function(.local_envir = parent.frame()) {
  sets <- function(iteration, infection, parents, latent, prob) {
    data.frame(iteration = iteration, infection = infection, parents = parents, latent = latent, prob = prob)
  }
  fake_store(
    file.path(withr::local_tempdir(.local_envir = .local_envir), "flow"),
    infections = c("1", "2", "3", "4"),
    parent_sets = rbind(
      sets(1, "1", "", TRUE, 1), sets(1, "2", "1", FALSE, 0.6), sets(1, "2", "", TRUE, 0.4),
      sets(1, "3", "1", TRUE, 0.7), sets(1, "3", "2", FALSE, 0.3), sets(1, "4", "3", FALSE, 1),
      sets(2, "1", "", TRUE, 1), sets(2, "2", "1", FALSE, 1), sets(2, "3", "", TRUE, 0.5),
      sets(2, "3", "1,2", FALSE, 0.5), sets(2, "4", "", TRUE, 1)
    ),
    observation_time = c("1" = 0, "2" = 2, "3" = 10),
    allowed_parents = list("2" = "1", "3" = c("1", "2"), "4" = c("3", "2"))
  )
}

flow_epi <- function() data.frame(infection = c("1", "2", "3"), district = c("x", "x", "y"))

read_flow_store <- function(path = local_flow_store(.local_envir = parent.frame()), epi = flow_epi()) {
  pt_read(path, epi = epi)
}

sorted <- function(x, keys) {
  x <- as.data.frame(dplyr::collect(x))
  x <- x[do.call(order, c(unname(as.list(x[keys])), method = "radix")), , drop = FALSE]
  rownames(x) <- NULL
  x
}

test_that("pt_attribution splits each recipient stratum's parent mass by source, Latent parent last", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  x <- pt_attribution(a, by = "district")
  expect_s3_class(x, "tbl_lazy")
  x <- sorted(x, c("run", "iteration", "recipient", "latent", "source"))
  expect_equal(names(x), c("run", "iteration", "recipient", "source", "latent", "share"))
  one <- x[x$iteration == 1, ]
  # Draw 1. x = {1, 2}: 2<-1 0.6 from x; Latent parent 1 + 0.4. y = {3}: 0.7 + 0.3 from x; Latent parent 0.7.
  # NA = {4}: 1 from y.
  expect_equal(one$recipient, rep(c("x", "y", NA), each = 4))
  expect_equal(one$source, rep(c("x", "y", NA, NA), 3))
  expect_equal(one$latent, rep(c(FALSE, FALSE, FALSE, TRUE), 3))
  expect_equal(one$share, c(0.3, 0, 0, 0.7, 1 / 1.7, 0, 0, 0.7 / 1.7, 0, 1, 0, 0), tolerance = 1e-12)
  # Draw 2. x: 2<-1 1, Latent parent 1. y: 0.5 + 0.5 from x, Latent parent 0.5. NA: Latent parent 1.
  two <- x[x$iteration == 2, ]
  expect_equal(two$share, c(0.5, 0, 0, 0.5, 1 / 1.5, 0, 0, 0.5 / 1.5, 0, 0, 0, 1), tolerance = 1e-12)
  sums <- tapply(x$share, paste(x$iteration, x$recipient), sum)
  expect_equal(as.vector(sums), rep(1, 6), tolerance = 1e-12)
})

test_that("pt_enrichment divides each cell's edge share by its candidate share, NA without candidates", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  x <- pt_enrichment(a, by = "district")
  expect_s3_class(x, "tbl_lazy")
  x <- sorted(x, c("run", "iteration", "source", "recipient"))
  expect_equal(names(x), c("run", "iteration", "source", "recipient", "enrichment"))
  # Candidates 2<-1 (x, x), 3<-1 and 3<-2 (x, y), 4<-3 (y, NA), 4<-2 (x, NA): shares 0.2, 0.4, 0.2, 0.2.
  # Draw 1 edge mass: (x, x) 0.6, (x, y) 1, (y, NA) 1, of 2.6. Draw 2: (x, x) 1, (x, y) 1, of 2.
  expect_equal(x$source, rep(rep(c("x", "y", NA), each = 3), 2))
  expect_equal(x$recipient, rep(c("x", "y", NA), 6))
  expect_equal(x$enrichment, c(0.6 / 2.6 / 0.2, 1 / 2.6 / 0.4, 0, NA, NA, 1 / 2.6 / 0.2, NA, NA, NA,
                               2.5, 1.25, 0, NA, NA, 0, NA, NA, NA), tolerance = 1e-12)
  # Candidates given by hand; unknown ids and repeats are dropped.
  candidates <- data.frame(parent = c("1", "1", "9", "3"), child = c("2", "2", "3", "4"))
  y <- sorted(pt_enrichment(a, by = "district", candidates = candidates), c("run", "iteration", "source", "recipient"))
  expect_equal(y$enrichment[y$iteration == 1], c(0.6 / 2.6 / 0.5, NA, NA, NA, NA, 1 / 2.6 / 0.5, NA, NA, NA),
               tolerance = 1e-12)
  expect_error(pt_enrichment(a, by = "district", candidates = data.frame(from = "1")), "parent")
})

test_that("pt_net_flow counts expected out-edges minus in-edges between different strata", {
  a <- read_flow_store()
  on.exit(pt_close(a))
  x <- pt_net_flow(a, by = "district")
  expect_s3_class(x, "tbl_lazy")
  x <- sorted(x, c("run", "iteration", "stratum"))
  expect_equal(names(x), c("run", "iteration", "stratum", "net_flow"))
  # Draw 1: x -> y 0.7 + 0.3, y -> NA 1 (2<-1 stays within x). Draw 2: x -> y 1.
  expect_equal(x$stratum, rep(c("x", "y", NA), 2))
  expect_equal(x$net_flow, c(1, 0, -1, 1, -1, 0), tolerance = 1e-12)
})
