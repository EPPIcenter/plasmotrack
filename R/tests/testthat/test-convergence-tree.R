# Convergence report (#125): tree diagnostics. Infection spread, top parents, disagreement classes, order pairs.

# Two Runs a and b of 2 Draws over Infections 1-4. `a` and `b` give each Run's Parent sets as a list: child id ->
# named numeric vector of set -> prob, the set written as observed ids joined by "," with "L" for the Latent parent
# (e.g. c("2" = 0.9, L = 0.1), c("2,L" = 1)). Every Draw of a Run has the same sets, so Run means are the probs.
# Infections not listed have the Latent parent alone.
local_tree_runs <- function(a, b, infections = c("1", "2", "3", "4"), .local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  sets <- function(spec) {
    for (id in setdiff(infections, names(spec))) spec[[id]] <- c(L = 1)
    rows <- lapply(names(spec), function(child) {
      p <- spec[[child]]
      members <- strsplit(names(p), ",", fixed = TRUE)
      do.call(rbind, lapply(1:2, function(it) data.frame(
        iteration = it, infection = child,
        parents = vapply(members, function(m) paste(setdiff(m, "L"), collapse = ","), character(1)),
        latent = vapply(members, function(m) "L" %in% m, logical(1)), prob = unname(p)
      )))
    })
    do.call(rbind, rows)
  }
  fake_store(file.path(parent, "a"), infections = infections, samples = 2, parent_sets = sets(a))
  fake_store(file.path(parent, "b"), infections = infections, samples = 2, parent_sets = sets(b))
  file.path(parent, c("a", "b"))
}

tree_row <- function(conv, infection) as.list(conv$infections[conv$infections$infection == infection, ])

class_of <- function(a, b) {
  x <- pt_read(local_tree_runs(a, b))
  on.exit(pt_close(x))
  tree_row(pt_convergence(x), "1")$class
}

test_that("each disagreement class is the first rule that matches", {
  # Both Runs' top parent is 2; weights 0.9 and 0.5.
  expect_identical(class_of(list("1" = c("2" = 0.9, L = 0.1)), list("1" = c("2" = 0.5, L = 0.5))), "same_top")
  # a: 1 has the Latent parent and infects 2 (0.9); b: 2 infects 1.
  expect_identical(class_of(list("2" = c("1" = 0.9, L = 0.1)), list("1" = c("2" = 1))), "root_flip")
  # a: 1 has the Latent parent, 2 does too; b: 2 infects 1.
  expect_identical(class_of(list(), list("1" = c("2" = 1))), "latent_vs_observed")
  # a: 2 infects 1, and 1 infects 3; b: 3 infects 1.
  expect_identical(class_of(list("1" = c("2" = 1), "3" = c("1" = 1)), list("1" = c("3" = 1))), "direction_flip")
  # a: 2 infects 1 and 3; b: 3 infects 1.
  expect_identical(class_of(list("1" = c("2" = 1), "3" = c("2" = 1)), list("1" = c("3" = 1))), "linked_parents")
  # a: 2 infects 1; b: 3 infects 1.
  expect_identical(class_of(list("1" = c("2" = 1)), list("1" = c("3" = 1))), "different_parent")
})

test_that("top-parent ties go to the smaller id, and an observed parent beats the Latent parent", {
  ids <- c("1", "10", "2", "3")
  x <- pt_read(local_tree_runs(
    a = list("1" = c("2" = 0.5, "10" = 0.5), "3" = c("2" = 0.5, L = 0.5)),
    b = list("1" = c("2" = 1), "3" = c(L = 1)), infections = ids
  ))
  on.exit(pt_close(x))
  conv <- pt_convergence(x)
  # Ids compare as strings: "10" < "2".
  expect_identical(tree_row(conv, "1")$top_parent[[1]], c(a = "10", b = "2"))
  expect_identical(tree_row(conv, "3")$top_parent[[1]], c(a = "2", b = "latent"))
})

test_that("an Infection id \"latent\" is refused: it would collide with the Latent parent's top_parent key", {
  x <- pt_read(local_tree_runs(list(`1` = c(latent = 1)), list(`1` = c(L = 1)), infections = c("1", "latent")))
  on.exit(pt_close(x))
  expect_error(pt_convergence(x), "latent")
  expect_error(pt_convergence(x), "Latent parent")
})

test_that("Infection spread uses zero-filled Edge marginals with the Latent parent as a key", {
  paths <- local_tree_runs(
    # 1: edge 1<-2 has mean 1 in both Runs; only the Latent parent differs (0.25 in a, 1 in b).
    # 2: 2<-3 only in a (0.6), 2<-4 only in b (0.7): the other Run counts 0.
    a = list("1" = c("2" = 0.75, "2,L" = 0.25), "2" = c("3" = 0.6, L = 0.4)),
    b = list("1" = c("2,L" = 1), "2" = c("4" = 0.7, L = 0.3))
  )
  x <- pt_read(paths)
  on.exit(pt_close(x))
  conv <- pt_convergence(x)
  inf <- conv$infections
  expect_identical(inf$infection, c("1", "2", "3", "4"))
  expect_equal(inf$spread, c(0.75, 0.7, 0, 0), tolerance = 1e-12)
  expect_identical(inf$flagged, c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(inf$class, c("same_top", "different_parent", NA, NA))
  expect_identical(pt_convergence(x, spread = 0.72)$infections$flagged, c(TRUE, FALSE, FALSE, FALSE))
  o <- oracle_read(paths)
  draws <- oracle_convergence_draws(o)
  q <- oracle_convergence_quantities(oracle_convergence_matrices(o, draws))
  expect_oracle_equal(inf, oracle_convergence_infections(o, draws, q), keys = "infection")
})

test_that("infections and pairs equal the oracle's (golden Runs)", {
  paths <- golden_store(c("single-chain", "replica-exchange"))
  a <- read_golden(paths)
  on.exit(pt_close(a))
  o <- oracle_read(paths)
  draws <- oracle_convergence_draws(o)
  m <- oracle_convergence_matrices(o, draws)
  q <- oracle_convergence_quantities(m)
  for (t in c(0.2, 0.05)) {
    conv <- pt_convergence(a, spread = t)
    expect_oracle_equal(conv$infections, oracle_convergence_infections(o, draws, q, t), keys = "infection")
    expect_oracle_equal(conv$pairs, oracle_convergence_pairs(m, t), keys = c("i", "j"))
  }
})

test_that("print() shows flagged Infections by class, the 10 worst Infections and the order pairs", {
  ids <- as.character(1:12)
  # 2-11: 1 infects them with weight 0.1 ... 1 in a, never in b. 12: 2 infects it in a, 3 in b.
  a <- lapply(stats::setNames(2:11, 2:11), function(i) {
    p <- (i - 1) / 10
    if (p == 1) c("1" = 1) else c("1" = p, L = 1 - p)
  })
  a[["12"]] <- c("2" = 1)
  x <- pt_read(local_tree_runs(a = a, b = list("12" = c("3" = 1)), infections = ids))
  on.exit(pt_close(x))
  conv <- pt_convergence(x)
  expect_identical(sum(conv$infections$flagged), 9L)
  expect_snapshot(print(conv))
})
