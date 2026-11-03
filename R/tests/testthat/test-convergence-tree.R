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
