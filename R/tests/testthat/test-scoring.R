# pt_score_edges() and pt_score_pairs() against hand-computed AP, PR-AUC and top-k.

# One Draw over Infections 1, 2, 3. Edge marginals (parent -> child): 1->2 0.8, 1->3 0.3, 3->1 0.3, 2->1 0.1; the
# other two ordered pairs score 0.
local_scored_store <- function(.local_envir = parent.frame()) {
  sets <- data.frame(
    iteration = 1,
    infection = c("2", "2", "3", "3", "1", "1", "1"),
    parents = c("1", "", "1", "", "3", "2", ""),
    latent = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE),
    prob = c(0.8, 0.2, 0.3, 0.7, 0.3, 0.1, 0.6)
  )
  a <- pt_read(c(r = fake_store(parent_sets = sets, infections = c("1", "2", "3"), .local_envir = .local_envir)))
  withr::defer(pt_close(a), envir = .local_envir)
  a
}

test_that("pt_score_edges() ranks direct edges, directed and undirected", {
  a <- local_scored_store()
  # Truth 1 -> 2, 1 -> 3, plus edges outside the Analysis (4 -> 1) that the directed universe ignores.
  truth <- pt_read_truth(write_truth(edges = list(c(1, 2), c(1, 3), c(4, 1)), nodes = lapply(1:4, truth_node)))
  got <- pt_score_edges(a, truth)
  # Directed, descending: 0.8 T | 0.3 T, 0.3 F | 0.1 F | 0, 0 F.
  #   curve (recall, precision): (1/2, 1), (1, 2/3), (1, 1/2), (1, 1/3)
  #   AP = 1/2 * 1 + 1/2 * 2/3 = 5/6; PR-AUC from (0, 1): 1/2 * 1 + 1/2 * (1 + 2/3) / 2 = 11/12
  #   top-k: k = 2, the 2nd score is 0.3, three pairs reach it, two true: 2/3
  # Undirected: {1,2} 0.9 T, {1,3} 0.6 T, {2,3} 0 F: AP 1, PR-AUC 1, top-k 1.
  expect_equal(as.data.frame(got), data.frame(
    direction = c("directed", "undirected"), n_edge_universe = c(6L, 3L), n_true_edges = c(2L, 2L),
    n_distinct_scores = c(4L, 3L), average_precision = c(5 / 6, 1), pr_auc = c(11 / 12, 1),
    top_k_precision = c(2 / 3, 1)
  ))
  curve <- pt_score_edges(a, truth, curve = TRUE)
  expect_equal(as.data.frame(curve), data.frame(
    direction = rep(c("directed", "undirected"), c(4, 3)),
    threshold = c(0.8, 0.3, 0.1, 0, 0.9, 0.6, 0),
    recall = c(1 / 2, 1, 1, 1, 1 / 2, 1, 1),
    precision = c(1, 2 / 3, 1 / 2, 1 / 3, 1, 1, 2 / 3)
  ))
})

test_that("transitive truth counts a path through unobserved cases; no true edge gives NA and no curve", {
  a <- local_scored_store()
  # 5, 6 and 7 are outside the Analysis. 1 -> 5 -> 2 and 2 -> 6 -> 3 are paths through unobserved cases only;
  # 1 -> 5 -> 2 -> 6 -> 3 passes through the observed 2, so 1 -> 3 is not true. No direct edge joins two Infections.
  truth <- pt_read_truth(write_truth(
    edges = list(c(1, 5), c(5, 2), c(2, 6), c(6, 3), c(3, 7)), nodes = lapply(1:7, truth_node)
  ))
  direct <- pt_score_edges(a, truth)
  expect_equal(direct$n_true_edges, c(0L, 0L))
  expect_equal(direct$n_edge_universe, c(6L, 3L))
  expect_equal(direct$n_distinct_scores, c(4L, 3L))
  expect_true(all(is.na(direct[c("average_precision", "pr_auc", "top_k_precision")])))
  empty <- pt_score_edges(a, truth, curve = TRUE)
  expect_equal(nrow(empty), 0)
  expect_equal(names(empty), c("direction", "threshold", "recall", "precision"))

  # Transitive truth 1 -> 2 and 2 -> 3. Directed: 0.8 T | 0.3 F, 0.3 F | 0.1 F | 0 T, 0 F.
  #   curve: (1/2, 1), (1/2, 1/3), (1/2, 1/4), (1, 1/3)
  #   AP = 1/2 + 1/2 * 1/3 = 2/3; PR-AUC = 1/2 + 1/2 * (1/4 + 1/3) / 2 = 31/48; top-k: 2nd score 0.3, 1 of 3
  # Undirected: {1,2} 0.9 T, {1,3} 0.6 F, {2,3} 0 T: AP = 1/2 + 1/2 * 2/3 = 5/6;
  #   PR-AUC = 1/2 + 1/2 * (1/2 + 2/3) / 2 = 19/24; top-k: 2nd score 0.6, 1 of 2.
  got <- pt_score_edges(a, truth, edges = "transitive")
  expect_equal(as.data.frame(got), data.frame(
    direction = c("directed", "undirected"), n_edge_universe = c(6L, 3L), n_true_edges = c(2L, 2L),
    n_distinct_scores = c(4L, 3L), average_precision = c(2 / 3, 5 / 6), pr_auc = c(31 / 48, 19 / 24),
    top_k_precision = c(1 / 3, 1 / 2)
  ))
  expect_error(pt_score_edges(a, truth, edges = "indirect"), "edges")
})

test_that("pt_score_edges() needs every Infection in the truth", {
  a <- local_scored_store()
  truth <- pt_read_truth(write_truth(edges = list(c(1, 2)), nodes = lapply(1:2, truth_node)))
  expect_error(pt_score_edges(a, truth), "not in the truth")
  expect_error(pt_score_edges(a, list()), "pt_truth")
})

test_that("pt_score_pairs() keeps the largest score per unordered pair and ignores outside ids", {
  truth <- pt_read_truth(write_truth(edges = list(c(1, 2), c(3, 4), c(4, 99)),
                                     nodes = lapply(c(1:4, 99), truth_node)))
  scores <- data.frame(
    a = c("1", "2", "1", "3", "2", "1", "1", "2", "3"),
    b = c("2", "1", "3", "4", "4", "4", "99", "2", "2"),
    score = c(0.1, 0.5, 0.1 + 0.2, 0.3, 1.4, -0.2, 0.9, 0.7, NA)
  )
  got <- pt_score_pairs(scores, truth, c("1", "2", "3", "4"))
  # Unordered pairs: {1,2} 0.5 T (not 0.1), {1,3} 0.3 F (0.1 + 0.2 rounds to tie 0.3), {3,4} 0.3 T, {2,4} 1 F
  # (clipped), {1,4} 0 F (clipped), {2,3} 0 F (NA ignored). The 1-99 and 2-2 rows are ignored.
  # Undirected: 1 F | 0.5 T | 0.3 F, 0.3 T | 0, 0: curve (0, 0), (1/2, 1/2), (1, 1/2), (1, 1/3)
  #   AP = 1/2 * 1/2 + 1/2 * 1/2 = 1/2; PR-AUC = 1/2 * (0 + 1/2) / 2 + 1/2 * 1/2 = 3/8; top-k: 2nd score 0.5, 1 of 2
  # Directed: each pair twice, true 1 -> 2 and 3 -> 4 only: curve (0, 0), (1/2, 1/4), (1, 1/4), (1, 1/6)
  #   AP = 1/2 * 1/4 + 1/2 * 1/4 = 1/4; PR-AUC = 1/2 * (0 + 1/4) / 2 + 1/2 * 1/4 = 3/16; top-k: 2nd score 1, 0 of 2
  expect_equal(as.data.frame(got), data.frame(
    direction = c("directed", "undirected"), n_edge_universe = c(12L, 6L), n_true_edges = c(2L, 2L),
    n_distinct_scores = c(4L, 4L), average_precision = c(1 / 4, 1 / 2), pr_auc = c(3 / 16, 3 / 8),
    top_k_precision = c(0, 1 / 2)
  ))
  expect_error(pt_score_pairs(data.frame(x = 1), truth, "1"), "score")
})

test_that("pt_score_edges() matches the oracle on both golden Stores", {
  paths <- c(one = golden_store(), four = golden_store("replica-exchange"))
  a <- read_golden(paths)
  on.exit(pt_close(a))
  # Golden Infections 1, 2, 3; 4 is a case outside the Analysis.
  truth <- pt_read_truth(write_truth(edges = list(c(1, 3), c(3, 2), c(4, 2)), nodes = lapply(1:4, truth_node)))
  expected <- oracle_score_edges(oracle_read(paths), data.frame(parent = c("1", "3"), child = c("3", "2")))
  expect_oracle_equal(pt_score_edges(a, truth), expected, "direction")
})
