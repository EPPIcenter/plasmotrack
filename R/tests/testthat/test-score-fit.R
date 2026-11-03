# Scoring Parent sets, the MAP graph, genotypes and COI against hand-written truth, values computed by hand.

# --- Genotypes and COI -------------------------------------------------------------------------------------------
#
# The Store is test-genotypes.R's: Infections 10, 2, 3 at Loci L1 (3 Alleles in the input), L2 (4), L3 (2); Infection
# 2 not genotyped at L3, Infection 3 not at L2. Sampling Draws 1 and 2.
score_genotype_rows <- function() {
  data.frame(
    iteration = rep(c(0, 1, 2), each = 7),
    infection = rep(c("10", "10", "10", "2", "2", "3", "3"), 3),
    locus = rep(c("L1", "L2", "L3", "L1", "L2", "L1", "L3"), 3),
    alleles = c(
      "100", "1000", "10", "100", "1000", "100", "10",
      "110", "0100", "01", "100", "1011", "001", "11",
      "110", "0010", "01", "100", "0011", "011", "10"
    )
  )
}

local_score_genotype_store <- function(.local_envir = parent.frame()) {
  path <- fake_store(infections = c("10", "2", "3"), loci = c(L1 = 3L, L2 = 4L, L3 = 2L), samples = 2,
                     tables = list(latent_genotypes = score_genotype_rows()), .local_envir = .local_envir)
  a <- pt_read(c(r = path))
  withr::defer(pt_close(a), envir = .local_envir)
  a
}

# The simulator has 5 Alleles at L1, 4 at L2 (the input kept them all), 3 at L3. Among the Analysis's Infections,
# L1 shows simulator Alleles 0, 2, 4 (input 0, 1, 2) and L3 shows 1, 2 (input 0, 1). Case 99 is not in the Analysis:
# its L1 Allele 1 does not count as observed.
score_truth <- function(l1_observed_by_3 = c(0, 0, 1, 0, 0), .local_envir = parent.frame()) {
  path <- write_truth(
    edges = list(c(99, 10)),
    nodes = list(
      truth_node(10, genotype = list(list(L1 = c(1, 0, 0, 0, 0), L2 = c(0, 1, 0, 0), L3 = c(0, 0, 1)),
                                     list(L1 = c(0, 0, 1, 0, 0), L2 = c(0, 1, 0, 0), L3 = c(0, 0, 1))),
                 observed = list(L1 = c(1, 0, 1, 0, 0), L2 = c(0, 1, 1, 0), L3 = c(0, 0, 1))),
      # Strain 3 carries simulator Allele 1 at L1, which no Infection of the Analysis shows: unrepresentable.
      truth_node(2, genotype = list(list(L1 = c(0, 0, 0, 0, 1), L2 = c(1, 0, 0, 0), L3 = c(1, 0, 0)),
                                    list(L1 = c(0, 0, 0, 0, 1), L2 = c(0, 0, 1, 0), L3 = c(1, 0, 0)),
                                    list(L1 = c(0, 1, 0, 0, 0), L2 = c(0, 0, 0, 1), L3 = c(1, 0, 0))),
                 observed = list(L1 = c(0, 0, 0, 0, 1), L2 = c(1, 0, 1, 1))),
      # Two strains that differ only at L2, where Infection 3 was not genotyped: identifiable COI 1.
      truth_node(3, genotype = list(list(L1 = c(0, 0, 1, 0, 0), L2 = c(1, 0, 0, 0), L3 = c(0, 1, 0)),
                                    list(L1 = c(0, 0, 1, 0, 0), L2 = c(0, 1, 0, 0), L3 = c(0, 1, 0))),
                 observed = list(L1 = l1_observed_by_3, L3 = c(0, 1, 0))),
      truth_node(99, genotype = list(list(L1 = c(0, 1, 0, 0, 0))), observed = list(L1 = c(0, 1, 0, 0, 0)))
    ),
    loci = list(L1 = 5, L2 = 4, L3 = 3),
    .local_envir = .local_envir
  )
  pt_read_truth(path)
}

test_that("pt_read_truth() reads each case's observed genotype in simulator numbering", {
  truth <- score_truth()
  observed <- truth$observed
  expect_equal(colnames(observed), c("infection", "locus", "allele"))
  got <- as.data.frame(observed[observed$infection %in% c("2", "99"), ])
  rownames(got) <- NULL
  expect_equal(got, data.frame(infection = c("2", "2", "2", "2", "99"), locus = c("L1", "L2", "L2", "L2", "L1"),
                               allele = c(4L, 0L, 2L, 3L, 1L)))
})

test_that("pt_score_genotypes() maps truth by the observed-Alleles rule and calls presence at 0.5", {
  a <- local_score_genotype_store()
  scores <- pt_score_genotypes(a, score_truth())
  expect_s3_class(scores, "tbl_df")
  expect_equal(colnames(scores), c("infection", "locus", "true_alleles", "called_alleles", "true_positives",
                                   "false_positives", "false_negatives", "unrepresentable", "exact"))
  expect_equal(scores$infection, c("10", "10", "10", "2", "2", "3", "3"))
  expect_equal(scores$locus, c("L1", "L2", "L3", "L1", "L2", "L1", "L3"))
  # In the input's numbering. Infection 10, L1: simulator 0, 2 -> 0, 1. L3: simulator 2 -> 1.
  expect_equal(scores$true_alleles, list(c(0L, 1L), 1L, 1L, 2L, c(0L, 2L, 3L), 1L, 0L))
  # Presence over Draws 1 and 2; 0.5 is called.
  expect_equal(scores$called_alleles, list(c(0L, 1L), c(1L, 2L), 1L, 0L, c(0L, 2L, 3L), c(1L, 2L), c(0L, 1L)))
  expect_equal(scores$true_positives, c(2L, 1L, 1L, 0L, 3L, 1L, 1L))
  expect_equal(scores$false_positives, c(0L, 1L, 0L, 1L, 0L, 1L, 1L))
  expect_equal(scores$false_negatives, c(0L, 0L, 0L, 1L, 0L, 0L, 0L))
  expect_equal(scores$unrepresentable, c(0L, 0L, 0L, 1L, 0L, 0L, 0L))
  expect_equal(scores$exact, c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE))
})

test_that("pt_score_genotypes() errors naming a Locus whose mapped Allele count differs from the input's", {
  a <- local_score_genotype_store()
  # Infection 3 no longer shows simulator Allele 2 at L1, but Infection 10 still does: L1 maps 3 Alleles still.
  expect_no_error(pt_score_genotypes(a, score_truth(l1_observed_by_3 = c(0, 0, 0, 0, 0))))
  # Infection 3 shows Allele 3 as well: 4 observed Alleles at L1 against the input's 3.
  expect_error(pt_score_genotypes(a, score_truth(l1_observed_by_3 = c(0, 0, 1, 1, 0))), "L1")
  # A truth with fewer Alleles than the input, or without the Locus.
  fewer <- score_truth()
  fewer$loci$n_alleles[fewer$loci$locus == "L2"] <- 3L
  expect_error(pt_score_genotypes(a, fewer), "L2")
  lacking <- score_truth()
  lacking$loci <- lacking$loci[lacking$loci$locus != "L3", ]
  expect_error(pt_score_genotypes(a, lacking), "L3")
})

test_that("pt_score_coi() compares the posterior COI median with the identifiable true COI", {
  a <- local_score_genotype_store()
  scores <- pt_score_coi(a, score_truth())
  expect_s3_class(scores, "tbl_df")
  expect_equal(as.data.frame(scores), data.frame(
    infection = c("10", "2", "3"),
    # Per-Draw COI: 10 has 2, 2; 2 has 3, 2; 3 has 2, 2.
    median = c(2, 2.5, 2),
    # Most distinct strain genotypes at a genotyped Locus: 10 has 2 at L1; 2 has 3 at L2 (L3 not genotyped);
    # 3's two strains differ only at L2, where it was not genotyped.
    true_coi = c(2L, 3L, 1L),
    strains = c(2L, 3L, 2L),
    error = c(0, -0.5, 1)
  ))
})
