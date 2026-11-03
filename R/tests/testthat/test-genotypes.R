# Latent genotypes and COI, on a fake Store (hand-computed) and golden Stores (oracle).

# Infections 10, 2, 3 at Loci L1 (3 Alleles), L2 (4), L3 (2). Infection 2 is not genotyped at L3, Infection 3 not at
# L2: those pairs have no rows. Sampling Draws 1 and 2; iteration 0 is burn-in.
genotype_rows <- function() {
  data.frame(
    iteration = rep(c(0, 1, 2), each = 7),
    infection = rep(c("10", "10", "10", "2", "2", "3", "3"), 3),
    locus = rep(c("L1", "L2", "L3", "L1", "L2", "L1", "L3"), 3),
    alleles = c(
      "100", "1000", "10", "100", "1000", "100", "10",   # iteration 0
      "110", "0100", "01", "100", "1011", "001", "11",   # iteration 1: COI 2, 3, 2
      "110", "0010", "01", "100", "0011", "011", "10"    # iteration 2: COI 2, 2, 2
    )
  )
}

local_genotype_store <- function(.local_envir = parent.frame()) {
  path <- fake_store(infections = c("10", "2", "3"), loci = c(L1 = 3L, L2 = 4L, L3 = 2L), samples = 2,
                     tables = list(latent_genotypes = genotype_rows()), .local_envir = .local_envir)
  a <- pt_read(c(r = path))
  withr::defer(pt_close(a), envir = .local_envir)
  a
}

test_that("pt_alleles() gives one row per Allele present; character k is Allele k", {
  a <- local_genotype_store()
  alleles <- pt_alleles(a)
  expect_s3_class(alleles, "tbl_lazy")
  expect_equal(colnames(alleles), c("run", "iteration", "infection", "locus", "allele"))
  got <- dplyr::collect(dplyr::filter(alleles, iteration == 1))
  got <- as.data.frame(got[order(got$infection, got$locus, got$allele, method = "radix"), ])
  rownames(got) <- NULL
  expect_equal(got, data.frame(
    run = "r", iteration = 1,
    infection = c("10", "10", "10", "10", "2", "2", "2", "2", "3", "3", "3"),
    locus = c("L1", "L1", "L2", "L3", "L1", "L2", "L2", "L2", "L1", "L3", "L3"),
    allele = c(0, 1, 1, 1, 0, 0, 2, 3, 2, 0, 1)
  ))
  burn_in <- dplyr::collect(pt_alleles(a, phase = "burn_in"))
  expect_equal(unique(burn_in$iteration), 0)
  expect_equal(nrow(burn_in), 7)
  expect_true(all(burn_in$allele == 0))
})

test_that("pt_coi_draws() is the most Alleles present at any genotyped Locus", {
  a <- local_genotype_store()
  draws <- pt_coi_draws(a)
  expect_s3_class(draws, "tbl_lazy")
  expect_equal(colnames(draws), c("run", "iteration", "infection", "coi"))
  got <- as.data.frame(dplyr::arrange(dplyr::collect(draws), iteration, infection))
  expect_equal(got, data.frame(run = "r", iteration = c(1, 1, 1, 2, 2, 2), infection = c("10", "2", "3"),
                               coi = c(2, 3, 2, 2, 2, 2)))
})

golden_genotype_paths <- function() c(one = golden_store(), four = golden_store("replica-exchange"))

test_that("per-Draw genotype tables match the oracle on golden and fake Stores", {
  a <- pt_read(golden_genotype_paths())
  on.exit(pt_close(a))
  o <- oracle_read(golden_genotype_paths())
  both <- c("burn_in", "sampling")
  expect_oracle_equal(pt_alleles(a), oracle_alleles(o), c("run", "iteration", "infection", "locus", "allele"))
  expect_oracle_equal(pt_alleles(a, runs = "four", phase = both), oracle_alleles(o, "four", both),
                      c("run", "iteration", "infection", "locus", "allele"))
  expect_oracle_equal(pt_coi_draws(a), oracle_coi_draws(o), c("run", "iteration", "infection"))
  expect_oracle_equal(pt_coi_draws(a, phase = "burn_in"), oracle_coi_draws(o, phase = "burn_in"),
                      c("run", "iteration", "infection"))
  expect_error(pt_alleles(a, runs = "five"), "five")
  expect_error(pt_coi_draws(a, phase = "warmup"), "phase")

  fake <- local_genotype_store()
  of <- oracle_read(c(r = fake$runs$path))
  expect_oracle_equal(pt_alleles(fake, phase = both), oracle_alleles(of, phase = both),
                      c("run", "iteration", "infection", "locus", "allele"))
  expect_oracle_equal(pt_coi_draws(fake), oracle_coi_draws(of), c("run", "iteration", "infection"))
})
