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

# Runs r (the Store above, 2 Draws) and s (1 Draw, where Infection 2 has all 4 Alleles at L2: COI 4).
local_two_genotype_runs <- function(.local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  loci <- c(L1 = 3L, L2 = 4L, L3 = 2L)
  fake_store(file.path(parent, "r"), infections = c("10", "2", "3"), loci = loci, samples = 2,
             tables = list(latent_genotypes = genotype_rows()))
  s <- genotype_rows()[genotype_rows()$iteration < 2, ]
  s$alleles[s$iteration == 1 & s$infection == "2" & s$locus == "L2"] <- "1111"
  fake_store(file.path(parent, "s"), infections = c("10", "2", "3"), loci = loci, samples = 1,
             tables = list(latent_genotypes = s))
  paths <- file.path(parent, c(r = "r", s = "s"))
  names(paths) <- c("r", "s")
  expect_warning(a <- pt_read(paths), "Draw counts differ")
  withr::defer(pt_close(a), envir = .local_envir)
  a
}

test_that("pt_coi() gives the type-7 median and the distribution over Draws", {
  a <- local_genotype_store()
  coi <- pt_coi(a)
  expect_s3_class(coi, "tbl_df")
  expect_equal(names(coi), c("infection", "median", "distribution"))
  expect_equal(coi$infection, c("10", "2", "3"))
  expect_equal(coi$median, c(2, 2.5, 2))
  expect_equal(lapply(coi$distribution, as.data.frame), list(
    data.frame(coi = 2, fraction = 1), data.frame(coi = c(2, 3), fraction = c(0.5, 0.5)),
    data.frame(coi = 2, fraction = 1)
  ))
})

test_that("pt_coi() pools Runs by their Draws, or gives one row per Run", {
  a <- local_two_genotype_runs()
  coi <- pt_coi(a)
  expect_equal(coi$median, c(2, 3, 2))
  expect_equal(as.data.frame(coi$distribution[[2]]), data.frame(coi = c(2, 3, 4), fraction = 1 / 3))
  by_run <- pt_coi(a, by_run = TRUE)
  expect_equal(names(by_run), c("run", "infection", "median", "distribution"))
  expect_equal(by_run$run, rep(c("r", "s"), each = 3))
  expect_equal(by_run$median, c(2, 2.5, 2, 2, 4, 2))
  expect_equal(as.data.frame(by_run$distribution[[5]]), data.frame(coi = 4, fraction = 1))
  expect_equal(pt_coi(a, runs = "s")$median, c(2, 4, 2))
  expect_error(pt_coi(a, by_run = NA), "by_run")
})

# The distribution list column, unnested, to compare with the oracle's long form.
coi_distribution_long <- function(coi) {
  keys <- intersect(c("run", "infection"), names(coi))
  rows <- lapply(seq_len(nrow(coi)), function(i) cbind(coi[i, keys], as.data.frame(coi$distribution[[i]])))
  do.call(rbind, rows)
}

test_that("pt_coi() matches the oracle on golden and fake Stores", {
  for (paths in list(golden_genotype_paths(), NULL)) {
    a <- if (is.null(paths)) local_two_genotype_runs() else pt_read(paths)
    o <- oracle_read(if (is.null(paths)) stats::setNames(a$runs$path, a$runs$run) else paths)
    for (by_run in c(FALSE, TRUE)) {
      keys <- c(if (by_run) "run", "infection")
      coi <- pt_coi(a, by_run = by_run)
      expected <- oracle_coi(o, by_run = by_run)
      expect_oracle_equal(coi[setdiff(names(coi), "distribution")], expected$median, keys)
      expect_oracle_equal(coi_distribution_long(coi), expected$distribution, c(keys, "coi"))
    }
    if (!is.null(paths)) pt_close(a)
  }
})
