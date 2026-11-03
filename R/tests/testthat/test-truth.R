# pt_read_truth() on hand-written simulator JSON, and the true Parent sets it implies for an Analysis.

test_that("pt_read_truth() reads edges, infections, genotypes and loci, gzipped or not", {
  path <- write_truth(
    edges = list(c(1, 2), c(2, 10)),
    nodes = list(
      truth_node(1, genotype = list(list(L1 = c(0, 1, 0)), list(L1 = c(1, 0, 0)))),
      truth_node(2, infection_time = 3.5, observation_time = 9, duration = 5.5, symptomatic = FALSE,
                 genotype = list(list(L1 = c(0, 0, 1)))),
      truth_node(10)
    ),
    loci = list(L1 = 3)
  )
  for (file in c(path, gzip_copy(path))) {
    truth <- pt_read_truth(file)
    expect_s3_class(truth, "pt_truth")
    expect_equal(as.data.frame(truth$edges), data.frame(parent = c("1", "2"), child = c("2", "10")))
    expect_equal(as.data.frame(truth$infections), data.frame(
      infection = c("1", "10", "2"), infection_time = c(0, 0, 3.5), observation_time = c(10, 10, 9),
      duration = c(10, 10, 5.5), symptomatic = c(TRUE, TRUE, FALSE)
    ))
    # Ids sort as strings. Allele k is position k of the simulator's 0/1 vector, counted from 0 as in the Store;
    # strains count from 1.
    expect_equal(as.data.frame(truth$genotypes), data.frame(
      infection = c("1", "1", "2"), strain = c(1L, 2L, 1L), locus = "L1", allele = c(1L, 0L, 2L)
    ))
    expect_equal(as.data.frame(truth$loci), data.frame(locus = "L1", n_alleles = 3L))
  }
  expect_snapshot(print(pt_read_truth(path)))
})

test_that("pt_read_truth() errors on a missing or malformed file", {
  expect_error(pt_read_truth(file.path(tempdir(), "no-such-truth.json")), "no-such-truth")
  bad <- withr::local_tempfile(fileext = ".json")
  writeLines('{"nodes": []}', bad)
  expect_error(pt_read_truth(bad), "network")
})

test_that("a true Parent set holds the Latent parent for a root or a parent outside the Analysis", {
  # Truth: 1 -> 2, 1 -> 3, 4 -> 3, 5 -> 6. The Analysis dropped 4 and 5 (a missing_50 fit).
  path <- write_truth(edges = list(c(1, 2), c(1, 3), c(4, 3), c(5, 6)), nodes = lapply(1:6, truth_node))
  truth <- pt_read_truth(path)
  store <- fake_store(infections = c("1", "2", "3", "6"))
  a <- pt_read(c(r = store))
  on.exit(pt_close(a))
  sets <- true_parent_sets(a, truth)
  expect_equal(sets$infection, c("1", "2", "3", "6"))
  expect_equal(sets$parents, list(character(), "1", "1", character()))
  expect_equal(sets$has_latent_parent, c(TRUE, FALSE, TRUE, TRUE))
})
