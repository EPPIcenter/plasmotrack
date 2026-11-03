# Keep DuckDB to one thread in tests: CI runners and the development machine are shared.
withr::local_options(plasmotrackr.threads = 1, .local_envir = testthat::teardown_env())
