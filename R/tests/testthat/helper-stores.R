# Store fixtures shared by every test file. Later tickets add helper-fake-store.R (fake_store()) and
# helper-oracle.R next to this one.

# Path of an installed golden Store: "single-chain" or "replica-exchange". Both share one input, so they open
# together as one Analysis. Read only; copy with local_golden_copy() before changing anything.
golden_store <- function(name = "single-chain") {
  system.file("extdata", "golden", name, package = "plasmotrackr", mustWork = TRUE)
}

# Copy a golden Store into a temporary directory removed when `.local_envir` exits. The copy keeps the Store's
# directory name, so its Run label is `name`. Returns the copy's path.
local_golden_copy <- function(name = "single-chain", .local_envir = parent.frame()) {
  parent <- withr::local_tempdir(.local_envir = .local_envir)
  file.copy(golden_store(name), parent, recursive = TRUE)
  file.path(parent, name)
}

# Overwrite fields of a Store's manifest.json, such as edit_manifest(path, null_model = TRUE).
edit_manifest <- function(path, ...) {
  file <- file.path(path, "manifest.json")
  manifest <- jsonlite::read_json(file)
  manifest[names(list(...))] <- list(...)
  jsonlite::write_json(manifest, file, auto_unbox = TRUE, digits = NA, pretty = TRUE)
  invisible(path)
}

# Make a Store look like a Run killed mid-commit: for every table, write a Window after the last closed one with no
# `_closed/` marker, and a `.in-flight.parquet` file. Rows of both are copies of the last closed Window's rows with
# iterations moved past the Run's end (and, in parent_set_members, set ids no real set has), so a reader that reads
# either returns rows the original Store does not have.
add_unclosed_window <- function(path) {
  markers <- list.files(file.path(path, "_closed"), full.names = TRUE)
  last <- max(vapply(markers, function(f) jsonlite::read_json(f)$window, numeric(1)))
  unmarked <- sprintf("window-%06d.parquet", last + 1)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  q <- function(x) as.character(DBI::dbQuoteString(con, x))
  for (table in setdiff(list.dirs(path, full.names = FALSE, recursive = FALSE), c("_closed", "structure"))) {
    files <- list.files(file.path(path, table), pattern = "^window-\\d+\\.parquet$", full.names = TRUE)
    source <- files[which.max(as.numeric(gsub("\\D", "", basename(files))))]
    columns <- DBI::dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM read_parquet(%s)", q(source)))$column_name
    replace <- c(
      iteration = "iteration + 1000", window = sprintf("%d::BIGINT", as.integer(last + 1)),
      set_id = "set_id + 1000"
    )
    replace <- replace[names(replace) %in% columns]
    select <- if (length(replace)) {
      sprintf("SELECT * REPLACE (%s)", paste(replace, "AS", DBI::dbQuoteIdentifier(con, names(replace)), collapse = ", "))
    } else {
      "SELECT *"
    }
    for (target in file.path(path, table, c(unmarked, ".in-flight.parquet"))) {
      DBI::dbExecute(con, sprintf(
        "COPY (%s FROM read_parquet(%s)) TO %s (FORMAT parquet)", select, q(source), q(target)
      ))
    }
  }
  invisible(path)
}
