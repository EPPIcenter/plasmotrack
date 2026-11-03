# One DuckDB connection per Analysis, opened with the engine options and closed by pt_close() or garbage collection.

engine_option <- function(value, option, default = NULL) {
  if (!is.null(value)) return(value)
  getOption(option, default)
}

# DuckDB 1.5.6's common_subplan optimizer returns wrong rows when one query left-joins two aggregations of the same
# joined subquery grouped by different keys (pt_net_flow()'s out- and in-edges): it is disabled.
open_engine <- function(memory_limit = NULL, temp_directory = NULL, threads = NULL) {
  config <- list(temp_directory = as.character(temp_directory), disabled_optimizers = "common_subplan")
  if (!is.null(memory_limit)) config$memory_limit <- as.character(memory_limit)
  if (!is.null(threads)) config$threads <- as.character(threads)
  drv <- duckdb::duckdb(shared_home = FALSE, config = config)
  engine <- new.env(parent = emptyenv())
  engine$drv <- drv
  engine$con <- DBI::dbConnect(drv)
  reg.finalizer(engine, close_engine, onexit = TRUE)
  engine$settings <- as.list(DBI::dbGetQuery(engine$con, paste(
    "SELECT current_setting('memory_limit') AS memory_limit,",
    "current_setting('temp_directory') AS temp_directory, current_setting('threads') AS threads"
  )))
  engine
}

close_engine <- function(engine) {
  if (!is.null(engine$con)) {
    DBI::dbDisconnect(engine$con)
    duckdb::duckdb_shutdown(engine$drv)
    engine$con <- NULL
    engine$drv <- NULL
  }
  invisible(NULL)
}

analysis_con <- function(a, call = rlang::caller_env()) {
  con <- a$engine$con
  if (is.null(con)) cli::cli_abort("This Analysis is closed. Open it again with {.fn pt_read}.", call = call)
  con
}

sql_string <- function(con, x) as.character(DBI::dbQuoteString(con, x))

sql_list <- function(con, x) paste0("[", paste(sql_string(con, x), collapse = ", "), "]")

# Rows of a query as a tibble, for small reads made eagerly (markers, structure, Draw counts).
query <- function(con, sql) tibble::as_tibble(DBI::dbGetQuery(con, sql))

# Text shown verbatim in a cli message: braces doubled so cli does not interpolate them.
cli_escape <- function(x) gsub("}", "}}", gsub("{", "{{", x, fixed = TRUE), fixed = TRUE)
