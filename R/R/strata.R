# Strata: `by =` names one epi column; each Infection's value is its stratum. Infections with an NA value, or missing
# from epi, form the NA stratum, which every output keeps.
#
#   check_by(a, by)                 NULL, or one epi column name (not `infection`).
#   infection_strata(a, by)         eager: infection, stratum for every Infection of the Analysis (factor -> character).
#   strata_keys(by)                 "stratum" with `by`, else nothing: the extra grouping key of a stratified table.
#   with_strata(values, a, by, column)  lazy: `values` with `stratum` joined on `column` (an Infection id column);
#                                   unchanged when `by` is NULL.
#   stratum_levels(stratum)         sorted distinct strata, NA last.

check_by <- function(a, by, call = rlang::caller_env()) {
  if (is.null(by)) return(invisible(NULL))
  if (!is.character(by) || length(by) != 1 || is.na(by)) {
    cli::cli_abort("{.arg by} must be one epi column name, as a string.", call = call)
  }
  if (is.null(a$epi)) {
    cli::cli_abort(c("{.arg by} needs epi data, and this Analysis has none.",
                     i = "Give {.code pt_read(epi = )} a data frame with an {.field infection} column."), call = call)
  }
  columns <- setdiff(names(a$epi), "infection")
  if (!by %in% columns) {
    cli::cli_abort(c("{.arg by} must name an epi column: {.val {cli_escape(by)}} is not one.",
                     i = "Epi columns: {.field {cli_escape(columns)}}."), call = call)
  }
  value <- a$epi[[by]]
  if (!(is.atomic(value) || is.factor(value)) || is.complex(value) || is.raw(value)) {
    cli::cli_abort("Epi column {.field {cli_escape(by)}} must hold strings, numbers, logicals or dates.", call = call)
  }
  invisible(NULL)
}

infection_strata <- function(a, by) {
  value <- a$epi[[by]]
  if (is.factor(value)) value <- as.character(value)
  tibble::tibble(infection = a$infections$infection, stratum = value[match(a$infections$infection, a$epi$infection)])
}

strata_keys <- function(by) if (is.null(by)) character() else "stratum"

with_strata <- function(values, a, by, column = "infection") {
  if (is.null(by)) return(values)
  strata <- dbplyr::copy_inline(analysis_con(a), infection_strata(a, by))
  dplyr::left_join(values, strata, by = stats::setNames("infection", column))
}

stratum_levels <- function(stratum) sort(unique(stratum), na.last = TRUE, method = "radix")
