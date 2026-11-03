pt_acceptance <- function(a, level = c("type", "sampler"), phase = "sampling", by_run = FALSE, per_window = FALSE) {
  level <- rlang::arg_match(level)
  phase <- check_phase(phase)
  keys <- c(block_keys(by_run, per_window), "type", if (level == "sampler") "sampler")
  rows <- pt_table(a, "acceptance", phase = phase) |>
    dplyr::mutate(type = split_part(.data$sampler, "/", 1L))
  sum_block(a, rows, keys, c("accepted", "proposed"))
}

# Sum a Block's count columns within `keys`, collect, and add `rate` = first count / second count (NA when the
# second is 0). Rows are ordered by the keys, Runs in the Analysis's order.
sum_block <- function(a, rows, keys, counts) {
  rows <- rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(counts), \(x) sum(x, na.rm = TRUE)), .groups = "drop") |>
    dplyr::collect()
  for (count in counts) rows[[count]] <- as.numeric(rows[[count]])
  rows$rate <- ifelse(rows[[counts[2]]] > 0, rows[[counts[1]]] / rows[[counts[2]]], NA_real_)
  sort_keys <- as.list(rows[keys])
  if ("run" %in% keys) sort_keys$run <- match(rows$run, a$runs$run)
  rows <- rows[do.call(order, unname(sort_keys)), c(keys, counts, "rate")]
  tibble::as_tibble(rows)
}

# The leading keys every Block summary groups by.
block_keys <- function(by_run, per_window) {
  c(if (by_run) "run", if (per_window) c("window", "phase"))
}
