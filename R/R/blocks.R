pt_acceptance <- function(a, level = c("type", "sampler"), phase = "sampling", by_run = FALSE, per_window = FALSE) {
  level <- rlang::arg_match(level)
  phase <- check_phase(phase)
  keys <- c("type", if (level == "sampler") "sampler")
  rows <- pt_table(a, "acceptance", phase = phase) |>
    dplyr::mutate(type = split_part(.data$sampler, "/", 1L))
  sum_block(rows, keys, c("accepted", "proposed"))
}

# Sum a Block's count columns within `keys`, collect, and add `rate` = first count / second count (NA when the
# second is 0). Rows are ordered by the keys.
sum_block <- function(rows, keys, counts) {
  rows <- rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(counts), \(x) sum(x, na.rm = TRUE)), .groups = "drop") |>
    dplyr::collect()
  for (count in counts) rows[[count]] <- as.numeric(rows[[count]])
  rows$rate <- ifelse(rows[[counts[2]]] > 0, rows[[counts[1]]] / rows[[counts[2]]], NA_real_)
  rows <- rows[do.call(order, unname(as.list(rows[keys]))), c(keys, counts, "rate")]
  tibble::as_tibble(rows)
}
