# Scoring against simulation truth: edge recovery ranked by Edge marginals or by pairwise scores.

#' Score edge recovery against the truth
#'
#' Ranks every pair of the Analysis's Infections by its pooled Edge marginal (see [pt_edge_marginals()]) and
#' measures how well the ranking recovers the true transmission edges.
#'
#' * **Truth.** With `edges = "direct"`, a pair is true when the simulator's network holds that transmission and both
#'   ends are Infections of the Analysis. With `edges = "transitive"`, a pair is also true when a path joins them
#'   through cases outside the Analysis only (unsampled or dropped cases), so a fit that links an Infection to its
#'   nearest sampled ancestor is not penalised.
#' * **Directed** ranks every ordered pair of distinct Infections, `n (n - 1)` of them, by its Edge marginal (0 when
#'   the edge never has weight). **Undirected** ranks every unordered pair, `n (n - 1) / 2`, by the sum of both
#'   directions, true when either direction is.
#' * **Scores** are clipped to \[0, 1\] and rounded to 10 decimal places before ranking, so sums that differ only by
#'   floating-point error tie.
#' * **PR curve.** One point per distinct score, from the highest down; tied pairs enter together. At a threshold,
#'   `recall` = TP / (true pairs) and `precision` = TP / (pairs scoring at least the threshold).
#' * `average_precision` \eqn{= \sum_k (r_k - r_{k-1}) p_k} with \eqn{r_0 = 0}. `pr_auc` is the trapezoid area under
#'   the curve with the point (recall 0, precision 1) prepended. `top_k_precision` takes k = the number of true pairs
#'   (at most the universe) and is the precision over every pair scoring at least the k-th highest score, ties
#'   included.
#'
#' With no true pair among the Analysis's Infections, the three metrics are `NA` and the curve is empty.
#'
#' @inheritParams pt_table
#' @param truth A `pt_truth` from [pt_read_truth()]. Every Infection of the Analysis must be a case in it.
#' @param edges `"direct"` or `"transitive"`: which pairs count as true.
#' @param curve `FALSE` for the metrics, `TRUE` for the PR curve.
#' @return With `curve = FALSE`, a tibble with one row per `direction` (`"directed"`, `"undirected"`):
#'   `n_edge_universe`, `n_true_edges`, `n_distinct_scores`, `average_precision`, `pr_auc`, `top_k_precision`. With
#'   `curve = TRUE`, a tibble of `direction`, `threshold`, `recall`, `precision`, thresholds descending.
#' @seealso [pt_score_pairs()] for pairwise scores such as IBD estimates.
#' @export
#' @examples
#' \dontrun{
#' a <- pt_read("runs/sim_0_s7")
#' truth <- pt_read_truth("sim_0_s7.json.gz")
#' pt_score_edges(a, truth)
#' pt_score_edges(a, truth, edges = "transitive", curve = TRUE)
#' }
pt_score_edges <- function(a, truth, edges = c("direct", "transitive"), curve = FALSE) {
  edges <- rlang::arg_match(edges)
  check_flag(curve)
  score_edges(a, truth, edges, curve, runs_with_draws(a))
}

# pt_score_edges() over `runs` (already from runs_with_draws()).
score_edges <- function(a, truth, edges, curve, runs, call = rlang::caller_env()) {
  ids <- truth_infection_ids(a, truth, call)
  marginals <- pt_edge_marginals(a, runs = runs)
  n <- length(ids)
  score <- matrix(0, n, n, dimnames = list(ids, ids))
  score[cbind(match(marginals$parent, ids), match(marginals$child, ids))] <- marginals$marginal
  true_edges <- if (edges == "direct") truth_direct_edges(truth, ids) else truth_transitive_edges(truth, ids)
  true <- matrix(FALSE, n, n, dimnames = list(ids, ids))
  true[cbind(match(true_edges$parent, ids), match(true_edges$child, ids))] <- TRUE
  score_matrices(score, true, curve)
}

#' Score pairwise relatedness against the truth
#'
#' Ranks pairs of Infections by an outside score, such as an IBD estimate from Dcifer, and measures edge recovery the
#' way [pt_score_edges()] does for Edge marginals: the baseline a transmission fit should beat.
#'
#' Each unordered pair keeps its largest score; pairs with no score score 0; rows naming an id outside `infections`,
#' pairing an id with itself, or with a missing score are ignored. A pair is true when the truth holds a direct
#' transmission between them, in either direction. Scores are clipped and rounded as in [pt_score_edges()].
#'
#' A pairwise score has no direction, so the directed row ranks every ordered pair with the score of its unordered
#' pair, each pair scoring the same both ways; the undirected row ranks the unordered pairs.
#'
#' @param scores A data frame with columns `a`, `b` (Infection ids) and `score`.
#' @param truth A `pt_truth` from [pt_read_truth()].
#' @param infections The Infection ids to score over, such as `pt_infections(a)$infection`.
#' @return A tibble like [pt_score_edges()]'s: one row per `direction` with `n_edge_universe`, `n_true_edges`,
#'   `n_distinct_scores`, `average_precision`, `pr_auc`, `top_k_precision`.
#' @export
#' @examples
#' \dontrun{
#' truth <- pt_read_truth("sim_0_s7.json.gz")
#' ibd <- read.csv("dcifer_ibd.csv", colClasses = c(a = "character", b = "character"))
#' pt_score_pairs(ibd, truth, pt_infections(a)$infection)
#' }
pt_score_pairs <- function(scores, truth, infections) {
  call <- rlang::current_env()
  check_truth(truth, call)
  if (!is.data.frame(scores) || !all(c("a", "b", "score") %in% names(scores))) {
    cli::cli_abort("{.arg scores} must be a data frame with columns {.field a}, {.field b} and {.field score}.",
                   call = call)
  }
  if (!is.numeric(scores$score)) cli::cli_abort("{.field score} must be numeric.", call = call)
  ids <- sort(unique(as.character(infections)), method = "radix")
  a <- as.character(scores$a)
  b <- as.character(scores$b)
  keep <- a %in% ids & b %in% ids & a != b & !is.na(scores$score)
  n <- length(ids)
  i <- match(a[keep], ids)
  j <- match(b[keep], ids)
  s <- scores$score[keep]
  # Largest score per unordered pair, placed both ways.
  score <- matrix(-Inf, n, n, dimnames = list(ids, ids))
  for (k in order(s)) {
    score[i[k], j[k]] <- s[k]
    score[j[k], i[k]] <- s[k]
  }
  score[score == -Inf] <- 0
  true_edges <- truth_direct_edges(truth, ids)
  true <- matrix(FALSE, n, n, dimnames = list(ids, ids))
  true[cbind(match(true_edges$parent, ids), match(true_edges$child, ids))] <- TRUE
  # The undirected score of a pair is its score, not the sum of both directions.
  score_matrices(score, true, curve = FALSE, undirected = score)
}

# Directed and undirected ranking metrics from an n x n score matrix (parent rows, child columns) and truth
# matrix. `undirected`: the unordered pairs' scores (upper triangle); default the sum of both directions.
score_matrices <- function(score, true, curve, undirected = score + t(score)) {
  off <- row(score) != col(score)
  upper <- row(score) < col(score)
  directed <- rank_scores(normalise_scores(score[off]), true[off])
  undirected <- rank_scores(normalise_scores(undirected[upper]), (true | t(true))[upper])
  if (curve) {
    return(dplyr::bind_rows(
      tibble::tibble(direction = rep("directed", nrow(directed$curve)), directed$curve),
      tibble::tibble(direction = rep("undirected", nrow(undirected$curve)), undirected$curve)
    ))
  }
  dplyr::bind_rows(
    tibble::tibble(direction = "directed", !!!directed$metrics),
    tibble::tibble(direction = "undirected", !!!undirected$metrics)
  )
}

# Clip to [0, 1], round to 10 places: near-equal sums tie.
normalise_scores <- function(s) round(pmin(1, pmax(0, s)), 10)

# Ranking metrics of scores against true flags: one PR point per distinct score, descending, ties entering together.
rank_scores <- function(score, true) {
  n_true <- sum(true)
  thresholds <- sort(unique(score), decreasing = TRUE)
  metrics <- list(n_edge_universe = length(score), n_true_edges = as.integer(n_true),
                  n_distinct_scores = length(thresholds), average_precision = NA_real_, pr_auc = NA_real_,
                  top_k_precision = NA_real_)
  empty <- tibble::tibble(threshold = numeric(), recall = numeric(), precision = numeric())
  if (n_true == 0) return(list(metrics = metrics, curve = empty))

  at <- match(score, thresholds)
  tp <- cumsum(tabulate(at[true], length(thresholds)))
  reached <- cumsum(tabulate(at, length(thresholds)))
  recall <- tp / n_true
  precision <- tp / reached
  metrics$average_precision <- sum(diff(c(0, recall)) * precision)
  metrics$pr_auc <- sum(diff(c(0, recall)) * (precision + c(1, precision[-length(precision)])) / 2)
  k <- min(n_true, length(score))
  kth <- sort(score, decreasing = TRUE)[[k]]
  metrics$top_k_precision <- sum(true[score >= kth]) / sum(score >= kth)
  list(metrics = metrics, curve = tibble::tibble(threshold = thresholds, recall = recall, precision = precision))
}

#' Score Latent genotypes against the truth
#'
#' Compares each Infection's true Alleles at each Locus it was genotyped at with its posterior genotype, called
#' Allele by Allele: an Allele is called present when its pooled presence (see [pt_genotype_presence()]) is at least
#' 0.5.
#'
#' **Allele numbering.** The simulator numbers a Locus's Alleles over the whole simulated population; an input built
#' from a subset of cases may keep only the Alleles its cases show and renumber them. Scoring maps the truth into the
#' input's numbering by this rule, the convention of the thesis inputs:
#' * A Locus whose Allele count in the input equals the simulator's keeps the simulator's numbering.
#' * Otherwise the input's Alleles at the Locus are the simulator's Alleles shown in the observed genotype (the
#'   truth's `observed`) of any of the Analysis's Infections, in simulator order: the smallest is Allele 0.
#' * When that count differs from the input's Allele count, the simulator has fewer Alleles than the input, or the
#'   truth lacks the Locus, scoring errors and names the Locus.
#'
#' A true Allele that no Infection of the Analysis shows has no number in the input: the model cannot represent it.
#' It is counted in `unrepresentable` and left out of every other column.
#'
#' @inheritParams pt_score_edges
#' @return A tibble with one row per (infection, locus) that has a Latent genotype, sorted by `infection` and
#'   `locus` as strings: `true_alleles` and `called_alleles` (list columns of integer Allele numbers in the input's
#'   numbering, 0-based as in [pt_alleles()]), `true_positives`, `false_positives`, `false_negatives`,
#'   `unrepresentable` (the true Alleles left out) and `exact` (whether the called Alleles are the true ones).
#' @seealso [pt_score_coi()]; [pt_read_truth()].
#' @export
#' @examples
#' \dontrun{
#' a <- pt_read("runs/sim_0_s7")
#' truth <- pt_read_truth("sim_0_s7.json.gz")
#' scores <- pt_score_genotypes(a, truth)
#' mean(scores$exact)
#' }
pt_score_genotypes <- function(a, truth) {
  ids <- truth_infection_ids(a, truth)
  map <- truth_allele_map(a, truth, ids)
  presence <- pt_genotype_presence(a)
  keys <- unique(presence[c("infection", "locus")])
  key <- paste(keys$infection, keys$locus, sep = "\r")

  called <- presence[presence$presence >= 0.5, ]
  called <- split(as.integer(called$allele), factor(paste(called$infection, called$locus, sep = "\r"), levels = key))

  true <- unique(truth$genotypes[truth$genotypes$infection %in% ids, c("infection", "locus", "allele")])
  true$key <- paste(true$infection, true$locus, sep = "\r")
  true <- true[true$key %in% key, ]
  true$mapped <- map$allele[match(paste(true$locus, true$allele, sep = "\r"), paste(map$locus, map$simulator, sep = "\r"))]
  unrepresentable <- as.vector(table(factor(true$key[is.na(true$mapped)], levels = key)))
  true <- true[!is.na(true$mapped), ]
  true <- split(true$mapped, factor(true$key, levels = key))
  true <- lapply(true, function(x) sort(as.integer(x)))

  hits <- mapply(function(t, c) sum(c %in% t), true, called)
  tibble::tibble(
    infection = keys$infection, locus = keys$locus,
    true_alleles = unname(true), called_alleles = unname(called),
    true_positives = as.integer(hits),
    false_positives = as.integer(lengths(called) - hits),
    false_negatives = as.integer(lengths(true) - hits),
    unrepresentable = as.integer(unrepresentable),
    exact = unname(mapply(identical, true, called))
  )
}

#' Score COI against the truth
#'
#' Compares each Infection's posterior COI median (from [pt_coi()]) with its identifiable true COI: the most distinct
#' true strain genotypes at any Locus the Infection was genotyped at. Strains that differ only at Loci the Infection
#' was not genotyped at cannot be told apart, so the identifiable COI can be below the number of strains. Strain
#' genotypes are compared in the simulator's numbering, so no Allele mapping is needed.
#'
#' @inheritParams pt_score_edges
#' @return A tibble with one row per Infection genotyped at some Locus, sorted by `infection` as strings: `median`
#'   (the posterior COI median, possibly halfway between two values), `true_coi` (identifiable), `strains` (the
#'   number of true strains) and `error` (`median - true_coi`).
#' @seealso [pt_score_genotypes()]; [pt_read_truth()].
#' @export
#' @examples
#' \dontrun{
#' a <- pt_read("runs/sim_0_s7")
#' truth <- pt_read_truth("sim_0_s7.json.gz")
#' scores <- pt_score_coi(a, truth)
#' mean(scores$error == 0)
#' }
pt_score_coi <- function(a, truth) {
  ids <- truth_infection_ids(a, truth)
  coi <- pt_coi(a)
  genotyped <- pt_table(a, "latent_genotypes", phase = "sampling") |>
    dplyr::distinct(.data$infection, .data$locus) |>
    dplyr::collect()
  g <- truth$genotypes[truth$genotypes$infection %in% ids, ]
  g <- g[order(g$infection, g$strain, g$locus, g$allele, method = "radix"), ]
  # One genotype string per (infection, strain, locus), then distinct strings per (infection, locus).
  strain_key <- paste(g$infection, g$strain, g$locus, sep = "\r")
  first <- !duplicated(strain_key)
  per_strain <- tibble::tibble(
    infection = g$infection[first], locus = g$locus[first],
    genotype = vapply(split(g$allele, factor(strain_key, levels = unique(strain_key))), paste, "", collapse = ",")
  )
  per_strain <- per_strain[paste(per_strain$infection, per_strain$locus, sep = "\r") %in%
                             paste(genotyped$infection, genotyped$locus, sep = "\r"), ]
  distinct <- unique(per_strain)
  at_locus <- tapply(distinct$genotype, paste(distinct$infection, distinct$locus, sep = "\r"), length)
  infection_of <- sub("\r.*$", "", names(at_locus))
  true_coi <- tapply(as.vector(at_locus), infection_of, max)[coi$infection]
  strains <- tapply(g$strain, g$infection, function(s) length(unique(s)))[coi$infection]
  out <- tibble::tibble(
    infection = coi$infection, median = coi$median,
    true_coi = as.integer(ifelse(is.na(true_coi), 0L, true_coi)),
    strains = as.integer(ifelse(is.na(strains), 0L, strains))
  )
  out$error <- out$median - out$true_coi
  out
}

#' Score Parent sets and the MAP graph against the truth
#'
#' Compares each Infection's pooled Parent sets (see [pt_parent_sets_pooled()]) with its true Parent set: its true
#' parents among the Analysis's Infections, plus the Latent parent when it is a root of the simulated outbreak or
#' when any of its true parents is outside the Analysis (an unsampled or dropped case counts as the Latent parent).
#' Two sets are equal when they hold the same observed Infections and both hold, or both lack, the Latent parent.
#'
#' * `pt_score_parent_sets()`: the Credible collection at `level`, which admits every set tied in weight with the
#'   last set needed to reach the level, so its size and coverage never depend on how ties are broken.
#'   `credible_size` counts its sets; `covered` is whether it holds the true set; `top_probability` is the weight of
#'   the MAP Parent set.
#' * `pt_score_map()`: whether the MAP Parent set (rank 1, ties broken by the members list as in
#'   [pt_parent_sets_pooled()]) is the true set. The share of Infections whose MAP Parent set is right is
#'   `mean(correct)`.
#'
#' Sets are pooled over the sampling Draws of every Run; Runs with no sampling Draws are skipped with a message.
#'
#' @inheritParams pt_score_edges
#' @param level Level of the Credible collection, in (0, 1].
#' @return A tibble with one row per Infection, sorted by `infection` as strings. `pt_score_parent_sets()`:
#'   `credible_size`, `covered`, `top_probability`. `pt_score_map()`: the MAP Parent set's `members` (a list of
#'   observed parent ids in string order), `has_latent_parent` and `weight`; the true set's `true_members` and
#'   `true_has_latent_parent`; and `correct`.
#' @seealso [pt_parent_sets_pooled()], [pt_map_graph()]; [pt_read_truth()].
#' @export
#' @examples
#' \dontrun{
#' a <- pt_read("runs/sim_0_s7")
#' truth <- pt_read_truth("sim_0_s7.json.gz")
#' sets <- pt_score_parent_sets(a, truth)
#' mean(sets$covered)
#' mean(pt_score_map(a, truth)$correct)
#' }
pt_score_parent_sets <- function(a, truth, level = 0.8) {
  check_credible_level(level)
  score_parent_sets(a, truth, level, runs_with_draws(a))
}

# pt_score_parent_sets() over `runs` (already from runs_with_draws()).
score_parent_sets <- function(a, truth, level, runs, call = rlang::caller_env()) {
  true <- true_parent_sets(a, truth, call)
  sets <- pooled_sets(a, runs, FALSE, level)
  credible <- sets[sets$credible, ]
  hit <- set_key(credible$infection, credible$members, credible$has_latent_parent) %in%
    set_key(true$infection, true$parents, true$has_latent_parent)
  ids <- sort(true$infection, method = "radix")
  map <- sets[sets$rank == 1L, ]
  top <- map$weight[match(ids, map$infection)]
  tibble::tibble(
    infection = ids,
    credible_size = as.integer(table(factor(credible$infection, levels = ids))),
    covered = ids %in% credible$infection[hit],
    top_probability = ifelse(is.na(top), 0, top)
  )
}

#' @rdname pt_score_parent_sets
#' @export
pt_score_map <- function(a, truth) {
  true <- true_parent_sets(a, truth)
  sets <- pooled_sets(a, runs_with_draws(a), FALSE, 0.8)
  map <- sets[sets$rank == 1L, ]
  ids <- sort(true$infection, method = "radix")
  at <- match(ids, map$infection)
  truth_at <- match(ids, true$infection)
  out <- tibble::tibble(
    infection = ids,
    members = map$members[at],
    has_latent_parent = map$has_latent_parent[at],
    weight = map$weight[at],
    true_members = true$parents[truth_at],
    true_has_latent_parent = true$has_latent_parent[truth_at]
  )
  out$correct <- !is.na(at) &
    set_key(out$infection, out$members, out$has_latent_parent) ==
      set_key(out$infection, out$true_members, out$true_has_latent_parent)
  out
}

# One string per Parent set: infection, observed members (in the order given), Latent parent flag.
set_key <- function(infection, members, has_latent_parent) {
  paste(infection, vapply(members, paste, "", collapse = "\x1f"), has_latent_parent, sep = "\r")
}
