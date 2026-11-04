# Simulation truth: the simulator's JSON read into tibbles, and the true Parent sets it implies for an Analysis.

#' Read a simulation's truth
#'
#' Reads the simulator's JSON (the `sim_*.json` or `sim_*.json.gz` next to the input it wrote), gzipped or not. The
#' truth covers every simulated case; it is restricted to an Analysis's Infections only when scoring, so one
#' `pt_truth` serves both the full and the partial fit of a simulation.
#'
#' @param path The truth JSON file, plain or gzipped.
#' @return A `pt_truth`: an S3 list of tibbles, with every id a string and rows sorted by id as strings.
#'   * `edges`: `parent`, `child`, one row per transmission in the simulator's `network`.
#'   * `infections`: `infection`, `infection_time`, `observation_time`, `duration`, `symptomatic`, one row per
#'     simulated case (`nodes`).
#'   * `genotypes`: `infection`, `strain`, `locus`, `allele`, one row per Allele each strain carries (long format,
#'     from `latent_genotype`). Alleles are in the simulator's numbering: Allele k is position k of the simulator's
#'     vector, counted from 0 as in the Store. An input may renumber them; scoring maps them back. Strains count
#'     from 1.
#'   * `loci`: `locus`, `n_alleles` (the simulator's Allele count).
#'   * `observed`: `infection`, `locus`, `allele`, one row per Allele each case's observed genotype shows
#'     (`observed_genotype`, with any sequencing errors the simulator added), in the simulator's numbering. Scoring
#'     reads the input's Allele numbering off it (see [pt_score_genotypes()]). Empty when the file has none.
#' @seealso [pt_score_edges()], [pt_score_pairs()], [pt_score_parent_sets()], [pt_score_map()],
#'   [pt_score_genotypes()], [pt_score_coi()]
#' @export
#' @examples
#' \dontrun{
#' truth <- pt_read_truth("sim_0_s7.json.gz")
#' truth$edges
#' }
pt_read_truth <- function(path) {
  call <- rlang::current_env()
  if (!is.character(path) || length(path) != 1 || is.na(path)) {
    cli::cli_abort("{.arg path} must be one file path.", call = call)
  }
  if (!file.exists(path)) {
    cli::cli_abort("Truth file {.file {cli_escape(path)}} does not exist.", call = call)
  }
  # file() reads gzipped files transparently.
  json <- tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) {
      cli::cli_abort(c("Can't parse truth file {.file {cli_escape(path)}}.", x = cli_escape(conditionMessage(e))),
                     call = call)
    }
  )
  missing <- setdiff(c("network", "nodes"), names(json))
  if (length(missing)) {
    cli::cli_abort("Truth file {.file {cli_escape(path)}} has no {.field {missing}}: is it a simulator truth JSON?",
                   call = call)
  }

  field <- function(xs, name, as) vapply(xs, function(x) if (is.null(x[[name]])) as(NA) else as(x[[name]]), as(NA))
  edges <- tibble::tibble(parent = field(json$network, "src", as.character),
                          child = field(json$network, "dest", as.character))
  nodes <- json$nodes
  infections <- tibble::tibble(
    infection = field(nodes, "id", as.character),
    infection_time = field(nodes, "infection_time", as.numeric),
    observation_time = field(nodes, "observation_time", as.numeric),
    duration = field(nodes, "infection_duration", as.numeric),
    symptomatic = field(nodes, "symptomatic", as.logical)
  )
  genotypes <- truth_genotypes(nodes)
  observed <- truth_observed(nodes)
  loci <- tibble::tibble(locus = field(json$loci, "locus", as.character),
                         n_alleles = field(json$loci, "num_alleles", as.integer))

  by <- function(df, ...) df[do.call(order, c(unname(list(...)), method = "radix")), ]
  structure(list(
    edges = by(edges, edges$parent, edges$child),
    infections = by(infections, infections$infection),
    genotypes = by(genotypes, genotypes$infection, genotypes$strain, genotypes$locus, genotypes$allele),
    loci = by(loci, loci$locus),
    observed = by(observed, observed$infection, observed$locus, observed$allele)
  ), class = "pt_truth")
}

# One row per (infection, locus, allele shown) from the nodes' observed_genotype: a bitstring per Locus, or a 0/1
# vector.
truth_observed <- function(nodes) {
  rows <- lapply(nodes, function(node) {
    loci <- node$observed_genotype
    alleles <- lapply(loci, function(g) {
      bits <- unlist(g$genotype)
      if (is.character(bits)) bits <- strsplit(paste(bits, collapse = ""), "", fixed = TRUE)[[1]]
      which(as.integer(bits) == 1L) - 1L
    })
    n <- lengths(alleles)
    list(infection = rep(as.character(node$id), sum(n)),
         locus = rep(vapply(loci, function(g) as.character(g$locus), character(1)), n),
         allele = unlist(alleles))
  })
  tibble::tibble(
    infection = as.character(unlist(lapply(rows, `[[`, "infection"))),
    locus = as.character(unlist(lapply(rows, `[[`, "locus"))),
    allele = as.integer(unlist(lapply(rows, `[[`, "allele")))
  )
}

# One row per (infection, strain, locus, allele present) from the nodes' latent_genotype.
truth_genotypes <- function(nodes) {
  rows <- lapply(nodes, function(node) {
    strains <- node$latent_genotype
    per_strain <- lapply(seq_along(strains), function(s) {
      alleles <- lapply(strains[[s]], function(g) which(as.integer(unlist(g$genotype)) == 1L) - 1L)
      n <- lengths(alleles)
      list(strain = rep(s, sum(n)),
           locus = rep(vapply(strains[[s]], function(g) as.character(g$locus), character(1)), n),
           allele = unlist(alleles))
    })
    n <- sum(vapply(per_strain, function(p) length(p$strain), integer(1)))
    list(infection = rep(as.character(node$id), n),
         strain = unlist(lapply(per_strain, `[[`, "strain")),
         locus = unlist(lapply(per_strain, `[[`, "locus")),
         allele = unlist(lapply(per_strain, `[[`, "allele")))
  })
  tibble::tibble(
    infection = as.character(unlist(lapply(rows, `[[`, "infection"))),
    strain = as.integer(unlist(lapply(rows, `[[`, "strain"))),
    locus = as.character(unlist(lapply(rows, `[[`, "locus"))),
    allele = as.integer(unlist(lapply(rows, `[[`, "allele")))
  )
}

#' @export
print.pt_truth <- function(x, ...) {
  roots <- setdiff(x$infections$infection, x$edges$child)
  cat(sprintf("<pt_truth> %d Infections (%d roots), %d edges, %d Loci, %d strain-Allele rows\n",
              nrow(x$infections), length(roots), nrow(x$edges), nrow(x$loci), nrow(x$genotypes)))
  invisible(x)
}

check_truth <- function(truth, call = rlang::caller_env()) {
  if (!inherits(truth, "pt_truth")) {
    cli::cli_abort("{.arg truth} must be a {.cls pt_truth} from {.fn pt_read_truth}.", call = call)
  }
}

# The Analysis's Infections, checked against the truth: each must be a simulated case.
truth_infection_ids <- function(a, truth, call = rlang::caller_env()) {
  check_truth(truth, call)
  ids <- a$infections$infection
  absent <- setdiff(ids, truth$infections$infection)
  if (length(absent)) {
    cli::cli_abort(c(
      "{length(absent)} Infection{?s} of the Analysis {?is/are} not in the truth: is it this simulation's truth?",
      rlang::set_names(cli_escape(utils::head(absent, 10)), "x")
    ), call = call)
  }
  ids
}

# True edges among `ids` (both ends among them, no self-edge): parent, child.
truth_direct_edges <- function(truth, ids) {
  e <- truth$edges
  e[e$parent %in% ids & e$child %in% ids & e$parent != e$child, c("parent", "child")]
}

# Observed pairs joined by a path whose inner cases are all outside `ids` (the direct edges included): parent,
# child. A depth-first walk from each observed case that stops at observed cases.
truth_transitive_edges <- function(truth, ids) {
  e <- truth$edges[truth$edges$parent != truth$edges$child, ]
  children <- split(e$child, e$parent)
  observed <- stats::setNames(rep(TRUE, length(ids)), ids)
  reach <- lapply(ids, function(from) {
    found <- character()
    seen <- character()
    stack <- children[[from]]
    while (length(stack)) {
      node <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      if (node %in% seen) next
      seen <- c(seen, node)
      if (isTRUE(observed[node])) {
        found <- c(found, node)
      } else {
        stack <- c(stack, children[[node]])
      }
    }
    setdiff(found, from)
  })
  tibble::tibble(parent = rep(ids, lengths(reach)), child = as.character(unlist(reach)))
}

# Per Infection of the Analysis: its true observed parents (in id string order) and whether the true Parent set
# holds the Latent parent: a root, or a true parent outside the Analysis (a dropped case counts as the Latent
# parent). Comparable with a fitted set's observed members and has_latent_parent.
true_parent_sets <- function(a, truth, call = rlang::caller_env()) {
  ids <- truth_infection_ids(a, truth, call)
  e <- truth$edges[truth$edges$parent != truth$edges$child, ]
  all_parents <- lapply(ids, function(id) e$parent[e$child == id])
  tibble::tibble(
    infection = ids,
    parents = lapply(all_parents, function(p) sort(unique(p[p %in% ids]), method = "radix")),
    has_latent_parent = vapply(all_parents, function(p) length(p) == 0 || any(!p %in% ids), logical(1))
  )
}

# The input's Allele numbering of the truth at each Locus of the Analysis: locus, simulator (the simulator's Allele),
# allele (the input's), both 0-based. A Locus whose input kept the simulator's Allele count is numbered as the
# simulator numbers it. A Locus with fewer Alleles in the input was renumbered: its Alleles are the simulator's
# Alleles shown in any observed genotype of the Analysis's Infections (`ids`), in simulator order. Errors naming
# every Locus where that count differs from the input's, where the simulator has fewer Alleles than the input, or
# that the truth lacks. Simulator Alleles with no row are
# unrepresentable.
truth_allele_map <- function(a, truth, ids, call = rlang::caller_env()) {
  loci <- a$loci
  simulator <- truth$loci$n_alleles[match(loci$locus, truth$loci$locus)]
  absent <- loci$locus[is.na(simulator)]
  if (length(absent)) {
    cli::cli_abort("{cli::qty(absent)}Locus{?/Loci} {.val {cli_escape(absent)}} of the Analysis {?is/are} not in the truth.",
                   call = call)
  }
  shown <- truth$observed[truth$observed$infection %in% ids, c("locus", "allele")]
  parts <- lapply(seq_len(nrow(loci)), function(i) {
    n <- loci$n_alleles[[i]]
    from <- if (simulator[[i]] == n) seq_len(n) - 1L else sort(unique(shown$allele[shown$locus == loci$locus[[i]]]))
    tibble::tibble(locus = rep(loci$locus[[i]], length(from)), simulator = as.integer(from),
                   allele = seq_along(from) - 1L)
  })
  mapped <- vapply(parts, nrow, 1L)
  wrong <- mapped != loci$n_alleles | simulator < loci$n_alleles
  if (any(wrong)) {
    detail <- sprintf("%s: the input has %d Allele%s, the simulator %d, of which the Analysis's Infections show %d",
                      loci$locus[wrong], loci$n_alleles[wrong], ifelse(loci$n_alleles[wrong] == 1, "", "s"),
                      simulator[wrong], mapped[wrong])
    cli::cli_abort(c(
      "Can't map the truth's Alleles to the input's at {cli::qty(sum(wrong))}Locus{?/Loci} {.val {cli_escape(loci$locus[wrong])}}.",
      rlang::set_names(cli_escape(detail), "x"),
      i = "A renumbered Locus must hold exactly the simulator's Alleles shown in the Analysis's Infections."
    ), call = call)
  }
  dplyr::bind_rows(parts)
}
