# Builds the `nodes` and `loci` structures the model consumes.
#
# The central invariant: for a given locus, bit position k of every genotype
# string, entry k of `allele_freqs`, and `num_alleles` all refer to the same
# allele. The original code established this ordering twice, independently
# (once via levels(as.factor()) for the bitstrings, once via dcifer's internal
# sort for the frequencies), and trusted the two to agree. Here a single
# ordering is computed once and both consumers are indexed by allele name.

build_locus_table <- function(dat) {
    dat |>
        dplyr::distinct(locus, allele) |>
        dplyr::arrange(locus, allele) |>
        dplyr::group_by(locus) |>
        dplyr::summarize(allele_levels = list(allele), .groups = "drop") |>
        dplyr::mutate(num_alleles = lengths(allele_levels))
}

# Attaches population allele frequencies, indexed by allele name so the vector
# order follows `allele_levels` no matter how dcifer ordered its output.
attach_allele_freqs <- function(locus_tbl, afreq) {
    missing <- setdiff(locus_tbl$locus, names(afreq))
    if (length(missing) > 0) {
        stop(sprintf(
            "Allele frequencies missing for locus/loci: %s",
            paste(utils::head(missing, 5), collapse = ", ")
        ), call. = FALSE)
    }

    locus_tbl$allele_freqs <- lapply(seq_len(nrow(locus_tbl)), function(k) {
        locus <- locus_tbl$locus[k]
        levels_k <- locus_tbl$allele_levels[[k]]
        freqs <- afreq[[locus]][levels_k]
        if (anyNA(freqs)) {
            stop(sprintf("Allele frequency lookup failed at locus %s.", locus), call. = FALSE)
        }
        # A zero frequency makes the model's log-likelihood -Inf, so nudge the
        # floor up and renormalise rather than emitting an unusable input.
        freqs <- pmax(as.numeric(freqs), .Machine$double.eps)
        unname(freqs / sum(freqs))
    })
    locus_tbl
}

# Returns a samples x loci character matrix of presence/absence bitstrings.
# Absent sample-locus combinations become all-zero strings, which the model
# reads as "not observed at this locus".
encode_genotypes <- function(dat, locus_tbl, sample_ids) {
    loci <- locus_tbl$locus
    by_locus <- split(dat, factor(dat$locus, levels = loci))

    out <- matrix("", nrow = length(sample_ids), ncol = length(loci),
        dimnames = list(sample_ids, loci)
    )

    for (k in seq_along(loci)) {
        levels_k <- locus_tbl$allele_levels[[k]]
        sub <- by_locus[[k]]
        bits <- matrix(0L, nrow = length(sample_ids), ncol = length(levels_k))

        if (!is.null(sub) && nrow(sub) > 0) {
            si <- match(sub$sample_id, sample_ids)
            ai <- match(sub$allele, levels_k)
            ok <- !is.na(si) & !is.na(ai)
            bits[cbind(si[ok], ai[ok])] <- 1L
        }

        # asplit + paste0 concatenates across columns for all rows at once,
        # avoiding a per-sample apply over ~166 loci.
        out[, k] <- do.call(paste0, asplit(bits, 2))
    }

    out
}

# Rebases times onto a common origin, keeping `lead_time` days of clearance
# below the earliest observation.
#
# The clearance is not cosmetic. The model does not use observation time
# directly: it orders infections by an inferred *infection* time, obtained by
# subtracting an infection duration drawn from the --symptomatic-idp /
# --asymptomatic-idp priors. Those priors put mass out to 43 days (symptomatic)
# and 364 days (asymptomatic), so an observation sitting at t = 0 forces the
# infection that caused it to a negative time. Rebasing to zero, which is what
# this function used to do, therefore breaks every node observed inside the
# first year of the series even though `observation_time` itself is
# non-negative and passes the binary's own check.
rebase_observation_time <- function(times, origin = "min", lead_time = 0) {
    if (identical(origin, "min")) {
        return(times - min(times, na.rm = TRUE) + lead_time)
    }
    base <- if (is.numeric(origin)) origin else as.numeric(as.Date(as.character(origin)))
    times - base
}

build_nodes <- function(dat, epi, locus_tbl, allowed_parents, origin = "min",
                        lead_time = 0) {
    sample_ids <- sort(unique(dat$sample_id))
    genotypes <- encode_genotypes(dat, locus_tbl, sample_ids)
    loci <- locus_tbl$locus

    observed_genotype <- lapply(seq_along(sample_ids), function(i) {
        data.frame(locus = loci, genotype = unname(genotypes[i, ]), stringsAsFactors = FALSE)
    })

    nodes <- tibble::tibble(id = sample_ids, observed_genotype = observed_genotype) |>
        dplyr::left_join(allowed_parents, by = c("id" = "sample_id")) |>
        dplyr::left_join(
            dplyr::select(epi, sample_id, symptomatic, observation_time),
            by = c("id" = "sample_id")
        )

    incomplete <- is.na(nodes$observation_time) | is.na(nodes$symptomatic)
    if (any(incomplete)) {
        stop(sprintf(
            "%d node(s) lack observation_time or symptomatic status, e.g. %s. Drop them before building nodes.",
            sum(incomplete), paste(utils::head(nodes$id[incomplete], 5), collapse = ", ")
        ), call. = FALSE)
    }

    nodes$observation_time <- rebase_observation_time(
        nodes$observation_time, origin, lead_time
    )
    nodes$symptomatic <- as.logical(nodes$symptomatic)

    # Candidate parents were chosen before any node was dropped, so prune any
    # that no longer correspond to a node in this stratum.
    nodes$allowed_parents <- lapply(nodes$allowed_parents, function(p) {
        if (is.null(p)) character(0) else intersect(as.character(p), sample_ids)
    })

    nodes[, c("id", "observed_genotype", "allowed_parents", "symptomatic", "observation_time")]
}

# Bit position k of each genotype string maps to a real allele label; without
# this table the model's per-locus output cannot be traced back to a sequence.
allele_index_table <- function(locus_tbl) {
    dplyr::bind_rows(lapply(seq_len(nrow(locus_tbl)), function(k) {
        data.frame(
            locus = locus_tbl$locus[k],
            bit_position = seq_along(locus_tbl$allele_levels[[k]]),
            allele = locus_tbl$allele_levels[[k]],
            allele_freq = locus_tbl$allele_freqs[[k]],
            stringsAsFactors = FALSE
        )
    }))
}
