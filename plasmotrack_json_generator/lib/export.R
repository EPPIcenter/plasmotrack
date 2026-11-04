# Serialisation of the model input, plus a validator that checks the same
# invariants the tx_model binary enforces at parse time.
#
# The binary is a C++ program using nlohmann::json; when an input violates its
# expectations it fails with terse messages ("Observation time for node X is
# negative.", "Locus X does not exist.", "type must be number, but is array")
# after the JSON has already been loaded. Checking here turns those into
# actionable errors that name the offending node or locus.

MODEL_MAX_ALLELES <- 128L

# Longest infection duration the supplied priors can draw, in days. The binary
# infers an infection time as observation_time minus a duration sampled from
# --symptomatic-idp / --asymptomatic-idp; the shipped asymptomatic prior runs to
# 364 days, the symptomatic one to 43. An observation closer to zero than this
# can be pushed to a negative infection time, which the binary's own
# "observation time is negative" check does not catch because it only inspects
# the observation.
MODEL_MAX_INFECTION_DURATION <- 364L

validate_model_input <- function(nodes, loci,
                                 max_infection_duration = MODEL_MAX_INFECTION_DURATION) {
    problems <- character()
    add <- function(...) problems <<- c(problems, sprintf(...))

    if (anyNA(nodes$id) || any(!nzchar(nodes$id))) add("Some node ids are missing or empty.")
    if (anyDuplicated(nodes$id) > 0) {
        dupes <- unique(nodes$id[duplicated(nodes$id)])
        add("Duplicate node ids: %s", paste(utils::head(dupes, 5), collapse = ", "))
    }

    if (!is.numeric(nodes$observation_time) || anyNA(nodes$observation_time)) {
        add("observation_time must be numeric and complete.")
    } else {
        negative <- nodes$id[nodes$observation_time < 0]
        if (length(negative) > 0) {
            add(
                "Negative observation_time for node(s): %s",
                paste(utils::head(negative, 5), collapse = ", ")
            )
        }
        if (any(!is.finite(nodes$observation_time))) add("observation_time contains non-finite values.")

        # Warned rather than refused: how much clearance is really needed
        # depends on which duration prior applies to each node, and a cohort
        # of purely symptomatic infections only needs 43 days.
        short <- sum(nodes$observation_time < max_infection_duration)
        if (short > 0) {
            warning(sprintf(
                paste0(
                    "%d of %d node(s) have observation_time < %d, the longest ",
                    "duration the infection-duration priors can draw. Those ",
                    "infections can be assigned a negative infection time. ",
                    "Raise epi.observation_time.lead_time if the model rejects ",
                    "this input or produces non-finite likelihoods."
                ),
                short, nrow(nodes), max_infection_duration
            ), call. = FALSE)
        }
    }

    if (!is.logical(nodes$symptomatic) || anyNA(nodes$symptomatic)) {
        add("symptomatic must be logical and complete.")
    }

    if (any(loci$num_alleles > MODEL_MAX_ALLELES)) {
        wide <- loci$locus[loci$num_alleles > MODEL_MAX_ALLELES]
        add(
            "Loci exceed the model's %d-allele bitset: %s",
            MODEL_MAX_ALLELES, paste(utils::head(wide, 5), collapse = ", ")
        )
    }

    freq_len <- lengths(loci$allele_freqs)
    if (any(freq_len != loci$num_alleles)) {
        bad <- loci$locus[freq_len != loci$num_alleles]
        add(
            "allele_freqs length does not match num_alleles at: %s",
            paste(utils::head(bad, 5), collapse = ", ")
        )
    }
    freq_sums <- vapply(loci$allele_freqs, sum, numeric(1))
    if (any(abs(freq_sums - 1) > 1e-6)) {
        bad <- loci$locus[abs(freq_sums - 1) > 1e-6]
        add(
            "allele_freqs do not sum to 1 at: %s",
            paste(utils::head(bad, 5), collapse = ", ")
        )
    }
    if (any(vapply(loci$allele_freqs, function(p) any(p <= 0), logical(1)))) {
        add("Some allele frequencies are zero or negative; the model takes their log.")
    }

    expected_width <- stats::setNames(loci$num_alleles, loci$locus)
    known_loci <- loci$locus

    for (i in seq_len(nrow(nodes))) {
        g <- nodes$observed_genotype[[i]]
        unknown <- setdiff(g$locus, known_loci)
        if (length(unknown) > 0) {
            add(
                "Node %s references unknown locus/loci: %s",
                nodes$id[i], paste(utils::head(unknown, 3), collapse = ", ")
            )
            next
        }
        if (!setequal(g$locus, known_loci)) {
            add("Node %s does not carry a genotype for every locus.", nodes$id[i])
            next
        }
        widths <- nchar(g$genotype)
        mismatch <- widths != expected_width[g$locus]
        if (any(mismatch)) {
            add(
                "Node %s genotype width mismatches num_alleles at locus %s (%d vs %d).",
                nodes$id[i], g$locus[which(mismatch)[1]],
                widths[which(mismatch)[1]], expected_width[[g$locus[which(mismatch)[1]]]]
            )
        }
        if (any(grepl("[^01]", g$genotype))) {
            add("Node %s has a genotype containing characters other than 0/1.", nodes$id[i])
        }
        # Only the first few failures are worth reporting; the rest will be the
        # same defect repeated across thousands of nodes.
        if (length(problems) > 20) break
    }

    all_ids <- nodes$id
    orphan_parents <- setdiff(unique(unlist(nodes$allowed_parents)), all_ids)
    if (length(orphan_parents) > 0) {
        add(
            "allowed_parents reference %d id(s) that are not nodes, e.g. %s",
            length(orphan_parents), paste(utils::head(orphan_parents, 5), collapse = ", ")
        )
    }
    self_parent <- nodes$id[mapply(function(id, p) id %in% p, nodes$id, nodes$allowed_parents)]
    if (length(self_parent) > 0) {
        add(
            "Node(s) listed as their own parent: %s",
            paste(utils::head(self_parent, 5), collapse = ", ")
        )
    }

    # Not a hard error: a node with no candidate parent asserts that nothing
    # observed could have infected it, which is legitimate for an importation.
    # It is worth surfacing though, because it is usually an unintended
    # consequence of raising the relatedness threshold, and it changes what the
    # model can do with that node.
    isolated <- sum(lengths(nodes$allowed_parents) == 0)
    if (isolated > 0) {
        warning(sprintf(
            paste0(
                "%d of %d node(s) have an empty allowed_parents set. ",
                "They can only be explained as importations. ",
                "Set relatedness.min_parents if that is not intended."
            ),
            isolated, nrow(nodes)
        ), call. = FALSE)
    }

    if (length(problems) > 0) {
        stop("Model input validation failed:\n  - ", paste(problems, collapse = "\n  - "),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

write_model_input <- function(nodes, loci, path, digits = NA) {
    validate_model_input(nodes, loci)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

    payload <- list(
        nodes = nodes,
        loci = loci[, c("locus", "allele_freqs", "num_alleles")]
    )

    # digits = NA keeps full double precision. jsonlite's default of 4 decimal
    # places materially degrades rare allele frequencies, which the model feeds
    # through log(); frequencies below 5e-5 would round to zero outright.
    jsonlite::write_json(payload, path, digits = digits, na = "null", null = "list")
    invisible(path)
}
