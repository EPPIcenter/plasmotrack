# Quality-control filters applied to long-form allele data.
#
# Every filter records how much it removed so a run leaves behind an auditable
# trail; `attr(dat, "filter_log")` holds the tally.

new_filter_log <- function() {
    tibble::tibble(
        step = character(), samples = integer(), loci = integer(),
        alleles = integer(), rows = integer()
    )
}

log_step <- function(log, dat, step) {
    dplyr::bind_rows(log, tibble::tibble(
        step = step,
        samples = dplyr::n_distinct(dat$sample_id),
        loci = dplyr::n_distinct(dat$locus),
        alleles = nrow(dplyr::distinct(dat, locus, allele)),
        rows = nrow(dat)
    ))
}

# Drops loci carrying more alleles than the model can represent. The model packs
# each locus into an AllelesBitSet<128>, so this is a correctness requirement
# rather than a tuning knob: keeping a wider locus corrupts the genotype.
drop_wide_loci <- function(dat, max_alleles) {
    counts <- dat |>
        dplyr::group_by(locus) |>
        dplyr::summarize(n_alleles = dplyr::n_distinct(allele), .groups = "drop")
    keep <- counts$locus[counts$n_alleles <= max_alleles]
    dat |> dplyr::filter(locus %in% keep)
}

# Expected heterozygosity per locus, from dcifer's COI-corrected population
# allele frequencies. Using the corrected frequencies matters here: naive
# counting inflates diversity in polyclonal infections, and ~70% of these
# infections are polyclonal at one or more loci.
locus_heterozygosity <- function(dat) {
    dsmp <- dcifer::formatDat(dat, svar = "sample_id", lvar = "locus", avar = "allele")
    coi <- dcifer::getCOI(dsmp, lrank = 2)
    afreq <- dcifer::calcAfreq(dsmp, coi, tol = 1e-5)
    vapply(afreq, function(p) 1 - sum(p^2), numeric(1))
}

# Drops loci that carry too little variation to inform relatedness, either by
# an absolute threshold or by keeping only the `top_n` most diverse. Both share
# one heterozygosity computation because it is not cheap.
#
# Selection is deliberately made within the stratum being prepared: a locus that
# is monomorphic in one country contributes nothing there even if it is diverse
# elsewhere, and the strata are modelled independently. The consequence is that
# different strata can end up with different locus panels, which matters if you
# intend to compare per-locus statistics across them.
filter_by_heterozygosity <- function(dat, min_he, top_n) {
    wants_threshold <- !is.null(min_he) && min_he > 0
    wants_top_n <- !is.null(top_n) && top_n < dplyr::n_distinct(dat$locus)
    if (!wants_threshold && !wants_top_n) {
        return(dat)
    }

    he <- locus_heterozygosity(dat)
    keep <- names(he)

    if (wants_threshold) {
        dropped <- sort(he[he < min_he])
        keep <- intersect(keep, names(he)[he >= min_he])
        if (length(dropped) > 0) {
            message(sprintf(
                "  dropped %d locus/loci with HE < %.3g: %s",
                length(dropped), min_he,
                paste(sprintf("%s (%.3f)", names(dropped), dropped), collapse = ", ")
            ))
        }
    }

    if (wants_top_n) {
        ranked <- names(sort(he[keep], decreasing = TRUE))
        keep <- utils::head(ranked, top_n)
    }

    if (length(keep) == 0) {
        stop("The heterozygosity filter removed every locus.", call. = FALSE)
    }
    dat |> dplyr::filter(locus %in% keep)
}

apply_filters <- function(dat, filters, quiet = FALSE) {
    log <- log_step(new_filter_log(), dat, "input")

    if (filters$min_reads_per_allele > 0) {
        dat <- dat |> dplyr::filter(reads >= filters$min_reads_per_allele)
        log <- log_step(log, dat, "min_reads_per_allele")
    }

    if (filters$min_allele_proportion > 0) {
        dat <- dat |>
            dplyr::group_by(sample_id, locus) |>
            dplyr::filter(reads / sum(reads) >= filters$min_allele_proportion) |>
            dplyr::ungroup()
        log <- log_step(log, dat, "min_allele_proportion")
    }

    if (filters$min_reads_per_sample_locus > 0) {
        dat <- dat |>
            dplyr::group_by(sample_id, locus) |>
            dplyr::filter(sum(reads) >= filters$min_reads_per_sample_locus) |>
            dplyr::ungroup()
        log <- log_step(log, dat, "min_reads_per_sample_locus")
    }

    if (filters$min_allele_contributors > 0) {
        dat <- dat |>
            dplyr::group_by(locus, allele) |>
            dplyr::filter(dplyr::n_distinct(sample_id) >= filters$min_allele_contributors) |>
            dplyr::ungroup()
        log <- log_step(log, dat, "min_allele_contributors")
    }

    if (!is.null(filters$drop_loci)) {
        dat <- dat |> dplyr::filter(!locus %in% filters$drop_loci)
        log <- log_step(log, dat, "drop_loci")
    }

    if (filters$min_samples_per_locus > 0) {
        dat <- dat |>
            dplyr::group_by(locus) |>
            dplyr::filter(dplyr::n_distinct(sample_id) >= filters$min_samples_per_locus) |>
            dplyr::ungroup()
        log <- log_step(log, dat, "min_samples_per_locus")
    }

    if (filters$min_loci_per_sample > 0) {
        dat <- dat |>
            dplyr::group_by(sample_id) |>
            dplyr::filter(dplyr::n_distinct(locus) >= filters$min_loci_per_sample) |>
            dplyr::ungroup()
        log <- log_step(log, dat, "min_loci_per_sample")
    }

    dat <- drop_wide_loci(dat, filters$max_alleles_per_locus)
    log <- log_step(log, dat, "max_alleles_per_locus")

    # Heterozygosity is judged last, so it is measured on data that has already
    # passed every coverage filter rather than on noise.
    dat <- filter_by_heterozygosity(
        dat, filters$min_expected_heterozygosity, filters$top_n_loci_by_he
    )
    log <- log_step(log, dat, "heterozygosity")

    if (nrow(dat) == 0) {
        stop("All rows were removed by QC filters; loosen the thresholds.", call. = FALSE)
    }

    if (!quiet) print(log, n = Inf)
    attr(dat, "filter_log") <- log
    dat
}
