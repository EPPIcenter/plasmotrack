# Data loading: turns arbitrarily-shaped source files into the pipeline's
# internal vocabulary. Every column name known to the rest of the pipeline is
# produced here, so downstream code never sees a dataset-specific name.

infer_delim <- function(path) {
    base <- sub("\\.gz$", "", tolower(path))
    if (grepl("\\.(tsv|tab|txt)$", base)) "\t" else ","
}

# Reads delimited text (transparently handling .gz), forcing the columns named
# in `as_character` to be read as text. IDs like "650001" would otherwise be
# guessed as numeric and stop joining against character IDs like "ZMC-238-35".
read_table_auto <- function(path, delim = NULL, as_character = character()) {
    if (!file.exists(path)) {
        stop(sprintf("Input file not found: %s", path), call. = FALSE)
    }
    delim <- if (is.null(delim)) infer_delim(path) else delim
    as_character <- unique(as_character[!is.na(as_character) & nzchar(as_character)])

    col_types <- readr::cols()
    if (length(as_character) > 0) {
        spec <- stats::setNames(
            replicate(length(as_character), readr::col_character(), simplify = FALSE),
            as_character
        )
        col_types <- do.call(readr::cols, spec)
    }

    readr::read_delim(
        path,
        delim = delim,
        col_types = col_types,
        progress = FALSE,
        show_col_types = FALSE
    )
}

# Normalises IDs so that joins survive round-tripping through numeric guessing,
# stray whitespace, and case differences between sources.
normalize_id <- function(x, ignore_case = FALSE) {
    x <- as.character(x)
    x <- trimws(x)
    x[!nzchar(x)] <- NA_character_
    if (ignore_case) x <- toupper(x)
    x
}

load_genetic_data <- function(cfg) {
    cols <- cfg$genetic$columns
    needed <- c(cols$sample_id, cols$locus, cols$reads, cols$allele, cfg$genetic$allele_from)
    needed <- needed[!vapply(needed, is.null, logical(1))]

    raw <- read_table_auto(
        cfg$genetic$path,
        delim = cfg$genetic$delim,
        as_character = c(cols$sample_id, cols$locus, cols$allele, cfg$genetic$allele_from)
    )

    missing <- setdiff(unlist(needed), names(raw))
    if (length(missing) > 0) {
        stop(sprintf(
            "Genetic file %s is missing column(s): %s",
            cfg$genetic$path, paste(missing, collapse = ", ")
        ), call. = FALSE)
    }

    dat <- tibble::tibble(
        sample_id = normalize_id(raw[[cols$sample_id]]),
        locus = as.character(raw[[cols$locus]]),
        reads = as.numeric(raw[[cols$reads]])
    )

    # An allele must be unique across the whole panel, not just within its
    # locus, so a bare variant label gets namespaced by its locus.
    dat$allele <- if (!is.null(cols$allele)) {
        as.character(raw[[cols$allele]])
    } else {
        paste0(dat$locus, ":", as.character(raw[[cfg$genetic$allele_from]]))
    }

    dat <- dat |>
        dplyr::filter(!is.na(sample_id), !is.na(locus), !is.na(allele))

    if (!is.null(cfg$genetic$locus_pattern)) {
        dat <- dat |> dplyr::filter(stringr::str_detect(locus, cfg$genetic$locus_pattern))
    }

    # Collapse technical replicates of the same allele call.
    dat |>
        dplyr::group_by(sample_id, locus, allele) |>
        dplyr::summarize(reads = sum(reads, na.rm = TRUE), .groups = "drop")
}

# Applies the configured `enrich` joins, each of which fills gaps in an epi
# column from a lookup table. Used, for example, to recover `country` for rows
# where it is blank but the facility is known.
apply_enrichment <- function(epi, enrich_specs) {
    for (spec in enrich_specs) {
        lookup <- read_table_auto(spec$path, delim = spec$delim)
        ignore_case <- spec$ignore_case %||% TRUE
        left_key <- normalize_id(epi[[spec$left_key]], ignore_case)

        for (target in names(spec$fill)) {
            source_col <- spec$fill[[target]]
            filled <- rep(NA_character_, nrow(epi))

            # Try each candidate key column in turn so a single field that
            # sometimes holds an ID and sometimes a name still resolves.
            for (right_key in spec$right_keys) {
                if (!right_key %in% names(lookup)) next
                map <- stats::setNames(
                    as.character(lookup[[source_col]]),
                    normalize_id(lookup[[right_key]], ignore_case)
                )
                map <- map[!is.na(names(map)) & !duplicated(names(map))]
                candidate <- unname(map[left_key])
                filled <- dplyr::coalesce(filled, candidate)
            }

            epi[[target]] <- if (target %in% names(epi)) {
                dplyr::coalesce(as.character(epi[[target]]), filled)
            } else {
                filled
            }
        }
    }
    epi
}

# Parses dates that may be written in several different layouts within the same
# column. The GenE8 "harmonized" export, for instance, mixes ISO-8601 with
# long-form text such as "Saturday, 09 December 2023" depending on the country.
parse_dates_multi <- function(x, formats, min_year, max_year) {
    x <- as.character(x)
    # Month and weekday names are matched in English regardless of the caller's
    # locale.
    old_locale <- Sys.getlocale("LC_TIME")
    on.exit(suppressWarnings(Sys.setlocale("LC_TIME", old_locale)), add = TRUE)
    suppressWarnings(Sys.setlocale("LC_TIME", "C"))

    parsed <- rep(NA_real_, length(x))
    for (fmt in formats) {
        remaining <- is.na(parsed)
        if (!any(remaining)) break
        candidate <- suppressWarnings(as.Date(x[remaining], format = fmt))
        parsed[remaining] <- as.numeric(candidate)
    }

    # A format that half-matches (e.g. "%Y/%m/%d" against "12/03/2023") yields a
    # syntactically valid but absurd date. Reject rather than propagate it.
    years <- as.numeric(format(as.Date(parsed, origin = "1970-01-01"), "%Y"))
    implausible <- !is.na(parsed) & (years < min_year | years > max_year)
    if (any(implausible)) {
        warning(sprintf(
            "Discarded %d date(s) outside %d-%d, e.g. %s",
            sum(implausible), min_year, max_year,
            paste(utils::head(unique(x[implausible]), 3), collapse = ", ")
        ), call. = FALSE)
        parsed[implausible] <- NA_real_
    }
    parsed
}

parse_observation_time <- function(epi, spec) {
    present <- intersect(spec$columns, names(epi))
    if (length(present) == 0) {
        stop(sprintf(
            "None of the observation_time columns (%s) are present in the epi file.",
            paste(spec$columns, collapse = ", ")
        ), call. = FALSE)
    }

    max_year <- if (is.null(spec$max_year)) {
        as.integer(format(Sys.Date(), "%Y")) + 1L
    } else {
        spec$max_year
    }

    parse_one <- function(x) {
        if (identical(spec$type, "numeric")) {
            return(suppressWarnings(as.numeric(x)))
        }
        parse_dates_multi(x, spec$date_formats, spec$min_year, max_year)
    }

    # Coalesce the candidate columns left to right, so the first is preferred
    # and later ones only cover its gaps.
    values <- rep(NA_real_, nrow(epi))
    for (col in present) {
        values <- dplyr::coalesce(values, parse_one(epi[[col]]))
    }
    values
}

resolve_symptomatic <- function(epi, spec) {
    if (!is.null(spec$constant)) {
        return(rep(as.logical(spec$constant), nrow(epi)))
    }
    if (!spec$column %in% names(epi)) {
        stop(sprintf("Symptomatic column `%s` not found in the epi file.", spec$column),
            call. = FALSE
        )
    }
    raw <- trimws(as.character(epi[[spec$column]]))
    out <- raw %in% as.character(spec$true_values)
    # A plain %in% test would fold missing values into FALSE, silently
    # relabelling unknown status as asymptomatic. Keep them NA so the runner
    # drops those samples instead.
    if (!is.null(spec$false_values)) {
        known <- raw %in% as.character(c(spec$true_values, spec$false_values))
        out[!known] <- NA
    } else {
        out[is.na(raw) | !nzchar(raw)] <- NA
    }
    out
}

# Produces one row per sample with the fields the model needs, plus any strata
# columns. `observation_time` is still on the raw scale here; it is rebased to a
# non-negative origin per stratum, once the cohort is known.
load_epi_data <- function(cfg) {
    ecfg <- cfg$epi
    raw <- read_table_auto(
        ecfg$path,
        delim = ecfg$delim,
        as_character = c(ecfg$columns$sample_id, ecfg$strata)
    )

    if (!ecfg$columns$sample_id %in% names(raw)) {
        stop(sprintf(
            "Epi file %s is missing the sample id column `%s`.",
            ecfg$path, ecfg$columns$sample_id
        ), call. = FALSE)
    }

    raw <- apply_enrichment(raw, ecfg$enrich)

    epi <- tibble::tibble(
        sample_id = normalize_id(raw[[ecfg$columns$sample_id]]),
        observation_time = parse_observation_time(raw, ecfg$observation_time),
        symptomatic = resolve_symptomatic(raw, ecfg$symptomatic)
    )

    for (col in ecfg$strata) {
        if (!col %in% names(raw)) {
            stop(sprintf("Strata column `%s` not found in the epi file.", col), call. = FALSE)
        }
        epi[[col]] <- as.character(raw[[col]])
    }

    # One record per sample: keep the earliest observation, matching the
    # "first episode wins" rule used for repeat cohort visits.
    epi |>
        dplyr::filter(!is.na(sample_id)) |>
        dplyr::arrange(sample_id, observation_time) |>
        dplyr::distinct(sample_id, .keep_all = TRUE)
}
