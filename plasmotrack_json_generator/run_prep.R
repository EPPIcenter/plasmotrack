#!/usr/bin/env Rscript
#
# Prepares amplicon + epidemiological data for the tx_model transmission
# network binary.
#
#   Rscript R/preprocessing/run_prep.R --config config/gene8.yaml
#   Rscript R/preprocessing/run_prep.R --config config/gene8.yaml --strata Namibia
#   Rscript R/preprocessing/run_prep.R --config config/gene8.yaml --dry-run
#
# Everything dataset-specific lives in the YAML config; this script and the
# helpers in lib/ make no assumptions about a particular study.

suppressPackageStartupMessages({
    library(dplyr)
})

# Path of this file, whether it was Rscript'd or source()'d. commandArgs()
# --file= is the outer script, so make_plasmotrack_json.R sourcing this file
# would otherwise look for lib/ next to itself.
this_file <- (function() {
    for (i in rev(seq_len(sys.nframe()))) {
        frame <- sys.frame(i)
        if (exists("ofile", envir = frame, inherits = FALSE)) {
            return(normalizePath(get("ofile", envir = frame)))
        }
    }
    args <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", args)
    if (length(hit) > 0) normalizePath(sub("^--file=", "", args[hit[1]])) else NULL
})()

lib_dir <- if (!is.null(this_file)) file.path(dirname(this_file), "lib") else "R/preprocessing/lib"
if (!dir.exists(lib_dir)) {
    stop("Cannot locate the pipeline library directory (expected at ", lib_dir, ").", call. = FALSE)
}
for (f in c("config.R", "io.R", "filters.R", "relatedness.R", "nodes.R", "export.R")) {
    source(file.path(lib_dir, f))
}

parse_args <- function(args) {
    out <- list(config = NULL, strata = NULL, dry_run = FALSE, cores = NULL)
    i <- 1
    while (i <= length(args)) {
        switch(args[i],
            "--config" = {
                out$config <- args[i + 1]
                i <- i + 1
            },
            "--strata" = {
                out$strata <- strsplit(args[i + 1], ",")[[1]]
                i <- i + 1
            },
            "--cores" = {
                out$cores <- as.integer(args[i + 1])
                i <- i + 1
            },
            "--dry-run" = out$dry_run <- TRUE,
            stop("Unknown argument: ", args[i], call. = FALSE)
        )
        i <- i + 1
    }
    if (is.null(out$config)) stop("--config is required.", call. = FALSE)
    out
}

# Turns the configured strata columns into a single label per sample. With no
# strata configured the whole cohort is one run.
assign_strata <- function(epi, strata_cols) {
    if (is.null(strata_cols) || length(strata_cols) == 0) {
        return(rep("all", nrow(epi)))
    }
    parts <- lapply(strata_cols, function(col) {
        value <- as.character(epi[[col]])
        value[is.na(value) | !nzchar(value)] <- "unknown"
        value
    })
    do.call(paste, c(parts, sep = "_"))
}

prepare_stratum <- function(stratum, dat, epi, cfg, output_dir) {
    message("\n=== ", stratum, " ===")
    message(sprintf(
        "%d samples, %d loci before QC",
        dplyr::n_distinct(dat$sample_id), dplyr::n_distinct(dat$locus)
    ))

    # Filters run inside the stratum so that allele counts, locus coverage, and
    # the 128-allele ceiling are all judged against the population actually
    # being modelled rather than the pooled multi-country dataset.
    dat <- apply_filters(dat, filters_for_stratum(cfg, stratum))
    filter_log <- attr(dat, "filter_log")

    epi <- epi |> dplyr::filter(sample_id %in% unique(dat$sample_id))
    dat <- dat |> dplyr::filter(sample_id %in% epi$sample_id)

    if (dplyr::n_distinct(dat$sample_id) < 2) {
        message("Skipping ", stratum, ": fewer than two usable samples.")
        return(NULL)
    }

    locus_tbl <- build_locus_table(dat)

    dsmp <- dcifer::formatDat(dat, svar = "sample_id", lvar = "locus", avar = "allele")
    coi <- dcifer::getCOI(dsmp, lrank = 2)
    afreq <- dcifer::calcAfreq(dsmp, coi, tol = 1e-5)
    locus_tbl <- attach_allele_freqs(locus_tbl, afreq)

    rel <- relatedness_for_stratum(cfg, stratum)
    ibd <- compute_pairwise_ibd_cached(
        dsmp, coi, afreq,
        cache_path = file.path(output_dir, "cache", paste0(stratum, "_ibd.rds")),
        enabled = rel$cache,
        cores = rel$cores,
        nr = rel$nr,
        confint = rel$confint
    )

    allowed_parents <- select_allowed_parents(
        ibd,
        sample_ids = sort(unique(dat$sample_id)),
        threshold = rel$ibd_threshold,
        max_parents = rel$max_parents,
        statistic = rel$statistic,
        min_parents = rel$min_parents
    )

    nodes <- build_nodes(
        dat, epi, locus_tbl, allowed_parents,
        origin = cfg$epi$observation_time$origin,
        lead_time = cfg$epi$observation_time$lead_time
    )

    json_path <- file.path(output_dir, paste0(stratum, ".json"))
    write_model_input(nodes, locus_tbl, json_path, digits = cfg$output$json_digits)

    readr::write_csv(filter_log, file.path(output_dir, paste0(stratum, "_filter_log.csv")))
    readr::write_csv(
        allele_index_table(locus_tbl),
        file.path(output_dir, paste0(stratum, "_allele_index.csv"))
    )
    readr::write_csv(ibd, file.path(output_dir, paste0(stratum, "_ibd.csv.gz")))

    parent_counts <- lengths(nodes$allowed_parents)
    message(sprintf(
        "Wrote %s: %d nodes, %d loci, median %d candidate parents (%d node(s) with none).",
        json_path, nrow(nodes), nrow(locus_tbl),
        as.integer(stats::median(parent_counts)), sum(parent_counts == 0)
    ))

    tibble::tibble(
        stratum = stratum, nodes = nrow(nodes), loci = nrow(locus_tbl),
        median_parents = stats::median(parent_counts),
        isolated_nodes = sum(parent_counts == 0), json = json_path
    )
}

run_pipeline <- function(cfg, only_strata = NULL, dry_run = FALSE) {
    output_dir <- cfg$output$dir
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    message("Loading genetic data from ", cfg$genetic$path)
    dat <- load_genetic_data(cfg)
    message(sprintf(
        "  %d rows, %d samples, %d loci",
        nrow(dat), dplyr::n_distinct(dat$sample_id), dplyr::n_distinct(dat$locus)
    ))

    message("Loading epi data from ", cfg$epi$path)
    epi <- load_epi_data(cfg)
    epi$stratum <- assign_strata(epi, cfg$epi$strata)

    sequenced <- unique(dat$sample_id)
    matched <- intersect(sequenced, epi$sample_id)
    message(sprintf(
        "  %d of %d sequenced samples matched to epi records (%.1f%%)",
        length(matched), length(sequenced), 100 * length(matched) / length(sequenced)
    ))

    usable <- epi |>
        dplyr::filter(sample_id %in% matched, !is.na(observation_time), !is.na(symptomatic))
    dropped <- length(matched) - nrow(usable)
    if (dropped > 0) {
        message(sprintf("  %d matched sample(s) dropped for missing time or status.", dropped))
    }

    dat <- dat |> dplyr::filter(sample_id %in% usable$sample_id)

    sizes <- table(usable$stratum)
    too_small <- names(sizes)[sizes < cfg$epi$min_stratum_size]
    excluded <- union(too_small, cfg$epi$drop_strata %||% character())
    if (length(excluded) > 0) {
        message(sprintf(
            "Excluding stratum/strata (%s): %s",
            "too small or blocklisted", paste(excluded, collapse = ", ")
        ))
        usable <- usable |> dplyr::filter(!stratum %in% excluded)
        dat <- dat |> dplyr::filter(sample_id %in% usable$sample_id)
    }

    strata <- sort(unique(usable$stratum))
    if (!is.null(only_strata)) {
        unknown <- setdiff(only_strata, strata)
        if (length(unknown) > 0) {
            stop("Unknown stratum/strata requested: ", paste(unknown, collapse = ", "),
                ". Available: ", paste(strata, collapse = ", "),
                call. = FALSE
            )
        }
        strata <- only_strata
    }

    counts <- usable |>
        dplyr::count(stratum, name = "samples") |>
        dplyr::arrange(dplyr::desc(samples))
    message("\nStrata:")
    print(as.data.frame(counts), row.names = FALSE)

    if (dry_run) {
        message("\nDry run: stopping before relatedness estimation.")
        return(invisible(NULL))
    }

    summaries <- lapply(strata, function(s) {
        ids <- usable$sample_id[usable$stratum == s]
        prepare_stratum(
            s,
            dat |> dplyr::filter(sample_id %in% ids),
            usable |> dplyr::filter(stratum == s),
            cfg, output_dir
        )
    })

    summary_tbl <- dplyr::bind_rows(summaries)
    if (nrow(summary_tbl) > 0) {
        summary_path <- file.path(output_dir, "run_summary.csv")
        # A --strata rerun should update those rows, not erase the rest.
        if (file.exists(summary_path) && !is.null(only_strata)) {
            existing <- readr::read_csv(summary_path, show_col_types = FALSE)
            existing <- existing |> dplyr::filter(!stratum %in% summary_tbl$stratum)
            summary_tbl <- dplyr::bind_rows(existing, summary_tbl) |>
                dplyr::arrange(stratum)
        }
        readr::write_csv(summary_tbl, summary_path)
        message("\nSummary:")
        print(as.data.frame(summary_tbl), row.names = FALSE)
    }
    invisible(summary_tbl)
}

main <- function() {
    opts <- parse_args(commandArgs(trailingOnly = TRUE))
    cfg <- load_config(opts$config)
    if (!is.null(opts$cores)) cfg$relatedness$cores <- opts$cores
    run_pipeline(cfg, only_strata = opts$strata, dry_run = opts$dry_run)
}

# Top-level frame means this file is the Rscript target. source() — from
# make_plasmotrack_json.R, or an interactive session — loads the helpers
# without starting a run.
if (sys.nframe() == 0L) {
    main()
}
