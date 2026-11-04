#!/usr/bin/env Rscript
#
# Builds the JSON files the plasmotrack transmission model reads.
#
# Send a collaborator this file together with, from the same directory:
#   run_prep.R
#   config.template.yaml
#   lib/
# legacy/ and networks/ are not needed.
#
#   Rscript make_plasmotrack_json.R --init my_study.yaml
#   # edit paths and column names in my_study.yaml
#   Rscript make_plasmotrack_json.R --config my_study.yaml --dry-run
#   Rscript make_plasmotrack_json.R --config my_study.yaml
#   Rscript make_plasmotrack_json.R --config my_study.yaml --strata SiteA --cores 8
#
# One JSON is written per stratum, named <stratum>.json under output.dir
# in the config. strata: [] produces all.json. Alongside each JSON the
# run also writes a filter log, an allele-index table (bit position to
# allele label), and the pairwise IBD estimates.

REQUIRED_PACKAGES <- c(
    "dplyr", "yaml", "jsonlite", "readr", "tibble", "stringr", "rlang", "dcifer"
)

usage <- function() {
    cat(
        "Usage:\n",
        "  Rscript make_plasmotrack_json.R --init [config.yaml]\n",
        "  Rscript make_plasmotrack_json.R --config config.yaml [--dry-run] [--strata A,B] [--cores N]\n\n",
        "Packages:\n",
        "  install.packages(c(\"dplyr\", \"yaml\", \"jsonlite\", \"readr\", \"tibble\", \"stringr\", \"rlang\", \"remotes\"))\n",
        "  remotes::install_github(\"EPPIcenter/dcifer\")\n",
        sep = ""
    )
}

this_file <- (function() {
    args <- commandArgs(trailingOnly = FALSE)
    hit <- grep("^--file=", args)
    if (length(hit) == 0) {
        return(NULL)
    }
    normalizePath(sub("^--file=", "", args[hit[1]]))
})()

parse_args <- function(args) {
    out <- list(
        config = NULL, strata = NULL, dry_run = FALSE,
        cores = NULL, init = NULL, help = FALSE
    )
    i <- 1
    while (i <= length(args)) {
        switch(args[i],
            "--config" = {
                out$config <- args[i + 1]
                i <- i + 1
            },
            "--strata" = {
                out$strata <- strsplit(args[i + 1], ",", fixed = TRUE)[[1]]
                i <- i + 1
            },
            "--cores" = {
                out$cores <- as.integer(args[i + 1])
                i <- i + 1
            },
            "--init" = {
                # Optional path; bare --init writes plasmotrack.yaml here.
                nxt <- args[i + 1]
                if (!is.na(nxt) && !startsWith(nxt, "--")) {
                    out$init <- nxt
                    i <- i + 1
                } else {
                    out$init <- "plasmotrack.yaml"
                }
            },
            "--dry-run" = out$dry_run <- TRUE,
            "--help" = out$help <- TRUE,
            "-h" = out$help <- TRUE,
            stop("Unknown argument: ", args[i], call. = FALSE)
        )
        i <- i + 1
    }
    out
}

check_packages <- function() {
    missing <- REQUIRED_PACKAGES[!vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) == 0) {
        return(invisible(TRUE))
    }
    stop(
        "Missing R package(s): ", paste(missing, collapse = ", "), ".\n",
        "  install.packages(c(\"dplyr\", \"yaml\", \"jsonlite\", \"readr\", \"tibble\", \"stringr\", \"rlang\", \"remotes\"))\n",
        "  remotes::install_github(\"EPPIcenter/dcifer\")",
        call. = FALSE
    )
}

write_template <- function(dest) {
    if (is.null(this_file)) {
        stop("Cannot find config.template.yaml; run this script with Rscript.", call. = FALSE)
    }
    src <- file.path(dirname(this_file), "config.template.yaml")
    if (!file.exists(src)) {
        stop("Template not found next to this script: ", src, call. = FALSE)
    }
    if (file.exists(dest)) {
        stop("Refusing to overwrite ", dest, ". Pass a different path.", call. = FALSE)
    }
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(src, dest)) {
        stop("Could not write ", dest, call. = FALSE)
    }
    message("Wrote ", normalizePath(dest))
    message("Edit the paths and column names, then rerun with --config ", dest)
}

main <- function() {
    opts <- parse_args(commandArgs(trailingOnly = TRUE))
    if (opts$help || (is.null(opts$config) && is.null(opts$init))) {
        usage()
        if (opts$help) {
            return(invisible(NULL))
        }
        stop("Pass --init or --config.", call. = FALSE)
    }

    if (!is.null(opts$init)) {
        write_template(opts$init)
        if (is.null(opts$config)) {
            return(invisible(NULL))
        }
    }

    check_packages()

    if (is.null(this_file)) {
        stop("Run this script with Rscript so it can find run_prep.R.", call. = FALSE)
    }
    prep <- file.path(dirname(this_file), "run_prep.R")
    if (!file.exists(prep)) {
        stop("run_prep.R is not next to this script (", prep, ").", call. = FALSE)
    }
    # Defines run_pipeline(). The sys.nframe() guard in that file keeps it
    # from starting its own run.
    source(prep, local = FALSE)

    cfg <- load_config(opts$config)
    if (!is.null(opts$cores)) cfg$relatedness$cores <- opts$cores
    run_pipeline(cfg, only_strata = opts$strata, dry_run = opts$dry_run)
}

if (sys.nframe() == 0L) {
    main()
}
