# Pairwise identity-by-descent estimation via dcifer.
#
# This replaces the original parallel_ibd.R, which materialised an
# nsmp^2-row expand.grid and then iterated it a row at a time with
# foreach/iterators. At this project's scale that is prohibitive: 2,712 samples
# is 3.7M pairs, and the grid alone is 7.4M rows before filtering. Here the
# pairs are enumerated implicitly, each worker writes into preallocated vectors,
# and results are cached so a rerun skips the whole computation.

# Splits the outer loop index across workers. Pair count falls linearly in `i`
# (index i is compared against nsmp - i partners), so a strided assignment
# gives each worker near-identical total work without any bookkeeping.
stride_partition <- function(n, workers) {
    if (n < 1) {
        return(list())
    }
    workers <- max(1L, min(workers, n))
    lapply(seq_len(workers), function(w) seq.int(w, n, by = workers))
}

# Mirrors the argument marshalling dcifer::ibdDat does internally, so that
# ibdPair can be called directly in the inner loop.
ibd_setup <- function(afreq, confint, nr, rnull, pval, reval = NULL) {
    mnewton <- !confint
    if (mnewton) {
        list(
            mnewton = TRUE, tol = 1 / nr, reval = NULL, logr = NULL,
            neval = NULL, inull = NULL, afreq = lapply(afreq, log)
        )
    } else {
        if (!inherits(reval, "matrix")) {
            reval <- dcifer::generateReval(M = 1, rval = reval, nr = nr)
        }
        neval <- ncol(reval)
        list(
            mnewton = FALSE, tol = NULL, reval = reval,
            logr = dcifer::logReval(reval, M = 1, neval = neval),
            neval = neval,
            inull = if (pval) which.min(abs(reval - rnull)) else NULL,
            afreq = lapply(afreq, log)
        )
    }
}

compute_pairwise_ibd <- function(dsmp, coi, afreq,
                                 cores = NULL, nr = 1000, confint = TRUE,
                                 pval = TRUE, rnull = 0, alpha = 0.05,
                                 quiet = FALSE) {
    nsmp <- length(dsmp)
    if (nsmp < 2) {
        stop("At least two samples are required to estimate relatedness.", call. = FALSE)
    }
    snames <- names(dsmp)
    setup <- ibd_setup(afreq, confint, nr, rnull, pval)
    nloc <- length(setup$afreq)
    total_pairs <- nsmp * (nsmp - 1) / 2

    cores <- if (is.null(cores)) max(1L, parallel::detectCores() - 1L) else cores
    # Forking is unavailable on Windows; fall back to a single process there.
    if (.Platform$OS.type != "unix") cores <- 1L
    cores <- max(1L, min(as.integer(cores), nsmp - 1L))

    if (!quiet) {
        message(sprintf(
            "Estimating IBD for %s pairs across %d samples on %d core(s)...",
            format(total_pairs, big.mark = ","), nsmp, cores
        ))
    }

    partitions <- stride_partition(nsmp - 1L, cores)

    run_partition <- function(idx) {
        starts <- partitions[[idx]]
        n_out <- sum(nsmp - starts)
        x <- integer(n_out)
        y <- integer(n_out)
        estimate <- numeric(n_out)
        p_value <- rep(NA_real_, n_out)
        ci_lower <- rep(NA_real_, n_out)
        ci_upper <- rep(NA_real_, n_out)

        pos <- 0L
        report_every <- max(1L, as.integer(n_out / 20))
        for (i in starts) {
            for (j in (i + 1L):nsmp) {
                res <- dcifer::ibdPair(
                    list(dsmp[[i]], dsmp[[j]]), c(coi[i], coi[j]), setup$afreq,
                    M = 1, pval = pval, confreg = confint,
                    rnull = rnull, alpha = alpha, mnewton = setup$mnewton,
                    freqlog = TRUE, reval = setup$reval, tol = setup$tol,
                    logr = setup$logr, neval = setup$neval, inull = setup$inull,
                    nloc = nloc
                )
                pos <- pos + 1L
                x[pos] <- i
                y[pos] <- j
                estimate[pos] <- res$rhat
                if (pval && !is.null(res$pval)) p_value[pos] <- res$pval
                # confreg is NULL unless confidence regions were requested;
                # range(NULL) would silently yield -Inf/Inf.
                if (confint && !is.null(res$confreg)) {
                    bounds <- range(res$confreg)
                    ci_lower[pos] <- bounds[1]
                    ci_upper[pos] <- bounds[2]
                }
                if (!quiet && idx == 1L && pos %% report_every == 0L) {
                    message(sprintf("  ...%d%%", as.integer(100 * pos / n_out)))
                }
            }
        }

        data.frame(
            x = snames[x], y = snames[y], estimate = estimate,
            p_value = p_value, CI_lower = ci_lower, CI_upper = ci_upper,
            stringsAsFactors = FALSE
        )
    }

    chunks <- if (cores == 1L) {
        lapply(seq_along(partitions), run_partition)
    } else {
        parallel::mclapply(seq_along(partitions), run_partition, mc.cores = cores)
    }

    failed <- vapply(chunks, inherits, logical(1), what = "try-error")
    if (any(failed)) {
        stop("IBD estimation failed in ", sum(failed), " worker(s): ",
            conditionMessage(attr(chunks[[which(failed)[1]]], "condition")),
            call. = FALSE
        )
    }

    dplyr::bind_rows(chunks)
}

# Wraps compute_pairwise_ibd with an on-disk cache keyed by the actual genotype
# data, so re-running a stratum after changing only downstream settings (parent
# thresholds, output paths) is instant, while changed input invalidates it.
compute_pairwise_ibd_cached <- function(dsmp, coi, afreq, cache_path, enabled = TRUE, ...) {
    if (!enabled) {
        return(compute_pairwise_ibd(dsmp, coi, afreq, ...))
    }
    key <- rlang::hash(list(dsmp, coi, afreq, list(...)))
    if (file.exists(cache_path)) {
        cached <- readRDS(cache_path)
        if (identical(cached$key, key)) {
            message("Reusing cached IBD estimates from ", cache_path)
            return(cached$value)
        }
        message("Cached IBD estimates are stale; recomputing.")
    }
    value <- compute_pairwise_ibd(dsmp, coi, afreq, ...)
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(key = key, value = value), cache_path)
    value
}

# Expands the upper-triangular estimates into both directions and, for each
# sample, keeps its most related partners as candidate transmission parents.
#
# `min_parents` guarantees each node retains at least that many candidates even
# when none clear `threshold`. Raising the threshold otherwise leaves some nodes
# with an empty candidate set, which is a meaningful modelling statement (the
# infection cannot be explained by anything observed) and not always what is
# intended. 0 leaves the behaviour purely threshold-driven.
select_allowed_parents <- function(ibd_df, sample_ids, threshold, max_parents,
                                   statistic = "estimate", min_parents = 0) {
    if (!statistic %in% names(ibd_df)) {
        stop(sprintf("IBD statistic `%s` not present in results.", statistic), call. = FALSE)
    }

    scored <- tibble::tibble(
        from = c(ibd_df$x, ibd_df$y),
        to = c(ibd_df$y, ibd_df$x),
        score = rep(ibd_df[[statistic]], 2)
    ) |>
        dplyr::filter(from != to, !is.na(score))

    parents <- scored |>
        dplyr::group_by(from) |>
        dplyr::arrange(dplyr::desc(score), .by_group = TRUE) |>
        dplyr::mutate(rank = dplyr::row_number()) |>
        # Keep everything above the threshold, plus enough top-ranked partners
        # to satisfy min_parents, then cap the total.
        dplyr::filter(score > threshold | rank <= min_parents) |>
        dplyr::slice_head(n = max_parents) |>
        dplyr::summarize(allowed_parents = list(to), .groups = "drop") |>
        dplyr::rename(sample_id = from)

    # Every node needs an entry, and it must be an empty character vector
    # rather than NULL so that it serialises as [] and not {}.
    tibble::tibble(sample_id = sample_ids) |>
        dplyr::left_join(parents, by = "sample_id") |>
        dplyr::mutate(allowed_parents = lapply(
            allowed_parents,
            function(p) if (is.null(p)) character(0) else as.character(p)
        ))
}
