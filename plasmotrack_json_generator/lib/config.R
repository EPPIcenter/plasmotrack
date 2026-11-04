# Configuration loading and validation for the transmission-network prep pipeline.
#
# A config is a YAML file describing where a dataset lives, how its columns map
# onto the pipeline's internal vocabulary (sample_id / locus / allele / reads),
# and which QC thresholds to apply. Everything dataset-specific belongs here so
# that the code under lib/ stays dataset-agnostic.

DEFAULT_CONFIG <- list(
    dataset = NULL,
    genetic = list(
        path = NULL,
        delim = NULL, # NULL => infer from file extension
        columns = list(
            sample_id = NULL,
            locus = NULL,
            allele = NULL,
            reads = NULL
        ),
        locus_pattern = NULL, # optional regex; keep only matching loci
        allele_from = NULL # optional: build allele as locus + sep + this column
    ),
    epi = list(
        path = NULL,
        delim = NULL,
        columns = list(sample_id = NULL),
        observation_time = list(
            columns = NULL, # coalesced left to right
            type = "date", # "date" or "numeric"
            # Tried in order; the first that parses a given value wins. Kept
            # deliberately narrow, because permissive formats such as %d/%m/%Y
            # and %m/%d/%Y silently accept each other's input.
            date_formats = c("%Y-%m-%d", "%A, %d %B %Y", "%d %B %Y", "%B %d, %Y"),
            # Anything outside this window is treated as a parse failure rather
            # than trusted; a mis-matched format typically yields a wild year.
            min_year = 1990,
            max_year = NULL, # NULL => current year + 1
            origin = "min", # "min" (per stratum) or a fixed date/number
            # Days of clearance left below the earliest observation when
            # origin = "min". The model subtracts a sampled infection duration
            # from each observation time to get an infection time, so the
            # series must start far enough above zero that this stays valid.
            # 365 covers the full support of both supplied duration priors
            # (43 days symptomatic, 364 asymptomatic); it is a constant offset,
            # so it costs nothing and changes no time difference or ordering.
            lead_time = 365
        ),
        symptomatic = list(
            constant = NULL, # TRUE/FALSE applied to every node
            column = NULL, # or read from a column
            true_values = NULL, # values in `column` meaning symptomatic
            # If given, anything outside true_values/false_values is treated as
            # unknown rather than assumed asymptomatic.
            false_values = NULL
        ),
        strata = NULL, # column(s) to split model runs on; NULL => single run
        min_stratum_size = 2, # skip strata too small to infer a network from
        drop_strata = NULL, # explicit stratum blocklist, e.g. ["unknown"]
        enrich = list() # list of joins used to fill in missing epi fields
    ),
    filters = list(
        min_reads_per_allele = 0,
        min_allele_proportion = 0, # within sample x locus
        min_allele_contributors = 0, # distinct samples carrying an allele
        min_reads_per_sample_locus = 0,
        min_loci_per_sample = 0,
        min_samples_per_locus = 0,
        max_alleles_per_locus = 128, # hard limit of the model's AllelesBitSet<128>
        # Locus informativeness. Both are evaluated within each stratum, so a
        # locus monomorphic in one population is dropped there but kept where it
        # is diverse. 0 / NULL disables.
        min_expected_heterozygosity = 0,
        top_n_loci_by_he = NULL, # keep the N most heterozygous loci
        drop_loci = NULL, # explicit locus blocklist
        by_stratum = NULL
    ),
    relatedness = list(
        ibd_threshold = 0.1,
        max_parents = 25,
        # Floor on candidate parents per node, applied regardless of threshold.
        # 0 means a node with nothing above the threshold gets an empty set.
        min_parents = 0,
        statistic = "estimate", # "estimate" or "CI_lower"
        cores = NULL, # NULL => parallel::detectCores() - 1
        nr = 1000,
        confint = TRUE,
        cache = TRUE,
        # Named list of per-stratum overlays, e.g. Eswatini: {max_parents: 10}.
        # Keys absent from an overlay keep the top-level relatedness value.
        by_stratum = NULL
    ),
    output = list(
        dir = NULL,
        json_digits = NA # NA => full double precision
    )
)

# Recursively overlay `override` onto `base`, keeping keys absent from
# `override` at their default. Leaf lists supplied by the user replace the
# default wholesale so that e.g. `strata: [country, district]` behaves.
# A key present but null in YAML (`foo: ~`) leaves the default in place rather
# than deleting it, since assigning NULL into a list removes the element.
merge_config <- function(base, override) {
    if (is.null(override)) {
        return(base)
    }
    for (key in names(override)) {
        value <- override[[key]]
        if (is.null(value)) {
            next
        }
        if (is.list(value) && !is.null(base[[key]]) && is.list(base[[key]]) &&
            !is.null(names(value))) {
            base[[key]] <- merge_config(base[[key]], value)
        } else {
            base[[key]] <- value
        }
    }
    base
}

require_field <- function(value, path) {
    empty <- is.null(value) ||
        length(value) == 0 ||
        (is.character(value) && all(!nzchar(value)))
    if (empty) {
        stop(sprintf("Config is missing required field `%s`.", path), call. = FALSE)
    }
    invisible(value)
}

validate_config <- function(cfg) {
    require_field(cfg$dataset, "dataset")
    require_field(cfg$genetic$path, "genetic.path")
    require_field(cfg$output$dir, "output.dir")

    for (col in c("sample_id", "locus", "reads")) {
        require_field(cfg$genetic$columns[[col]], paste0("genetic.columns.", col))
    }
    if (is.null(cfg$genetic$columns$allele) && is.null(cfg$genetic$allele_from)) {
        stop("Config must set either `genetic.columns.allele` or `genetic.allele_from`.",
            call. = FALSE
        )
    }

    if (!is.null(cfg$epi$path)) {
        require_field(cfg$epi$columns$sample_id, "epi.columns.sample_id")
        require_field(cfg$epi$observation_time$columns, "epi.observation_time.columns")
    }

    sym <- cfg$epi$symptomatic
    if (is.null(sym$constant) && is.null(sym$column)) {
        stop("Config must set either `epi.symptomatic.constant` or `epi.symptomatic.column`.",
            call. = FALSE
        )
    }
    if (!is.null(sym$column) && is.null(sym$true_values)) {
        stop("`epi.symptomatic.column` requires `epi.symptomatic.true_values`.", call. = FALSE)
    }

    # The model packs each locus into an AllelesBitSet<128>; anything wider is
    # silently truncated or rejected by the binary, so refuse it here instead.
    if (cfg$filters$max_alleles_per_locus > 128) {
        stop("`filters.max_alleles_per_locus` cannot exceed 128 (model bitset width).",
            call. = FALSE
        )
    }

    if (!cfg$relatedness$statistic %in% c("estimate", "CI_lower")) {
        stop("`relatedness.statistic` must be \"estimate\" or \"CI_lower\".", call. = FALSE)
    }

    invisible(cfg)
}

# Top-level section values, then any `by_stratum.<name>` overlay. Keys
# absent from the overlay keep the shared value.
overlay_for_stratum <- function(section, stratum) {
    overlay <- section$by_stratum[[stratum]]
    if (is.null(overlay)) {
        return(section)
    }
    merge_config(section, overlay)
}

relatedness_for_stratum <- function(cfg, stratum) {
    overlay_for_stratum(cfg$relatedness, stratum)
}

filters_for_stratum <- function(cfg, stratum) {
    overlay_for_stratum(cfg$filters, stratum)
}

load_config <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("Config file not found: %s", path), call. = FALSE)
    }
    cfg <- merge_config(DEFAULT_CONFIG, yaml::read_yaml(path))
    cfg$.config_path <- path
    validate_config(cfg)
    cfg
}
