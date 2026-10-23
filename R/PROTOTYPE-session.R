# PROTOTYPE, throwaway. Not the package.
#
# An imagined plasmotrackr session, written to react to, for the ticket "How should a plasmotrackr session look?"
# (https://github.com/EPPIcenter/plasmotrack/issues/102). No function here exists. Output after `#>` is invented, with
# Eswatini-sized numbers (457 Infections, 4 Runs of 500 Draws). Lines marked CHOICE are the decisions the ticket must
# settle; ALT names the option the sketch didn't take.

library(plasmotrackr)
library(dplyr)

# ---------------------------------------------------------------------------------------------------------------------
# 1. Engine settings
# ---------------------------------------------------------------------------------------------------------------------
# CHOICE: the DuckDB memory cap and spill directory are R options, read when pt_read() opens its connection, and
# overridable per call. One connection per Analysis object, owned by it, closed when it is garbage-collected or by
# pt_close().
# ALT: a pt_engine(memory_limit =, temp_directory =, threads =) setter that changes a package-global connection.

options(
  plasmotrackr.memory_limit   = "8GB",
  plasmotrackr.temp_directory = "/scratch/mmurphy/duckdb",  # default: tempdir()
  plasmotrackr.threads        = 4                           # default: DuckDB's own
)

# ---------------------------------------------------------------------------------------------------------------------
# 2. Opening an Analysis
# ---------------------------------------------------------------------------------------------------------------------
# CHOICE: one entry point. A path ending in .json is an Analysis manifest (the one pt_analysis reads); anything else
# is one or more Run directories. Bare Runs are labelled by directory name unless named.

esw <- pt_read("~/Workspace/plasmotrack-runs/eswatini-prod/analysis.json")

esw2 <- pt_read(c(
  a = "~/Workspace/plasmotrack-runs/eswatini-prod/run1",
  b = "~/Workspace/plasmotrack-runs/eswatini-prod/run2"
))

# Per-call engine override.
big <- pt_read("he01/analysis.json", memory_limit = "24GB")

# CHOICE: the object is an S3 "pt_analysis" list: the DuckDB connection, the Run table, the Infection and Locus
# lookup tables, and the manifest's optional paths (truth, epi, IBD scores, null-model Analysis). Nothing is read from
# the Traces at pt_read() time except the Window markers and the Draw count per Run.
# ALT: an R6 object with methods (esw$edges()). Rejected in the sketch as un-tidyverse.

esw
#> <pt_analysis> eswatini-prod
#> Input      eswatini.json: 457 Infections, 26 Loci
#> Runs       4, sampling Phase, 500 Draws each
#>   run   draws  iterations        windows (burn-in + sampling)  status
#>   run1    500  50,100-100,000    10 + 10                       finished
#>   run2    500  50,100-100,000    10 + 10                       finished
#>   run3    500  50,100-100,000    10 + 10                       finished
#>   run4    480  50,100- 98,000    10 + 10                       running, window 20 open
#> Truth      none
#> Epi        eswatini-epi.csv (district, season, travel)
#> Engine     DuckDB, memory_limit 8GB, temp /scratch/mmurphy/duckdb
#> ! Runs have unequal Draws (480-500). Diagnostics truncate to 480; pooled estimates use every Draw.

# CHOICE: the Run table is a tibble, and it carries D_r, which every per-Run mean divides by.
pt_runs(esw)
#> # A tibble: 4 x 8
#>   run   label path                      draws first_iteration last_iteration windows status
#>   <chr> <chr> <chr>                     <int>           <int>          <int>   <int> <chr>
#> 1 run1  run1  ~/.../eswatini-prod/run1    500           50100         100000      20 finished
#> ...

pt_infections(esw)   # the input's Infections: infection, observation_time, symptomatic, latent_parent
pt_loci(esw)         # locus, n_alleles

# Subsetting Runs gives another pt_analysis on the same connection.
esw_123 <- pt_runs_keep(esw, c("run1", "run2", "run3"))
# ALT: a runs = argument on every function (see section 3), instead of or as well as this.

# ---------------------------------------------------------------------------------------------------------------------
# 3. Raw tables: lazy
# ---------------------------------------------------------------------------------------------------------------------
# DECIDED: one accessor for every Store table, pt_table(), taking the Store's table name unchanged (doc/store-layout.md).
# Every named pt_ function is something the package computes; raw tables don't take names in that namespace. Its help
# page lists the tables and links the layout doc; an unknown name errors with the list.
# Returns a dbplyr lazy table over the closed Windows of every Run, unioned, with `run` as the first column, keeping
# the Store's own columns, including `phase` and `window`.
# CHOICE: arguments: runs = NULL (all), phase = "sampling" (or "burn_in", or both).

ps <- pt_table(esw, "parent_sets")
ps
#> # Source:   SQL [?? x 7]
#> # Database: DuckDB 1.4.1 [plasmotrackr: eswatini-prod]
#>   run   iteration phase    window infection set_id     prob
#>   <chr>     <int> <chr>     <int> <chr>      <int>    <dbl>
#> 1 run1      50100 sampling     10 E0001          0 0.912
#> 2 run1      50100 sampling     10 E0001          3 0.0871
#> ...

pt_table(esw, "run_scalars", phase = c("burn_in", "sampling")) |>
  mutate(log_posterior = prior + likelihood) |>
  select(run, iteration, phase, log_posterior) |>
  collect()

pt_table(esw, "parent_set_members")   # always both Phases: a set's members are written where it first appears
pt_table(esw, "infections")           # the Trace; pt_infections() is the input's Infections
pt_table(esw, "acceptance")           # Block: one row per sampler per Window, no iteration column
pt_table(esw, "parent_set")
#> Error in `pt_table()`:
#> ! No Store table "parent_set".
#> i Tables: run_scalars, infections, ordering, allele_frequencies, latent_genotypes, parent_sets,
#>   parent_set_members, chains, acceptance, swaps, events. See doc/store-layout.md.

# ---------------------------------------------------------------------------------------------------------------------
# 4. Per-Draw derived tables: still lazy
# ---------------------------------------------------------------------------------------------------------------------
# CHOICE: the lazy/eager line is "per Draw is lazy, over Draws is eager". A derived table with a row per Draw is
# another lazy query the user can filter before collecting. Anything that averages or ranks over Draws is computed
# and returned as a tibble.
#
# CHOICE: a Draw is identified by (run, iteration) in every table; there is no separate draw index column. Where
# posterior needs a 1..n index, pt_as_draws() makes it.

# Edge marginals per Draw: observed parents only.
pt_edge_draws(esw)
#> # Source:   SQL [?? x 5]
#>   run   iteration child parent  prob
#>   <chr>     <int> <chr> <chr>  <dbl>
#> 1 run1      50100 E0001 E0417  0.0871
#> ...

# CHOICE: the Latent parent never appears as a `parent`. CONTEXT.md says it is not a network edge, so the edge
# tables hold Infection-to-Infection rows only, and its share of an Infection's Parent sets is the column
# `p_latent_parent` (pt_analysis's q, "unobserved-source inclusion") in the per-Infection tables.
# ALT: keep pt_analysis's rows, with parent = "E0001_copy" and a logical `latent` column, and filter by default.
# Glossary: add "Latent parent probability", the probability that an Infection's Parent set includes its Latent parent.

pt_infection_draws(esw)
#>   run   iteration infection p_latent_parent expected_observed_parents entropy effective_size
#>   <chr>     <int> <chr>             <dbl>                     <dbl>   <dbl>          <dbl>
#>   ... then one column per composition class: unobserved_only, observed_1, observed_1_unobserved, ...

pt_out_degree(esw)            # run, iteration, infection, out_degree: expected children among the Infections
                              # (functionals over Infections are in section 6)

# Bitstring work: package functions, because gsub and str_count don't translate to DuckDB SQL.
pt_coi_draws(esw)             # run, iteration, infection, coi
pt_table(esw, "latent_genotypes") |>
  mutate(n_alleles = pt_sql_count_alleles(alleles)) # a dbplyr-translatable helper
# ALT: only the named derived tables, no SQL helpers exported.

# Order pairs are O(n^2): lazy, and the user filters to the pairs they care about first.
pt_order_pair_draws(esw) |> filter(i == "E0001")

# ---------------------------------------------------------------------------------------------------------------------
# 5. Summaries over Draws: eager tibbles, pooled by default
# ---------------------------------------------------------------------------------------------------------------------
# CHOICE: every summary takes by_run = FALSE. FALSE gives pooled values (Draw-weighted over Runs); TRUE gives one row
# per Run per key with `run` first, and the same columns otherwise. Never both in one tibble.
# ALT: always both, with run = NA for the pooled row (pt_analysis's convergence table does this with chain = null).

pt_edge_marginals(esw)
#> # A tibble: 18,212 x 4
#>   child parent marginal lag
#>   <chr> <chr>     <dbl> <dbl>
#> 1 E0102 E0099     0.981    12
#> 2 E0233 E0231     0.974     4
#> ...

pt_edge_marginals(esw, by_run = TRUE)
#> # A tibble: 61,807 x 5
#>   run   child parent marginal switches
#>   <chr> <chr> <chr>     <dbl>    <int>
#> ...

pt_parent_sets_pooled(esw)             # infection, members (list<chr>), n_observed, has_latent_parent, class,
                                       # weight, rank, cumulative_mass, credible
pt_parent_sets_pooled(esw, level = 0.95)
pt_map_graph(esw)                      # infection, members, class, weight: rank 1 of the above
pt_infection_summary(esw)              # the Infections table, P17 of the pt_analysis inventory
pt_genotype_presence(esw)              # infection, locus, allele, presence
pt_genotype_mode(esw)                  # infection, locus, alleles, fraction, agreement
pt_coi(esw)                            # infection, median, distribution (list)
pt_infection_times(esw)                # infection, lower, median, upper

# ---------------------------------------------------------------------------------------------------------------------
# 6. Functionals with intervals
# ---------------------------------------------------------------------------------------------------------------------
# CHOICE: one function per functional, named for what it estimates, no bundle. Two shapes, and the shape says which:
#  - a mean over Infections of a per-Draw value is a lazy per-Draw table (run, iteration, value). Its interval is an
#    explicit second step, pt_interval(), which works on any per-Draw table, including one the user builds;
#  - a point summary of the pooled Parent sets is a tibble with no interval columns.
# Each function's docs give the manuscript symbol (psi_edge, psi_MAP, psi_lo/psi_hi).

# Mean out-degree among the Infections, three estimates. Per stratum, all three count the stratum's Infections as
# sources (onward transmission).
pt_mean_out_degree(esw)                # lazy: run, iteration, mean_out_degree. Was psi_edge
#> # Source:   SQL [?? x 3]
#>   run   iteration mean_out_degree
#>   <chr>     <int>           <dbl>
#> 1 run1      50100           0.608
#> ...

pt_mean_out_degree(esw) |> pt_interval()
#> # A tibble: 1 x 5
#>   estimate lower upper level draws
#>      <dbl> <dbl> <dbl> <dbl> <int>
#> 1    0.612 0.581 0.644  0.95  1980

pt_map_out_degree(esw)                 # tibble: map_out_degree, infections. Was psi_MAP
pt_out_degree_range(esw, level = 0.8)  # tibble: lo, hi, level. Was psi_lo, psi_hi; per stratum counts sources,
                                       # so it always contains pt_map_out_degree() (pt_analysis changes to match)

# Recipient-side functionals: means over Infections of what each Infection received.
pt_mean_p_latent_parent(esw)           # lazy: run, iteration, mean_p_latent_parent
pt_composition(esw)                    # lazy, long: run, iteration, class, share
pt_composition(esw) |> pt_interval(by = class, level = 0.9)
pt_multi_source(esw)                   # lazy: run, iteration, augmented, observed_only
pt_entropy(esw)                        # lazy: run, iteration, mean_entropy, mean_effective_size

# Attribution (the recipient's view of who supplied it) is not an out-degree and is named apart. Open on the map
# under epi strata.
pt_attribution(esw, by = "travel")

# Strata: an argument on each functional, from the manifest's epi CSV. Adds a `stratum` column.
pt_mean_out_degree(esw, by = "district") |> pt_interval(by = stratum)

# ---------------------------------------------------------------------------------------------------------------------
# 7. Convergence
# ---------------------------------------------------------------------------------------------------------------------
# One row per quantity, pooled across Runs. The group names follow pt_analysis's convergence table.
conv <- pt_convergence(esw)
conv
#> # A tibble: 21,604 x 6
#>   group            name            ess  rhat  n_runs draws_per_run
#>   <chr>            <chr>         <dbl> <dbl>   <int>         <int>
#> 1 log_posterior    log_posterior  212.  1.01       4           480
#> 2 edge             E0102<-E0099    88.  1.13       4           480
#> 3 p_latent_parent  E0017          NA    NA         4           480   # constant: every set contains it
#> ...

pt_convergence(esw, groups = c("log_posterior", "functional"), by_run = TRUE)

summary(conv)                          # the pt_analysis diagnostics table: per group, median and min ESS,
                                       # share with R-hat > 1.1, max R-hat
pt_tree_diagnostics(esw)               # edge spread, parent-set TVD, order-pair switches: settled in the ticket
                                       # "Which tree convergence diagnostics ship, and how are they defined?"

# Escape hatch to posterior: a draws_array truncated to equal-length Runs, with a message saying so.
pt_as_draws(esw, groups = "functional") |> posterior::summarise_draws()
#> i Truncated Runs to 480 Draws each (run4 is shortest).

# Plots, ggplot2 in Suggests.
pt_plot_trace(esw, "log_posterior")
pt_mean_out_degree(esw) |> pt_plot_rank()     # plots take any per-Draw table

# ---------------------------------------------------------------------------------------------------------------------
# 8. Networks out
# ---------------------------------------------------------------------------------------------------------------------
# CHOICE: export, don't draw. Vertices are every Infection with pt_infection_summary() columns as attributes; edges
# are pooled Edge marginals at or above `floor`, child and parent as tidygraph's to and from.
g  <- pt_as_igraph(esw, floor = 0.05)
tg <- pt_as_tbl_graph(esw, floor = 0.05)
pt_as_igraph(esw, graph = "map")       # the MAP graph instead of thresholded marginals

pt_close(esw)
