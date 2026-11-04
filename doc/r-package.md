# plasmotrackr: spec

Status: accepted spec, 2026-11-03. Nothing here is built. Every decision below was settled on the map
[Find the way to an R analysis package spec](https://github.com/EPPIcenter/plasmotrack/issues/98); each section links
the ticket that holds the reasoning.

`plasmotrackr` is an R package that reads plasmotrack Run stores (`doc/store-layout.md`), checks convergence across
Runs, and hands users tidy data to analyze however they like. It is native R: it reads the Parquet Stores itself
through DuckDB and computes every quantity from the raw tables. The Infections' observation times and candidate
parents come from the input JSON. New Stores embed it (`input.json` or `input.json.gz`) and the IDPs
(`doc/store-layout.md`, "Run inputs"); old Stores don't. It needs no Python, never reads `derived/`, and writes
nothing into a Store.

## Destination

An installable package, `plasmotrackr`, in this repo's `R/` directory, for R and tidyverse users outside this codebase.
It replaces the 2021 scripts now in `R/`, which are deleted when the package is built. It also produces the
thesis-rerun fit summary, replacing `pt-analysis summarize` and `cohort`. It is installed from GitHub. CRAN is not a
goal, but nothing here rules it out.

## Terms

Language follows `CONTEXT.md`. The package relies on these entries there: Run, Analysis, Draw, Phase, Window,
Trace, Block, Store, Parent set, Latent parent, Latent parent probability, Edge marginal, Composition class, MAP Parent
set, Credible collection, Functional, Spread, Disputed Infection, Direction flip. In particular:

- One Run of an Analysis is a **Run**, never a chain. The column is `run`.
- A **Draw** is (`run`, `iteration`). No table has a separate draw index.
- The Latent parent is never an edge. It never appears as a `parent`, and its `<id>_copy` id never appears in output.
  Its share of an Infection's Parent sets is `p_latent_parent`, and on a set it is `has_latent_parent`.
- Each function's docs give the manuscript symbol it estimates (ψ_edge, ψ_MAP, ψ_lo and ψ_hi, q_i). The docs also map
  the manuscript's "unobserved source" names to the Latent-parent names used here.

## Shape of the package

From [How should a plasmotrackr session look?](https://github.com/EPPIcenter/plasmotrack/issues/102)
([sketch](https://github.com/EPPIcenter/plasmotrack/blob/prototype/plasmotrackr-session/R/PROTOTYPE-session.R), never
merged) and [Which derived quantities ship in v1, with which signatures?](https://github.com/EPPIcenter/plasmotrack/issues/104).

- **Per Draw is lazy, over Draws is eager.** A table with a row per Draw is a dbplyr lazy table the user can filter
  before `collect()`. Anything averaged or ranked over Draws is computed and returned as a tibble. `pt_interval()` is
  the explicit step from lazy to eager. Nothing is cached: users reuse a lazy table with `dplyr::compute()`.
- **Naming.**
  - Every function starts with `pt_`.
  - A lazy table with one row per Draw per key ends in `_draws`.
  - A functional (one value per Draw, a mean over Infections) has no suffix.
  - Raw Store tables come only through `pt_table()`, so every other named function is something the package computes.
  - There is one function per functional and no bundle function.
- **Columns.** `run` is the first column of every table that has it, and `iteration` the second.
- **Common arguments.**
  - `runs = NULL` is the only way to narrow Runs.
  - `phase = "sampling"` is taken only by `pt_table()`, the lazy per-Draw tables and the Block summaries
    ([How does plasmotrackr treat replica-exchange, resumed, killed and --start-from Runs?](https://github.com/EPPIcenter/plasmotrack/issues/110)).
    Pooled summaries, scoring and convergence always use sampling Draws.
  - `by_run = FALSE` on summaries gives pooled values. `TRUE` gives one row per Run per key with the same columns.
    One tibble never holds both.
  - `by =` is reserved for epi strata.

## Reading an Analysis

An Analysis is its Runs. There is no manifest: the user's script is the manifest, and a cohort is a loop over a
directory layout.

```r
esw <- pt_read("~/runs/eswatini-prod", epi = epi)        # finds the Run subdirectories
esw <- pt_read(c(a = "run1", b = "run2"), epi = epi)     # names become Run labels
esw <- pt_read("~/runs/eswatini-prod", input = "eswatini.json.gz")  # the input JSON, when pt_read() can't find it
pt_close(esw)
```

- **`pt_read(paths, epi = NULL, input = NULL, memory_limit =, temp_directory =, threads =)`** takes Run directories,
  or one parent directory whose Run subdirectories it finds. `epi` is a data frame keyed by Infection id.
  - It reads only the Window markers, the manifests, `structure/`, events and the Draw count per Run. It reads no
    Traces.
  - It checks that the tables and columns the package needs are present, and errors naming what is missing. There is
    no format-version check.
  - It errors when Runs differ in `input_sha256` or `null_model`. Differences in `thin`, `burn_in`, `samples`, `chains`
    or engine version are allowed.
  - It warns when Draw counts differ, when a Run is not finished, and when a Run has `non_finite` > 0.
  - **The input JSON** holds what the Store doesn't: observation times and each Infection's `allowed_parents`
    ([#139](https://github.com/EPPIcenter/plasmotrack/issues/139)). `pt_read()` looks for it in this order, gzipped
    or not:
    1. `input =`. Its SHA-256 must equal every Run manifest's `input_sha256.input`; a mismatch errors.
    2. The first Run's own copy, `input.json` or `input.json.gz` in the Store root, which newer engines write
       ([#140](https://github.com/EPPIcenter/plasmotrack/issues/140)).
    3. The input path in the first Run's manifest `cli_arguments` (`-i`, `--input` or the first positional). A
       relative path is resolved from the directory the engine was started in, when the Run directory's path ends
       with the recorded output directory, then from R's working directory.

    Steps 2 and 3 use a file only if it exists and its hash matches. When none is found, observation times are NA,
    `pt_candidates()` errors, and `pt_read()` gives one message listing what it tried and how to supply the input.
  - It prints one line per Run, plus the epi match count, the input's path (or that it wasn't found) and the engine
    settings.
- **The object** is an S3 `pt_analysis` list. It holds the DuckDB connection, the Run table, the Infections, the Loci,
  epi, and the input's path and candidate pairs.
- **Engine.** One DuckDB connection per Analysis, opened with `duckdb(shared_home = FALSE)` and closed by `pt_close()`
  or garbage collection.
  - Memory cap, spill directory and threads are R options (`plasmotrackr.memory_limit`,
    `plasmotrackr.temp_directory`, `plasmotrackr.threads`), overridable on `pt_read()`. `temp_directory` defaults to
    `tempdir()`.
  - The cap bounds each Analysis, not the session, and the docs say so.
  - No shared connection: comparisons across Analyses join small pooled tibbles in R.
  - DuckDB through `dbplyr`, never `arrow` or `duckplyr`
    ([Which R engine reads a Run store lazily within a memory cap?](https://github.com/EPPIcenter/plasmotrack/issues/99)).
- **Reader rule.** Each table is read from a file list built from the `_closed/` markers, never a glob, because a DuckDB
  glob reads `.in-flight.parquet` files. `parent_set_members` is read over every closed Window in both Phases.
  `window` must be quoted in literal SQL.
- **`pt_runs(a)`**, the Run table, has these columns:
  - run, label, path, draws, first_iteration, last_iteration, windows
  - status: `finished` when the last event is `run_ended` with `completed`, `interrupted` for `interrupted`, otherwise
    `incomplete`
  - resumes: the count of `run_resumed` events
  - started_from: the source `run_id`, or NA
  - thin, burn_in
  - sampling_seconds and iterations_per_second, from sampling-Window `window_closed` events
  - non_finite: the count of `non_finite_value` events
- **`pt_infections(a)`** returns infection, observation_time, symptomatic and latent_parent. observation_time is a
  double from the input JSON, NA when the input wasn't found. **`pt_loci(a)`** returns locus and n_alleles.
- **`pt_candidates(a)`** returns the candidate pairs, parent and child, from the input's `allowed_parents`, sorted by
  child then parent. It errors when `pt_read()` found no input.

**Runs that aren't one uninterrupted chain** (resolved in the edge-cases ticket above):
- **Replica exchange:** Draws are Hot-chain rows only. Tempered chains are reachable only through
  `pt_table(a, "chains")`.
- **Resumed Runs:** a resumed Run is one Run with continuous Draws.
- **Killed or running Runs:** their closed Windows are read.
- **Runs with 0 sampling Draws:** they stay in the Run table, and pooled summaries and diagnostics skip them with a
  message.
- **`--start-from` Runs:** an independent Run with its own burn-in. `pt_read()` never follows `started_from` and never
  joins Stores into one Run.

## Raw tables

`pt_table(a, name, runs = NULL, phase = "sampling")` returns a lazy table over the closed Windows of every Run, unioned,
with `run` first and the Store's own columns otherwise. Names are the Store's table names, unchanged. An unknown name
errors and lists the valid ones. The help page links `doc/store-layout.md`. `phase` may be `"burn_in"`, `"sampling"`
or both, and burn-in includes the iteration-0 starting State.

Allele frequencies, error rates and durations get no dedicated functions:
`pt_table(a, "allele_frequencies") |> pt_interval(by = c(locus, allele))` covers them.

## Per-Draw tables

All lazy, with `run` and `iteration` first, and taking `runs` and `phase`.

| Function | Columns after run, iteration |
|---|---|
| `pt_edge_draws` | child, parent, prob. Observed parents only |
| `pt_infection_draws` | infection, p_latent_parent, expected_observed_parents, entropy, effective_size, multi_source_augmented, multi_source_observed |
| `pt_class_draws` | infection, class, n_observed, has_latent_parent, share. Every class is zero-filled for every Infection |
| `pt_out_degree_draws` | infection, out_degree. A full zero-filled grid of Infections as sources |
| `pt_coi_draws` | infection, coi: the maximum over genotyped Loci of the Alleles present |
| `pt_alleles` | infection, locus, allele: one row per Allele present. Bitstrings are expanded in SQL |
| `pt_order_pair_draws` | i, j, before. i < j in string order. O(n²): filter first |

Per Draw and Infection, with `w` the stored `prob` of each of its Parent sets:
- `p_latent_parent = Σ w·[has_latent_parent]`
- `expected_observed_parents = Σ w·n_observed`
- `entropy = −Σ w ln w` over w > 0
- `effective_size = exp(entropy)`
- `multi_source_augmented = Σ w·[n_observed + has_latent_parent > 1]`
- `multi_source_observed = Σ w·[n_observed > 1]`
- `share` of a class `= Σ w·[class]`

## Functionals and intervals

Functionals are lazy per-Draw tables: a mean over Infections of a per-Draw value.

| Function | Value columns | Manuscript |
|---|---|---|
| `pt_mean_out_degree` | mean_out_degree | ψ_edge |
| `pt_mean_p_latent_parent` | mean_p_latent_parent | mean q |
| `pt_composition` | class, share (long) | composition classes |
| `pt_multi_source` | augmented, observed_only | multi-source probability |
| `pt_entropy` | mean_entropy, mean_effective_size (mean of exp(H)) | parent-set entropy |

**`pt_interval(x, by = NULL, value = NULL, level = 0.95)`** takes any per-Draw table with `run`, `iteration` and
numeric value columns, including tables the user builds. It returns the `by` keys followed by estimate, lower, upper,
level and draws. The estimate is the mean pooled over every Draw of every Run, and the interval is equal-tailed type 7.

```r
pt_mean_out_degree(esw) |> pt_interval()
pt_composition(esw) |> pt_interval(by = class, level = 0.9)
pt_mean_out_degree(esw, by = "district") |> pt_interval(by = stratum)
```

## Pooled summaries

All eager tibbles. `by_run = TRUE` gives one row per Run per key.

- **`pt_edge_marginals(a, floor = 0)`**: child, parent, marginal, lag. Pooled Edge marginals over every Draw.
  - Lag is the child's observation time minus the parent's.
  - With `by_run`, an edge seen in any Run gets a 0 row in each Run that lacks it.
- **`pt_parent_sets_pooled(a, level = 0.8)`**: infection, members, n_observed, has_latent_parent, class, weight, rank,
  cumulative_mass, credible.
  - `members` is a list of observed parent ids.
  - Sets are matched across Runs by (observed members, has_latent_parent).
  - `weight` is pooled over Draws.
- **`pt_map_graph(a)`**: the rank-1 rows of `pt_parent_sets_pooled()`: each Infection's MAP Parent set.
- **`pt_map_out_degree(a, by = NULL)`** (ψ_MAP): map_out_degree and infections.
- **`pt_out_degree_range(a, level = 0.8, by = NULL)`** (ψ_lo, ψ_hi): lo, hi and level.
  - ψ_lo(S) = n_S⁻¹ Σ_j min over U in j's Credible collection of |U_obs ∩ S|, and ψ_hi the same with max.
  - It reduces to the all-Infections definition and always contains ψ_MAP(S).
- **`pt_infection_summary(a)`**, one row per Infection:
  - observation time and p_latent_parent
  - expected observed parents and pooled class shares
  - mean_entropy and mean_effective_size (per Draw), posterior_entropy and posterior_effective_size (of the pooled
    weights; the docs explain the difference)
  - MAP set: members, class, weight
  - top observed parent and its marginal (ties to the smaller id)
  - expected out-degree and COI median

  It has no R-hat or ESS. Those are in `pt_convergence()$infections`, and the docs show the join.
- **`pt_coi(a)`**: infection, median and distribution (a list column of COI value and fraction).
- **`pt_genotype_presence(a)`**: infection, locus, allele, presence. Pooled, with 0 for Alleles that never appear.
- **`pt_genotype_mode(a)`**: infection, locus, alleles (the bitstring), fraction, agreement.
  - The mode is the bitstring with the most Draws, with ties going to the smaller bitstring.
  - Count bitstrings before expanding them: joining every row to its Alleles first cost 25 GB on Eswatini.
- **`pt_infection_times(a, level = 0.95)`**: infection, lower, median, upper. Observation time minus the type-7
  quantiles of `infection_duration`.
- **`pt_order_pairs(a, infections = NULL)`**: i, j (i < j in string order), p_before and draws.
- **`pt_co_inclusion(a, infections = NULL)`**: child, parent_a, parent_b (a < b) and probability: the pooled weight of
  the child's sets that hold both. Per recipient, from sets with at least two observed parents.

## Out-degree posteriors

From the 2026-11-03 amendment to the derived-quantities ticket. `pt_interval(pt_mean_out_degree(x))` covers only the
variance between Draws. Given a Draw, children's Parent sets are independent, so the out-degree's distribution given a
Draw is exact.

- **`pt_out_degree_dist(a, infections = NULL, level = 0.95, pmf = FALSE)`** gives one row per Infection: infection,
  estimate, lower, upper, level, var_within, var_between, p_zero and draws. It is Poisson-binomial given a Draw, over
  the Infection's candidate children. With `pmf = TRUE` the output is long: infection, k, probability.
- **`pt_realized_out_degree(a, by = NULL, level = 0.95)`**: stratum, n, estimate, lower, upper, level, var_within,
  var_between, draws, on the mean scale. It counts sources as `pt_mean_out_degree` does, and reads each Run once
  whatever the number of strata.
- **`pt_sample_networks(a, n)`**: sample, run, iteration, infection, members, has_latent_parent.
  - Each network picks one Draw uniformly, then one Parent set per child from its weights.
  - It uses R's RNG, so `set.seed()` reproduces it.
  - It serves functionals that aren't sums over children.

Conventions for these three:
- **Mixture:** the posterior is the Draw-weighted mean of the per-Draw distributions.
- **Weights:** used as stored. A child's count is 0 with the weight the other counts leave, so `estimate` equals
  `pt_mean_out_degree`'s pooled mean. Sampling normalises a child's stored weights.
- **Intervals:** equal-tailed quantiles of the discrete mixture: the smallest count with cumulative probability
  ≥ p − 1e-12. These are not type 7, and an interval need not contain the estimate.
- **Variance:** `var_within + var_between` is the posterior variance. var_between is what `pt_interval()` sees.
- **Support:** the convolution is capped at 12 conditional sd plus 12 above the largest per-Draw mean, with the mass
  beyond kept in the last entry.
- **Latent parent:** never counted.

## Strata and flows

From [How does plasmotrackr define attribution, flows and strata without an epi value?](https://github.com/EPPIcenter/plasmotrack/issues/109).

- `by =` takes one epi column and adds a `stratum` column. Cross-strata come from a combined column the user adds to
  `epi`.
- Infections with an NA value, or missing from `epi`, form an NA stratum. It is kept in every output, so per-stratum
  values weighted by `n_S` average back to the overall value.
- **Which functions take `by =`:** `pt_mean_out_degree`, `pt_map_out_degree`, `pt_out_degree_range`,
  `pt_mean_p_latent_parent`, `pt_composition`, `pt_multi_source`, `pt_entropy`, `pt_realized_out_degree`, and the flow
  functions below.
- **Who is counted:** out-degree estimates count a stratum's Infections as sources. Recipient-side functionals count
  recipients.

Flow functions, all lazy and per Draw:

- **`pt_attribution(a, by)`**: recipient, source, latent, share.
  - Per recipient stratum, its Infections' attributed parent mass (expected observed parents plus Latent parent
    probability) is split by source stratum and normalised within the recipient at each Draw.
  - The Latent parent row has `latent = TRUE` and `source = NA`.
- **`pt_enrichment(a, by, candidates = pt_candidates(a))`**: source, recipient, enrichment.
  - Enrichment is the cell's share of edge mass over its share of candidate pairs, and NA where the cell has no
    candidate pair.
  - `candidates` is a `parent, child` data frame. It defaults to `pt_candidates(a)`, the input's `allowed_parents`,
    from the Store's embedded copy of the input or, for old Stores, the input file.
- **`pt_net_flow(a, by)`**: stratum, net_flow: expected out-edges minus in-edges, counting only edges between different
  strata.
- **`pt_reassign_long_lags(a, max_lag, by = NULL)`**: mean_out_degree and mean_p_latent_parent after every edge with
  lag > `max_lag` is moved to the Latent parent.
  - Each Infection's Latent parent probability is clamped at 1 at each Draw.
  - Edges with a missing lag are never moved.

The lag distribution (a ggplot of `pt_edge_marginals()` with `weight = marginal`) and spatial scale (a join with epi)
get no functions. The vignette shows both.

## Convergence

From [Which tree convergence diagnostics ship, and how are they defined?](https://github.com/EPPIcenter/plasmotrack/issues/103),
with facts from [What does posterior offer for multi-Run convergence diagnostics?](https://github.com/EPPIcenter/plasmotrack/issues/100).

**`pt_convergence(a, runs = NULL, iterations = NULL, edge_floor = 0.01, rhat = 1.01, ess = NULL, spread = 0.2)`**
returns an S3 `pt_convergence`, a list of eager tibbles:

- `runs`: run, Draws used, Run-mean log posterior, and its gap to the best Run.
- `quantities`: group, name, status, rhat, ess_bulk, ess_tail, flagged.
- `infections`: infection, spread, class, top parent per Run (a named list column), worst R-hat, flagged.
- `pairs`: i, j, P(i before j) per Run, spread, flagged.

**Draw set.**
- Every diagnostic uses each Run's first *n* sampling Draws, where *n* is the shortest Run's count.
- `pt_convergence()` errors when Runs have different thinning strides.
- `iterations = c(first, last)` restricts every diagnostic to a range of iterations.
- Pooled estimates elsewhere still use every Draw.

**Estimators.**
- `posterior::rhat`, `ess_bulk` and `ess_tail` on every quantity.
- `ess_tail` is NA for 0/1 quantities, and the docs say why.

**Status.** Each quantity is classified before `posterior` is called:
- `constant`: one value over every Draw of every Run. R-hat is NA and the quantity is never flagged.
- `stuck`: every Run is constant but the Runs differ. R-hat is `Inf` and the quantity is always flagged.
- `mixing`: everything else. R-hat comes from `rhat()`, falling back to `rhat_basic()` when `rhat()` is NA (the 50/50
  fold case).

**Quantity groups.** log posterior (prior + likelihood), global scalars, error rates, durations, Allele frequencies,
edges, `p_latent_parent`, and the functionals.
- Edges and order pairs are included only when some Run's mean is ≥ `edge_floor`.
- Order pairs are further limited to pairs joined by such an edge in either direction.
- Quantities below the floor still count toward spread.

**Tree diagnostics.** These use zero-filled Edge marginals, with the Latent parent as one more parent key.
- *Run mean*: a quantity's mean over one Run's truncated Draws.
- *Spread*: max − min of a quantity's Run means.
- *Infection spread*: the max of spread over the Infection's parent keys.
- *Order-pair spread*: spread of the indicator "i before j".
- *Log-posterior gap*: spread of the log posterior.
- *Top parent* in a Run: the key with the largest Run mean. Ties go to the smaller id, and an observed parent wins a
  tie against the Latent parent.
- `top_parent` labels the Latent parent `"latent"`. Any string can be an Infection id, so no label is safe; instead
  `pt_convergence()` errors when an Infection's id is `"latent"`.

**Disagreement classes.** Each flagged Infection gets the first matching rule, with *t* the spread threshold:
1. `same_top`: every Run has the same top parent.
2. `root_flip`: some Run's top is the Latent parent, and in that Run the child has edge mean > *t* into another Run's
   top parent.
3. `latent_vs_observed`: one Run's top is the Latent parent, another's is observed, and rule 2 doesn't apply.
4. `direction_flip`: in some Run, another Run's top parent has the child as a parent with mean > *t*.
5. `linked_parents`: the top parents differ, and in some Run one of them infects another with mean > *t*.
6. `different_parent`: everything else.

**Flags.** Each is a single-tier `flagged` column:
- R-hat > 1.01
- `ess_bulk` < 100 × the number of Runs (the default when `ess = NULL`)
- Infection spread and order-pair spread > 0.2
- `stuck`: always flagged

The Log-posterior gap is reported but not flagged on its own.

**`print()`** shows the Runs and Draws used, the gap, flagged counts per group, flagged Infections by class, and the 10
worst Infections.

**`pt_as_draws(a, groups =)`** builds a `posterior` draws array, truncated to the shortest Run with a message. It is the
escape hatch: `posterior::rstar()` on it is shown in the docs, since R* needs `caret` and isn't in v1.

**Dropped:** parent-set TVD, edge switches, a separate top-parent agreement, R*, and matching `pt_analysis`'s ESS.

## Sampler Blocks

From [Do sampler Blocks get plasmotrackr functions in v1?](https://github.com/EPPIcenter/plasmotrack/issues/111).

- **`pt_acceptance(a, level = c("type", "sampler"), phase = "sampling", by_run = FALSE, per_window = FALSE)`**: type
  (and sampler), accepted, proposed, rate. Rate is NA when nothing was proposed. Counts are summed over Windows, and
  over Runs unless `by_run`.
- **`pt_swaps(a, phase = "sampling", by_run = TRUE, per_window = FALSE)`**: temperature_pair, accepted, attempted, rate.
  Empty for one-chain Runs.
- No acceptance band, sampler plots, ladder or proposal-variance views, or events parser. `pt_table(a, "events")`
  returns payloads as JSON strings.

## Truth and scoring

From [How does plasmotrackr read simulation truth and score against it?](https://github.com/EPPIcenter/plasmotrack/issues/105).
Metric definitions are entries M1–M16 of the
[pt_analysis inventory](https://github.com/EPPIcenter/plasmotrack/blob/research/pt-analysis-inventory/docs/research/pt-analysis-inventory.md)
unless stated here.

**`pt_read_truth(path)`** reads the simulator JSON, gzipped or not, into an S3 `pt_truth` made of tibbles:
- `edges`: parent, child
- `infections`: infection, infection_time, observation_time, duration, symptomatic
- `genotypes`: infection, strain, locus, allele, long, in simulator numbering
- `loci`
- `observed`: infection, locus, allele, each case's observed genotype in simulator numbering (for the Allele mapping)

Ids are strings. Truth is restricted to the Analysis's Infections only when scoring, so one `pt_truth` serves both the
full and the partial fit of a simulation.

**Allele mapping.**
- **The rule:** an input's Alleles at a Locus are the simulator's Alleles observed in any of that input's Infections,
  in simulator order. This held for every renumbered Locus checked: 683 in full inputs and 185 in partial ones.
- **Check:** scoring errors, naming the Locus, when the mapped count differs from the Locus's `n_alleles`.
- A Locus whose input kept the simulator's Allele count keeps the simulator's numbering.
- **Unrepresentable Alleles:** true Alleles that no Infection shows are counted as `unrepresentable` and excluded from
  genotype accuracy.

**True Parent set.** An Infection's true parents among the Analysis's Infections, plus the Latent parent when the
Infection is a root or any true parent is outside the Analysis. In a `missing_50` fit, a dropped parent counts as the
Latent parent. Comparison uses observed members plus `has_latent_parent`.

**Scoring functions.**
- **`pt_score_edges(a, truth, edges = c("direct", "transitive"), curve = FALSE)`**: AP, PR-AUC and top-k precision, directed and
  undirected, from zero-filled pooled Edge marginals.
  - Directed universe: every ordered pair.
  - Undirected universe: the sum of both directions, clipped to [0, 1].
  - Scores are rounded to 10 places before ranking, and tied pairs enter the curve together.
  - `transitive` counts an observed pair joined by a path through only unobserved cases as true.
  - It returns one row per direction (directed, undirected): direction, n_edge_universe, n_true_edges,
    n_distinct_scores, average_precision, pr_auc, top_k_precision, as the fit summary names them. With `curve = TRUE`
    it returns the PR curve instead: direction, threshold, recall, precision.
- **`pt_score_parent_sets(a, truth, level = 0.8)`**: one row per Infection: infection, credible_size, covered,
  top_probability. The fit summary's credible size, coverage and top-probability columns are its means and median.
- **`pt_score_map(a, truth)`**: per Infection, whether the MAP Parent set equals the true Parent set (`correct`); the
  share that do is `mean(correct)`.
- **`pt_score_genotypes(a, truth)`**: per (infection, locus), true Alleles against posterior presence called at 0.5,
  plus the unrepresentable count.
- **`pt_score_coi(a, truth)`**: posterior COI median against the identifiable true COI, meaning the most distinct true
  genotypes at any observed Locus.
- **`pt_score_pairs(scores, truth, infections)`**: the IBD baseline, from a data frame of `a, b, score`.
  - Each unordered pair keeps its largest score, and ids outside `infections` are ignored.
  - It returns the directed ranking columns (each pair scores the same both ways) and the undirected ones.

**Fit summary.** `pt_fit_summary(a, truth, fit_type = "posterior", missingness = "full", factors = NULL, scores =
NULL)` returns one row with the 70 columns of `simulation_fit_summary.rds`, names unchanged, so the existing table and
figure scripts read it as is. The committed `fit-summary-columns.csv` holds the layout.
- `factors` is a one-row data frame bound in front of the other columns. With `NULL`, the 12 thesis factor columns are
  `NA`. Its `missingness_condition` and `fit_type` columns stand in for those arguments when they are not given.
- Parsing the thesis directory layout into the 14 factor columns is `tools/plasmotrackr/thesis-factors.R`, not part of
  the package.
- The null model is its own row: `pt_fit_summary()` on the null Analysis with `fit_type = "null_prior"`.
- With `fit_type = "dcifer"`, the `pt_score_pairs()` columns of `scores` over the Analysis's Infections fill the layout
  with NA functionals.
- The functional columns are built from the functions above. The observed-only and mixed composition columns are
  rollups of `pt_composition`, and `mean_no_parent_probability` is always 0.
- There is no cohort function: users `bind_rows` the per-fit rows.
- On a SLURM cluster it runs as an array, one simulation per task, each writing a one-row file.

## Plots

From [Which plots ship in v1, with which signatures?](https://github.com/EPPIcenter/plasmotrack/issues/107).

**Contract.**
- Every plot returns a ggplot object, with no printing and no file output.
- Plots use ggplot2's default theme and scales, with colour mapped to `run`.
- A plot given a lazy table collects only the columns it needs.
- ggplot2 is checked with `rlang::check_installed()`.
- No bayesplot.

**Functions.**
- **`pt_plot_trace(x, value =, by = NULL)`**: one line per Run over iteration, faceted by `by`. Burn-in is shaded when
  the table carries it.
- **`pt_plot_rank(x, value =, by = NULL, bins = 20)`**: rank histograms per Run.
- Both error past `max_series = 25` keys, with a message telling the user to filter first.
- **`pt_plot_interval(x, truth = NULL)`**: pointranges from any table with estimate, lower and upper. That covers
  `pt_interval()` and the out-degree posteriors.
- **`autoplot.pt_convergence(object, type = c("runs", "rhat", "spread", "classes"))`**: registered lazily for ggplot2's
  generic. `plot()` calls it with `"runs"`.
- **`pt_plot_run_pairs(x, value =)`**: each pair of Runs' values from a `by_run = TRUE` summary, with the y = x line.

Cohort figures from bound fit-summary rows stay in the scripts that make them. COI, genotype and network plots are
vignette snippets.

## Network export

Networks are exported, not drawn. `pt_as_igraph(a, floor = 0.05, graph = c("marginal", "map"))` and
`pt_as_tbl_graph()`, with the same arguments, build a directed graph:
- **Vertices:** every Infection, with `pt_infection_summary()` columns as attributes.
- **Edges, `"marginal"`:** pooled Edge marginals at or above `floor`, parent → child (tidygraph's `from` and `to`),
  with the marginal and lag as attributes.
- **Edges, `"map"`:** the MAP graph's observed-parent edges.

igraph and tidygraph are in Suggests.

## Conventions

These apply everywhere unless a section says otherwise.

1. **Draw count.** `D_r` is the number of distinct sampling iterations in Run r's `run_scalars`. It is the denominator
   of every per-Run mean.
2. **Zero-filling.** A set, edge or class with no row at a Draw has weight 0 there. Every mean over Draws divides by
   `D_r`.
3. **Pooling.** A pooled value sums over every Draw of every Run and divides by `Σ_r D_r`, so each Run weighs by its
   Draws.
4. **Weights as stored.** The Store drops sets below `parent_set_floor` without renormalising, and so does the package.
   An Infection's weights may sum to slightly less than 1. The docs bound the error.
5. **Quantiles.** Type 7 everywhere, the COI median included. The out-degree posteriors are the one exception.
6. **Ranking tie-break.** Weights tie when they agree to 12 significant digits (`signif(w, 12)`): pooled weights
   are sums, and summation order changes their last bits. Ties in weight go to the members list compared element by
   element. The list is the observed parents in id string order, with the Latent parent last. This applies to the MAP
   Parent set and `rank`.
7. **Credible collection.** The highest-ranked sets whose cumulative weight reaches ≥ level − 1e-12. It also admits
   every set tied in weight (as in 6) with the last set needed, and holds every set when the weights never reach the
   level. Default level 0.8.
8. **Levels.** Credible collection and size range default to 0.8. Intervals and infection times default to 0.95.
9. **Ids** are strings and sort as strings (`"10" < "2"`).
10. **Bitstrings.** Character k of an `alleles` bitstring is Allele k.

## Testing

From [How is plasmotrackr tested?](https://github.com/EPPIcenter/plasmotrack/issues/106). `pt_analysis` is not a
reference: agreement with it is not a goal.

- **Oracle.** `helper-oracle.R` reads raw Parquet with `nanoparquet` (no DuckDB or dbplyr) and implements the formulas
  as plain base-R loops.
  - It runs live, so no expected numbers are committed for real Stores.
  - Comparisons use relative 1e-12, and quantiles must match exactly.
- **Edge cases.** `fake_store()` in `helper-fake-store.R` writes a valid Run store from small R data frames, with Parquet
  written through DuckDB `COPY`. Each case gets a store small enough to check by hand, with hand-computed expected
  values next to the oracle check:
  - sets dropped by the floor
  - weights that never reach the level
  - exact ties, in the MAP rule and at the credible cutoff
  - COI above 1
  - Ordering flips
  - unequal Runs, and more than two Runs
  - ids where `"10" < "2"`
  - ungenotyped Loci
  - constant and stuck quantities
  - different thinning strides
- **Golden stores.** `test/data/golden/` is copied into `R/inst/extdata/golden/` (about 300 KB). A CI job fails when the
  two copies differ, which is how an engine format change reaches the R tests.
- **Reader rule.** A helper copies a golden store to `tempdir()` and adds an unmarked Window and a `.in-flight.parquet`
  file. The test checks that neither is read.
- **Fit summary.** Its layout (names, order, types) is checked against a committed `fit-summary-columns.csv` made from
  the rds. Its values are checked by the oracle on a golden store with a hand-written truth JSON.
- **Out-degree posteriors.** The oracle enumerates every combination of one Parent set per child at every Draw.
- **Scale.** CI runs no scale tests. A unit test checks that `memory_limit`, `temp_directory` and `threads` are applied.
  Manual tests in `tests/manual/` run when `PLASMOTRACKR_SCALE_STORE` is set: every function completes under a 4 GB
  cap, and timings are reported.
- **Layout.** testthat 3e, one file per family:
  - `test-read`, `test-reader-rule`, `test-tables`
  - `test-edges`, `test-parent-sets`, `test-functionals`, `test-interval`, `test-strata`
  - `test-genotypes`, `test-order`, `test-convergence`, `test-truth`, `test-scoring`, `test-fit-summary`

  `print()` methods use snapshot tests.
- **CI.** GitHub Actions `r-lib/actions` `check-r-package`:
  - **Platforms:** macOS, Ubuntu and Windows on R release, plus Ubuntu on oldrel-1.
  - **Trigger:** changes to `R/**` or `test/data/golden/**`.
  - **Sync job:** a separate job runs the golden-copy check.

  CI never builds the C++ engine.

## Documentation and example data

From [Which vignettes and example data ship with plasmotrackr?](https://github.com/EPPIcenter/plasmotrack/issues/112).

**Vignettes.**
1. **Getting started**, a package vignette on the golden stores: `pt_read`, the Run table, `pt_table`, lazy versus
   eager, functionals with `pt_interval()`, and network export.
2. **Convergence**: `pt_convergence()`, its plots, disputed Infections, and what to do when Runs disagree.
3. **Scoring a simulation**: truth, scoring, `pt_fit_summary()`, and a cohort loop as a SLURM
   array.
4. **Strata and flows**: `by =`, attribution, enrichment, net flow, long-lag reassignment, the lag histogram and spatial
   scale.

Vignettes 2–4 are pkgdown articles in `vignettes/articles/`, Rbuildignored. They evaluate on downloaded data when the
site builds, so `R CMD check` and installs never download anything. The vignettes also carry these snippets:
- COI and genotype plots
- the Edge-marginal ggplot
- `posterior::rstar()` on `pt_as_draws()`
- event payloads read with `jsonlite`
- low acceptance by sampler type

**Example data.** **`pt_example()`** downloads a GitHub release asset of about 16–20 MB into
`tools::R_user_dir("plasmotrackr", "cache")`. It holds:
- one thesis simulation's four full-model Runs (200 Infections, MOI 2, medium diversity, without checkpoints)
- its null-model Run
- the truth JSON and the input JSON
- a synthetic epi data frame with district and travel columns

No real-cohort data (Eswatini, Zanzibar, He01) is shipped or downloaded.

**pkgdown.**
- **Hosting:** GitHub Pages, built by a workflow on pushes to `master` that touch `R/`.
- **Reference index:** grouped as reading, raw tables, per-Draw tables, functionals, summaries, out-degree posteriors,
  strata and flows, convergence, sampler Blocks, truth and scoring, plots, and export.

## Package metadata and install

From [What are plasmotrackr's Imports, Suggests and supported platforms?](https://github.com/EPPIcenter/plasmotrack/issues/113).

```
Package: plasmotrackr
Version: 0.1.0
Depends: R (>= 4.2.0)
Imports: DBI, duckdb (>= 1.5.6), dbplyr (>= 2.6.0), dplyr (>= 1.1.0), tibble, rlang, cli, jsonlite, posterior
Suggests: ggplot2, igraph, tidygraph, testthat (>= 3.0.0), nanoparquet, withr, knitr, rmarkdown
License: MIT + file LICENSE
```

- **Size:** 33 packages are installed recursively.
- **People:** Maxwell Murphy is author and maintainer.
- **Links:** `URL` lists the repo and the pkgdown site, and `BugReports` the repo's issues.
- **License file:** the package has its own R-style `LICENSE`, because the repo-root `LICENSE` file is empty in the
  working tree.
- **Install:** `pak::pak("EPPIcenter/plasmotrack/R")`, or `remotes::install_github("EPPIcenter/plasmotrack",
  subdir = "R")`.
- **Supported platforms:** macOS (arm64 and x86_64), Linux x86_64 and Windows, all in CI. On Linux clusters without
  internet on compute nodes, install from Posit Package Manager's binaries on a login node into a user library.

## Build order

Each step ends with tests passing and can become one or more implementation tickets.

1. Package skeleton, DESCRIPTION, CI (check plus golden sync), `fake_store()`, the oracle helper and the golden copy.
2. `pt_read` and its engine options, the reader rule, the Run table, `pt_table`, `pt_infections`, `pt_loci` and
   `pt_close`.
3. Per-Draw tables and functionals, plus `pt_interval`.
4. Pooled summaries: Edge marginals, Parent sets, MAP, size range, infection summary, genotypes, COI, times, order
   pairs and co-inclusion.
5. `pt_convergence`, `pt_as_draws`, and the convergence `print` method.
6. Truth, scoring and `pt_fit_summary`, then the thesis factor script in `tools/` and a dry run on a few
   simulations.
7. Strata and flows, and the out-degree posteriors.
8. Plots, network export and sampler Blocks.
9. Vignettes, `pt_example()` and its release asset, and the pkgdown site. Delete the 2021 `R/` scripts.

## Decisions

Every decision lives on its ticket. In map order:
- [What does pt_analysis compute that plasmotrackr must match?](https://github.com/EPPIcenter/plasmotrack/issues/101)
- [What does posterior offer for multi-Run convergence diagnostics?](https://github.com/EPPIcenter/plasmotrack/issues/100)
- [Which R engine reads a Run store lazily within a memory cap?](https://github.com/EPPIcenter/plasmotrack/issues/99)
- [How should a plasmotrackr session look?](https://github.com/EPPIcenter/plasmotrack/issues/102)
- [Which tree convergence diagnostics ship, and how are they defined?](https://github.com/EPPIcenter/plasmotrack/issues/103)
- [Which derived quantities ship in v1, with which signatures?](https://github.com/EPPIcenter/plasmotrack/issues/104),
  amended twice
- [How does plasmotrackr read simulation truth and score against it?](https://github.com/EPPIcenter/plasmotrack/issues/105)
- [How is plasmotrackr tested?](https://github.com/EPPIcenter/plasmotrack/issues/106)
- [Which plots ship in v1, with which signatures?](https://github.com/EPPIcenter/plasmotrack/issues/107)
- [How does plasmotrackr define attribution, flows and strata without an epi value?](https://github.com/EPPIcenter/plasmotrack/issues/109)
- [How does plasmotrackr treat replica-exchange, resumed, killed and --start-from Runs?](https://github.com/EPPIcenter/plasmotrack/issues/110)
- [Do sampler Blocks get plasmotrackr functions in v1?](https://github.com/EPPIcenter/plasmotrack/issues/111)
- [Which vignettes and example data ship with plasmotrackr?](https://github.com/EPPIcenter/plasmotrack/issues/112)
- [What are plasmotrackr's Imports, Suggests and supported platforms?](https://github.com/EPPIcenter/plasmotrack/issues/113)

## Open questions

None. The last three (network export, the shape `pt_score_edges()` returns, per-Infection Parent-set scoring)
were settled when the spec was accepted on 2026-11-03.

## Out of scope

- CRAN submission. Nothing rules it out later.
- Old-format (gzip CSV) output, including the original thesis fits.
- The Python/React analysis GUI: it is neither replaced nor fed by this package.
- A `reticulate` wrapper around `pt_analysis`.
- Writing anything into a Store, including `derived/`, and any cache.
- Launching, resuming or tuning Runs.
