# pt-analysis

Reads plasmotrack Run stores (`doc/store-layout.md`) and Analysis manifests (`doc/analysis-gui.md`). It owns the
closed-Window reader rule: a Window exists only if its marker under `_closed/` exists, so table files of an unmarked
Window and in-flight files are never read. Each table is a DuckDB relation over the Parquet files in place. Traces are
read over the sampling Phase by default, so every row is a Draw; `parent_set_members` is read over both Phases.

## Install

Needs [uv](https://docs.astral.sh/uv/). From this directory:

```sh
uv sync
```

This creates `.venv/` with Python 3.13, DuckDB, pyarrow, Starlette, uvicorn and pytest, and installs the `pt-analysis`
command into it.

## Run the tests

```sh
./test.sh          # every test: pytest, then the frontend's type check and vitest
uv run pytest      # the Python tests only
```

The tests read the golden stores in `test/data/golden/` and check the benchmark scripts in `tools/benchmark/`, which
read stores through this package. `test.sh` needs pnpm for the frontend (`web/README.md`).

## Describe an Analysis

An Analysis is a set of Runs of one input, described by a JSON manifest. Paths are relative to the manifest file.

```json
{
  "label": "example, two Chains",
  "input": "example_input.json",
  "symptomatic_idp": "symptomatic_idp.txt",
  "asymptomatic_idp": "asymptomatic_idp.txt",
  "runs": {"chain1": "runs/chain1", "chain2": "runs/chain2"},
  "truth": "sim_0_s0.json.gz",
  "epi": "epi.csv",
  "pairwise_scores": "ibd.csv",
  "null_model": "../null/analysis.json",
  "timings": "timings.json"
}
```

`input`, the two IDP paths and `runs` are required; every Run directory must exist. The rest is optional; `label`
names the Analysis in the GUI and defaults to the manifest's directory name.

```sh
uv run pt-analysis info analysis.json
```

prints, for each Run, its closed Windows per Phase and its Draw count.

## Derive tables

```sh
uv run pt-analysis derive RUN_DIR          # one Run
uv run pt-analysis derive analysis.json    # every Run of an Analysis, then the Analysis's convergence table
```

writes Parquet tables under `RUN_DIR/derived/`. `edges`, `set_summary`, `infection_draws` and `functional_draws`
come from the sampling Phase, so every weight is a Draw's; `set_classes` comes from `parent_set_members` over both
Phases.

| Table | Columns | Rows |
|---|---|---|
| `set_classes` | infection, set_id, n_observed, has_unobserved, class | one per Parent set, from `parent_set_members` |
| `edges` | iteration, child, parent, prob | per Draw and child, one per parent: summed weight of the child's sets containing it |
| `set_summary` | infection, set_id, weight, rank, cumulative_mass | one per Parent set with weight in some Draw |
| `infection_draws` | iteration, infection, q, expected_observed_parents, entropy, effective_size, one column per class, multi_source_augmented, multi_source_observed | one per Draw and Infection |
| `functional_draws` | iteration, psi_edge, mean_q, mean_entropy, `mean_<class>` per class | one per Draw |

- `n_observed` counts a set's observed parents (members that are Infection ids). `has_unobserved` says whether the
  unobserved source, the Infection's Latent parent (`<id>_copy`), is a member.
- `class` is the composition class: `unobserved_only` for a set with no observed parent (its only member is the
  Latent parent; there is never an empty set), otherwise `observed_<n>` with `_unobserved` appended when the Latent
  parent is a member too. A Run allowing at most two observed parents has the five labels `observed_1`,
  `observed_1_unobserved`, `observed_2`, `observed_2_unobserved` and `unobserved_only`. Labels do not change between
  versions. They sort as listed (by count, each count's sets without the Latent parent first, `unobserved_only`
  last) while counts stay below 10; sort on `n_observed, has_unobserved` to order any count.
- In `edges` the unobserved source appears as a parent under the child's Latent parent id. A parent with no row at a
  Draw has weight 0 there: the store drops sets below the Run's floor. These are the per-Draw entries of the
  conditional transmission matrix.
- `weight` in `set_summary` is the pooled weight: the set's weight summed over Draws and divided by the number of
  Draws, the manuscript's estimator, not a mean over the Draws where the set has a row. Sets are ranked from 1 by
  weight descending, then set_id; `cumulative_mass` is the running sum down the ranking.
- `infection_draws` holds functions of one Infection's Parent-set weights at one Draw: `q`, the unobserved-source
  inclusion (weight of sets containing the Latent parent); `expected_observed_parents` (weight times `n_observed`,
  summed); `entropy` (minus the sum of weight times its log) and `effective_size` (exp of the entropy); the weight of
  each composition class present in the Run, in a column named by its `set_classes` label (0 where the Infection
  has no set of that class); and the multi-source probabilities, the weight of sets with more than one member
  (`multi_source_augmented`, the Latent parent counting as a member) and with more than one observed parent
  (`multi_source_observed`). The class columns sum to the Infection's weight, which falls short of 1 by at most the
  floor for each set the store dropped.
- `functional_draws` averages `infection_draws` over Infections per Draw. `psi_edge`, the edge-sum, is the mean of
  `expected_observed_parents`; `mean_q`, `mean_entropy` and `mean_<class>` are the means of the others. The 95%
  interval of a functional is the 2.5% and 97.5% quantiles of its column over Draws.

`derived/manifest.json` records the closed Windows per Phase the tables came from and the package version that wrote
them. `derive` does nothing when both match the store and this package, and replaces every table otherwise, for
example once the Run has closed another Window. The manifest is written last, so an interrupted derive is redone.

### The convergence table

An Analysis's own tables go beside its manifest in `<manifest stem>.derived/`, for example `analysis.derived/` for
`analysis.json`. Its `manifest.json` records each Run's directory and closed Windows and the package version; the
tables are rewritten when any of those changes.

`convergence.parquet` has columns group, name, chain, ess, rhat: for every quantity, one row per Chain (`chain` is
the Run's label) with that Chain's ESS, then a pooled row (`chain` null) with the multi-chain ESS and split R-hat over
the Chains with Draws. `rhat` is null on Chain rows. ESS and R-hat are null for a quantity that is constant, or
constant up to rounding, and below 4 Draws; split R-hat also needs 2 Chains.

| group | name | quantity per Draw |
|---|---|---|
| log posterior | `log posterior` | prior plus likelihood, from `run_scalars` |
| global scalars | `mean_coi`, `mean_strains_tx`, `parent_set_size_prob` | from `run_scalars` |
| false positive rate | Infection id | `expected_false_positives`, from `infections` |
| false negative rate | Infection id | `expected_false_negatives`, from `infections` |
| infection duration | Infection id | `infection_duration`, from `infections` |
| allele frequency | `<locus>[<allele>]` | from `allele_frequencies` |
| edge | `<child><-<parent>` | Edge marginal from an observed parent, from `edges`; 0 at Draws with no row |
| unobserved-source inclusion | Infection id | `q`, from `infection_draws` |
| functionals | `psi_edge`, `mean_q`, `mean_entropy` | from `functional_draws` |

Groups are in this order and names sorted within a group. The statistics, in `pt_analysis.diagnostics`, are those
the benchmark scripts used:

- `multichain_ess`: chains truncated to the shortest; within-chain variance W from autocovariances with divisor n
  scaled by n/(n-1), between-chain variance B, pooled variance (n-1)/n W + B/n; autocorrelations against the pooled
  variance, summed in pairs with Geyer's initial monotone sequence; ESS = m n / tau (Vehtari et al. 2021 without
  rank normalisation).
- `split_rhat`: each chain's two halves, every half truncated to the shortest half's length n; sqrt(((n-1)/n W +
  B/n) / W) over the halves; infinite when every half is constant but the halves differ.

## Serve the analysis GUI

```sh
(cd web && pnpm install && pnpm build)            # once, and after the frontend changes
uv run pt-analysis serve analysis.json [--port 8765] [--no-browser]
```

derives what is not current, then serves the page in `web/dist/` and a JSON API on `127.0.0.1` and opens the
browser. The page opens on the convergence table: one row per quantity with its pooled ESS, split R-hat and each
Chain's ESS, sorted by clicking the ESS or R-hat header and filtered by group and by name. Clicking a row opens the
trace panel: one line per Chain over Sampling iterations, burn-in greyed on request, a rank-normalised histogram per
Chain, and ESS and R-hat for the Draws in the range brushed on the trace.

| Endpoint | Returns |
|---|---|
| `GET /api/analysis` | label, Draws per Chain and in all, Infections, Loci |
| `GET /api/convergence` | `{rows}`: the convergence table's rows |
| `GET /api/trace?group=&name=[&burn_in=true]` | per Chain, `iteration`, `phase` and `value` lists of a quantity |
| `GET /api/range?group=&name=[&first=][&last=]` | ESS per Chain, pooled ESS and split R-hat over the Draws with Sampling iteration in [first, last] |

A quantity is a group and name of the convergence table. Burn-in rows come only on request and only for the groups
read from the store's Traces (log posterior to allele frequency); edges, unobserved-source inclusion and functionals
come from the derived tables, which hold Draws only. The trace and range queries are in `pt_analysis.traces`. The
manifest's optional `label` names the Analysis on the page (default: the manifest's directory name).

### Accuracy against truth (U6)

With a `truth` file in the manifest, the **Ranking** page draws the directed and undirected precision-recall curves
of the posterior Edge marginals against the true network, with the null-model Analysis and the pairwise scores
overlaid when the manifest names them, and a table of their metrics. The **Parent sets** page lists, per Infection,
the effective parent-set size, the credible collection's size at a level (0.8 by default), whether the true set is in
it, and the modal set's probability, beside the null-model fit, with a summary line per fit. It works without truth,
less the coverage column.

| Endpoint | Returns |
|---|---|
| `GET /api/ranking` | realised truth, and per fit (posterior, null, IBD) its directed and undirected metrics and curve; 404 without truth |
| `GET /api/parent-sets[?level=0.8]` | per fit its summary columns, and one row per Infection, `null_` columns for the null fit |

The definitions follow the thesis pipeline's per-fit summary (`remote/R/final_analysis/R`), with weights pooled as in
the rest of this package:

- **Truth** is the simulation's JSON (plain or gzipped): `network: [{src, dest}]` and `nodes: [{id, ...}]`. True
  edges are those with both ends among the Analysis's Infections. An Infection's true Parent set is its true parents
  among the Infections, plus the unobserved source when it is a root or has a true parent outside them; it matches a
  fitted set whose Latent parent `<id>_copy` is read as `<id>`. `n_true_cases` counts the truth's nodes,
  `n_observed_cases` the Infections, `true_average_outdegree` is true edges per Infection and
  `true_unobserved_source_proportion` the share of Infections whose true set has the unobserved source.
- **Edge scores** are the posterior transmission matrix's observed entries: weight summed over every Draw of every
  Chain, divided by the number of Draws (not the R pipeline's mean over Draws with a row). Scores are clipped to [0, 1]
  and rounded to 10 decimals. The directed universe is every ordered pair of distinct Infections, a missing score
  counting 0; the undirected universe is every unordered pair, scored by the sum of its two directions clipped to 1,
  and true when a true edge joins it either way.
- **Ranking metrics**, directed and undirected: the precision-recall curve at each distinct score (tied pairs enter
  together); `average_precision` = sum of the recall gained at each breakpoint times the precision there; `pr_auc` =
  trapezoid area with (recall 0, precision 1) prepended; `top_k_precision` = precision of every pair scoring at least
  the k-th highest score, k the number of true edges (ties included); `n_edge_universe`, `n_true_edges`,
  `n_distinct_scores`.
- **Baselines**: the null-model Analysis (`null_model`) is scored the same way against the same truth. The pairwise
  scores CSV (`pairwise_scores`) has columns `a, b, score`, one row per pair of Infection ids in either order (the
  largest score of a pair's rows counts; pairs without a row score 0; ids outside the Analysis are ignored), and is
  ranked undirected only.
- **Parent-set metrics**: an Infection's sets are pooled over Chains (matched by members), ranked by weight. Its
  effective size is exp of the entropy of the pooled weights, its credible collection the smallest prefix with mass at
  least the level (as `functionals.credible_collection`), its top probability the modal set's weight. The summary
  columns are `mean_credible_parent_set_size_80`, `median_credible_parent_set_size_80`,
  `credible_parent_set_coverage_80` and `mean_top_parent_set_probability`; a level other than 0.8 has the suffix
  `_a<percent>`, for example `_a095`.

`web/README.md` says how to develop the frontend and how to add a view.

## Summarise a Run

These read the derived tables, deriving them first when they are not current.

```sh
uv run pt-analysis edges RUN_DIR [--top N] [--unobserved]   # pooled Edge marginals, largest first
uv run pt-analysis map RUN_DIR [--level 0.8]                # MAP graph, psi_MAP, parent-set size range
uv run pt-analysis credible RUN_DIR INFECTION [--level 0.8] # an Infection's credible collection
```

- The posterior transmission matrix entry of an edge is the mean over Draws of its `edges` weight. `edges` leaves
  out edges from the unobserved source unless given `--unobserved`.
- The MAP graph gives every Infection its modal set (rank 1). psi_MAP is the mean over Infections of the modal set's
  `n_observed`.
- An Infection's credible collection at a level is its smallest prefix of sets by rank whose mass is at least the
  level (within 1e-12), or all its sets if their mass never reaches it.
- The parent-set size range [psi_lo, psi_hi] at a level is the min and the max `n_observed` in each Infection's
  credible collection, each averaged over Infections.

## Summarise an Analysis's Chains

```sh
uv run pt-analysis summary analysis.json [--pairs A:B ...] [--show N]
```

derives what is not current, then prints three sections:

- The diagnostics table: per group of the convergence table, the number of quantities with an ESS, the median and
  smallest pooled ESS, the share of those quantities with R-hat above 1.1, and the largest R-hat.
- With `--pairs`, a line per Chain giving, for each pair, the Draws at which A precedes B in the Ordering out of all
  Draws, and how many times that order changed between adjacent Draws. Read from the `ordering` Trace.
- Agreement on parents: per Infection, its observed parent with the largest pooled Edge marginal (the mean over
  Chains of each Chain's pooled marginal), that probability, the pooled `unobserved_only` probability, and the
  largest difference between Chains in any of those probabilities. The `N` Infections (default 15) with the largest
  spread are listed, then the median, 90th percentile and largest spread over all Infections.

`tools/benchmark/chain_summary.py RUN_DIR` writes `RUN_DIR/chain_summary.json`, an Analysis of the stores
`RUN_DIR/chain<K>`, and prints this summary of it.

## Use from Python

```python
from pt_analysis.analysis import Analysis
from pt_analysis.store import Store

analysis = Analysis.open("analysis.json")
store = analysis.runs["chain1"]                  # a Store
store.closed_windows()                           # [Window(number, first_iteration, last_iteration, phase), ...]
store.table("parent_sets")                       # DuckDB relation, sampling Phase
store.table("run_scalars", phases=("burn_in", "sampling"))
store.rows("infections")                         # list of dicts in Sampling-iteration order
store.structure("infections")                    # lookup table written at Run start
store.draws()

from pt_analysis.derive import derive, derived
from pt_analysis import functionals

derive(store)                                    # True if the tables were written
derived(store, "edges")                          # DuckDB relation over derived/edges.parquet

from pt_analysis.derive import derive_analysis, derived_analysis
from pt_analysis.diagnostics import multichain_ess, split_rhat

derive_analysis(analysis)                        # every Run, then the convergence table
derived_analysis(analysis, "convergence")        # DuckDB relation over analysis.derived/convergence.parquet
multichain_ess([chain1_series, chain2_series])
split_rhat([chain1_series, chain2_series])
functionals.top_edges(store, n=10)               # [Edge(parent, child, weight), ...]
functionals.map_graph(store)                     # {infection: RankedSet} of modal sets
functionals.psi_map(store)
functionals.size_range(store, level=0.8)         # (psi_lo, psi_hi)
functionals.credible_collection(store, "3", level=0.8)   # [RankedSet(set_id, members, weight, rank, ...), ...]

from pt_analysis import metrics

truth = metrics.analysis_truth(analysis)         # Truth: edges, parent_sets, n_true_cases, ...; None without truth
metrics.edge_scores(analysis)                    # {(parent, child): pooled Edge marginal}
metrics.edge_rankings(analysis, truth)           # {"directed": Ranking, "undirected": Ranking}
metrics.ranking_fits(analysis, truth)            # posterior, null and IBD Rankings, as the manifest names them
metrics.rank(scores, truth_flags)                # Ranking of any scores: curve, average_precision, pr_auc, ...
metrics.parent_set_rows(analysis, truth, 0.8)    # [InfectionParentSets(infection, effective_size, credible_size,
                                                 #   covered, top_probability), ...]
metrics.fit_summary_metrics(analysis, truth)     # {fit-summary column name: value} for the ranking and parent sets
```

The benchmark scripts run in this environment:

```sh
uv run python ../benchmark/chain_summary.py RUN_DIR
```
