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

This creates `.venv/` with Python 3.13, DuckDB, pyarrow and pytest, and installs the `pt-analysis` command into it.

## Run the tests

```sh
uv run pytest
```

The tests read the golden stores in `test/data/golden/` and check the benchmark scripts in `tools/benchmark/`, which
read stores through this package.

## Describe an Analysis

An Analysis is a set of Runs of one input, described by a JSON manifest. Paths are relative to the manifest file.

```json
{
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

`input`, the two IDP paths and `runs` are required; every Run directory must exist. The rest is optional.

```sh
uv run pt-analysis info analysis.json
```

prints, for each Run, its closed Windows per Phase and its Draw count.

## Derive tables

```sh
uv run pt-analysis derive RUN_DIR          # one Run
uv run pt-analysis derive analysis.json    # every Run of an Analysis
```

writes Parquet tables under `RUN_DIR/derived/`. `edges` and `set_summary` come from the sampling Phase, so every
weight is a Draw's; `set_classes` comes from `parent_set_members` over both Phases.

| Table | Columns | Rows |
|---|---|---|
| `set_classes` | infection, set_id, n_observed, has_unobserved, class | one per Parent set, from `parent_set_members` |
| `edges` | iteration, child, parent, prob | per Draw and child, one per parent: summed weight of the child's sets containing it |
| `set_summary` | infection, set_id, weight, rank, cumulative_mass | one per Parent set with weight in some Draw |

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

`derived/manifest.json` records the closed Windows per Phase the tables came from and the package version that wrote
them. `derive` does nothing when both match the store and this package, and replaces every table otherwise, for
example once the Run has closed another Window. The manifest is written last, so an interrupted derive is redone.

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
functionals.top_edges(store, n=10)               # [Edge(parent, child, weight), ...]
functionals.map_graph(store)                     # {infection: RankedSet} of modal sets
functionals.psi_map(store)
functionals.size_range(store, level=0.8)         # (psi_lo, psi_hi)
functionals.credible_collection(store, "3", level=0.8)   # [RankedSet(set_id, members, weight, rank, ...), ...]
```

The benchmark scripts run in this environment:

```sh
uv run python ../benchmark/chain_summary.py RUN_DIR
```
