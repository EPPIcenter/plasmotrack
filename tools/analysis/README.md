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
```

The benchmark scripts run in this environment:

```sh
uv run python ../benchmark/chain_summary.py RUN_DIR
```
