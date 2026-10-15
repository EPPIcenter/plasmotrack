# Run store layout

A Run writes its output to one directory, the store: Parquet tables committed one Window at a time, a Run manifest,
lookup tables, and checkpoints. This page describes what the engine writes, so that a Run can be read without reading
the C++. `test/data/golden/single-chain` and `test/data/golden/replica-exchange` are small, real stores to try it on.

The writer is `src/impl/model/Model/RunRecorder.cpp` (which rows exist and when) and
`src/core/io/ParquetStoreSink.cpp` (files and commits). Where this page and the code disagree, the code wins, and this
page needs fixing.

## Terms

- **Sampling iteration**: one pass in which every parameter gets an update opportunity.
- **Phase**: `burn_in` or `sampling`. Only sampling-Phase rows are posterior samples.
- **Trace**: a table with one row (or one row per Infection, Locus, and so on) per written Sampling iteration.
- **Block**: a table with one row per Window, holding counts accumulated over every Sampling iteration of that
  Window.
- **Window**: a contiguous range of Sampling iterations, committed to disk whole. A killed Run loses at most the
  Window in flight.
- **Hot chain**: the replica at temperature 1 (temperature slot 0). Every Trace except `chains` describes the Hot
  chain.

## Directory tree

```
<run>/
├── manifest.json                     Run manifest, written once at Run start
├── structure/
│   ├── loci.parquet                  Allele positions per Locus
│   ├── infections.parquet            each Infection's Latent parent and symptomatic flag
│   └── samplers.parquet              each sampler id's type and target
├── run_scalars/window-NNNNNN.parquet one file per table per closed Window
├── infections/window-NNNNNN.parquet
├── ordering/…
├── allele_frequencies/…
├── latent_genotypes/…
├── latent_parent_genotypes/…
├── parent_sets/…
├── parent_set_members/…
├── chains/…
├── acceptance/…
├── swaps/…
├── events/…
├── _closed/window-NNNNNN             commit marker, one per closed Window
├── checkpoint.cbor                   latest checkpoint, replaced after every Window
└── checkpoint-burn-in.cbor           checkpoint taken when burn-in ended
```

- `NNNNNN` is the Window number, zero-padded to at least six digits: `window-1000000` follows `window-999999`, so
  sort by number, not by name. Windows are numbered 0, 1, 2, … across both Phases, and a resumed Run carries on the
  numbering.
- A table with no rows in a Window has no file for that Window. A table with no rows in the whole Run has no directory.
  For example, a one-chain Run has no `swaps/`, and `parent_set_members/` has files only for Windows in which a new
  Parent set appeared.
- While a Window is being written, its files are named `.in-flight.parquet`. They are renamed into place at commit.
  `arrow` and `pyarrow` skip files whose names start with `.` or `_`.
- The `window` value is a column in every row. It is not encoded in a directory name, so the tables aren't
  Hive-partitioned.

## Reader rule

**A Window exists only if `_closed/window-NNNNNN` exists.** Read only rows whose `window` is in the set of markers.

A Window commits in this order: each table is written under a hidden name, synced, renamed into place, and its
directory synced; the marker is written last, the same way; then `checkpoint.cbor` is replaced, the same way. A Run killed mid-commit can leave table files of a Window that has
no marker. Those rows must be ignored. `--resume` deletes them.

A marker is a JSON object:

```json
{"first_iteration":21,"last_iteration":45,"phase":"sampling","window":1}
```

| Field | Meaning |
|---|---|
| `window` | Window number |
| `first_iteration`, `last_iteration` | Sampling iterations the Window covers, inclusive |
| `phase` | Phase of every row in the Window |

Every Window covers at least one Sampling iteration, so `first_iteration` ≤ `last_iteration`.

## Iterations, thinning and Windows

- `iteration` counts Sampling iterations from the start of the Run, across both Phases. Iteration 0 is the starting
  State. Burn-in is iterations 1 to `burn_in`. Sampling is iterations `burn_in + 1` to `burn_in + samples`. A Run
  started with `--start-from` numbers its iterations from 0 again.
- Every Trace has a row at iteration 0 (Phase `burn_in`), the starting State. Blocks have no row for it.
- Trace rows are written every `thin` (`-t`) Sampling iterations, counted from the start of each Phase: burn-in rows
  at iterations `t, 2t, …`, and sampling rows at `burn_in + t, burn_in + 2t, …`. The same grid applies to every Trace,
  so tables join on `iteration`.
- A Window holds `window` (`--window`, default 100) Sampling iterations, also counted from the start of each Phase.
  `--window` must be a multiple of `-t`, so a full Window holds `window / thin` rows per Trace. Burn-in end, an
  interrupt and the end of the Run end the in-flight Window early. Every Window belongs to exactly one Phase.
- A Window that is full, or that burn-in end is recorded in, stays open until the next Sampling iteration is recorded,
  and is committed then. Run events raised after its last Sampling iteration (a temperature adaptation, burn-in end,
  an interrupt, the end of the Run) go into it, so no Window holds Run events alone.
- An interrupt writes no extra Trace row.

## Tables

All ids are strings: Infection ids are the input's node `id`, Latent parent ids are `<infection id>_copy`, and Locus
ids are the input's `locus` labels. Probabilities and parameter values are `double`. Counts and indices are `int64`.
`prior` and `likelihood` values are natural logs.

Every Trace begins with these three columns:

| Column | Type | Meaning |
|---|---|---|
| `iteration` | int64 | Sampling iteration |
| `phase` | string | `burn_in` or `sampling` |
| `window` | int64 | Window the row was committed in |

### `run_scalars` (Trace, one row per written iteration)

| Column | Type | Meaning |
|---|---|---|
| `mean_coi` | double | Mean of the zero-truncated Poisson COI distribution |
| `mean_strains_tx` | double | Mean number of strains transmitted |
| `parent_set_size_prob` | double | Parameter of the zero-truncated geometric distribution of Parent set size |
| `prior` | double | Log prior of the Hot chain |
| `likelihood` | double | Log likelihood of the Hot chain. The log posterior is `prior + likelihood` |
| `observation` | double | Log likelihood of the Observation process, summed over Infections |
| `transmission` | double | Log likelihood of the Transmission process, summed over Infections |
| `source` | double | Always 0. Kept for the table's layout |

`likelihood = observation + transmission`. Latent parents are summed out of the Transmission process (each Parent set
with a Latent parent scores its expectation over the Source), so there is no separate Source term.

### `infections` (Trace, one row per Infection per written iteration)

| Column | Type | Meaning |
|---|---|---|
| `infection` | string | Infection id |
| `expected_false_positives` | double | Expected number of false-positive Allele calls over the Infection's Observed genotypes. The per-Allele rate is this divided by the Infection's total Alleles |
| `expected_false_negatives` | double | Expected number of false-negative Allele calls, as above |
| `infection_duration` | double | Time from infection to observation, in the time unit of the IDP files (the IDP is indexed by the rounded duration) |

### `ordering` (Trace, one row per Infection per written iteration)

| Column | Type | Meaning |
|---|---|---|
| `infection` | string | Infection id |
| `rank` | int64 | Position in the Ordering, 0 = acquired first |

### `allele_frequencies` (Trace, one row per Locus × Allele per written iteration)

| Column | Type | Meaning |
|---|---|---|
| `locus` | string | Locus id |
| `allele` | int64 | Allele index, 0-based, as in `structure/loci` |
| `frequency` | double | Allele frequency in the Source population. Sums to 1 over a Locus |

### `latent_genotypes` and `latent_parent_genotypes` (Trace, one row per Infection or Latent parent × Locus per written iteration)

| Column | Type | Meaning |
|---|---|---|
| `infection` / `latent_parent` | string | Infection id, or Latent parent id |
| `locus` | string | Locus id |
| `alleles` | string | Latent genotype as an Allele bitstring, such as `00101000`. Character `k` (0-based, from the left) is Allele `k`, and `1` means present. Same convention as the input's `genotype` strings |

Only Loci the Infection was genotyped at are written. `latent_parent_genotypes` has no rows: Latent parents are summed
out of the model and have no genotypes to record. A reader should treat the table as possibly absent.

### `parent_sets` (Trace, one row per Infection × Parent set per written iteration)

| Column | Type | Meaning |
|---|---|---|
| `infection` | string | The Infection (child) |
| `set_id` | int64 | Parent set id, see `parent_set_members` |
| `prob` | double | Probability of this Parent set given the current State |

- Probabilities are normalised over every Parent set, and then rows with `prob` below the Run's `parent_set_floor`
  (manifest, default 1e-9) are dropped without renormalising. So the rows of one Infection and iteration may sum to
  slightly less than 1.
- Rows whose `prob` is exactly 0 (the Parent set is impossible, or its probability underflowed) are dropped whatever
  the floor. A floor of 0 keeps every non-zero row.
- There is no separate row for "Source only". The Parent set whose only member is the Infection's Latent parent is the
  probability that the Infection has no sampled parent.
- The marginal probability that Infection `j` infected Infection `i` is the sum of `prob` over `i`'s Parent sets that
  contain `j`.
- `set_id` is written without dictionary encoding, and this table's data pages are up to 8 MiB, so a page holds about
  2^20 values of an 8-byte column. Arrow 22's C++ writer cuts pages by bytes and has no rows-per-page setting. Other
  tables keep the writer's defaults. Every table is compressed with zstd level 3.

### `parent_set_members` (lookup table, one row per Parent set member)

| Column | Type | Meaning |
|---|---|---|
| `infection` | string | The Infection the Parent set belongs to |
| `set_id` | int64 | Parent set id |
| `member` | string | An Infection id or the Infection's Latent parent id |

- `set_id`s are numbered per Infection from 0, in order of first appearance, and are stable for the whole Run,
  including across `--resume`.
- A Parent set's members are written once, in the Window in which the set first appears, which may be a burn-in
  Window. Read this table over every closed Window, in both Phases, before joining it to sampling-Phase `parent_sets`.
- This table has no `iteration` or `phase` column.

### `chains` (Trace, one row per temperature slot per written iteration)

| Column | Type | Meaning |
|---|---|---|
| `temperature_slot` | int64 | 0 = Hot chain; temperature falls as the slot number rises |
| `replica` | int64 | Which replica occupies the slot. Replicas move between slots when they swap |
| `temperature` | double | Temperature of the slot |
| `prior` | double | Log prior of that replica |
| `likelihood` | double | Log likelihood of that replica, untempered. The slot's target is `temperature × likelihood + prior` |

### `acceptance` (Block, one row per sampler per Window)

| Column | Type | Meaning |
|---|---|---|
| `window` | int64 | Window |
| `phase` | string | Phase of the Window |
| `sampler` | string | Sampler id, `<type>/<target>`, as in `structure/samplers` |
| `accepted` | int64 | Proposals accepted in the Window |
| `proposed` | int64 | Proposals made in the Window |

- Blocks count every Sampling iteration of the Window, not only the thinned ones.
- Only the Hot chain is counted. At each Sampling iteration the counts come from whichever replica occupied
  temperature slot 0 during it.
- Summed over Windows, the rows give Run totals for each Phase.

### `swaps` (Block, one row per neighbouring pair of temperature slots per Window)

| Column | Type | Meaning |
|---|---|---|
| `window` | int64 | Window |
| `phase` | string | Phase of the Window |
| `temperature_pair` | int64 | `k` = swaps between slots `k` and `k + 1` |
| `accepted` | int64 | Swaps accepted in the Window |
| `attempted` | int64 | Swaps attempted in the Window |

Only Runs with more than one chain have this table.

### `events` (Run events, one row per event)

| Column | Type | Meaning |
|---|---|---|
| `iteration` | int64 | The last Sampling iteration recorded when the event happened |
| `phase` | string | Phase at that moment |
| `window` | int64 | Window the event was committed with |
| `wall_time` | string | UTC time, ISO 8601 to the millisecond: `2026-10-11T12:00:00.000Z` |
| `type` | string | Event type, below |
| `payload` | string | JSON object, below |

| `type` | When | `payload` |
|---|---|---|
| `run_started` | Start of a new Run, including one started with `--start-from` | `started_at` (UTC string), `burn_in`, `samples` (Sampling iterations planned) |
| `run_resumed` | Start of a `--resume` | `iteration` (last Sampling iteration the checkpoint covers), `window` (last closed Window), `checkpoint` (path), `samples` (sampling iterations now planned, which `-s` may have raised) |
| `temperatures_adapted` | Every 5th burn-in iteration after the 10th, when there is more than one chain | `applied` (false if the new ladder had NaNs and the old one was kept), `temperatures` (ladder by slot), `barriers` (swap barrier per neighbouring pair) |
| `burn_in_ended` | Burn-in is over and adaptation stops | `proposal_variances` (object: sampler id → array of final proposal variances, one per proposal scale; SALT samplers have one per Allele; samplers without a scale are left out), `temperatures` (final ladder) |
| `window_closed` | Each Window commit | `sampling_iterations` (Sampling iterations in the Window, the starting State not counted), `seconds` (wall time from the previous Window's last Sampling iteration to this Window's last, or to the end of the Run for the Window the Run's end commits), `iterations_per_second` (null if `seconds` is 0) |
| `non_finite_value` | A temperature slot's prior or likelihood becomes NaN or ±inf, checked every Sampling iteration, thinned or not. Raised again only after it has been finite | `temperature_slot`, `term` (`prior` or `likelihood`), `value` (`"nan"`, `"inf"` or `"-inf"`) |
| `interrupted` | First and second interrupt signal | `signal` (such as `"SIGINT"`), `second` (bool) |
| `run_ended` | End of the Run | `reason` (`completed` or `interrupted`), `seconds` (wall time since this process started or resumed the Run) |

A hard kill writes no event. The Run's events stop at its last closed Window.

## Lookup tables (`structure/`)

Written once at Run start. They are not part of any Window and have no marker.

| Table | Columns | Meaning |
|---|---|---|
| `loci` | `locus` string, `allele` int64 | One row per Allele of each Locus, in bitstring order: `allele` k is character k of an `alleles` bitstring |
| `infections` | `infection` string, `latent_parent` string, `symptomatic` bool | Each Infection's Latent parent id, and whether it was symptomatic (which IDP it uses) |
| `samplers` | `sampler` string, `type` string, `target` string | Each sampler id split at its first `/`. For example `SALTSampler/allele_frequencies/AS1` has type `SALTSampler` and target `allele_frequencies/AS1` |

## Run manifest (`manifest.json`)

Written once at Run start. A resumed Run keeps its manifest, and its `run_resumed` event records the new `-s`.

| Field | Meaning |
|---|---|
| `run_id` | 128-bit random hex id of the Run |
| `started_at` | UTC wall time the Run started |
| `seed` | Seed used. If `--seed` was -1, the seed that was drawn |
| `cli_arguments` | Every command-line argument, `argv[0]` included |
| `engine_version` | Engine version, as `--version` prints it |
| `git_commit` | Commit the engine was built from |
| `input_sha256` | Object: `input`, `symptomatic_idp`, `asymptomatic_idp` → SHA-256 of each file |
| `chains` | Number of chains (`-n`) |
| `gradient` | `-g`. Absent for a Run started with `--start-from`, whose temperature ladder comes from the checkpoint; see `started_from` |
| `burn_in` | `-b` |
| `samples` | `-s`, as planned at Run start |
| `null_model` | Whether the Run uses the null model (`--null-model`, which may be abbreviated). `--resume` reads it from here |
| `thin` | `-t` |
| `window` | `--window` |
| `parent_set_floor` | `--parent-set-floor` |
| `started_from` | Only for `--start-from`: `run_id` and `engine_version` of the Run the checkpoint came from, and `temperatures`, the checkpoint's ladder by temperature slot, which the new Run starts with |

## Checkpoints and the two CLI options

`checkpoint.cbor` is written after every Window's marker. It holds the State at the Window's last Sampling iteration,
with any temperature adaptation after it, even though a full Window is committed only once the next Sampling iteration
has run. It is written under a temporary name, synced, and renamed over the old one, so a crash leaves either the old
checkpoint or the new one. `checkpoint-burn-in.cbor` is the same checkpoint for the Window burn-in end is recorded in,
written when that Window commits. A second interrupt closes the in-flight Window but writes no checkpoint, since the
State may be mid-iteration.

Both files are CBOR (`format_version` 5). An engine refuses a checkpoint of another format version. Python `cbor2` or
R `cbor` can read them for debugging:

| Key | Contents |
|---|---|
| `format_version` | 5 (4 did not save `burn_in_iterations`; 3 saved the scalars and per-Infection values under fixed keys; 2 saved Latent parent genotypes, which are now summed out of the model) |
| `identity` | `run_id`, `input_hash` (SHA-256 of the input), `engine_version` |
| `run_position` | `iteration`, `phase` (Phase of the next Sampling iteration), `last_closed_window`, `burn_in_remaining`, `sampling_remaining` |
| `parent_set_ids` | Per Infection: `sets` (array of `members`, `set_id`) and `next_id` |
| `run_rng` | State of the top-level random engine, which draws swaps |
| `replicas` | Per replica: `state` (`scalars` by Quantity name: `mean_coi`, `mean_strains_tx`, `parent_set_size_prob`; `per_infection` by Quantity name, each an array by Infection: `expected_false_positives`, `expected_false_negatives`, `infection_duration`; `allele_frequencies` by Locus; `latent_genotypes` packed bits; `ordering`), `rng`, `samplers` (by id: `counts`, `adaptation`), `scheduler_iterations` |
| `exchange` | Temperature slots: `replicas` (replica in each slot), `temperatures`, `swap_barriers`, `barrier_iterations`, `burn_in_iterations` (burn-in iterations run, which decide when the ladder next adapts), `second_pair_set_next`, `swap_counts` |

Model caches are not saved; they are rebuilt from the State.

**`--resume <run>`** continues the same Run in its own store. Give the same input and IDPs and no output directory.

- Before continuing, it deletes every Window newer than `checkpoint.cbor`, markers first, and any in-flight file.
- It refuses a different input, IDP, engine version, `-n`, `-t`, `--window`, `-b`, `-g`, `--parent-set-floor`,
  `--seed` or `--null-model`. Options left out take the Run's values.
- `-s` may grow, which extends a finished Run. It may not shrink.
- Iteration and Window numbering continue, no starting-State row is written, and a `run_resumed` event is recorded.
- A resumed Run that ends before it records a Sampling iteration (for example, interrupted at once) commits nothing:
  its `run_resumed` event is not kept, and the store stays as its checkpoint left it.
- Resume is statistically valid, not bitwise identical to an uninterrupted Run.

**`--start-from <checkpoint>`** starts a new Run, with its own store, manifest and seed, from a checkpoint's State and
adaptation state. The checkpoint may be `checkpoint-burn-in.cbor` or `checkpoint.cbor` of any Run.

- It needs the same input and the same chain count. `-n` defaults to the checkpoint's.
- The engine version may differ; it is recorded in `started_from`.
- The temperature ladder is the checkpoint's, recorded in `started_from`; `-g` is ignored, and the manifest has no
  `gradient`. Resuming such a Run refuses `-g`.
- It runs its own `-b` burn-in from the saved State (`-b 0` is allowed) and numbers iterations from 0. Iteration 0 is
  the saved State. The ladder adapts on the iterations a fresh Run's does, from a clean slate: swap barriers the
  source Run had pending are not carried over.

## Reading a store

### R (`arrow`)

```r
library(arrow)
library(dplyr)

run <- "path/to/run"
closed <- vapply(list.files(file.path(run, "_closed"), pattern = "^window-", full.names = TRUE),
                 function(f) jsonlite::fromJSON(f)$window, numeric(1))

# The files of one table in closed Windows: the reader rule.
closed_files <- function(table) {
  files <- file.path(run, table, sprintf("window-%06d.parquet", closed))
  files[file.exists(files)]
}

# Hot-chain log posterior over the sampling Phase.
scalars <- open_dataset(closed_files("run_scalars")) |>
  filter(phase == "sampling") |>
  mutate(log_posterior = prior + likelihood) |>
  arrange(iteration) |>
  collect()

# Marginal probability of each parent of each Infection, per Sampling iteration.
members <- open_dataset(closed_files("parent_set_members")) |> collect()
edges <- open_dataset(closed_files("parent_sets")) |>
  filter(phase == "sampling") |>
  collect() |>
  inner_join(members, by = c("infection", "set_id"), relationship = "many-to-many") |>
  group_by(iteration, infection, member) |>
  summarise(prob = sum(prob), .groups = "drop")
```

For a table with a `window` column, `open_dataset(file.path(run, "run_scalars")) |> filter(window %in% closed)` is
equivalent. `parent_set_members` has no `window` column, so select its files by name.

### Python (`pyarrow`)

```python
import json
from pathlib import Path

import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

run = Path("path/to/run")
closed = sorted(json.loads(p.read_text())["window"] for p in (run / "_closed").glob("window-*"))

def closed_files(table):
    return [str(f) for w in closed if (f := run / table / f"window-{w:06d}.parquet").exists()]

# Hot-chain log posterior over the sampling Phase.
scalars = ds.dataset(closed_files("run_scalars")).to_table(filter=ds.field("phase") == "sampling")
scalars = scalars.sort_by("iteration")
log_posterior = pc.add(scalars["prior"], scalars["likelihood"])

# Marginal probability of each parent of each Infection, per Sampling iteration.
members = ds.dataset(closed_files("parent_set_members")).to_table()
sets = ds.dataset(closed_files("parent_sets")).to_table(filter=ds.field("phase") == "sampling")
edges = (sets.join(members, ["infection", "set_id"])
             .group_by(["iteration", "infection", "member"])
             .aggregate([("prob", "sum")]))
```

Selecting files by closed Window, as above, applies the reader rule to every table, including `parent_set_members`.
`tools/benchmark/mixing_benchmark.py` reads stores the same way.
