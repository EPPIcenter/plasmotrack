# Which R engine reads a Run store lazily within a memory cap?

Research for [Which R engine reads a Run store lazily within a memory cap?](https://github.com/EPPIcenter/plasmotrack/issues/99),
part of [Find the way to an R analysis package spec](https://github.com/EPPIcenter/plasmotrack/issues/98).

Domain terms follow `CONTEXT.md`: Run, Draw, Window, Window marker, Store, Parent set, Trace. A Draw is one written
sampling-Phase Sampling iteration.

## Answer

Use **DuckDB through `duckdb` and `dbplyr`**, on a connection the package owns, with `memory_limit` and
`temp_directory` set. Don't depend on `arrow` and don't use `duckplyr` as the user-facing interface.

- DuckDB is the only one of the three that caps memory and spills to disk. On a synthetic store with 234M sampling
  `parent_sets` rows, per-Draw edge marginals ran in 3.1 GB RSS under a 2 GB `memory_limit`, spilling 5.5 GB. The
  `arrow` (Acero) version of the same query passed 10 GB and was killed.
- `dbplyr` translated every verb the derivations need, lazily: joins, grouped summaries, exact median and quantiles,
  grouped window functions (`row_number`, `lag`, `cumsum`), grouped `filter`, `slice_max`, `pivot_longer`,
  `pivot_wider`, `complete` and `cross_join`. `arrow` errored on every window function and every tidyr verb.
- `duckdb` has one hard dependency (`DBI`), a 56 MB install, CRAN binaries for macOS and Windows, Posit Package
  Manager binaries for RHEL 9, and OK on all 13 CRAN check flavors. `arrow` is 108 MB, has 14 dependencies, and its
  Linux source install downloads `libarrow` or needs `cmake` and a C++20 compiler.
- Neither engine applies the reader rule by itself. With a glob, DuckDB reads hidden `.in-flight.parquet` files.
  Opening a directory with `arrow` skips those but reads unmarked Windows. The package should build each table's file
  list from the `_closed/` markers and never glob a table directory. Both engines read that list correctly.

Two findings matter beyond the engine choice:

- **Query shape matters as much as the engine.** Pooled edge marginals computed directly, as one sum per (child,
  parent) with no per-Draw grouping, ran in 0.6 GB (DuckDB) and 0.8 GB (`arrow`) on the 234M-row store. Only per-Draw
  tables, which R-hat and ESS need, require spilling.
- **The 25 GB `derive` peak on Eswatini came from `genotype_summary`, not from edges.** That query joins every
  `latent_genotypes` row to every Allele of its Locus before counting. Reached 12 GB unbounded (killed), and spilled
  24 GB over 61 s under a 2 GB cap. Counting distinct bitstrings first and expanding afterwards gives an identical table
  in 0.2 s and about 2 GB. Edges for the same Run need 0.3 GB.

## Method

All experiments ran on a 48 GB Apple M-series laptop (14 threads) under R 4.5.2 with duckdb 1.5.6 (DuckDB v1.5.6),
dbplyr 2.6.0, dplyr 1.2.1, tidyr 1.3.2, arrow 25.0.1 and duckplyr 1.2.2, all CRAN macOS arm64 binaries. Each query
ran in a fresh `Rscript` process under a watchdog that sampled RSS every 0.1 s, killed the process above 10 GB, and
sampled the size of the DuckDB spill directory. Peak RSS includes the R process itself, about 0.2 GB with the
packages loaded. RSS and spill figures are GiB. DuckDB reads `memory_limit` values as decimal units: it reports
`6GB` as 5.5 GiB, `2GB` as 1.86 GiB and `1GB` as 0.93 GiB. Scratch scripts lived in `/tmp/rengine/`, not in any store, and are not committed.

Stores:

| Store | Sampling `parent_sets` rows | Draws | Notes |
|---|---|---|---|
| `test/data/golden/single-chain` | 90 | 10 | 3 Windows |
| `test/data/golden/replica-exchange` | 90 | 10 | has `chains`, `swaps` |
| Eswatini `runs/chain1` (`~/Workspace/plasmotrack-runs/eswatini-prod`) | 1,460,868 | 500 | 457 Infections, 15 Windows, 54 MB; the largest Parquet store on disk |
| Eswatini x40 (synthetic) | 58,434,720 | 20,000 | chain1's sampling Windows repeated 40 times at new iterations and Window numbers, with markers; 444 Windows, 170 MB |
| Eswatini x160 (synthetic) | 233,738,880 | 80,000 | as above, 160 times; 671 MB |

He01 (about 7k Parent-set rows per Infection per iteration) exists only in the old gzip CSV format, so the synthetic
stores stand in for He01-scale volume. They repeat chain1's 1,918 (child, parent) pairs, so they stress row count
and Draw count, not the number of distinct edges.

The heavy query is the one the ticket names: sampling-Phase `parent_sets` joined to `parent_set_members` (read over
both Phases), summed per (iteration, child, parent) into per-Draw edge marginals, then pooled over Draws as
`sum / number of Draws`, zero-filled as in `pt_analysis`. "perdraw" writes the per-Draw table to Parquet; "pooled"
collects the pooled table. On every store, every engine returned the same pooled checksum, and on chain1 it matched
`pt_analysis`'s own `derived/edges.parquet` exactly (sum of pooled means 653.198386 over 1,918 pairs).

## Reader rule

Test store: a copy of `single-chain` (closed Windows 0, 1, 2) plus, in `parent_sets`, `run_scalars` and
`parent_set_members`, an unmarked `window-000003.parquet` (a Run killed mid-commit) and a `.in-flight.parquet` (Window
4). The `parent_set_members` extras carry marker members `GHOST` and `INFLIGHT`.

| How the table is opened | `parent_sets` Windows read | Bad `parent_set_members` rows | Correct? |
|---|---|---|---|
| DuckDB `read_parquet('<table>/*.parquet')` | 0 1 2 3 **4** | GHOST, **INFLIGHT** | no |
| DuckDB `read_parquet([marker file list])` | 0 1 2 | none | yes |
| `arrow::open_dataset('<table>/')` | 0 1 2 3 | GHOST | no |
| `open_dataset('<table>/') \|> filter(window %in% closed)` | 0 1 2 | n/a: no `window` column | only for tables with `window` |
| `arrow::open_dataset(marker file list)` | 0 1 2 | none | yes |

- **DuckDB globs match hidden files.** `*.parquet` matched `.in-flight.parquet`. DuckDB's Parquet docs describe globs
  and file lists but say nothing about hidden files ([DuckDB 1.5 Parquet overview](https://duckdb.org/docs/1.5/data/parquet/overview.html)).
  `doc/store-layout.md` says "`arrow` and `pyarrow` skip files whose names start with `.` or `_`". That holds, but a
  DuckDB reader can't rely on it.
- **`arrow` skips the in-flight file when it discovers a directory**, but `selector_ignore_prefixes` is "Not valid
  when providing a vector of file paths" ([arrow `open_dataset`](https://arrow.apache.org/docs/r/reference/open_dataset.html)),
  and directory discovery still reads unmarked Windows.
- **The marker-derived file list is correct in both engines and is the only approach that covers
  `parent_set_members`**, which has no `window` column. It's what `pt_analysis.store.Store.table` does
  (`tools/analysis/src/pt_analysis/store.py` on `analysis-gui-react`). `read_parquet` takes a list of paths
  ([DuckDB 1.5 Parquet overview](https://duckdb.org/docs/1.5/data/parquet/overview.html)). A 444-file list caused no
  trouble.
- Filtering to sampling Windows in the file list, not with `filter(phase == "sampling")`, means burn-in files are never
  opened.
- **`window` is a reserved word in DuckDB SQL.** `SELECT ... AS window` is a parser error. `dbplyr` quotes identifiers,
  so lazy-table code is fine, but any hand-written SQL in the package must write `"window"`.

## Memory cap and spill

DuckDB's `memory_limit` defaults to 80% of RAM and `temp_directory` to `.tmp` in in-memory mode;
`max_temp_directory_size` defaults to 90% of free disk ([DuckDB 1.5 configuration](https://duckdb.org/docs/1.5/configuration/overview.html)).
DuckDB "supports larger-than-memory processing" for grouping, joins, sorting and windowing. But it "may still throw an
out-of-memory exception" when several blocking operators share a query, and `list()`, `string_agg()` and `PIVOT`
can't spill ([DuckDB 1.5: tuning workloads](https://duckdb.org/docs/1.5/guides/performance/how_to_tune_workloads.html)).
The R client docs add that the limit "is only applied to the memory DuckDB uses"
([DuckDB 1.5 R client](https://duckdb.org/docs/1.5/clients/r.html)).

`arrow`'s query engine, Acero, has no memory limit or spill setting. Its docs describe neither
([arrow `acero`](https://arrow.apache.org/docs/r/reference/acero.html),
[Acero overview, Arrow 26](https://arrow.apache.org/docs/cpp/acero/overview.html)). Arrow's own tracker says the hash
join holds its hash table and inputs in memory, with spilling still unresolved
([ARROW-16389](https://issues.apache.org/jira/browse/ARROW-16389)).

Measured peak RSS (GB) and peak spill (GB). "killed" means the watchdog stopped it at 10 GB.

| Store | Query | DuckDB no limit | DuckDB 6 GB | DuckDB 2 GB | DuckDB 1 GB | DuckDB 500 MB | `arrow` | duckplyr 1 GB |
|---|---|---|---|---|---|---|---|---|
| golden, both | pooled | | 0.24 | | | | 0.22 | |
| Eswatini chain1 | pooled | | 0.24 | | | | 0.54 | |
| x40 | pooled | 2.60 | 2.59 | 2.27 | 1.53 | 1.12 (3.6 s) | 5.86 | 1.57 |
| x40 | perdraw | 2.92 | 2.90 | 2.50 | 1.72 | 1.47 (4.8 s) | 5.91 | 1.68 |
| x160 | perdraw | | 6.13 / spill 3.2 | 3.12 / spill 5.5 | OOM error | | killed | |
| x160 | pooled | | 6.31 / spill 3.3 | 3.10 / spill 5.6 | OOM error | | killed | |
| x160 | pooled directly | | | | 0.62 | 0.60 | 0.79 | |

- **DuckDB holds its cap and spills.** On x160, the per-Draw query (102.6M output rows) ran in 13 s at both 6 GB and
  2 GB. Lowering the cap moved memory to disk.
- **RSS runs above `memory_limit`:** on the spilling queries, about 0.6–0.8 GiB over at `6GB` and 1.2–1.3 GiB over at
  `2GB`. On x40 at `1GB` it was 0.6–0.8 GiB over. That covers R itself, Parquet reader buffers and anything else the
  buffer manager doesn't track. To hold a process under N GB, set `memory_limit` to roughly N − 1.5 GB.
- **Too small a cap fails cleanly.** At 1 GB, with 14 or 4 threads, the x160 queries raised `Out of Memory Error: could
  not allocate block of size 256.0 KiB` as an R error, after spilling 5–7 GB. The R session survived. The package
  should catch `duckdb_error` of type `OUT_OF_MEMORY` and tell the user how to raise the limit.
- **`arrow` memory grows with the data.** 5.9 GB for a 170 MB store, then past 10 GB for x160.
- **Pooling directly needs little memory in either engine.** A sum over Draws of per-Draw sums is one sum per (child,
  parent). Grouping by (child, parent) alone, with no per-Draw intermediate, keeps the aggregate at 1,918 groups.
  Both engines then stream the 234M rows.
- **`arrow` → `to_duckdb()` is worse than DuckDB alone.** The same x40 pooled query failed with OOM under a 1 GB limit,
  where DuckDB reading the Parquet itself succeeded at 500 MB.

### The Eswatini `derive` peak

`eswatini-prod/derive.log` records 25–26 GB max RSS for `pt-analysis derive` on each of chain2–4 and 25 GB for the
Analysis, with no `memory_limit` set (`derive.sh`), so DuckDB's default of 80% of 48 GB applied. The same SQL as
`derive.py` on `analysis-gui-react`, run one table at a time on chain1 from R:

| Derived table | No limit | 2 GB limit |
|---|---|---|
| `edges` | 0.29 GB | 0.26 GB |
| `coi` | 0.17 GB | 0.18 GB |
| `genotype_mode` | 1.82 GB | 2.27 GB, spill 2.0 GB |
| `order_pairs` | 3.95 GB | 2.38 GB, spill 0.8 GB |
| `genotype_summary` | **killed at about 12 GB** | 2.48 GB, **spill 24.1 GB, 61 s** |
| `genotype_summary`, counted first | 2.07 GB, 0.2 s | 2.28 GB, spill 2.0 GB, 0.7 s |

`genotype_summary` joins each of chain1's sampling `latent_genotypes` rows (Infection × Locus × Draw) to every Allele
of the Locus in `structure/loci` before it counts. The rewrite counts each (Infection, Locus, bitstring) first, as
`genotype_mode` already does, then joins the counts to the Alleles and sums. It produced an identical table: 629,746
rows, summed presence count 49,311,083. `plasmotrackr` should define `genotype_summary` the counted way, and
`pt_analysis` can take the same fix.

## dplyr translation

Each operation ran on lazy tables over the golden single-chain store and was then collected. duckplyr ran with
`prudence = "stingy"`, so an operation it can't run in DuckDB errors instead of silently materializing.

| Operation | DuckDB + dbplyr | duckplyr | `arrow` |
|---|---|---|---|
| `inner_join`, `left_join`, `anti_join`, many-to-many join | lazy | lazy | lazy |
| `group_by() \|> summarise(mean, sum, n())` | lazy | error: `group_by` unsupported | lazy |
| `n_distinct` | lazy (needs dbplyr ≥ 2.6.0) | error | lazy |
| `median`, `quantile` | lazy, exact | error | lazy, approximate (t-digest, warns) |
| `summarise(.by =)` | lazy | lazy | lazy |
| grouped `mutate(row_number())`, `lag`, `cumsum` | lazy (`OVER (PARTITION BY ... ORDER BY ...)`) | error, except ungrouped `lag` | **error: "Expression not supported in Arrow"** |
| grouped `filter(prob == max(prob))` | lazy | error | **error** |
| grouped `slice_max` | lazy | error | **error: "Slicing grouped data not supported"** |
| `pivot_longer`, `pivot_wider` | lazy | error | **error: no method** |
| `complete(..., fill =)` | lazy | error | **error: no method** |
| `cross_join` | lazy | error | **error: no method** |
| `case_when`, `if_else`, `substr` | lazy | error | lazy |
| COI via `nchar(gsub("0", "", alleles))` | **error: no `gsub` in DuckDB** | error | lazy |
| COI via `nchar(replace(alleles, "0", ""))` | lazy (unknown functions pass through as SQL) | | |
| `stringr::str_count` | **error: not available in this SQL variant** | error | lazy |
| `union_all`, `distinct`, `arrange`, `count` | lazy | lazy except `count` | lazy |

- `arrow`'s docs list "window functions not currently supported" under `summarise()`
  ([arrow `acero`](https://arrow.apache.org/docs/r/reference/acero.html)). Every grouped window the derivations use
  failed: `set_summary`'s rank and cumulative mass, `order_pairs`'s `lag`, and the per-Draw ranks behind diagnostics.
- `duckplyr` doesn't run `group_by()` in DuckDB, only `.by`. With its default prudence it falls back to dplyr in R
  memory, which defeats a memory cap. Ordinary tidyverse users write `group_by()`, so `duckplyr` is a poor fit for
  lazy tables handed to them.
- `duckdb` 1.5.6 warns at load that its `dbplyr` backend "needs dbplyr 2.6.0 or later", and without it translations
  such as `n_distinct()` fail. Declare `dbplyr (>= 2.6.0)`.
- Bitstring helpers (COI, Allele presence) should be written with DuckDB functions that `dbplyr` passes through
  (`replace`, `substr`, `length`), or wrapped in the package's own functions, not left to users' `gsub`/`stringr`.

## Many Runs

Each Run's lazy edge table carries a `run` column and its own Draw count, and the Runs are combined with
`union_all()` into one query that pools per Run.

| Runs pooled | Rows | DuckDB 2 GB | `arrow` |
|---|---|---|---|
| 4 real Eswatini Runs (chain1–4) | 5.8M | 0.27 GB, 0.7 s | 0.93 GB, 0.5 s |
| 4 Runs of x40 (234M rows) | 234M | 3.26 GB, spill 2.6 GB, 8 s | **killed** |
| 8 Runs of x40 | 467M | 3.16 GB, spill 6.8 GB, 18 s | not run |
| 4 Runs of x40, per-Draw to Parquet | 234M | 3.27 GB, spill 2.5 GB, 20 s | not run |

Pooling Runs in one DuckDB query stays inside the cap, and memory doesn't grow with the number of Runs. The SQL is a
`UNION ALL` of per-Run subqueries. The 4-Run real Eswatini query read 60 files with no problems.

## Install and platforms

| | `duckdb` 1.5.6 | `arrow` 25.0.1 |
|---|---|---|
| Hard dependencies (recursive, non-base) | 1 (`DBI`) | 14 |
| Installed size, macOS arm64 | 56 MB | 108 MB |
| CRAN checks | OK on all 13 flavors ([results](https://cran.r-project.org/web/checks/check_results_duckdb.html)) | OK on all 13 flavors ([results](https://cran.r-project.org/web/checks/check_results_arrow.html)) |
| CRAN Linux check install time (source build) | 2,700–4,926 s | 2,469–4,271 s |
| macOS, Windows | CRAN binaries | CRAN binaries, `libarrow` bundled ([arrow install](https://arrow.apache.org/docs/r/articles/install.html)) |
| Linux binaries | Posit Package Manager (listed for RHEL 9), conda-forge `r-duckdb` 1.5.5, r-universe | Posit Package Manager (listed for RHEL 9), conda-forge `r-arrow` 25.0.0 |
| Linux source install | Self-contained C++, no system libraries; `SystemRequirements: xz` | Downloads prebuilt `libarrow` (x86_64, allow-listed distros), else full build; needs C++20 (gcc ≥ 10) and `cmake`; offline installs build from source ([arrow install](https://arrow.apache.org/docs/r/articles/install.html)) |

- `dbplyr` 2.6.0 has 22 recursive dependencies, almost all of them already pulled in by `dplyr` and `tidyr`.
- **corehpc has no R.** No `R` on the path and no R module. It runs Rocky Linux 9.8 on x86_64 with gcc 11.5, and the
  login node reaches CRAN and Posit Package Manager. Running `plasmotrackr` there first needs an R install (conda-forge
  `r-base` with `r-duckdb`, or R plus Posit Package Manager's RHEL 9 binaries). A source build of `duckdb` takes about
  an hour on CRAN's Linux check machines.
- Windows was not tested on a machine. The CRAN check results are the evidence.
- DuckDB prints a note on every connection about where it keeps extensions. `duckdb(shared_home = FALSE)` silences
  it. The package needs no extensions: Parquet support is built in.

## Recommendation for the spec

1. **Engine:** `Imports: DBI, duckdb, dbplyr (>= 2.6.0), dplyr`. No `arrow`, no `duckplyr`, no Python.
2. **Connection:** one package-owned in-memory DuckDB connection, created lazily, so users' DuckDB connections are
   untouched. Set `memory_limit` from an option, such as `plasmotrackr.memory_limit`. The default should be well
   below DuckDB's 80% of RAM (for example 4 GB, or half of RAM if that's smaller), since RSS runs about 1–1.5 GB over
   the limit. Set `temp_directory` from an option defaulting under `tempdir()`, set `preserve_insertion_order = false`,
   and expose `threads`.
3. **Reader rule:** each raw table is `tbl(con, sql("SELECT * FROM read_parquet([...])"))` over the file list built
   from `_closed/` markers. Traces and Blocks default to sampling Windows, and `parent_set_members` reads both Phases.
   Never glob a table directory. Quote `"window"` in any literal SQL. A Run with no file for a table yields an empty
   lazy table with the right columns, as `derive.table` does.
4. **Derivations:** write pooled quantities so they aggregate straight to their final groups, and build per-Draw
   tables only where a diagnostic needs the per-Draw series. Define `genotype_summary` counted first.
5. **Errors:** turn DuckDB's out-of-memory error into a message naming the option to raise and the spill directory.

## What this changes on the map

- **Reusing `derived/`** (Not yet specified): memory doesn't force reuse. Under a cap, R computes edges for 234M
  Parent-set rows in about 3 GB, and the real Eswatini Runs need under 0.3 GB for edges. The 25 GB was
  `genotype_summary`'s join-then-count plus the uncapped default. Reuse or caching can be decided on speed and on
  agreement with `pt_analysis`. The `genotype_summary` rewrite applies to `pt_analysis` too.
- **Dependencies and install** (Not yet specified): the engine is `duckdb` + `dbplyr (>= 2.6.0)`. `arrow` is not
  needed. corehpc needs an R install before any `plasmotrackr` check can run there.
- [How should a plasmotrackr session look?](https://github.com/EPPIcenter/plasmotrack/issues/102): the lazy tables
  are `dbplyr` `tbl_duckdb_connection` objects. Users can call any `dbplyr`-translatable verb, including
  `pivot_wider` and grouped window functions, but `gsub` and `stringr::str_count` don't translate. That argues for
  package helpers for bitstring work and for a `pt_` setter (or option) for the memory cap.
- [How is plasmotrackr tested against pt_analysis?](https://github.com/EPPIcenter/plasmotrack/issues/106): both
  readers use DuckDB. The pooled edges matched `pt_analysis`'s `derived/edges.parquet` exactly on Eswatini chain1. A
  reader-rule fixture (a golden store copy plus an unmarked Window and a `.in-flight.parquet`, as built here) catches
  glob-based regressions on both sides.
