# Storage candidate survey: SQLite vs Parquet for concurrent-read MCMC output

Research notes for [issue #15](https://github.com/EPPIcenter/plasmotrack/issues/15).
**This document does not pick a winner.** It gathers sourced facts and flags what must be
measured. The choice belongs to the separate benchmark-driven decision ticket.

Every factual claim below carries a URL. Where the honest answer is "measure it," that is
stated explicitly and marked **[MEASURE]**.

## Workload being sized against

- 200+ genetic loci, 5–25 alleles per locus, ~3,000 infections, 1,000–5,000 sampling iterations.
- Streams: per-iteration scalars (tens of columns), per-iteration medium-width vectors,
  per-(infection, locus, allele) accumulated counters (~9 million counters), diagnostics.
- Single process, C++20, OpenMP over tempered replica-exchange chains.
- A **separate monitoring process must read while the run writes**.
- Local workstation or HPC **node-local scratch**. Never a networked filesystem.
- Downstream consumer is **R**. Jobs get SIGKILLed.
- Conan 2 + CMake. Existing deps: Boost, Eigen, fmt, nlohmann_json, zlib, OpenMP.

---

## 1. Concurrent reading of a live writer

This is the hardest requirement and it separates the candidates most sharply.

### SQLite in WAL mode — designed for exactly this

- **Readers do not block the writer and the writer does not block readers.** "WAL provides more
  concurrency as readers do not block writers and a writer does not block readers. Reading and
  writing can proceed concurrently." ([sqlite.org/wal.html §1](https://sqlite.org/wal.html))
- **Snapshot isolation.** "In WAL mode, SQLite exhibits 'snapshot isolation'. When a read
  transaction starts, that reader continues to see an unchanging 'snapshot' of the database file
  as it existed at the moment in time when the read transaction started."
  ([sqlite.org/isolation.html](https://sqlite.org/isolation.html)) Mechanically, a reader records
  the last valid commit record in the WAL as its "end mark" and that end mark is fixed for the
  duration of its transaction ([wal.html §2.2](https://sqlite.org/wal.html)).
- **Reader never sees a partial transaction.** "the reader is only able to see complete committed
  transactions from the writer. Partial changes by the writer that have not been committed are
  invisible to the reader… regardless of whether the two database connections are in the same
  thread, in different threads of the same process, or in different processes."
  ([isolation.html](https://sqlite.org/isolation.html))
- **A live monitor sees new data by ending and restarting its read transaction.** "If X wants to
  see the changes that Y made, then X must end its read transaction and start a new one."
  ([isolation.html](https://sqlite.org/isolation.html)) So a polling monitor loop is the correct
  pattern, not a long-lived transaction.

Caveats that matter here, all from the official WAL page:

| Caveat | Detail | Relevance |
| --- | --- | --- |
| Same-host only | "All processes using a database must be on the same host computer; WAL does not work over a network filesystem." The wal-index lives in shared memory backed by a mmapped `-shm` file next to the database. | **Not a problem** — node-local scratch is explicitly the target. |
| Extra files | An extra `-wal` and `-shm` file exist alongside the database. "The WAL file is part of the persistent state of the database and should be kept with the database if the database is copied or moved." | Copying results off a node must move all three files, or the database must be cleanly closed first. |
| Checkpoint starvation | "if a database has many concurrent overlapping readers and there is always at least one active reader, then no checkpoints will be able to complete and hence the WAL file will grow without bound." | A monitor that polls in a tight loop with no gaps can make the WAL grow unboundedly. Mitigation named by the docs: ensure "reader gaps," or run manual `SQLITE_CHECKPOINT_RESTART`/`TRUNCATE` checkpoints (which do block readers briefly). |
| Read performance degrades with WAL size | "read performance deteriorates as the WAL file grows in size since each reader must check the WAL file for the content and the time needed to check the WAL file is proportional to the size of the WAL file." | Argues for keeping `wal_autocheckpoint` enabled. |
| `wal_autocheckpoint` default | Checkpoint runs automatically when the WAL reaches 1000 pages (~4 MB at the 4 KiB default page size). Setting the interval to zero or negative disables autocheckpointing. ([pragma.html#pragma_wal_autocheckpoint](https://sqlite.org/pragma.html)) | Default is probably fine; disabling it is the documented way to get an unbounded WAL. |
| Occasional `SQLITE_BUSY` for readers | Possible if another connection holds `locking_mode=EXCLUSIVE`, during last-connection cleanup, or while a crashed database is being recovered. ([wal.html §9](https://sqlite.org/wal.html)) | The monitor must retry on `SQLITE_BUSY`; it is not a hard failure. |
| Known corruption bug | The "WAL-reset bug" is present in 3.7.0 through 3.51.2, fixed in 3.51.3 (2026-03-13). It "only affects databases in WAL mode when there are two or more database connections open on the same file, in separate threads or processes, and when those two connections attempt to write or checkpoint at the same instant." ([wal.html §11](https://sqlite.org/wal.html)) | Our monitor is read-only, so the two-writer precondition does not hold — but pin ≥ 3.51.3 anyway. ConanCenter has 3.51.3 and 3.53.4. ([conan.io/center/recipes/sqlite3](https://conan.io/center/recipes/sqlite3)) |

Note also: writes are serialised. "since there is only one WAL file, there can only be one writer
at a time" ([wal.html §2.2](https://sqlite.org/wal.html)). Our engine is one process, so this only
constrains internal design (one writer thread, or a mutex).

### Parquet — a reader can see nothing until a file is closed

- **The footer is written last and the file is only interpretable from the footer.** The layout is
  magic `PAR1`, then column chunks, then File Metadata, then a 4-byte metadata length, then `PAR1`
  again. "File metadata is written after the data to allow for single pass writing. Readers are
  expected to first read the file metadata to find all the column chunks they are interested in."
  ([parquet.apache.org/docs/file-format/](https://parquet.apache.org/docs/file-format/))
  Consequence: an in-progress file has no footer, so a reader cannot locate any row group — not
  even the completed ones. Row-group boundaries do not make data visible; closing the file does.
- **Arrow's own docs say concurrent read/write is unsupported.** "The Dataset API offers no
  transaction support or any ACID guarantees. This affects both reading and writing. Concurrent
  reads are fine. Concurrent writes or writes concurring with reads may have unexpected behavior."
  ([Arrow R dataset article](https://arrow.apache.org/docs/r/articles/dataset.html))
- **Conventional workaround: chunked files / a dataset directory.** Arrow's dataset layer reads a
  directory of Parquet files, optionally Hive-partitioned (`key=value` folder names), and
  `open_dataset()` infers a unified schema across them
  ([Arrow R dataset article](https://arrow.apache.org/docs/r/articles/dataset.html)). So the
  incremental pattern is: write chunk *k* to a temp name, close it (footer written), atomically
  rename into the dataset directory; the monitor re-scans the directory each poll. Arrow's docs
  recommend exactly this class of mitigation: "using a unique basename template for each writer, a
  temporary directory for new files, or separate storage of the file list instead of relying on
  directory discovery" (same page).
- **But that scheme has its own sizing limits.** Arrow's guidance: "Avoid files smaller than 20MB
  and larger than 2GB. Avoid partitioning layouts with more than 10,000 distinct partitions," and
  "each additional file adds a little overhead… It also increases the overall Dataset size since
  each file has some shared metadata"
  ([Arrow R dataset article](https://arrow.apache.org/docs/r/articles/dataset.html)). A monitor
  wanting fresh data every few seconds over a multi-hour run pushes toward many small files, which
  collides with the 20 MB floor. The freshness/file-count tradeoff is a real design constraint, not
  a detail.

**Bottom line for criterion 1:** SQLite WAL gives live querying as a documented, first-class
guarantee. Parquet cannot do it at all within a file; it can only approximate it by finalising
files frequently, and the finalisation cadence directly sets both monitor latency and file count.

---

## 2. C++ write path: maturity, ergonomics, build weight

### SQLite

- **Amalgamation**: "Over 100 separate source files are concatenated into a single large file of
  C-code named `sqlite3.c`… The amalgamation contains everything an application needs to embed
  SQLite." As of 2021-12-29 it was "about 238K lines of code (or 145K if you omit blank lines and
  comments) and is over 8.4 megabytes in size."
  ([sqlite.org/amalgamation.html](https://sqlite.org/amalgamation.html)) Single translation unit,
  no transitive dependencies.
- **ConanCenter**: `sqlite3/3.53.4`, with prebuilt binaries for Linux x86_64 (gcc 13), macOS
  arm64/x86_64 (apple-clang 17), and Windows. No listed third-party requirements.
  ([conan.io/center/recipes/sqlite3](https://conan.io/center/recipes/sqlite3))
- **API**: plain C — `sqlite3_open`, `sqlite3_prepare_v2`, `sqlite3_bind_*`, `sqlite3_step`. No
  C++ wrapper needed, but also no type safety; the project would want a thin RAII layer.
- **Longevity/format stability**: SQLite is a US Library of Congress Recommended Storage Format for
  datasets ([sqlite.org/locrsf.html](https://sqlite.org/locrsf.html)) — relevant if run output is
  archived alongside a paper.
- **Licensing**: public domain, no license required. Paid options (a $6,000 warranty-of-title
  license) exist but are optional. ([sqlite.org/prosupport.html](https://sqlite.org/prosupport.html))

### Parquet via Arrow C++

- **ConanCenter**: `arrow/24.0.0`, with prebuilt binaries listed for Linux, Windows, macOS
  (Intel and Apple Silicon). ([conan.io/center/recipes/arrow](https://conan.io/center/recipes/arrow))
- **Parquet is not free-standing.** The recipe validates: `parquet=True` requires `with_thrift=True`,
  and `with_thrift=True` requires `with_boost=True` and `with_zlib=True`. Recipe defaults already set
  `parquet=True, with_thrift=True, with_boost=True, with_zlib=True, with_mimalloc=True`.
  ([recipes/arrow/all/conanfile.py](https://github.com/conan-io/conan-center-index/blob/master/recipes/arrow/all/conanfile.py))
  So a minimal Parquet-only Arrow still pulls **thrift + boost + zlib + mimalloc**, and from
  Arrow 21.0.0 onward `parquet=True` also pulls **rapidjson** (same file, `_requires_rapidjson`).
- **Compiler floor**: Arrow ≥ 23.0.0 requires C++20; Arrow ≥ 10.0.0 requires C++17 (same recipe
  file). We are C++20, so no problem.
- **Arrow's Conan recipe has a history of option/dependency breakage** — e.g. `parquet=True`
  failing on a missing Thrift dependency
  ([conan-center-index#21627](https://github.com/conan-io/conan-center-index/issues/21627)) and
  user options not being respected after the Conan 2 migration
  ([conan-center-index#22506](https://github.com/conan-io/conan-center-index/issues/22506)). Both
  are reported fixed, but it is a maintenance signal for a project whose only current deps are
  header-mostly libraries.
- **Write APIs are good.** Three levels: `parquet::arrow::WriteTable()` for a whole table;
  `parquet::arrow::FileWriter` to "write each batch as a row_group" then `writer->Close()` to
  "write file footer and close"; and `parquet::StreamWriter` with C++ `operator<<` and an explicit
  `EndRowGroup` modifier. ([Arrow C++ Parquet docs](https://arrow.apache.org/docs/cpp/parquet.html))
  `StreamWriter` is the most ergonomic fit for row-at-a-time MCMC emission, at a documented cost:
  its read-side twin is noted as slower "due to the type checking and the fact that column values
  are processed one at a time" (same page).

### **[MEASURE]** Build weight

I could not find authoritative published figures for Arrow C++ build time or installed size with a
Parquet-only option set, nor for the compiled size of the SQLite amalgamation in our configuration.
Do not take numbers from blog posts. Measure directly in this repo:

1. `conan install` + build time and `~/.conan2/p` on-disk size for
   `arrow/*:parquet=True` with everything else off, on a cold cache, vs. the same for `sqlite3`.
2. Incremental rebuild time of plasmotrack after touching a header that includes the storage layer
   (Arrow headers are heavy; SQLite is one C header).
3. Final linked binary size delta for each.

---

## 3. Write throughput and amortisation

### SQLite

The controlling factors are all documented, even though absolute rates are not:

- **Transaction batching is the dominant factor.** SQLite's own benchmark harness notes that the
  `--multitrans` option, which puts each row in its own transaction, "makes SQLite much slower, and
  uncompetitive with direct disk I/O. This option proves, yet again, that to get the most
  performance out of SQLite, you should group as much database interaction as possible within a
  single transaction." ([sqlite.org/fasterthanfs.html §2.2](https://sqlite.org/fasterthanfs.html))
- **`synchronous` sets the fsync cadence.** With `synchronous=NORMAL` in WAL mode, "the WAL file is
  synchronized before each checkpoint and the database file is synchronized after each completed
  checkpoint… but no sync operations occur during most transactions." With `FULL`, an extra WAL sync
  happens after every commit. The docs' own recommendation: "The `synchronous=NORMAL` setting
  provides the best balance between performance and safety for most applications running in WAL
  mode." ([pragma.html#pragma_synchronous](https://sqlite.org/pragma.html))
- **WAL avoids double-writing.** "Write transactions are very fast since they only involve writing
  the content once (versus twice for rollback-journal transactions) and because the writes are all
  sequential." ([wal.html §2.3](https://sqlite.org/wal.html))
- **Very large transactions are the failure mode to avoid.** "WAL works best with smaller
  transactions… For transactions larger than about 100 megabytes, traditional rollback journal
  modes will likely be faster. For transactions in excess of a gigabyte, WAL mode may fail with an
  I/O or disk-full error." ([wal.html §1](https://sqlite.org/wal.html)) Note the same section adds
  that from 3.11.0 "WAL mode works as efficiently with large transactions as does rollback mode" —
  the two statements sit uneasily together on the same page, so treat >100 MB transactions as
  untrusted territory and keep commits well below that.
- **Checkpoint cost is real and can be moved off the hot path.** "Checkpointing does require sync
  operations… These factors combine to make checkpoints slower than write transactions." By default
  the committing thread runs the checkpoint, "causing most COMMIT operations to be very fast but an
  occasional COMMIT (those that trigger a checkpoint) to be much slower. If that effect is
  undesirable, then the application can disable automatic checkpointing and run the periodic
  checkpoints in a separate thread, or separate process." ([wal.html §2.3](https://sqlite.org/wal.html))
- **Prepared statements**: bind parameters rather than embedding literals — "It also runs faster
  since the large string does not need to be parsed or copied as much."
  ([sqlite.org/limits.html](https://sqlite.org/limits.html))

**A sane flush cadence, derived from the above rather than measured:** one prepared statement per
output stream, reused; one transaction spanning one iteration's worth of output across all streams;
`journal_mode=WAL`, `synchronous=NORMAL`; keep `wal_autocheckpoint` at its 1000-page default, or
disable it and checkpoint from a dedicated thread if commit-latency jitter shows up in profiling.
The per-iteration transaction keeps commits far below the 100 MB WAL danger zone and bounds data
loss on SIGKILL to at most one iteration.

**[MEASURE]** There is no official SQLite insert-rate figure to cite. Rates quoted online
(e.g. "50k inserts/sec") are not from sqlite.org and should not be used. Measure: rows/second for
the ~9 million-counter stream at commit granularities of 1 iteration, 10 iterations, and 100
iterations, at `synchronous=NORMAL` and `OFF`, on both a workstation SSD and node-local NVMe.

### Parquet

- **Official row-group guidance is for HDFS-era analytics, not for us**: "We recommend large row
  groups (512MB - 1GB)… An optimized read setup would be: 1GB row groups, 1GB HDFS block size, 1
  HDFS block per HDFS file." Data pages: "We recommend 8KB for page sizes."
  ([Parquet Configurations](https://parquet.apache.org/docs/file-format/configurations/))
- **Large row groups cost memory in the writer**: "Larger groups also require more buffering in the
  write path (or a two pass write)" (same page). Given the current design's ~10 GB of in-process
  buffers is itself a stated problem, GB-scale row-group buffering is a direct conflict.
- **Arrow lets you cap it**: `WriterProperties::Builder().max_row_group_length(...)` "sets an upper
  bound on the number of rows per row group that takes precedent over the `chunk_size` passed in the
  write methods." ([Arrow C++ Parquet docs](https://arrow.apache.org/docs/cpp/parquet.html))
- **Arrow's counter-guidance on small groups**: "Row groups can provide parallelism when reading and
  allow data skipping based on statistics, but very small groups can cause metadata to be a
  significant portion of file size. Arrow's file writer provides sensible defaults for group sizing
  in most cases." ([Arrow R dataset article](https://arrow.apache.org/docs/r/articles/dataset.html))

---

## 4. On-disk size for highly repetitive data

Our data is dominated by values that barely change between iterations, which is the case Parquet's
encodings were built for.

### Parquet's advantages are structural, not incidental

- **Dictionary encoding is on by default in Arrow C++**: "By default, the writer will attempt to
  dictionary encode all supported columns, unless the dictionary grows too large."
  ([Arrow C++ Parquet docs](https://arrow.apache.org/docs/cpp/parquet.html))
- **Dictionary indices are stored with RLE/bit-packing hybrid**, so a column that holds one repeated
  value collapses to a run header plus one value: "The dictionary will be stored in a dictionary page
  per column chunk. The values are stored as integers using the RLE/Bit-Packing Hybrid encoding."
  ([Parquet Encodings](https://parquet.apache.org/docs/file-format/data-pages/encodings/))
- **`DELTA_BINARY_PACKED`** for integers degenerates to almost nothing for constant or
  linearly-increasing sequences — "a block containing all the same values will be bit packed to a
  zero bit width thus being only a header" (same page). Directly relevant to iteration indices and
  slowly-drifting counters.
- **`BYTE_STREAM_SPLIT`** for FLOAT/DOUBLE "does not reduce the size of the data but can lead to a
  significantly better compression ratio and speed when a compression algorithm is used afterwards"
  (same page). Relevant to per-iteration scalar and vector streams of doubles.
- **Block compression is available but off by default in C++**: "Column compression is off by default
  in C++," with SNAPPY, GZIP, BROTLI, LZ4 and ZSTD supported, selectable per column
  ([Arrow C++ Parquet docs](https://arrow.apache.org/docs/cpp/parquet.html)). This must be
  explicitly enabled — a default-configured Parquet writer is *uncompressed*.

### SQLite is row-oriented with no built-in compression

- **No page compression in core SQLite.** Compression requires an extension:
  - **ZIPVFS** — official, from the SQLite developers, but **proprietary and paid: $4,000 one-time
    fee**, supplied as source only. ([sqlite.org/prosupport.html](https://sqlite.org/prosupport.html),
    [ZIPVFS home](https://sqlite.org/zipvfs/doc/trunk/www/index.wiki))
  - **`sqlite_zstd_vfs`** (mlin) — third-party, open source, but the README is blunt about its
    limits: "**USE AT YOUR OWN RISK**… it's young and unlikely to have zero bugs. This project is
    not associated with the SQLite developers," plus "EXCLUSIVE locking mode always applies (a
    writer excludes all other connections for the lifetime of its own connection)" and
    "**WAL mode isn't yet supported; do not touch.**"
    ([README](https://github.com/mlin/sqlite_zstd_vfs))
    **This is disqualifying for our concurrent-read requirement**: exclusive locking plus no WAL
    means no live monitor. Its reported ratios (1182 MiB → 433 MiB tuned on a TPC-H database, same
    README) are therefore of academic interest only.
- **SQLite is not space-hostile in absolute terms**, though: for many small blobs, "a single SQLite
  database holding 10-kilobyte blobs uses about 20% less disk space than storing the blobs in
  individual files," because files are padded to filesystem block size
  ([sqlite.org/fasterthanfs.html](https://sqlite.org/fasterthanfs.html)). That is a comparison
  against the *current* one-file-per-quantity design, not against Parquet.
- **Size limits are not a concern**: max database size is ~281 TB at the maximum 64 KiB page size,
  ~17.5 TB at the 4 KiB default. ([sqlite.org/limits.html](https://sqlite.org/limits.html))

**[MEASURE]** The size gap is the single most decision-relevant unknown and it cannot be read out of
documentation, because it depends entirely on how repetitive our specific columns are. Measure: take
one real (or faithfully simulated) run's output and write it four ways — (a) SQLite with default
page size, (b) SQLite with `page_size=65536`, (c) Parquet with dictionary+RLE and no block
compression, (d) Parquet with dictionary+RLE plus ZSTD — and report bytes on disk plus wall-clock
write time for each. Also record what the current gzip-CSV design produces, as the baseline being
replaced.

---

## 5. R read path

### SQLite → `RSQLite` / `DBI`

- `RSQLite` "Embeds the SQLite database engine in R, providing a DBI-compliant interface"
  ([RSQLite README](https://cran.r-project.org/web/packages/RSQLite/readme/README.html)).
- **Chunked fetching avoids materialising a huge result.** The documented pattern is
  `dbSendQuery()` then `dbFetch(res, n = ...)` in a loop while `!dbHasCompleted(res)` (same page).
  So "read a billion rows into a data.frame" is avoidable, but the analyst has to opt into it —
  `dbReadTable()` / `dbGetQuery()` materialise everything.
- **Push aggregation into SQL.** Because SQLite is a query engine, `SELECT ... GROUP BY` returns
  summaries instead of raw rows, which is the natural fit for convergence plots.
- **Can R safely query a live writer?** Yes, subject to the same WAL rules as any other reader: it
  is a separate database connection in a separate process on the same host, and therefore gets
  committed-only snapshot isolation ([isolation.html](https://sqlite.org/isolation.html)). Two
  practical notes: the R side must be prepared for occasional `SQLITE_BUSY`
  ([wal.html §9](https://sqlite.org/wal.html)), and a monitor that never closes its read transaction
  will starve the writer's checkpointer ([wal.html §6](https://sqlite.org/wal.html)).
- **dplyr backend**: `dbplyr` provides a lazy `tbl()` interface over any DBI connection; this is the
  standard R ecosystem path (documented in DuckDB's own R page, which demonstrates the same
  `tbl(con, ...) |> group_by() |> summarise() |> collect()` pattern,
  [duckdb.org/docs/stable/clients/r](https://duckdb.org/docs/stable/clients/r)).

### Parquet → `arrow`, `nanoparquet`, or `duckdb`

- **`arrow`**: `open_dataset()` reads a directory of Parquet files without loading them. "the data
  values are not loaded into memory. Instead, Arrow scans the data directory… and reads headers of
  the data files to construct a Schema." Queries are lazy: "arrow adopts a lazy evaluation approach
  to queries: when dplyr verbs are called on the Dataset, they record their actions but do not
  evaluate those actions on the data until you run `collect()`", and "all work is pushed down to the
  individual data files… you don't have to load the whole data set in memory to slice from it."
  ([Arrow R dataset article](https://arrow.apache.org/docs/r/articles/dataset.html))
  Caveat from the same page: unsupported dplyr verbs raise an error on Datasets rather than silently
  falling back, so analyst scripts have to stay inside the supported verb set until `collect()`.
- **`nanoparquet`**: MIT-licensed, "Self-sufficient reader and writer for flat 'Parquet' files",
  `Depends: R (>= 4.0.0)` with everything else in `Suggests` — i.e. no heavyweight linked
  dependency, unlike `arrow`.
  ([CRAN: nanoparquet](https://cran.r-project.org/web/packages/nanoparquet/index.html))
  It has no dataset/lazy layer, so it is a "read this one file" tool, not a query engine.
- **`duckdb`**: reads Parquet globs lazily through dbplyr —
  `tbl(con, "read_parquet('dataset/**/*.parquet', hive_partitioning = true)") |> ... |> collect()`,
  described as a way "to avoid reading 12 Parquet files into R's memory"
  ([duckdb.org/docs/stable/clients/r](https://duckdb.org/docs/stable/clients/r)).
  DuckDB can also read SQLite files via its SQLite extension, so choosing SQLite for the writer does
  not preclude a DuckDB-based analysis path.

**Materialisation summary:** all three routes have a lazy/chunked path, so "R will blow up on a
large result" is avoidable in every case. What differs is the default: `dbReadTable()` and
`nanoparquet::read_parquet()` are eager; `arrow::open_dataset()` and `dbplyr::tbl()` are lazy.

---

## 6. Crash and kill robustness (SIGKILL mid-write)

This is where the two options differ most and most simply.

### SQLite: committed transactions survive an application kill, unconditionally

- **The explicit guarantee**: "Transactions are durable across application crashes regardless of the
  synchronous setting or journal mode." ([pragma.html#pragma_synchronous](https://sqlite.org/pragma.html))
  A SIGKILL is an application crash, not an OS crash or power loss — the data is already in the OS
  page cache and the kernel will write it out. So with `synchronous=NORMAL` (or even `OFF`), every
  transaction that returned from `COMMIT` before the kill is recoverable.
- **What is left on disk**: the database plus a non-empty `-wal` and `-shm`. "if the last process to
  have the database open exits without cleanly shutting down the database connection… then the WAL
  file might be retained on disk after all connections to the database have been closed."
  ([wal.html §4](https://sqlite.org/wal.html))
- **Recovery is automatic on next open**: "If the last connection to a database crashed, then the
  first new connection to open the database will start a recovery process."
  ([wal.html §9](https://sqlite.org/wal.html)) The uncommitted tail of the WAL is discarded; the
  reader sees everything up to the last commit record.
- **Do not separate the files**: "If a database file is separated from its WAL file, then
  transactions that were previously committed to the database might be lost, or the database file
  might become corrupted." ([wal.html §4](https://sqlite.org/wal.html)) So the job-teardown copy
  step must move `*.db`, `*.db-wal`, `*.db-shm` together — or open-and-close the database once to
  fold the WAL in ("The only way to safely remove a WAL file is to open the database file using one
  of the `sqlite3_open()` interfaces then immediately close the database using `sqlite3_close()`",
  same section).
- **What is lost**: only the in-flight transaction. With a per-iteration commit, that is one
  iteration.

### Parquet: an unfinalised file has no footer and is not a Parquet file

- Arrow states the consequence directly: "Unexpectedly killing the process while a write is in
  progress can leave the system in an inconsistent state."
  ([Arrow R dataset article](https://arrow.apache.org/docs/r/articles/dataset.html))
- The one piece of good news, same page: "Most file formats have magic numbers which are written at
  the end. This means a partial file write can safely be detected and discarded." Parquet's trailing
  `PAR1` magic ([file-format](https://parquet.apache.org/docs/file-format/)) is exactly this — so a
  torn file is *detectable and discardable*, not silently corrupt. Everything in it is lost, though;
  there is no partial recovery of completed row groups without a footer.
- **Therefore what survives a SIGKILL is exactly the set of files already closed.** Data loss on kill
  equals the file-finalisation interval. This is the same knob as monitor freshness in §1: the two
  requirements pull in the same direction (finalise often), and both fight the 20 MB-minimum
  file-size guidance.
- A recovery tool could in principle rebuild a footer by scanning page headers, but nothing in the
  Parquet spec or Arrow C++ docs offers one; writing that is not a reasonable ask for this project.

---

## 7. Node-local scratch on HPC

- **Node-local NVMe is opt-in and ephemeral.** On OLCF Frontier, each compute node has 2× 1.92 TB
  NVMe SSDs ("Burst Buffer") with peak sequential 5500 MB/s read / 2000 MB/s write; "The NVMes on
  Frontier are local to each node." Access requires `-C nvme` on `sbatch`/`salloc`/`srun`, mounts at
  `/mnt/bb/<user>`, and — critically — "**Users are responsible for moving data to/from the NVMe
  before/after their jobs.**"
  ([Frontier User Guide](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html))
  **Implication:** if the job is SIGKILLed on a scheduler timeout, the copy-back step in the batch
  script may never run and node-local output can vanish. Either budget time for a copy-back and use
  a Slurm signal (`--signal`) to trigger it early, or periodically mirror to the parallel filesystem.
  Whichever store is chosen, this risk is identical — it is a job-script problem, not a format
  problem.
- **Shared scratch is purged on a timer, so it is not an archive either.** NERSC Perlmutter scratch:
  "Any file not accessed within 8 weeks is subject to deletion." OLCF Orion: "files older than 90
  days are purged."
  ([NERSC Perlmutter Scratch](https://docs.nersc.gov/filesystems/perlmutter-scratch/),
  [Frontier User Guide](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html))
- **Networked-filesystem locking is out of scope by assumption**, and that removes SQLite's single
  biggest documented disqualifier: WAL "does not work over a network filesystem" because the
  wal-index needs shared memory ([wal.html §1, §2.2](https://sqlite.org/wal.html)). On node-local
  storage this constraint is satisfied by construction.
- **`fsync` cost**: the tunable is `synchronous`. With `NORMAL` in WAL mode "no sync operations occur
  during most transactions" and the checkpoint carries the syncs
  ([pragma.html#pragma_synchronous](https://sqlite.org/pragma.html)). Combined with "if an
  application therefore runs checkpoint in a separate thread or process, the main thread or process
  that is doing database queries and updates will never block on a sync operation"
  ([wal.html §2.3](https://sqlite.org/wal.html)), fsync cost is controllable.
- **File-size limits**: SQLite's own ceiling (17.5 TB at 4 KiB pages) is far above anything the
  filesystem will allow ([sqlite.org/limits.html](https://sqlite.org/limits.html)); the binding
  constraint will be the node's NVMe capacity or the project quota, not the format.
- **tmpfs vs local SSD**: **[MEASURE]** — I found no authoritative general statement about tmpfs
  suitability across sites, and it is site-specific (tmpfs consumes node RAM, which competes with a
  memory-hungry MCMC). Check the target site's docs before assuming `/tmp` is disk-backed; on
  Frontier, `/tmp` and `/mnt/bb/$USER` are named as distinct node-local locations
  ([Frontier User Guide](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html)).

---

## 8. Third candidates, and why they are or are not worth naming

### DuckDB — worth naming, but the concurrency model is wrong for us

- **Available in ConanCenter** as `duckdb/1.4.3` with prebuilt binaries for Linux/macOS/Windows
  ([conan.io/center/recipes/duckdb](https://conan.io/center/recipes/duckdb)), and it has an
  excellent R client with dbplyr support
  ([duckdb.org/docs/stable/clients/r](https://duckdb.org/docs/stable/clients/r)). Columnar, so it
  would get compression benefits similar to Parquet.
- **But cross-process concurrent read of a live writer is not supported in the stable product.** The
  two in-process options are "one process can both read and write" or "multiple processes can read
  from the database, but no processes can write (`access_mode = 'READ_ONLY'`)." Multi-process write
  access requires either the Quack remote protocol — "in beta stage as of DuckDB v1.5.2, and is
  expected to become mature by DuckDB v2.0 in fall 2026" — or DuckLake with a PostgreSQL catalog.
  ([DuckDB Concurrency](https://duckdb.org/docs/stable/connect/concurrency))
  The R docs restate it plainly: "Read-only mode is required if multiple R processes want to access
  the same database file at the same time."
  ([duckdb.org/docs/stable/clients/r](https://duckdb.org/docs/stable/clients/r))
- **Verdict**: this is precisely the requirement in §1, and DuckDB's own documentation says the
  stable answer is no. Introducing a Postgres catalog on an HPC compute node to work around it is
  disproportionate. DuckDB remains highly relevant as a *read-side* tool over whichever store we
  pick — it reads both Parquet and SQLite.

### HDF5 — worth naming because SWMR exists, but the constraints bite

- **SWMR does what we want, in principle**: "The Single-Writer / Multiple-Reader (SWMR) feature
  enables multiple processes to read an HDF5 file while it is being written to (by a single process)
  without using locks or requiring communication between processes."
  ([HDF5 Intro to SWMR](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html))
- **Available in ConanCenter** as `hdf5/1.14.6` with prebuilt binaries
  ([conan.io/center/recipes/hdf5](https://conan.io/center/recipes/hdf5)).
- **The limitations are the problem**, all from the same official page:
  - The schema is frozen once SWMR starts: "The writer process is only allowed to modify raw data of
    existing datasets by: Appending data along any unlimited dimension. Modifying existing data."
    The writer "cannot add" or delete objects. Every dataset must be created before writing begins —
    workable for a fixed MCMC schema, but rigid.
  - "The writer must periodically flush data" and "The SWMR reader must continually poll for new
    data" by "checking the size of the dataset" — visibility is manual, not transactional.
  - "File space recycling is not allowed. As a result the size of a file [written by a SWMR writer
    may grow]."
  - Requires the latest file format via `H5Pset_libver_bounds`, and requires POSIX write semantics.
  - Startup ordering is fragile: "If a reader process opens a file before a [writer] process, the
    writer will not be able to modify the file due to the lock the [reader] puts on the file… It is
    the user's responsibility to start the writer and readers in the correct order."
    ([SWMR User's Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf))
  - On non-atomic systems the library "might possibly read inconsistent metadata" and relies on
    bounded read-retries via `H5Pset_metadata_read_attempts()`
    ([HDF5 SWMR API group](https://docs.hdfgroup.org/documentation/hdf5/latest/group___s_w_m_r.html)).
- **Verdict**: HDF5 gives no transactional commit boundary, no query language for the R monitor to
  aggregate with, and a reader/writer startup-ordering hazard. It buys us concurrent read at the
  price of hand-rolling everything SQLite gives for free. Name it, then set it aside unless the
  benchmark shows a decisive size or throughput win.

### Not worth naming

- **Zarr, TileDB, LMDB, Feather/Arrow IPC**: none were requested and none clear the bar of "queryable
  by a live R monitor *and* already packaged for Conan 2 *and* first-class in R." Arrow IPC in
  particular shares Parquet's finalisation problem while giving up its compression.

---

## 9. Consolidated open questions for the benchmark ticket

Each of these is a "measure, don't read" item. They are the deliverables the decision ticket needs.

1. **On-disk size**, one realistic run, four ways: SQLite (default page size), SQLite
   (`page_size=65536`), Parquet (dictionary+RLE, no block compression), Parquet (+ZSTD). Include the
   current gzip-CSV baseline. (§4)
2. **Write throughput and commit-cadence sensitivity** for SQLite: rows/sec for the ~9M-counter
   stream at 1 / 10 / 100-iteration commit granularity, `synchronous=NORMAL` vs `OFF`, workstation
   SSD vs node-local NVMe. (§3)
3. **Parquet writer memory** at row-group sizes small enough to give a live monitor useful latency —
   does it stay under budget, given the existing ~10 GB buffer problem? (§3)
4. **Monitor freshness vs. cost**: for SQLite, query latency for a typical convergence-plot query
   against a live writer at several WAL sizes, and whether checkpoint starvation appears under a
   realistic polling loop. For Parquet, the file-count / staleness curve at several finalisation
   intervals, checked against Arrow's 20 MB-minimum guidance. (§1)
5. **Build weight**: cold-cache Conan install+build time and on-disk size for `arrow` with
   `parquet=True` and everything else off, vs `sqlite3`; plus incremental rebuild time and linked
   binary size delta for each. (§2)
6. **Kill test**: SIGKILL a running job at a random point, then verify what a fresh reader recovers,
   for each candidate. Confirm SQLite's "durable across application crashes" claim empirically on
   the target filesystem, and confirm that a torn Parquet file is cleanly rejected rather than
   silently truncated. (§6)
7. **Target-site scratch behaviour**: whether `/tmp` is tmpfs or disk on the actual HPC systems in
   use, and whether the batch scripts survive a timeout kill with output intact. (§7)

---

## Source index

SQLite (all sqlite.org, primary):
[WAL](https://sqlite.org/wal.html) ·
[Isolation](https://sqlite.org/isolation.html) ·
[PRAGMA](https://sqlite.org/pragma.html) ·
[Amalgamation](https://sqlite.org/amalgamation.html) ·
[35% Faster Than The Filesystem](https://sqlite.org/fasterthanfs.html) ·
[Limits](https://sqlite.org/limits.html) ·
[Pro Support / licensing](https://sqlite.org/prosupport.html) ·
[ZIPVFS](https://sqlite.org/zipvfs/doc/trunk/www/index.wiki) ·
[LoC Recommended Storage Format](https://sqlite.org/locrsf.html)

Parquet / Arrow (primary):
[Parquet file format](https://parquet.apache.org/docs/file-format/) ·
[Parquet configurations](https://parquet.apache.org/docs/file-format/configurations/) ·
[Parquet encodings](https://parquet.apache.org/docs/file-format/data-pages/encodings/) ·
[Arrow C++ Parquet](https://arrow.apache.org/docs/cpp/parquet.html) ·
[Arrow R datasets](https://arrow.apache.org/docs/r/articles/dataset.html)

Packaging (primary):
[ConanCenter sqlite3](https://conan.io/center/recipes/sqlite3) ·
[ConanCenter arrow](https://conan.io/center/recipes/arrow) ·
[ConanCenter duckdb](https://conan.io/center/recipes/duckdb) ·
[ConanCenter hdf5](https://conan.io/center/recipes/hdf5) ·
[arrow conanfile.py](https://github.com/conan-io/conan-center-index/blob/master/recipes/arrow/all/conanfile.py) ·
[CCI #21627](https://github.com/conan-io/conan-center-index/issues/21627) ·
[CCI #22506](https://github.com/conan-io/conan-center-index/issues/22506)

R (primary):
[RSQLite README](https://cran.r-project.org/web/packages/RSQLite/readme/README.html) ·
[CRAN nanoparquet](https://cran.r-project.org/web/packages/nanoparquet/index.html) ·
[DuckDB R client](https://duckdb.org/docs/stable/clients/r)

Other engines (primary):
[DuckDB concurrency](https://duckdb.org/docs/stable/connect/concurrency) ·
[HDF5 SWMR intro](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html) ·
[HDF5 SWMR API](https://docs.hdfgroup.org/documentation/hdf5/latest/group___s_w_m_r.html) ·
[HDF5 SWMR User's Guide (PDF)](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf) ·
[sqlite_zstd_vfs](https://github.com/mlin/sqlite_zstd_vfs)

HPC (primary):
[OLCF Frontier User Guide](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html) ·
[NERSC Perlmutter Scratch](https://docs.nersc.gov/filesystems/perlmutter-scratch/)
