# Which store formats meet the output constraints?

Research for [Which store formats meet the output constraints?](https://github.com/EPPIcenter/plasmotrack/issues/27), part of [Find the way to a unified output spec](https://github.com/EPPIcenter/plasmotrack/issues/26).

This note is facts against the locked constraints, not a pick.

Domain terms follow `CONTEXT.md`: Trace, Block, Window, Run, Sampling iteration, Run event, Run manifest. Streams are never called logs.

## Question

Which store formats can satisfy all of the locked constraints, and how do they compare?

## Locked constraints

All of these must hold:

1. C++ can write during a Run (online), not only as a post-hoc convert.
2. A Window can be committed as the durability boundary; a crash loses at most the in-flight Window.
3. Another process can read a consistent prefix of closed Windows while the Run is still writing.
4. An existing native R package can read the result (no Python, no new R package). Named packages to check: Rarr, rhdf5, arrow.
5. Local POSIX must work. Shared filesystem (NFS/Lustre) should be characterized. Object store is a plus, not required.
6. One store per Run, not millions of files per quantity.
7. Streams include dense numeric Traces and ragged / variable-length streams (parent-set distributions, genotypes, Orderings).

## Method

Claims below are taken from official specifications, first-party library documentation, or first-party library source. Secondary write-ups are not used as evidence.

## Verdict at a glance

| Candidate | Survives all constraints? | Short reason |
| --- | --- | --- |
| Zarr via TensorStore (C++) | Yes, with layout rules | Online C++ write; chunk/key commit; Rarr reads v2 and v3; POSIX + object stores. File-count and ragged encodings are layout choices, not format kills. |
| Zarr via other C++ (z5, xtensor-zarr) | Same format, weaker writer | Same on-disk format as TensorStore Zarr, so Rarr still reads. No first-party transaction / `fsync` story comparable to TensorStore. |
| HDF5 (libhdf5 / HighFive) + rhdf5 | No, if R is the mid-run reader | Library SWMR matches the Window/prefix story, but rhdf5 does not support SWMR. SWMR also forbids variable-length datatypes and does not work on NFS. |
| TileDB (C++) + tiledb (CRAN) | Yes, with layout rules | Fragment + commit file is a Window boundary; concurrent readers see only committed fragments. R package exists. Numeric variable-length lists are not a first-class R read path. |
| Arrow / Parquet, one growing file | No | Footer is written after the data; a reader cannot open an incomplete file; a crash before close loses the whole file. |
| Arrow / Parquet, one closed file per Window | Yes, with a file-count budget | C++ can write; R `arrow` can `open_dataset` a directory of complete files. File count scales with Windows × streams. |
| Arrow IPC streaming format | Yes for prefix reads; R path is thinner | Stream messages have no footer. R `open_dataset(..., format = "ipc")` is the random-access / Feather-v2 file format, not a documented mid-stream reader. |
| SQLite (WAL) + RSQLite | Local POSIX only | Transaction = Window. WAL lets readers see committed data. Official docs: WAL does not work over a network filesystem. |

## How a Window maps onto a store

A Window is the contiguous range of Sampling iterations a Block covers, and the unit of loss when a Run is killed (`CONTEXT.md`). The store must make “this Window is closed” a durable, reader-visible event.

Three mechanisms appear in the sources:

- **Atomic key / file replace.** Zarr stores one encoded chunk under one key. TensorStore’s `file` kvstore writes via a temporary file and rename, and can `fsync` so a committed key survives a crash ([file kvstore](https://google.github.io/tensorstore/kvstore/file/index.html)).
- **Flush of appended dataset extent.** HDF5 SWMR appends along an unlimited dimension and calls `H5Dflush`; readers call `H5Drefresh` and poll extent ([SWMR technical note](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)).
- **Commit sidecar.** TileDB writes an immutable fragment, then a commit file as the last step; readers ignore fragments without a commit ([Commits](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/commits/)).
- **Closed footer / closed transaction.** Parquet and Arrow IPC *file* format become readable only after the footer is written ([Parquet file format](https://parquet.apache.org/docs/file-format/); [Arrow columnar / IPC](https://arrow.apache.org/docs/format/Columnar.html)). SQLite COMMIT appends a commit record to the WAL ([WAL](https://www.sqlite.org/wal.html)).

A consistent prefix of closed Windows means: a reader must never observe a Window that is only half written, and must be able to ignore the in-flight Window.

## Ragged streams (cross-cutting)

The constraint is that streams include both dense numeric Traces and ragged / variable-length streams. Native variable-length types are not required if a layout of two dense arrays (offsets + payload) can be written and read by the same R package.

That layout is available in every surviving candidate. Native variable-length support is uneven and is called out per candidate, because it is a failure mode if the spec later insists on a native vlen type rather than an offsets+payload encoding.

## Filesystems (cross-cutting)

| Store | Local POSIX | NFS | Lustre / GPFS | Object store |
| --- | --- | --- | --- | --- |
| Zarr + TensorStore `file` | Supported; `fsync` optional ([file kvstore](https://google.github.io/tensorstore/kvstore/file/index.html)) | “Concurrent access from multiple machines to a network filesystem is also safe” if shared locking is supported ([file kvstore](https://google.github.io/tensorstore/kvstore/file/index.html)). Locking uses `flock` by default. | Not separately documented. Same locking caveat as any POSIX network FS. | First-party `gcs` and `s3` kvstores ([kvstore index](https://google.github.io/tensorstore/kvstore/index.html)) |
| HDF5 SWMR | Required: POSIX `write()` semantics ([SWMR note](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)) | Officially does not work: write ordering may not be preserved ([SWMR User’s Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)) | Tested on Lustre and GPFS ([SWMR User’s Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)) | rhdf5 can open `http(s)` via the HDF5 read-only S3 VFD ([H5Fopen](https://huber-group-embl.github.io/rhdf5/reference/H5Fopen.html)); not a writer target |
| TileDB | Supported | Lock-free for ordinary read/write. Vacuuming uses file locks; Lustre needs `-o flock` or `-o localflock` ([Concurrency](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/concurrency/)) | Same locking note as above | First-class S3 / GCS / Azure ([CRAN tiledb](https://cran.r-project.org/web/packages/tiledb/index.html)) |
| Arrow / Parquet files | Ordinary POSIX files | Write-once complete files; no SWMR-style in-place append | Same | R `open_dataset` accepts a FileSystem / URI ([open_dataset](https://arrow.apache.org/docs/r/reference/open_dataset.html)) |
| SQLite WAL | Supported | Officially does not work: WAL requires shared memory on one host ([WAL](https://www.sqlite.org/wal.html); [use over a network](https://www.sqlite.org/useovernet.html)) | Same class of problem as NFS | Not a store |

The map’s note that shared filesystem should be characterized if that is where Runs write is material: HDF5 SWMR and SQLite WAL are documented as NFS-unsafe. Zarr-on-files and TileDB fragment writes do not depend on POSIX write-order of a single growing file.

---

## 1. Zarr written via TensorStore

### Writer

TensorStore is a C++ library (Python bindings sit on top) for reading and writing large arrays, including Zarr, with transactions and optimistic concurrency ([TensorStore](https://google.github.io/tensorstore/)).

- C++ open of a Zarr array on the local filesystem is a first-party API: `tensorstore::Open` with `{"driver", "zarr"}` and `{"kvstore", {{"driver", "file"}, {"path", ...}}}` ([tensorstore::Open](https://google.github.io/tensorstore/cpp/api/tensorstore.Open.html)).
- The `zarr2` / `"zarr"` driver reads, writes, creates, and **resizes** Zarr v2 arrays ([zarr2 driver](https://google.github.io/tensorstore/driver/zarr/index.html)).
- The `zarr3` driver does the same for Zarr v3, and will insert a `sharding_indexed` codec when read and write chunk shapes differ ([zarr3 driver](https://google.github.io/tensorstore/driver/zarr3/index.html)).
- `tensorstore::Resize` grows or shrinks bounds ([tensorstore::Resize](https://google.github.io/tensorstore/cpp/api/tensorstore.Resize.html)).

Online write during a Run is therefore a supported C++ path, not a convert-after-the-fact path.

### R reader, versions, codecs

[Rarr](https://github.com/Huber-group-EMBL/Rarr) is a native R/C package (Bioconductor; no Python). The supported-features vignette states ([Rarr features](https://huber-group-embl.github.io/Rarr/articles/features.html)):

- Zarr version 2 and Zarr version 3.
- Stores: `FilesystemStore` read/write; `S3Store` read only.
- Compression codecs, all read and write: zlib/gzip, bzip2, blosc, LZMA, LZ4, Zstd.
- Other codecs: `endian` read only; `transpose` read/write; `VLen-UTF8` read/write (write only via `update_zarr_array`); `sharding_indexed` **read only**.
- Numeric types through float64 are readable; `void *` is not; V2/V3 structured types are read-only.

TensorStore’s default Zarr v2 compressor is `{"id": "blosc"}`; also documented: zlib, bz2, zstd ([zarr2 driver](https://google.github.io/tensorstore/driver/zarr/index.html)). Those ids are in Rarr’s read set.

TensorStore Zarr v2 currently does **not** support filters (`filters` must be `null`) ([zarr2 driver](https://google.github.io/tensorstore/driver/zarr/index.html)). Zarr v2 object / vlen codecs (`vlen-utf8`, `vlen-bytes`, `vlen-array`) are filters in the v2 model ([Zarr Python data types](https://zarr.readthedocs.io/en/3.1.x/user-guide/data_types.html)). So TensorStore cannot write native Zarr vlen object arrays. Ragged streams have to be encoded as dense arrays (offsets + payload, or a fixed-width struct), which Rarr can read as ordinary numeric / string arrays.

### Window commit

Two complementary mechanisms:

1. **Per-chunk durability.** The `file` kvstore, with default `file_io_sync: true`, “ensures all writes are durable” via `fsync`. Writes go through a temporary file and rename except in `non_atomic` mode. “Partial writes will not be observed; unconditional writes are still atomic” even with locking `none` ([file kvstore](https://google.github.io/tensorstore/kvstore/file/index.html)).
2. **Multi-key isolation.** Isolated transactions: “no modifications made are visible or persist outside the transaction until the transaction is committed.” Atomic-isolated transactions additionally require a single atomic commit and fail if the store cannot provide it ([tensorstore.Transaction](https://google.github.io/tensorstore/python/api/tensorstore.Transaction.html); [TransactionMode](https://google.github.io/tensorstore/cpp/api/tensorstore.TransactionMode.html)). Isolated (not atomic-isolated) is the mode that works on a directory of files: each key is atomic, the set of keys is not one rename.

A Window commit that matches the durability constraint is therefore:

- buffer the Window in memory;
- write every chunk that belongs to that Window (and any Run-event / Run-manifest keys that close with it) inside an isolated transaction, or write chunks then a single “closed-through-Window-N” marker key last;
- `fsync` on (default);
- readers treat missing chunks as fill ([Zarr v2 spec](https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html); [Zarr v3 core](https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html)) and trust only the marker / resized shape for “how far the Run has closed.”

A crash during the transaction loses that Window’s uncommitted keys. Already-committed Windows remain.

### Mid-run read

Zarr is a key/value of independent chunks. A reader that opens the same directory sees whatever keys exist. Unwritten chunks are fill ([Zarr v2](https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html)). TensorStore documents that isolated writes are invisible to other readers until commit ([Transaction](https://google.github.io/tensorstore/python/api/tensorstore.Transaction.html)).

Rarr’s `read_zarr_array` reads a path and an index; `zarr_overview` reads `.zarray` / `zarr.json` ([Rarr manual](https://bioc.r-universe.dev/Rarr/doc/manual.html)). There is no SWMR protocol. Consistency is “complete objects on disk.” That is enough if the writer never replaces a chunk in place and publishes a Window only after its objects exist.

TensorStore `recheck_cached_data` defaults to revalidate on every read; `recheck_cached_metadata` defaults to `"open"` ([zarr2 driver](https://google.github.io/tensorstore/driver/zarr/index.html)). A mid-run Rarr process that keeps metadata cached across Window closes would need to re-read `.zarray` / `zarr.json` / the marker key.

### File-count behavior

Without sharding, each chunk is one storage object ([Zarr v2 chunks](https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html); [Zarr v3](https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html)). One store per Run is one hierarchy directory. File count is `chunks × arrays`, not “one file per Sampling iteration,” unless the chunk shape is one iteration.

The v3 `sharding_indexed` codec packs many inner chunks into one shard object specifically because “the file block size and maximum inode number restrict the usage of numerous small files” ([sharding codec](https://zarr-specs.readthedocs.io/en/latest/v3/codecs/sharding-indexed/)). TensorStore inserts that codec when write chunks (shards) are larger than read chunks ([zarr3 driver](https://google.github.io/tensorstore/driver/zarr3/index.html)). Rarr can **read** `sharding_indexed` but cannot write it ([Rarr features](https://huber-group-embl.github.io/Rarr/articles/features.html)). That is acceptable: C++ writes, R reads.

TensorStore also has an `ocdbt` kvstore adapter that packs many keys into a B+tree on a base store ([ocdbt](https://google.github.io/tensorstore/kvstore/ocdbt/index.html)). Rarr talks to a filesystem or S3 *Zarr* store, not OCDBT. Using OCDBT would hide keys from Rarr and fail constraint 4.

### Failure modes that would kill this spec

- Chunking at one Sampling iteration without sharding: millions of files per Trace.
- Using OCDBT or any non-Zarr kvstore adapter under TensorStore: Rarr cannot open it.
- Relying on TensorStore v2 filters / vlen object dtypes: not implemented; ragged must be dense encodings.
- `atomic_isolated` transactions across many files: the file kvstore cannot commit many keys as one atomic rename; the transaction fails rather than lying ([TransactionMode](https://google.github.io/tensorstore/cpp/api/tensorstore.TransactionMode.html)). Isolated + last-key marker is the workable pattern.
- `file_io_sync: false`: committed data may be lost on crash ([file kvstore](https://google.github.io/tensorstore/kvstore/file/index.html)).
- Network FS without working shared locks, if the writer uses read-modify-write on a shard. Unconditional whole-chunk writes remain atomic by rename.

---

## 2. Zarr written some other C++ way

First-party C++ writers besides TensorStore:

- **z5** — C++ library with filesystem and S3; Zarr v2 and v3; `createDataset` / `writeSubarray` ([z5](https://github.com/constantinpape/z5)).
- **xtensor-zarr** — C++ Zarr v2/v3 hierarchy on a filesystem or GCS store ([xtensor-zarr](https://xtensor-zarr.readthedocs.io/en/latest/); [basic usage](https://xtensor-zarr.readthedocs.io/en/latest/basic_usage.html)).

The on-disk format is still Zarr, so the Rarr story, chunk-as-object story, and fill-value mid-run story are the same as §1.

What these writers do **not** document, and TensorStore does:

- isolated transactions with invisible-until-commit writes;
- `fsync` as a first-class durability switch;
- automatic v3 sharding insertion.

xtensor-zarr’s own docs describe a v3-oriented stack; maintenance status is not a primary-source claim used here. For this spec they are viable only if the application itself does temp-file+rename (or equivalent) per chunk and publishes a Window marker last. That is implementable, but it is not a library guarantee.

Failure mode that would kill a non-TensorStore C++ Zarr writer: writing a chunk by opening the final path and streaming into it, so a crash or a concurrent Rarr read observes a partial object. The Zarr spec does not define atomic writes; atomicity is a store implementation property ([Zarr v2](https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html)).

---

## 3. HDF5

### Writer

The HDF5 C library writes during a process lifetime. SWMR is the documented online-append path ([SWMR technical note](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)). HighFive is a header-only C++ wrapper and exposes `File::ReadSWMR`, `File::WriteSWMR`, and `startSWMRWrite()` for HDF5 ≥ 1.10 ([HighFive File](https://highfive-devs.github.io/highfive/class_high_five_1_1_file.html); [H5File_misc.hpp](https://highfive-devs.github.io/highfive/_h5_file__misc_8hpp_source.html)).

### R reader, versions, codecs

[rhdf5](https://huber-group-embl.github.io/rhdf5/) is the named native R package. It wraps the C API and ships compression via Rhdf5lib / rhdf5filters ([h5createDataset](https://huber-group-embl.github.io/rhdf5/reference/h5createDataset.html)):

- GZIP / ZLIB / DEFLATE, SZIP, BZIP2, BLOSC_{BLOSCLZ,LZ4,LZ4HC,SNAPPY,ZLIB,ZSTD}, LZF, NONE.

`H5Fopen` flags are only `H5F_ACC_RDWR` and `H5F_ACC_RDONLY`. The man page states: “HDF5's Single Write Multiple Reader (SWMR) mode is not currently supported via rhdf5” ([H5Fopen](https://huber-group-embl.github.io/rhdf5/reference/H5Fopen.html)).

[hdf5r](https://cran.r-project.org/package=hdf5r) is another existing native R package (CRAN). It is a near-complete R6 wrapper of the C API and loads constants from the linked libhdf5 ([H5constants.R](https://raw.githubusercontent.com/hhoeflin/hdf5r/master/R/H5constants.R)). The package manual does not document a SWMR workflow; whether `H5F_ACC_SWMR_READ` is present depends on the linked HDF5 ≥ 1.10. That is not a documented, supported mid-run path.

### Window commit

SWMR writer: create the file with `H5F_LIBVER_LATEST`, create all groups/datasets (no new objects after SWMR starts), `H5Fstart_swmr_write`, append, `H5Dflush` ([SWMR note](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)). `H5Dflush` “flushes the specified dataset’s metadata from the metadata cache to the file and all raw data buffers associated with the dataset to the file. If the dataset is chunked, raw data chunks are written to the file when this call is issued” ([SWMR User’s Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)).

A Window commit is: extend the unlimited dimension by the Window’s rows, write, `H5Dflush`. The official interrupted-writer example leaves a valid file containing data up to the interrupt ([SWMR note](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)).

### Mid-run read

Library-level: reader opens with `H5F_ACC_RDONLY | H5F_ACC_SWMR_READ`, polls extent, `H5Drefresh` ([SWMR note](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)). No IPC, no locks between writer and readers. File locking exists only to keep a second writer out ([HDF5 file locking](https://support.hdfgroup.org/documentation/hdf5/latest/_file_lock.html)).

rhdf5 cannot open in SWMR-read mode ([H5Fopen](https://huber-group-embl.github.io/rhdf5/reference/H5Fopen.html)). Opening a file that is being written, without SWMR, is the pre-1.10 unsupported case: “concurrent access to an HDF5 file by multiple processes, when one or more processes is a writer, was not supported” ([file locking](https://support.hdfgroup.org/documentation/hdf5/latest/_file_lock.html)). Analysis is in R. If the mid-run reader is R via rhdf5, constraint 3 fails.

### File-count behavior

One HDF5 file holds the whole Run: groups for Traces, Blocks, Run events, Run manifest. This is the strongest “one store per Run” shape among the candidates.

### Ragged streams

SWMR limitations ([SWMR note](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html); [SWMR User’s Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)):

- writer cannot add or delete objects after SWMR starts;
- writer cannot modify or append variable-length, string, or region-reference datatypes;
- no file-space recycling.

Native HDF5 vlen for genotypes / parent-set distributions / Orderings is incompatible with SWMR. An offsets+payload pair of numeric datasets, created before `H5Fstart_swmr_write` and appended as fixed-size data, is allowed. rhdf5 can create those datasets and read them after the Run; it still cannot SWMR-read them mid-Run.

### Failure modes that would kill this spec

- **rhdf5 is the R reader and mid-run read is required.** rhdf5 has no SWMR.
- **NFS (or SMB) is where Runs write.** “The SWMR feature does not work on network file systems (NFS, SMB/Windows file shares, and so on) since write ordering may not be preserved” ([SWMR User’s Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)). `flock` also does not lock over NFS; HDF5 did not consider that a problem *because* SWMR is unsupported there ([file locking](https://support.hdfgroup.org/documentation/hdf5/latest/_file_lock.html)).
- **Native vlen / variable-length strings as the ragged encoding under SWMR.** Forbidden.
- Creating new datasets after the Run has entered SWMR (for example a late-discovered stream). Forbidden.
- Using HDF5 1.8 readers on a 1.10 SWMR file: new chunk index; 1.8 tools cannot read it ([SWMR User’s Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)).

HDF5 remains a strong *library* match for Window flush and prefix reads on local POSIX or Lustre/GPFS. It does not survive the conjunction of “R reads mid-run” and “rhdf5 is the reader,” and it does not survive NFS.

---

## 4. TileDB

### Writer

TileDB Embedded is a C++ library with a documented C++ API ([TileDB C++ API](https://tiledb-inc-tiledb.readthedocs-hosted.com/en/stable/c++-api.html)). Writes are online: each write produces a fragment.

### R reader, versions, codecs

[tiledb](https://cran.r-project.org/web/packages/tiledb/index.html) on CRAN (0.33.0 as of 2025-10-01) is an existing native R package: Rcpp bindings to TileDB Embedded, R ≥ 4.3. It is not one of the three named packages, but it satisfies “existing native R package” if that package is acceptable.

The package description lists chunked arrays, multiple compression / encryption / checksum filters, S3 / GCS / Azure, and time travel ([CRAN tiledb](https://cran.r-project.org/web/packages/tiledb/index.html)). Variable-length string attributes are a documented R read/write path. Variable-length *lists of numerics* can be written from R via low-level offset buffers; the official tutorial then says reading those lists “is currently supported with the Python API. Please use the Python API instead” ([Variable-length attributes](https://documentation.cloud.tiledb.com/academy/structure/arrays/tutorials/basics/var-length-attributes/)). That is a ragged-stream gap unless the layout is offsets+payload as ordinary attributes, or the values are stored as variable-length strings / bytes.

### Window commit

A fragment is a timestamped write. It is invisible until a commit file is written in `__commits` as the last step, and only after all fragment files exist ([Fragments](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/fragments/); [Commits](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/commits/)). “This guarantees write atomicity and consistency, even if a fragment write fails with partially written files.”

A Window commit is one write query that becomes one fragment plus one commit file. A crash before the commit file leaves the fragment invisible. Readers see a consistent prefix of committed Windows.

### Mid-run read

“You can mix concurrent reads and writes. Fragments are not visible unless the write query completes successfully (and a `.ok` file appears).” Lock-free multiple reader/writer ([Concurrency](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/concurrency/)). Readers that opened earlier see the fragment set at open time; a new open is required to pick up later commits (multiversion / open-time view; [ACID](https://documentation.cloud.tiledb.com/academy/structure/tables/foundation/key-concepts/acid/)).

### File-count behavior

One array URI is one store. Each write creates a fragment directory with per-attribute files plus metadata ([Fragments](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/fragments/)). File count grows with Windows unless fragments are consolidated. Consolidation is process-safe and can run while readers continue; vacuuming deletes old fragments and needs file locks ([Concurrency](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/concurrency/); [Consolidation](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/consolidation/index.html)).

One fragment per Window, then periodic consolidate+vacuum, keeps file count at “one store per Run” rather than “one file per Sampling iteration.” One fragment per Sampling iteration would be the millions-of-files failure.

### Failure modes that would kill this spec

- Treating `tiledb` as unacceptable because it is not Rarr/rhdf5/arrow. Then there is no named-package reader.
- Writing a fragment per Sampling iteration and never consolidating.
- Vacuuming on a lock-less object store while a reader is time-traveling to the fragments being deleted ([Concurrency](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/concurrency/)).
- Relying on R to read native numeric variable-length lists. Official docs send that read to Python.
- Creating the array from two processes at once: array creation is neither thread-safe nor process-safe ([Concurrency](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/concurrency/)).

---

## 5. Arrow / Parquet

### Writer

Apache Arrow C++ writes Parquet via `parquet::ParquetFileWriter` / `StreamWriter` ([C++ Parquet](https://arrow.apache.org/docs/dev/cpp/parquet.html); [file_writer.h](https://github.com/apache/arrow/blob/master/cpp/src/parquet/file_writer.h)). Row groups can be appended *while the writer is still open*. Closing the writer writes the footer.

### R reader, versions, codecs

The CRAN [arrow](https://cran.r-project.org/web/packages/arrow/index.html) package (25.0.0 as of 2026-07-16) binds Arrow C++. It reads and writes Parquet, Arrow IPC / Feather v2, and multi-file datasets ([arrow R docs](https://arrow.apache.org/docs/r/)). `open_dataset` defaults to Parquet; `format` may be `"parquet"`, `"ipc"` / `"arrow"` / `"feather"` (Feather v2 only), CSV, TSV, JSON ([open_dataset](https://arrow.apache.org/docs/r/reference/open_dataset.html)).

Writer-side compression is configured on `WriterProperties` (Snappy, Zstd, Brotli, …) ([C++ Parquet](https://arrow.apache.org/docs/dev/cpp/parquet.html)). The R package reads whatever the linked Arrow C++ build decodes.

Ragged streams are first-class: Arrow list / variable-size types are in the columnar spec; IPC converts their offset buffers on endian swap ([Arrow IPC C++](https://arrow.apache.org/docs/cpp/api/ipc.html)).

### Window commit — single growing file (dies)

Parquet layout ([file format](https://parquet.apache.org/docs/file-format/)):

```
PAR1 | row groups… | File Metadata | metadata length | PAR1
```

“File metadata is written after the data to allow for single pass writing. Readers are expected to first read the file metadata.” The format repo notes that if an error happens while writing the file metadata, all the data written will be unreadable ([parquet-format](https://github.com/apache/parquet-format)).

A mid-run reader of a file the C++ process is still appending cannot find the footer. A crash before `Close()` loses every Window in that file, not only the in-flight Window. That fails constraints 2 and 3.

### Window commit — one closed Parquet file per Window (survives)

Write Window *k* to `run/<stream>/window=k/part-0.parquet` (or any Hive-style path), close the writer (footer present), then start the next file. `open_dataset` on the parent directory discovers complete files ([open_dataset](https://arrow.apache.org/docs/r/reference/open_dataset.html); [dataset article](https://arrow.apache.org/docs/r/articles/dataset.html)). The in-flight file is either not created yet or is excluded (`exclude_invalid_files` / prefix ignore).

A crash loses only the unclosed file = the in-flight Window.

### Mid-run read

For the directory-of-files layout: a new `open_dataset` sees every file that already has a footer. There is no need for SWMR. The reader must not open the writer’s current temp/partial path.

### File-count behavior

One directory per Run is one dataset. File count is `Windows × streams` (or fewer if several streams share a file). That is one store per Run, but it is not one file. A long Run with a small Window produces many files. That is the constraint-6 risk: not millions per quantity if Window is coarse, millions if Window is one Sampling iteration and every quantity is its own tree.

### Failure modes that would kill this spec

- One Parquet file opened at Run start and closed at Run end.
- Window = one Sampling iteration and one file per quantity: file explosion.
- Treating an in-flight `.parquet` as readable.
- Expecting to append to a closed Parquet file from C++ as the mid-run protocol. The format is single-pass; the footer is the catalog.

---

## 6. Arrow IPC (streaming vs file)

Arrow defines two IPC binaries ([Columnar format](https://arrow.apache.org/docs/format/Columnar.html)):

- **Streaming format.** Schema, then a sequence of encapsulated messages. End of stream is a continuation + zero-length or closing the stream. No footer. A reader can consume complete record batches from the start without waiting for the writer to finish.
- **File / random-access format.** Magic `ARROW1`, an embedded stream, then a footer with offsets. Readers seek to the end for the footer.

C++ writers: `MakeStreamWriter` / `MakeFileWriter` → `RecordBatchWriter` ([C++ IPC](https://arrow.apache.org/docs/dev/cpp/ipc.html)).

A Window as one record batch (or one closed stream file) matches constraint 2 for the *streaming* format: a crash loses the incomplete last message, not prior batches. Mid-run prefix read is what the stream format is for.

R `open_dataset(..., format = "ipc")` is documented as Feather v2 / Arrow *file* format ([open_dataset](https://arrow.apache.org/docs/r/reference/open_dataset.html)). That path needs a footer. A growing stream file is not a documented R dataset source. Constraint 4 therefore holds for **closed IPC files per Window** (same as Parquet), and is not documented for **one growing stream** that R opens while C++ is still writing.

Failure modes: using the file format as if it were a stream (no footer until close); assuming R `arrow` can tail a live stream.

---

## 7. SQLite (extra candidate)

### Writer

The SQLite C library (`sqlite3`) writes from C++. A transaction is the durability unit.

### R reader

[RSQLite](https://cran.r-project.org/package=RSQLite) is an existing native R package. Not in the named three; same “or equivalent” status as tiledb.

### Window commit and mid-run read

WAL: “A COMMIT occurs when a special record indicating a commit is appended to the WAL.” Readers continue from the original database plus committed WAL frames; “readers do not block writers and a writer does not block readers” ([WAL](https://www.sqlite.org/wal.html)). A crash rolls back the open transaction = the in-flight Window.

### File-count

One `.db` (plus `-wal` / `-shm` in WAL mode) per Run.

### Filesystems

“All processes using a database must be on the same host computer; WAL does not work over a network filesystem” ([WAL](https://www.sqlite.org/wal.html)). Network filesystems do not support simultaneous reads and writes while keeping the database consistent ([use over a network](https://www.sqlite.org/useovernet.html)). If Runs write on NFS/Lustre, this candidate dies. Object store is not a SQLite store.

### Ragged

BLOBs or child tables. Fine for genotypes / Orderings / parent-set distributions as relational rows.

### Failure modes that would kill this spec

- Shared filesystem between writer host and reader host under WAL.
- Very large transactions (WAL docs warn above tens of megabytes historically; large transactions are now supported but still a design smell) ([WAL](https://www.sqlite.org/wal.html)).
- Requiring object-store output.

---

## Comparison of survivors

Assuming the mid-run reader may be R, and local POSIX is required:

| | Zarr + TensorStore + Rarr | TileDB + tiledb | Parquet/IPC files per Window + arrow | HDF5 + rhdf5 | SQLite + RSQLite |
| --- | --- | --- | --- | --- | --- |
| C++ online write | Yes | Yes | Yes | Yes | Yes |
| Window = durable commit | Chunk keys + marker / isolated txn | Fragment commit file | Close file (footer) | `H5Dflush` | `COMMIT` |
| Mid-run consistent prefix | Complete keys / fill | Committed fragments only | Complete files only | Library SWMR yes; **rhdf5 no** | WAL on one host |
| Named R package | Rarr | tiledb (not named) | arrow | rhdf5 (no SWMR) | RSQLite (not named) |
| Dense Traces | Yes | Yes | Yes | Yes | Yes |
| Ragged native | No in TensorStore v2; encode dense | Strings yes; numeric lists not a supported R read | Arrow lists | SWMR forbids vlen | Tables / BLOB |
| One store / not millions of files | Hierarchy; shard or coarse chunks | Array URI; consolidate | Directory; file count = Windows | **One file** | **One file** |
| Local POSIX | Yes | Yes | Yes | Yes | Yes |
| NFS | Locking-dependent; no single-file write-order | Lock-free R/W; lock for vacuum | Write-once files | **No (SWMR)** | **No (WAL)** |
| Lustre | Undocumented as such | Documented flock mounts | Write-once files | Tested | Same as NFS |
| Object store | gcs, s3 | S3, GCS, Azure | URI / FileSystem | S3 VFD read-only | No |

## What dies, cleanly

- **A single growing Parquet or Arrow-IPC-file.** Footer-after-data. Fails Window durability and mid-run prefix.
- **HDF5 SWMR + rhdf5 as the mid-run R reader.** rhdf5 does not support SWMR.
- **HDF5 SWMR on NFS.** Officially unsupported.
- **HDF5 SWMR + native vlen ragged streams.** Officially forbidden.
- **SQLite WAL on NFS/Lustre or across hosts.** Officially unsupported.
- **TensorStore OCDBT (or zip) under Zarr that Rarr must open.** Rarr is a Zarr filesystem/S3 reader, not an OCDBT client.
- **Zarr without sharding or coarse chunks, chunked at one Sampling iteration.** File-count constraint.

## Biggest risks among survivors

1. **File count vs Window size.** Zarr shards and TileDB consolidation exist to keep object counts down. Parquet-per-Window has no equivalent except “make Windows large.”
2. **Ragged encoding.** TensorStore and HDF5 SWMR both push ragged streams into dense offsets+payload. TileDB’s R reader does not document numeric vlen-list reads. Arrow lists are the only named-package native ragged type.
3. **Mid-run R vs mid-run any-process.** HDF5’s library story is excellent and rhdf5’s is not. If “another process” may be C++ or `h5watch` and R only reads after the Run, HDF5 re-enters. The map says analysis is in R and mid-run readers see closed Windows; that pairing is the HDF5 kill.
4. **Shared filesystem.** If Runs actually write on NFS, HDF5 SWMR and SQLite are out; Zarr-on-files and TileDB fragments remain, with locking caveats for shard RMW and vacuum.

## Sources

- [Zarr v2 storage specification](https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html)
- [Zarr v3 core specification](https://zarr-specs.readthedocs.io/en/latest/v3/core/index.html)
- [Zarr v3 sharding_indexed codec](https://zarr-specs.readthedocs.io/en/latest/v3/codecs/sharding-indexed/)
- [Zarr Python data types (v2 object codecs)](https://zarr.readthedocs.io/en/3.1.x/user-guide/data_types.html)
- [TensorStore](https://google.github.io/tensorstore/)
- [TensorStore zarr2 driver](https://google.github.io/tensorstore/driver/zarr/index.html)
- [TensorStore zarr3 driver](https://google.github.io/tensorstore/driver/zarr3/index.html)
- [TensorStore file kvstore](https://google.github.io/tensorstore/kvstore/file/index.html)
- [TensorStore kvstore index](https://google.github.io/tensorstore/kvstore/index.html)
- [TensorStore ocdbt kvstore](https://google.github.io/tensorstore/kvstore/ocdbt/index.html)
- [tensorstore::Open](https://google.github.io/tensorstore/cpp/api/tensorstore.Open.html)
- [tensorstore::Resize](https://google.github.io/tensorstore/cpp/api/tensorstore.Resize.html)
- [tensorstore.Transaction](https://google.github.io/tensorstore/python/api/tensorstore.Transaction.html)
- [tensorstore::TransactionMode](https://google.github.io/tensorstore/cpp/api/tensorstore.TransactionMode.html)
- [Rarr supported features](https://huber-group-embl.github.io/Rarr/articles/features.html)
- [Rarr GitHub](https://github.com/Huber-group-EMBL/Rarr)
- [z5](https://github.com/constantinpape/z5)
- [xtensor-zarr](https://xtensor-zarr.readthedocs.io/en/latest/)
- [HDF5 SWMR technical note](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)
- [HDF5 SWMR User’s Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)
- [HDF5 file locking](https://support.hdfgroup.org/documentation/hdf5/latest/_file_lock.html)
- [HighFive File](https://highfive-devs.github.io/highfive/class_high_five_1_1_file.html)
- [rhdf5 H5Fopen](https://huber-group-embl.github.io/rhdf5/reference/H5Fopen.html)
- [rhdf5 h5createDataset](https://huber-group-embl.github.io/rhdf5/reference/h5createDataset.html)
- [hdf5r CRAN](https://cran.r-project.org/package=hdf5r)
- [TileDB fragments](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/fragments/)
- [TileDB commits](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/commits/)
- [TileDB concurrency](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/concurrency/)
- [TileDB ACID](https://documentation.cloud.tiledb.com/academy/structure/tables/foundation/key-concepts/acid/)
- [TileDB variable-length attributes](https://documentation.cloud.tiledb.com/academy/structure/arrays/tutorials/basics/var-length-attributes/)
- [TileDB C++ API](https://tiledb-inc-tiledb.readthedocs-hosted.com/en/stable/c++-api.html)
- [CRAN tiledb](https://cran.r-project.org/web/packages/tiledb/index.html)
- [Parquet file format](https://parquet.apache.org/docs/file-format/)
- [parquet-format](https://github.com/apache/parquet-format)
- [Arrow C++ Parquet](https://arrow.apache.org/docs/dev/cpp/parquet.html)
- [Arrow columnar / IPC formats](https://arrow.apache.org/docs/format/Columnar.html)
- [Arrow C++ IPC](https://arrow.apache.org/docs/dev/cpp/ipc.html)
- [arrow R open_dataset](https://arrow.apache.org/docs/r/reference/open_dataset.html)
- [arrow R datasets article](https://arrow.apache.org/docs/r/articles/dataset.html)
- [CRAN arrow](https://cran.r-project.org/web/packages/arrow/index.html)
- [SQLite WAL](https://www.sqlite.org/wal.html)
- [SQLite over a network](https://www.sqlite.org/useovernet.html)
