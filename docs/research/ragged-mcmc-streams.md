# How do candidate stores hold ragged MCMC streams?

Research for [How do candidate stores hold ragged MCMC streams?](https://github.com/EPPIcenter/plasmotrack/issues/28), a child of [Find the way to a unified output spec](https://github.com/EPPIcenter/plasmotrack/issues/26).

Facts only. This note does not pick a store. It records how each named candidate can represent the non-rectangular streams a plasmotrack Run already writes, and what would force those streams into a side format.

Primary sources: official specs and first-party library documentation, plus this repo's writer and datatype source.

## Question

How do zarr, HDF5, TileDB, Arrow/Parquet (and any other clearly relevant store) represent the ragged and non-rectangular streams a Run already writes, well enough that a spec can choose a format without discovering this later?

Streams that are not a dense scalar-per-iteration:

- Parent-set distribution: many (set, probability) pairs per Sampling iteration; cardinality changes
- Latent genotype: per Infection × Locus, a variable-width Allele bitstring / set
- Latent parent genotype: same shape
- Ordering: a permutation of Infection IDs, length = number of Infections
- Allele frequency: a simplex per Locus, length = number of Alleles at that Locus

Today parent-set files and Ordering are text; genotypes are binary strings.

For each store: idiomatic representation, what breaks mid-run append / Window commit / native R read, and what would force a side format for the ragged streams.

## What a Run already writes

These encodings are the current C++ writers, not a proposed spec.

**Parent-set distribution.** `ParentSetDistLogger` writes one gzip CSV per Infection under `parent_sets/`, header `parent_set,prob,iter`. Each Sampling iteration emits one row per Parent set in `calcParentSetDist()`, plus an optional source row `{S}` ([`src/core/io/loggers/ParentSetDistLogger.h`](src/core/io/loggers/ParentSetDistLogger.h), [`src/impl/model/Model/ModelLogger.cpp`](src/impl/model/Model/ModelLogger.cpp)). A Parent set serializes as `{id;id;…}` ([`src/core/io/serialize.h`](src/core/io/serialize.h)). Cardinality is not fixed: it is the number of Parent sets with mass at that iteration, not a dense vector over all combinatorially possible sets. Parent-set size is bounded in the model (`MAX_PARENTS = 2`, `MAX_PARENT_SET_SIZE = 3` in [`src/impl/model/Model/config.h`](src/impl/model/Model/config.h)), but the *number of pairs written per iteration* still changes.

**Latent genotype and Latent parent genotype.** One gzip CSV per Infection × Locus under `parameters/genotypes/` and `parameters/latent_parents/` ([`src/impl/model/Model/StateLogger.cpp`](src/impl/model/Model/StateLogger.cpp)). The production genetics type is `AllelesBitSet<MAX_ALLELES>` with `MAX_ALLELES = 128` ([`src/impl/model/Model/config.h`](src/impl/model/Model/config.h)). `serialize()` writes `allelesStr()`: a `'0'`/`'1'` string of length `total_alleles_`, not the full 128-bit capacity ([`src/core/datatypes/Alleles.h`](src/core/datatypes/Alleles.h)). Width therefore varies **across Loci**, not across Sampling iterations of one Locus. Resume reads the last line back as a string (`hotloadString` in [`src/impl/model/Model/State.cpp`](src/impl/model/Model/State.cpp)).

**Ordering.** One gzip CSV `parameters/infection_order.csv.gz`. `ValueLogger` serializes `Ordering::value()`, a `std::vector<std::shared_ptr<Infection>>`, as comma-separated Infection IDs ([`src/core/io/loggers/ValueLogger.h`](src/core/io/loggers/ValueLogger.h), [`src/core/io/serialize.h`](src/core/io/serialize.h), [`src/core/parameters/Ordering.h`](src/core/parameters/Ordering.h), [`src/core/containers/Infection.h`](src/core/containers/Infection.h) `serialize()` returns `id_`). Length equals the number of Infections in the Run. Infections are parsed from the input at construction and are not added mid-Run ([`src/impl/model/Model/State.cpp`](src/impl/model/Model/State.cpp)).

**Allele frequency.** One gzip CSV per Locus under `parameters/allele_frequencies/`. `Simplex::serialize()` writes `frequencies()` as a comma-separated vector of doubles ([`src/core/datatypes/Simplex.cpp`](src/core/datatypes/Simplex.cpp)). Internally the coefficients live in `std::array<double, 128>` with `total_elements_` = Alleles at that Locus ([`src/core/datatypes/Simplex.h`](src/core/datatypes/Simplex.h)). Width varies **across Loci**, not across Sampling iterations of one Locus. Loci and their Allele sets come from the input JSON ([`src/impl/model/Model/State.cpp`](src/impl/model/Model/State.cpp)).

## Which of these are actually ragged

Two different non-rectangularities are easy to conflate.

| Stream | Varies across Sampling iterations? | Varies across entities in one Run? | Fixed bound in this codebase? |
| --- | --- | --- | --- |
| Parent-set distribution | Yes — pair count changes | Yes — per Infection | Pair count unbounded by a small constant; each set has size ≤ `MAX_PARENTS` |
| Latent genotype / Latent parent genotype | No — width is the Locus Allele count | Yes — per Locus | `MAX_ALLELES = 128` bits |
| Ordering | No — length is Infection count | No | Infection count is fixed for the Run |
| Allele frequency | No — width is the Locus Allele count | Yes — per Locus | 128 coefficients allocated; live width is the Locus Allele count |

So:

- **Ordering** is a dense 2-D Trace: Sampling iteration × Infection rank. The current text line is a serialization choice, not a ragged type.
- **Allele frequency** is a dense 2-D Trace per Locus: Sampling iteration × Allele. Across Loci the second axis differs. That is a *collection of rectangles*, not a per-iteration ragged list.
- **Latent genotype** (and Latent parent genotype) is a dense Trace per Infection × Locus of a bitset whose live width is the Locus Allele count. Across Loci the width differs. The in-memory type already pads to 128 bits.
- **Parent-set distribution** is the only stream whose *cell* is a variable-cardinality list of records. It is the stream that actually exercises native ragged types.

A spec can still choose a ragged encoding for the per-Locus streams (one array with a variable-width cell per Locus). It does not have to. Padding to `MAX_ALLELES`, or one rectangular array per Locus, removes raggedness from every stream except the Parent-set distribution.

## Shared representation patterns

Every candidate can hold these streams without a side format if the spec picks one of the rectangular encodings. Native ragged types are optional convenience, not a requirement.

1. **Long table (fully rectangular columns).** One row per atom: `(sampling_iteration, infection_id, parent_set_id…, probability)` for Parent-set; `(sampling_iteration, infection_rank, infection_id)` for Ordering; `(sampling_iteration, allele_index, frequency)` for Allele frequency; `(sampling_iteration, infection_id, locus_id, bit_i)` or a packed integer for genotypes. Works in every store that can hold numeric/string columns. Mid-run append is ordinary row/chunk/fragment growth. Native R readers that handle dense numerics can read it.

2. **Offset + values (CSR / two arrays).** One dense values array plus a dense offsets (or length) array of length `n_iterations + 1`. This is the structure the Zarr v3 core spec names as the likely random-access encoding for variable-length elements, without committing to it as a core type ([Zarr v3 core](https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html)). TileDB stores variable-length attributes this way on disk ([TileDB attributes](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/attributes/)). Arrow List is this layout in memory ([Arrow Columnar Format](https://arrow.apache.org/docs/format/Columnar.html)).

3. **Native variable-length cell.** One cell per Sampling iteration holding a list / VLEN / var-length attribute. First-class in HDF5 (`H5Tvlen_create`), TileDB (variable-length attributes), Arrow (`List` / `LargeList` / `Map`), NetCDF-4 (`nc_def_vlen`). Not a Zarr v3 *core* type ([Zarr v3 core](https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html)).

4. **Pad to a known max.** Genotypes to 128 bits; Allele frequency to 128 doubles; Ordering is already fixed-width. Parent-set padding is a poor fit: the pair count is the number of sets with mass, not `MAX_PARENTS`.

5. **Keep today's text.** Store the existing `{id;id}`, bitstring, or comma-separated line as a string cell. This is a side *encoding*, not a side *file format*, if the store's string type is readable by the native R package.

What forces a **side format** (a second file type beside the chosen store) is choosing a native ragged type that the C++ writer cannot emit, or that the named R reader cannot decode, *and* refusing the rectangular fallbacks above.

## Zarr (including tensorstore)

### Type system

Zarr v3 core data types are fixed-size only: bool, integers, floats, complex, and optional raw `r*` bit blobs that are a multiple of 8 bits ([Zarr v3 data types](https://zarr-specs.readthedocs.io/en/latest/v3/data-types/index.html), [Zarr v3 core](https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html)). The spec says additional types MAY be extensions with variable-sized elements, and explicitly does not standardize a variable-length layout: “the most likely contender is to have two arrays, one with the actual variable length data and one with fixed size (pointer + length)” ([Zarr v3 core](https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html), citing [zarr-specs#62](https://github.com/zarr-developers/zarr-specs/issues/62)).

Zarr-Python implements variable-length UTF-8 and bytes as `VariableLengthUTF8` / `VariableLengthBytes`, using object-codec id `vlen-utf8` / `vlen-bytes`. The implementation notes that variable-length UTF-8 “does not have a Zarr V3 specification” ([zarr-python data types](https://zarr.readthedocs.io/en/3.1.x/user-guide/data_types.html), [zarr.core.dtype.npy.string](https://zarr.readthedocs.io/en/v3.1.2/_modules/zarr/core/dtype/npy/string.html)). In v2, object arrays are disambiguated by the filter/codec: `"vlen-utf8"` vs `"vlen-bytes"` ([zarr-python data types](https://zarr.readthedocs.io/en/3.1.x/user-guide/data_types.html)).

Numcodecs, the first-party codec library, also ships `VLenArray(dtype)` (`codec_id = 'vlen-array'`) for a NumPy object array of 1-D numeric arrays packed “into a parquet-style byte array” ([numcodecs vlen](https://numcodecs.readthedocs.io/en/stable/other/vlen.html)). That is the closest zarr-native encoding of a jagged numeric list (Ordering IDs, Allele-frequency vectors, Parent-set ID lists). It is a Python/numcodecs convention, not a Zarr core type.

### Idiomatic representation of these streams

| Stream | Idiomatic zarr (portable, tensorstore-writable) | Zarr-Python-only alternative |
| --- | --- | --- |
| Parent-set distribution | Long table, or two arrays (offsets + compound/struct values) | Object array + `vlen-array` / `vlen-utf8` of today's text |
| Latent genotype | `uint8`/`r16` (16 bytes) or `uint64[2]` padded to 128 bits; or one `uint8` array per Locus of width = Allele count | `vlen-utf8` / `vlen-bytes` of the `'01…'` string |
| Ordering | 2-D `uint32`/`int32` of Infection ranks or dictionary-encoded IDs, shape `[iteration, n_infections]` | `vlen-array` of ID vectors, or `vlen-utf8` of the CSV line |
| Allele frequency | One 2-D `float64` per Locus, or 3-D padded to 128 | `vlen-array` of per-Locus vectors |

Tensorstore's zarr2 driver supports scalar and structured (multi-field) dtypes and can resize arrays. It states “Currently, filters are not supported” ([tensorstore zarr2 driver](https://google.github.io/tensorstore/driver/zarr/index.html)). V2 object/vlen codecs are filters. Tensorstore's zarr3 driver lists codecs `bytes`, `transpose`, `sharding_indexed`, `gzip`, `blosc`, `zstd`, `crc32c` — not `vlen-utf8` / `vlen-bytes` / `vlen-array` ([tensorstore zarr3 driver](https://google.github.io/tensorstore/driver/zarr3/index.html)). Tensorstore *has* `string` / `ustring` dtypes in its own type system ([tensorstore.string](https://google.github.io/tensorstore/python/api/tensorstore.string.html)); they do not become a portable zarr vlen array through these drivers.

A C++ tensorstore writer can therefore hold every stream as dense numeric (or structured) arrays, including an offsets+values pair. It cannot idiomatically emit zarr object/vlen chunks.

### Mid-run append / Window commit

Tensorstore `resize` persistently changes stored bounds; expanding fills new elements with the fill value ([tensorstore.TensorStore.resize](https://google.github.io/tensorstore/python/api/tensorstore.TensorStore.resize.html), [tensorstore::Resize](https://google.github.io/tensorstore/cpp/api/tensorstore.Resize.html)). A Window commit is: write the new chunk(s) along the Sampling-iteration axis, then `resize` the published shape (or write chunks that already lie inside a pre-reserved shape and publish shape in metadata).

For an offsets+values pair, a Window must update **both** arrays and the published shapes together. Zarr metadata is a JSON document (`.zarray` / zarr.json). There is no single-file atomic commit in the core spec. A crash between the two resizes leaves a reader with mismatched offsets and values. That is the ragged-specific Window hazard for zarr. A long table has one array (or several independent columns) and avoids the pair.

Mid-run readers see whatever shape is in the current metadata. Tensorstore's `recheck_cached_metadata` defaults to `"open"` (not revalidated per read) ([tensorstore zarr2 driver](https://google.github.io/tensorstore/driver/zarr/index.html)); a concurrent R reader must re-read metadata to observe a newly closed Window.

### Native R read

Rarr (Bioconductor) reads zarr v2 and v3 numeric types, fixed-length strings, and `VLen-UTF8` (read and write; write of VLen-UTF8 only via `update_zarr_array()`). It does **not** list `vlen-array`, `vlen-bytes`, or `void *` as supported. Structured/struct types are read-only. Writing is limited to R-native `"integer"` / `"double"` / `"string"` ([Rarr supported features](https://huber-group-embl.github.io/Rarr/articles/features.html), [Rarr vignette](https://bioconductor.org/packages/release/bioc/vignettes/Rarr/inst/doc/Rarr.html)).

So Rarr can read:

- padded / per-Locus rectangular numeric arrays
- an offsets+values pair of numerics
- today's text if stored as `vlen-utf8` (or fixed-length string)

Rarr cannot be assumed to read a `vlen-array` object array of integer lists. That encoding would force either a Python decode step (ruled out by the map) or a side format.

### What would force a side format

- Specifying zarr object arrays with `vlen-array` (or pickle/JSON object codecs) as the Parent-set / Ordering encoding, while the writer is tensorstore and the reader is Rarr.
- Specifying a v3 extension dtype that neither tensorstore nor Rarr implements.
- Using two arrays for offsets+values *without* a documented commit order, if the spec also requires a crash to lose at most the in-flight Window **and** a mid-run reader to see a consistent prefix — unless Window commit is defined as “both arrays plus metadata”, which is application protocol, not a zarr feature.

Keeping genotypes as `vlen-utf8` bitstrings does **not** force a side format for Rarr; it does force a non-tensorstore writer for those arrays, or a second C++ stack.

## HDF5

### Type system

HDF5 datatypes include atomic, compound, array, enumeration, bitfield, and variable-length. `H5Tvlen_create(base_id)` creates a one-dimensional variable-length sequence of any base type, including compound and nested VLEN. This is **not** how variable-length strings are created; those use `H5T_VARIABLE` on a string type ([HDF5 Variable-length Sequence Datatypes](https://docs.hdfgroup.org/documentation/hdf5/latest/group___v_l_e_n.html), [HDF5 Datatypes](https://portal.hdfgroup.org/documentation/hdf5/latest/_h5_t__u_g.html), [HDF5 Datatype Basics](https://support.hdfgroup.org/documentation/hdf5/latest/_l_b_datatypes.html)).

A VLEN element is an `hvl_t` `{size, p}` in memory. On disk the dataset cells hold pointers into a global heap; the payload is not in the chunk ([HDF5 VLTypes technical note](https://docs.hdfgroup.org/archive/support/HDF5/doc/TechNotes/VLTypes.html), [HDF5 Datatype Basics](https://support.hdfgroup.org/releases/hdf5/v2_0/v2_0_0/documentation/doxygen/_l_b_datatypes.html)).

Compound types can nest VLEN members ([HDF5 Datatypes](https://portal.hdfgroup.org/documentation/hdf5/latest/_h5_t__u_g.html)). That is the native encoding of a list of `{parent_set, probability}` records.

### Idiomatic representation of these streams

| Stream | Idiomatic HDF5 |
| --- | --- |
| Parent-set distribution | Dataset of VLEN compound `{set: VLEN int or fixed int[2], prob: float64}`, one cell per (Infection, Sampling iteration); **or** a 1-D extendible table of compound rows `(iter, infection, parent0, parent1, is_source, prob)` |
| Latent genotype | `H5T_BITFIELD` or opaque/uint8[16] padded to 128 bits; or fixed-length string of length = Allele count per Locus dataset |
| Ordering | 2-D chunked dataset `[unlimited, n_infections]` of integer IDs or enum |
| Allele frequency | 2-D `[unlimited, n_alleles]` `float64` per Locus |

The long-table / per-Locus rectangle encodings use only atomic types and unlimited dimensions. They avoid every VLEN restriction below.

### Mid-run append / Window commit

Chunked datasets with an unlimited dimension grow via `H5Dset_extent` / append along that axis. SWMR is the official mid-run read story: the writer appends along an unlimited dimension and flushes; readers poll `H5Dget_space` / `H5Drefresh` ([HDF5 SWMR introduction](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)).

SWMR **forbids** the writer from modifying or appending data with variable-length, string, or region-reference datatypes ([HDF5 SWMR introduction](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html), [HDF5 SWMR User's Guide](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)). The heap that holds VLEN payloads is not SWMR-safe. If a Window must be visible to another process while the Run is still writing, **HDF5 VLEN cannot be the Parent-set (or text Ordering / genotype) cell type**.

Further VLEN costs, official:

- Chunk filters, including compression, operate on the `(pointer, count)` cells, not the heap payload. “Chunking and filters, including compression, are not available for heaps.” ([HDF5 Datatype Basics](https://support.hdfgroup.org/releases/hdf5/v2_0/v2_0_0/documentation/doxygen/_l_b_datatypes.html), [HDF5 VLTypes](https://docs.hdfgroup.org/archive/support/HDF5/doc/TechNotes/VLTypes.html)).
- SZIP cannot be applied to VLEN (or compound, array, enum) ([HDF5 Compressed Datasets](https://support.hdfgroup.org/documentation/hdf5/latest/_l_b_com_dset.html)).
- Parallel MPI file drivers refuse VLEN object creation ([HDF5 VLTypes](https://docs.hdfgroup.org/archive/support/HDF5/doc/TechNotes/VLTypes.html)).
- SWMR also forbids adding new objects after `H5Fstart_swmr_write`. All datasets must exist before the writer starts appending ([HDF5 SWMR introduction](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)). Per-Infection groups created mid-Run would break SWMR; they must be created up front.

A Window commit for rectangular datasets is `H5Dset_extent` + `H5Dflush` / `H5Fflush`. For VLEN, there is no SWMR-safe equivalent.

### Native R read

rhdf5 reads and writes variable-length **strings** (`size = NULL` on character datasets). It documents that a VLEN-string dataset “will not be compressed by HDF5” whereas fixed-length strings can be ([rhdf5 manual](https://bioconductor.statistik.uni-dortmund.de/packages/3.24/bioc/manuals/rhdf5/man/rhdf5.pdf)). `h5read(..., compoundAsDataFrame = TRUE)` coerces a compound dataset to a data.frame if the dataset is not multi-dimensional ([rhdf5 h5read](https://huber-group-embl.github.io/rhdf5/reference/h5_read.html)).

rhdf5's public high-level API is documented around arrays, strings, and compounds. A nested `H5T_VLEN` of compound (list of `{set, prob}` per cell) is not a documented high-level target. A 1-D table of compound rows *is*. If the spec needs mid-run R reads of closed Windows, SWMR + VLEN is already disallowed on the writer side, so the R reader never sees a growing VLEN dataset under SWMR.

### What would force a side format

- Using `H5T_VLEN` (or VLEN strings) for Parent-set / Ordering / genotypes **and** requiring SWMR mid-run reads of those datasets.
- Using nested VLEN-of-compound as the only Parent-set encoding if rhdf5 cannot round-trip it to an R list/data.frame without a custom low-level walk — that is a reader-capability risk, not a file-format impossibility.
- Creating per-Infection datasets after SWMR starts.

A compound long table, or per-Locus rectangular `float64` / bitfield datasets with one unlimited dimension, does not force a side format and is SWMR-legal.

## TileDB

### Type system

A TileDB attribute is either fixed-length (one value, or a schema-fixed number of values, per cell) or variable-length (a list of values of one base type per cell). Variable-length attributes are stored as two files per fragment: serialized values, and starting byte offsets. Example: cells `"a"`, `"bb"`, `"ccc"` become values `abbccc` and offsets `0, 1, 3` ([TileDB attributes](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/attributes/)). Filters, including compression, can be set independently on the values file and the offsets file; positive-delta encoding is documented as particularly useful for offsets ([TileDB tile filters](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/tile-filters/)).

C++ `Attribute::create<std::string>` and `Attribute::create<std::vector<int>>` are the typed constructors for var-length string and var-length numeric lists ([TileDB C++ variable-length example](https://github.com/TileDB-Inc/TileDB/blob/1.6.3/examples/cpp_api/variable_length.cc), current conceptual docs at [TileDB attributes](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/attributes/)).

There is no first-class nested list-of-structs attribute. A Parent-set cell that is a list of `{set, prob}` pairs is two var-length attributes on the same cell (e.g. packed set IDs + packed probabilities) or a sparse/long-table layout.

### Idiomatic representation of these streams

| Stream | Idiomatic TileDB |
| --- | --- |
| Parent-set distribution | Sparse array with dimensions `(sampling_iteration, infection, pair_index)` and attributes `parent0`, `parent1`, `is_source`, `prob`; **or** dense `(iteration, infection)` with two var-length attributes (IDs, probs) |
| Latent genotype | Fixed-length `TILEDB_UINT8` with `ncells = 16`, or `TILEDB_BLOB`; or var-length `UINT8` of live width |
| Ordering | Dense `(iteration, rank)` of `INT32` IDs; or one var-length `INT32` attribute per iteration |
| Allele frequency | Dense `(iteration, allele)` per Locus, or var-length `FLOAT64` per Locus cell |

### Mid-run append / Window commit

Each write creates one or more **fragments**: timestamped, immutable subfolders under `__fragments`. A fragment is considered by readers only after its commit file exists; TileDB documents fragment writes as atomic and consistent for concurrent readers and writers ([TileDB fragments](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/fragments/), [TileDB writes](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/writes/)). That is a natural Window: one fragment per closed Window. Numerous fragments hurt read performance; consolidation compacts them ([TileDB fragments](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/fragments/)).

Var-length attributes do not change this model: the values file and offsets file are part of the same fragment, so they commit together. That is the opposite of zarr's two independently resized arrays.

A dense write must supply a single-range subarray and all cell values in that subarray ([TileDB writes](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/writes/)). Growing a dense domain mid-Run is a schema/domain concern; sparse writes name coordinates explicitly and do not require a pre-sized iteration axis. For a Trace that grows every Window, sparse `(iteration, …)` or a dense domain with a large unused iteration extent are the two documented fits.

### Native R read

The `tiledb` R package (TileDB-R) constructs attributes with `ncells = NA` for variable length (`tiledb_attr(..., ncells = NA)`, `tiledb_attribute_set_cell_val_num(attr, NA)`). The official API vignette writes and reads var-length `CHAR` via offset+data buffers ([TileDB-R documentation vignette](https://tiledb-inc.github.io/TileDB-R/articles/documentation.html), [tiledb_attr](https://rdrr.io/cran/tiledb/man/tiledb_attr.html)). High-level `tiledb_array` as `data.frame` is the usual analysis path; var-length numeric attributes are part of the same C API the R package wraps.

This is a first-party R reader for the same engine the C++ writer uses. Ragged cells do not, by themselves, force a side format.

### What would force a side format

- Encoding Parent-set as a nested list-of-structs with no TileDB attribute mapping (TileDB has var-length *scalars of one type*, not Arrow-style `List<Struct>`).
- Relying on a dense domain that cannot grow, then overflowing the reserved Sampling-iteration extent.
- Choosing a datatype the R package does not expose (none of the streams require that if they stay on integer/float/ASCII/blob).

Keeping today's text as `TILEDB_STRING_ASCII` / `TILEDB_STRING_UTF8` is supported in both C++ and R ([TileDB attributes](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/attributes/)). That is a side *encoding*, not a side format.

## Arrow / Parquet

### Type system

The Arrow columnar format defines, among others ([Arrow Columnar Format](https://arrow.apache.org/docs/format/Columnar.html)):

- `List<T>` / `LargeList<T>`: validity + offsets + child array (32- or 64-bit offsets)
- `ListView` / `LargeListView`: offsets + explicit sizes (out-of-order child segments allowed)
- `FixedSizeList<T>[N]`
- `Utf8` / `LargeUtf8` / `Binary` / `LargeBinary`: variable-size binary
- `FixedSizeBinary[N]`
- `Struct`
- `Map`: “Variable-size List of Structs” (keys + items)

Parquet's `LIST` logical type must be a 3-level group (`optional/required group (LIST) { repeated group list { element } }`). `MAP` is the corresponding nested encoding for key/value ([Parquet Logical Types](https://parquet.apache.org/docs/file-format/types/logicaltypes/), [parquet-format LogicalTypes.md](https://github.com/apache/parquet-format/blob/master/LogicalTypes.md)).

This is the only named candidate whose *native cell* is a nested record list. A Parent-set distribution is idiomatically `List<Struct<parent_ids: List<int32> | FixedSizeList<int32>[2], is_source: bool, probability: float64>>` or `Map<parent_set_key, float64>` per (Infection, Sampling iteration).

### Idiomatic representation of these streams

| Stream | Idiomatic Arrow / Parquet |
| --- | --- |
| Parent-set distribution | Column `list_of(struct(parent0, parent1, is_source, prob))`, or a long table (one row per pair). `map_of` if the set is a single key |
| Latent genotype | `fixed_size_binary(16)` or `fixed_size_list_of(boolean(), 128)` or `uint8`[16]; or `binary`/`utf8` of the live bitstring |
| Ordering | `fixed_size_list_of(utf8() or int32(), n_infections)` — width known at Run start — or a long table |
| Allele frequency | `list_of(float64())` per Locus column, or `fixed_size_list_of(float64(), n_alleles)` per Locus, or a long table |

The Arrow R package exposes `list_of`, `large_list_of`, `fixed_size_list_of`, `map_of`, `struct`, `fixed_size_binary` ([Arrow R data-type](https://arrow.apache.org/docs/23.0/r/reference/data-type.html)). Nested columns survive Parquet round-trip through the same package (`read_parquet` / `write_parquet`) ([Arrow R read/write](https://arrow.apache.org/docs/r/articles/read_write.html)).

### Mid-run append / Window commit

A Parquet file is a closed unit of row groups + footer. The Dataset API writes **new** files (`basename_template = "part{i}.parquet"`) under a directory, optionally Hive-partitioned ([Arrow C++ Datasets](https://arrow.apache.org/docs/cpp/dataset.html)). `ExistingDataBehavior` overwrites or deletes matching partitions; it is not an in-file append ([Arrow C++ Datasets](https://arrow.apache.org/docs/cpp/dataset.html)).

A Window is therefore one (or more) new Parquet file(s) whose directory listing becomes visible after the file is fully written and the footer is present. Mid-run readers use `open_dataset()` on the directory and see completed files only, provided the writer publishes via atomic rename into the dataset directory (application protocol; Arrow writes to the destination stream it is given).

In-file append of row groups to an existing Parquet file exists on `parquet::arrow::FileWriter` but is not the Dataset write path and is a poor fit for a crash-safe Window (a torn footer loses the file). Arrow IPC / Feather can append record batches more naturally; it is the same type system, a different container.

List/Map columns do not change the commit story: they live inside the same record batches / row groups as the dense columns.

C++ write is first-party (`arrow::dataset::FileSystemDataset::Write`, `parquet::arrow::WriteTable`). R read is first-party (`arrow::open_dataset`, `read_parquet`).

### What would force a side format

- Choosing a nested depth or extension type the installed `arrow` R package cannot convert to an R vector (unusual for List/Struct/Map of primitives).
- Treating “one Parquet file per Run that is appended in place” as the durability model — that is a container problem, not a ragged-type problem. The ragged streams do not make this worse than dense Traces.
- Storing genotypes as Arrow `Boolean` bitmaps of live (variable) length *without* a List wrapper — that is not a valid Array. `List<Boolean>` or `FixedSizeList<Boolean>[128]` or `FixedSizeBinary[16]` are valid.

A long table of primitive columns is the conservative encoding and is universally readable. Nested List/Struct is the idiomatic one and is in-spec for both C++ and R.

## Other store that is clearly relevant: NetCDF-4

NetCDF-4 (HDF5-backed) adds user-defined types, including “variable length array type: used to store ragged arrays” via `nc_def_vlen` ([NetCDF user-defined types](https://docs.unidata.ucar.edu/netcdf-c/current/group__user__types.html), [NetCDF Users Guide: types](https://docs.unidata.ucar.edu/nug/current/md_types.html)). Compound + VLEN nesting is allowed. This is the same HDF5 heap VLEN under a different API, so the SWMR / filter / heap caveats above still apply.

The native R package `ncdf4` supports “compression, chunking, groups, vlen strings, and multiple unlimited dimensions.” It states: “User-defined types and vlen data arrays (variable-length arrays) are not supported yet.” ([ncdf4 reference](https://cran.r-project.org/web/packages/ncdf4/refman/ncdf4.html)). So NetCDF-4's first-class ragged type is **not** readable by the obvious native R reader. Rectangular unlimited dimensions are. Using NetCDF-4 VLEN for Parent-set would force a side format or a different R package (`rhdf5` on the underlying HDF5, with the SWMR/VLEN issues unchanged).

NetCDF-4 is relevant as “HDF5 VLEN with a worse R story,” not as a new ragged capability.

## Cross-store comparison (ragged-specific)

| Concern | Zarr + tensorstore | HDF5 | TileDB | Arrow / Parquet | NetCDF-4 |
| --- | --- | --- | --- | --- | --- |
| Native list-of-records cell | No (core). Two arrays or long table. Python `vlen-array` is non-portable | Yes: VLEN of compound | Var-length of one base type; nest by pairing attributes or long table | Yes: `List<Struct>` / `Map` | Yes: `nc_def_vlen` (+ compound) |
| C++ writer for native ragged | Not via tensorstore filters | HDF5 C/C++ | TileDB C++ | Arrow C++ | netCDF-C |
| Native R reader for native ragged | Rarr: `vlen-utf8` only, not `vlen-array` | rhdf5: VLEN **strings** documented; nested VLEN-of-compound not a high-level feature | tiledb: `ncells = NA` | arrow: `list_of` / `map_of` / `struct` | ncdf4: VLEN arrays **unsupported** |
| Window commit of ragged + dense together | Application: resize multiple arrays + metadata | `H5Dflush`; VLEN illegal under SWMR | One fragment (values + offsets atomic) | One new Parquet/IPC file | Same as HDF5 |
| Mid-run consistent prefix | Re-read zarr metadata; two-array race | SWMR on **fixed-size** types only | Readers ignore uncommitted fragments | Readers see finished files only | SWMR/VLEN as HDF5 |
| Compression of ragged payload | Chunk codecs on numeric arrays; object-codec payload is a separate path | Heap payload not filtered | Yes, on values and offsets separately | Yes, Parquet page/column compression | As HDF5 |
| What forces a side format | Object/`vlen-array` as the only encoding | VLEN + required SWMR | Nested List-of-Struct with no mapping | In-place single-file append as the Window model (container, not type) | VLEN arrays (ncdf4 cannot read them) |

## Implications for a later layout spec (not a pick)

These are constraints a layout ticket will hit, not a format decision.

1. **Only the Parent-set distribution is ragged in the iteration axis.** Ordering, Allele frequency, and genotypes are rectangular once the spec fixes an Infection index and a per-Locus Allele index (or pads to 128). A spec that says “the store must have native ragged types” is justifying that requirement with one stream.

2. **A long table of atoms fits every named candidate and every named R reader** (Rarr, rhdf5, arrow, tiledb, ncdf4). It is the encoding that cannot force a side format. Parent-set `{id;id}` text can become two integer columns plus a source flag; genotype `'01…'` can become 16 bytes or 128 bits; Ordering can become `n_infections` integer columns or a 2-D array.

3. **Native ragged types diverge on the constraint that matters for Windows:** mid-run readers of a consistent prefix. TileDB fragments and Arrow-dataset files commit offsets and values together. HDF5 VLEN cannot be appended under SWMR. Zarr two-array + metadata is an application-level transaction.

4. **tensorstore + Rarr, as a pair, push the spec toward numeric rectangles (or offsets+values), not zarr object arrays.** Rarr's VLen-UTF8 support would let the spec keep today's text *inside* zarr, but tensorstore cannot write those codecs.

5. **Padding genotypes to 128 bits and Allele frequencies to the Locus Allele count (or 128) is already how the C++ types are laid out.** `AllelesBitSet<128>` and `Simplex`'s `std::array<double, 128>` are fixed-capacity ([`config.h`](src/impl/model/Model/config.h), [`Simplex.h`](src/core/datatypes/Simplex.h)). The raggedness in the current files is the text serializer (`allelesStr()` of live width, comma-separated `frequencies()`).

6. **Exact encodings** (bit-endian of the 128-bit genotype, integer vs string Infection IDs, how `{S}` is distinguished from an empty Parent set) remain a layout-ticket concern, as the map already notes. This ticket only needs: no named candidate is unable to hold these streams as numeric rectangles; the candidates that offer native ragged cells disagree on C++ write, R read, and Window atomicity.

## Sources

### This repository

- [`src/impl/model/Model/config.h`](src/impl/model/Model/config.h) — `MAX_ALLELES`, `MAX_PARENTS`, `GeneticsImpl`
- [`src/impl/model/Model/StateLogger.cpp`](src/impl/model/Model/StateLogger.cpp) — genotype, Allele-frequency, Ordering writers
- [`src/impl/model/Model/ModelLogger.cpp`](src/impl/model/Model/ModelLogger.cpp) — Parent-set writer
- [`src/impl/model/Model/State.cpp`](src/impl/model/Model/State.cpp) — input-fixed Infections/Loci; resume via `hotloadString` / `hotloadVector`
- [`src/core/io/loggers/ParentSetDistLogger.h`](src/core/io/loggers/ParentSetDistLogger.h)
- [`src/core/io/loggers/ValueLogger.h`](src/core/io/loggers/ValueLogger.h)
- [`src/core/io/serialize.h`](src/core/io/serialize.h)
- [`src/core/datatypes/Alleles.h`](src/core/datatypes/Alleles.h)
- [`src/core/datatypes/Simplex.h`](src/core/datatypes/Simplex.h), [`src/core/datatypes/Simplex.cpp`](src/core/datatypes/Simplex.cpp)
- [`src/core/parameters/Ordering.h`](src/core/parameters/Ordering.h)
- [`src/core/containers/Infection.h`](src/core/containers/Infection.h)

### Zarr / tensorstore / numcodecs / Rarr

- [Zarr v3 core specification](https://zarr-specs.readthedocs.io/en/latest/v3/core/v3.0.html)
- [Zarr v3 data types](https://zarr-specs.readthedocs.io/en/latest/v3/data-types/index.html)
- [Zarr v2 specification](https://zarr-specs.readthedocs.io/en/latest/v2/v2.0.html)
- [zarr-python data types](https://zarr.readthedocs.io/en/3.1.x/user-guide/data_types.html)
- [numcodecs variable-length codecs](https://numcodecs.readthedocs.io/en/stable/other/vlen.html)
- [tensorstore zarr2 driver](https://google.github.io/tensorstore/driver/zarr/index.html)
- [tensorstore zarr3 driver](https://google.github.io/tensorstore/driver/zarr3/index.html)
- [tensorstore.TensorStore.resize](https://google.github.io/tensorstore/python/api/tensorstore.TensorStore.resize.html)
- [tensorstore::Resize](https://google.github.io/tensorstore/cpp/api/tensorstore.Resize.html)
- [tensorstore.string](https://google.github.io/tensorstore/python/api/tensorstore.string.html)
- [Rarr supported Zarr features](https://huber-group-embl.github.io/Rarr/articles/features.html)
- [Rarr vignette](https://bioconductor.org/packages/release/bioc/vignettes/Rarr/inst/doc/Rarr.html)

### HDF5 / rhdf5

- [HDF5 variable-length sequence datatypes](https://docs.hdfgroup.org/documentation/hdf5/latest/group___v_l_e_n.html)
- [HDF5 Datatypes (user guide)](https://portal.hdfgroup.org/documentation/hdf5/latest/_h5_t__u_g.html)
- [HDF5 Datatype Basics](https://support.hdfgroup.org/documentation/hdf5/latest/_l_b_datatypes.html)
- [HDF5 VLTypes technical note](https://docs.hdfgroup.org/archive/support/HDF5/doc/TechNotes/VLTypes.html)
- [HDF5 Compressed Datasets](https://support.hdfgroup.org/documentation/hdf5/latest/_l_b_com_dset.html)
- [HDF5 SWMR introduction](https://docs.hdfgroup.org/documentation/hdf5/latest/_s_w_m_r_t_n.html)
- [HDF5 SWMR User's Guide (PDF)](https://support.hdfgroup.org/releases/hdf5/documentation/features/SWMR/HDF5_SWMR_Users_Guide.pdf)
- [rhdf5 manual](https://bioconductor.statistik.uni-dortmund.de/packages/3.24/bioc/manuals/rhdf5/man/rhdf5.pdf)
- [rhdf5 h5read](https://huber-group-embl.github.io/rhdf5/reference/h5_read.html)

### TileDB / TileDB-R

- [TileDB attributes (variable-length)](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/attributes/)
- [TileDB fragments](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/fragments/)
- [TileDB writes](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/compute/writes/)
- [TileDB tile filters](https://documentation.cloud.tiledb.com/academy/structure/arrays/foundation/key-concepts/storage/tile-filters/)
- [TileDB-R API vignette](https://tiledb-inc.github.io/TileDB-R/articles/documentation.html)
- [tiledb_attr](https://rdrr.io/cran/tiledb/man/tiledb_attr.html)

### Arrow / Parquet

- [Arrow Columnar Format](https://arrow.apache.org/docs/format/Columnar.html)
- [Parquet Logical Types](https://parquet.apache.org/docs/file-format/types/logicaltypes/)
- [parquet-format LogicalTypes.md](https://github.com/apache/parquet-format/blob/master/LogicalTypes.md)
- [Arrow C++ Tabular Datasets](https://arrow.apache.org/docs/cpp/dataset.html)
- [Arrow R data-type](https://arrow.apache.org/docs/23.0/r/reference/data-type.html)
- [Arrow R read/write](https://arrow.apache.org/docs/r/articles/read_write.html)

### NetCDF-4 / ncdf4

- [NetCDF user-defined types (VLEN)](https://docs.unidata.ucar.edu/netcdf-c/current/group__user__types.html)
- [NetCDF Users Guide: types](https://docs.unidata.ucar.edu/nug/current/md_types.html)
- [ncdf4 reference manual](https://cran.r-project.org/web/packages/ncdf4/refman/ncdf4.html)
