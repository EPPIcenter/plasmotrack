"""A cohort of fit summaries (U7): `summarize` rows of many Analyses concatenated to one Parquet file, with their
per-edge rows in a sibling file, and read back for the cohort view.

A row file `X.parquet` (`fit_summary.SCHEMA`) has its per-edge rows, when written, in `X.edges.parquet`
(`fit_summary.EDGE_SCHEMA`); a cohort `C.parquet` likewise has `C.edges.parquet`. The cohort adds one column to the
70, `manifest`: the summarised Analysis's manifest path from the row file's metadata, relative to the cohort file's
directory (null when the row file does not name one), so the cohort view can open a fit's Analysis. R reads the cohort
like the RDS; the extra column is last.
"""

import glob
import os
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from pt_analysis.fit_summary import EDGE_SCHEMA, FIT_KEY, MANIFEST_KEY, SCHEMA

EDGES_SUFFIX = ".edges.parquet"
COHORT_SCHEMA = SCHEMA.append(pa.field("manifest", pa.string()))
BINS = 10


class CohortError(Exception):
    """Row files cannot be concatenated, or a cohort file cannot be read."""


def edges_path(path):
    """The per-edge file beside a row or cohort file: `X.parquet` -> `X.edges.parquet`."""
    path = Path(path)
    return path.with_name(path.name.removesuffix(".parquet") + EDGES_SUFFIX)


def row_files(sources):
    """The row files that `sources` name, each a directory (its `*.parquet`) or a file or glob, sorted, leaving out
    per-edge files."""
    found = set()
    for source in sources:
        source = str(source)
        paths = Path(source).glob("*.parquet") if Path(source).is_dir() else map(Path, glob.glob(source))
        found.update(p.resolve() for p in paths if not p.name.endswith(EDGES_SUFFIX))
    return sorted(found)


def with_manifest(path, out_dir):
    """A row file's table with the `manifest` column."""
    table = pq.read_table(path)
    if table.schema.names != SCHEMA.names:
        raise CohortError(f"{path}: not a fit-summary row file; its columns differ from the fit summary's")
    named = (table.schema.metadata or {}).get(MANIFEST_KEY)
    manifest = os.path.relpath(Path(named.decode()).resolve(), out_dir) if named else None
    table = table.replace_schema_metadata(None).cast(SCHEMA)
    return table.append_column("manifest", pa.array([manifest] * table.num_rows, pa.string()))


def concatenate(sources, out):
    """Concatenate the row files `sources` name into `out`, and their per-edge files into `edges_path(out)` when any
    has one. Returns the number of rows."""
    files = row_files(sources)
    if not files:
        raise CohortError(f"no fit-summary rows in {', '.join(map(str, sources))}")
    out = Path(out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tables = [with_manifest(f, out.parent) for f in files]
    pq.write_table(pa.concat_tables(tables), out)
    edge_files = [edges_path(f) for f in files if edges_path(f).is_file()]
    if edge_files:
        with pq.ParquetWriter(edges_path(out), EDGE_SCHEMA) as writer:
            for f in edge_files:
                writer.write_table(pq.read_table(f).cast(EDGE_SCHEMA))
    elif edges_path(out).exists():
        edges_path(out).unlink()
    return sum(t.num_rows for t in tables)


def kind(t):
    if pa.types.is_string(t):
        return "string"
    return "integer" if pa.types.is_integer(t) else "number"


def calibration(cohort_path, bins=BINS):
    """Per fit (cohort row index) and score bin of width 1/bins (the last closed), the per-edge rows' count, true
    count and summed score, as columns; None without a per-edge file. Edges join their fit by `FIT_KEY`."""
    edges = edges_path(cohort_path)
    if not edges.is_file():
        return None
    on = " AND ".join(f"e.{k} IS NOT DISTINCT FROM r.{k}" for k in FIT_KEY)
    connection = duckdb.connect()  # its own: the cohort is not an Analysis and has no Store
    rows = connection.execute(f"""
        WITH r AS (SELECT row_number() OVER () - 1 AS row, {', '.join(FIT_KEY)} FROM read_parquet($cohort)),
             e AS (SELECT *, least(floor(score * {bins}), {bins - 1})::INTEGER AS bin FROM read_parquet($edges))
        SELECT r.row, e.bin, count(*) AS n, count(*) FILTER (WHERE e.truth) AS n_true, sum(e.score) AS score_sum
        FROM e JOIN r ON {on} GROUP BY ALL ORDER BY r.row, e.bin""",
                              {"cohort": str(cohort_path), "edges": str(edges)}).fetchall()
    connection.close()
    names = ("row", "bin", "n", "n_true", "score_sum")
    return {"bins": bins, **{name: [r[i] for r in rows] for i, name in enumerate(names)}}


def load(path, bins=BINS):
    """The cohort view's payload: columns ({name, kind}: string, integer or number), the data by column, and the
    calibration counts per fit and bin."""
    path = Path(path)
    if not path.is_file():
        raise CohortError(f"no cohort file at {path}")
    try:
        table = pq.read_table(path)
    except (pa.ArrowInvalid, OSError) as e:
        raise CohortError(f"{path}: not a Parquet file ({e})") from None
    return {"file": str(path), "n": table.num_rows,
            "columns": [{"name": f.name, "kind": kind(f.type)} for f in table.schema],
            "data": table.to_pydict(), "calibration": calibration(path, bins)}


def manifests(path):
    """The cohort's `manifest` column, relative to its directory; Nones without the column."""
    table = pq.read_table(path, columns=["manifest"]) if "manifest" in pq.read_schema(path).names else None
    return table.column("manifest").to_pylist() if table is not None else [None] * pq.read_metadata(path).num_rows
