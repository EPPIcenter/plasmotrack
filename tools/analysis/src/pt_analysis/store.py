"""One Run's store, read under the closed-Window rule (doc/store-layout.md).

A Window exists for a reader only if its marker `_closed/window-NNNNNN` exists. Table files of a Window without a
marker (left by a Run killed mid-commit) and in-flight files are never read. Each table is exposed as a DuckDB
relation over its Parquet files in place.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import duckdb

SAMPLING = ("sampling",)
BOTH_PHASES = ("burn_in", "sampling")

# Tables read over both Phases unless the caller says otherwise. `parent_set_members` is a lookup table: a set's
# members are written in the Window where the set first appears, which may be in burn-in, and every later Draw's
# `parent_sets` rows refer to them by set_id.
DEFAULT_PHASES = {"parent_set_members": BOTH_PHASES}


@dataclass(frozen=True, order=True)
class Window:
    """A closed Window: its number, the Sampling iterations it covers (inclusive) and the Phase of every row in it."""

    number: int
    first_iteration: int
    last_iteration: int
    phase: str


class Store:
    """A Run's store directory."""

    def __init__(self, path, connection=None):
        self.path = Path(path)
        self.connection = connection or duckdb.default_connection()

    def closed_windows(self):
        """The Windows a reader may use, in Window order: those whose marker exists."""
        closed = self.path / "_closed"
        if not closed.is_dir():
            return []
        windows = []
        for marker in closed.glob("window-*"):
            m = json.loads(marker.read_text())
            windows.append(Window(m["window"], m["first_iteration"], m["last_iteration"], m["phase"]))
        return sorted(windows)

    def windows(self, phases=SAMPLING):
        """The closed Windows whose Phase is one of `phases`."""
        return [w for w in self.closed_windows() if w.phase in phases]

    def table(self, name, phases=None):
        """Table `name` over the closed Windows of `phases`, as a DuckDB relation over its Parquet files in place.
        Traces and Blocks default to the sampling Phase, so every row is a Draw; `parent_set_members` defaults to both
        Phases. None if the table has no rows in those Windows: a table with no rows in a Window has no file for it."""
        if phases is None:
            phases = DEFAULT_PHASES.get(name, SAMPLING)
        files = [self.path / name / f"window-{w.number:06d}.parquet" for w in self.windows(phases)]
        files = [str(f) for f in files if f.is_file()]
        if not files:
            return None
        return self.connection.read_parquet(files)

    def rows(self, name, phases=None):
        """`table(name, phases)` as a list of dicts, in Sampling-iteration order when the table has one; [] when the
        table has no rows in those Windows. For code that walks rows in Python."""
        relation = self.table(name, phases)
        if relation is None:
            return []
        if "iteration" in relation.columns:
            relation = relation.order("iteration")
        return relation.to_arrow_table().to_pylist()

    def structure(self, name):
        """One of the lookup tables written at Run start: `infections`, `loci` or `samplers`."""
        return self.connection.read_parquet(str(self.path / "structure" / f"{name}.parquet"))

    @property
    def manifest(self):
        """The Run manifest."""
        return json.loads((self.path / "manifest.json").read_text())

    def draws(self):
        """The number of Draws: written sampling-Phase Sampling iterations, counted on `run_scalars`."""
        scalars = self.table("run_scalars")
        if scalars is None:
            return 0
        return scalars.aggregate("count(distinct iteration)").fetchone()[0]

    def sampling_time(self):
        """Wall time and Sampling iterations of the sampling Phase, summed over its `window_closed` events.
        `sampling_seconds` is None when no sampling Window has closed."""
        closes = [json.loads(e["payload"]) for e in self.rows("events") if e["type"] == "window_closed"]
        iterations = sum(c["sampling_iterations"] for c in closes)
        seconds = sum(c["seconds"] for c in closes) if iterations else None
        return {"sampling_seconds": seconds, "samples": iterations}  # the keys mixing_benchmark's timings.json uses
