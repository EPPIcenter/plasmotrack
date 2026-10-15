"""Epi data: the CSV an Analysis manifest names under `epi`, one row per Infection, any columns.

The key is the column named `infection`, or the first column when none is. Keys are read as text, so `01` and `1`
are different ids. Every other column is typed by DuckDB's CSV inspection: `date` when it parses as dates or
timestamps, `numeric` when it parses as numbers, `categorical` otherwise (booleans included). A value is the CSV's
own text, except that a date is written as an ISO date; an empty cell has no value.

Keys are checked against the Infections of the Analysis's first Run (`structure/infections`), which are the input's:
`unmatched` lists keys that are no Infection, `missing` lists Infections without a row.
"""

from dataclasses import dataclass
from datetime import datetime

import duckdb
import pyarrow as pa

KEY = "infection"

NUMERIC_TYPES = ("TINYINT", "SMALLINT", "INTEGER", "BIGINT", "HUGEINT", "UTINYINT", "USMALLINT", "UINTEGER",
                 "UBIGINT", "UHUGEINT", "FLOAT", "DOUBLE")


class EpiError(Exception):
    """The epi CSV cannot be read, or a column it does not have was asked for."""


def column_type(duckdb_type):
    """"date", "numeric" or "categorical" for the type DuckDB's CSV inspection gave a column."""
    name = str(duckdb_type).upper()
    if name == "DATE" or name.startswith("TIMESTAMP"):
        return "date"
    if name in NUMERIC_TYPES or name.startswith("DECIMAL"):
        return "numeric"
    return "categorical"


@dataclass
class Epi:
    path: str
    key: str
    columns: dict[str, str]  # column -> "date", "numeric" or "categorical", in CSV order, the key left out
    infections: list[str]  # the Analysis's Infections, sorted
    unmatched: list[str]  # keys that are no Infection of the Analysis, sorted
    missing: list[str]  # Infections without a row, sorted
    typed: pa.Table  # the CSV as DuckDB typed it, the key as text
    text: pa.Table  # the CSV as text

    def check(self, column):
        if column not in self.columns:
            raise EpiError(f"the epi data {self.path} has no column '{column}'; its columns are "
                           f"{', '.join(self.columns)}")

    def _cells(self, column):
        """(key, value as text, typed value) for every row of the CSV."""
        keys = self.text[self.key].to_pylist()
        typed = self.typed[column].to_pylist()
        if self.columns[column] == "date":
            texts = [None if v is None else v.isoformat(sep=" ") if isinstance(v, datetime) else v.isoformat()
                     for v in typed]
        else:
            texts = [None if v is None or v.strip() == "" else v.strip() for v in self.text[column].to_pylist()]
        return zip(keys, texts, typed)

    def values(self, column):
        """{Infection: its value of `column` as text, None where it has none}, for every Infection of the Analysis."""
        self.check(column)
        by_id = {key: text for key, text, _ in self._cells(column)}
        return {infection: by_id.get(infection) for infection in self.infections}

    def stratum_values(self, column):
        """The distinct values of `column` among the Analysis's Infections, as text, in the column's own order
        (dates and numbers ascending, text alphabetically), then None when some Infection has no value."""
        self.check(column)
        known = set(self.infections)
        cells = {(text, typed) for key, text, typed in self._cells(column) if key in known and text is not None}
        ordered = [text for text, _ in sorted(cells, key=lambda c: (c[1], c[0]))]
        unique = list(dict.fromkeys(ordered))
        return unique + ([None] if None in self.values(column).values() else [])


def load_epi(analysis):
    """The Analysis's epi data, or None when its manifest names none."""
    if analysis.epi is None:
        return None
    path = analysis.epi
    if not path.is_file():
        raise EpiError(f"epi data not found: {path}")
    connection = duckdb.connect()
    try:
        text = connection.read_csv(str(path), all_varchar=True).to_arrow_table()
        names = text.column_names
        key = KEY if KEY in names else names[0]
        typed = connection.read_csv(str(path), dtype={key: "VARCHAR"}).to_arrow_table()
    except duckdb.Error as e:
        raise EpiError(f"cannot read the epi data {path}: {e}") from None
    columns = {f.name: column_type(_duckdb_type(f.type)) for f in typed.schema if f.name != key}

    ids = text[key].to_pylist()
    repeated = sorted({i for i in ids if ids.count(i) > 1})
    if repeated:
        raise EpiError(f"the epi data {path} has more than one row for "
                       f"{', '.join(repr(i) for i in repeated)}; it needs one row per Infection")
    store = next(iter(analysis.runs.values()))
    infections = sorted(i for (i,) in store.structure("infections").select("infection").fetchall())
    known, present = set(infections), set(ids)
    return Epi(path=str(path), key=key, columns=columns, infections=infections,
               unmatched=sorted(i for i in present if i is not None and i not in known),
               missing=[i for i in infections if i not in present], typed=typed, text=text)


def _duckdb_type(arrow_type):
    """The DuckDB type family of an Arrow column type, for `column_type`."""
    if pa.types.is_date(arrow_type):
        return "DATE"
    if pa.types.is_timestamp(arrow_type):
        return "TIMESTAMP"
    if pa.types.is_integer(arrow_type):
        return "BIGINT"
    if pa.types.is_floating(arrow_type) or pa.types.is_decimal(arrow_type):
        return "DOUBLE"
    return "VARCHAR"
