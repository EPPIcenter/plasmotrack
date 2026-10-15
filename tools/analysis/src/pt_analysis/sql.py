"""SQL over DuckDB relations without leaving names behind on the connection.

A query that joins relations held in Python names each one through a `Scope`: the relation is registered on the
connection under a fresh name (`pt_<n>`, unique in the process) for as long as the scope is open, and unregistered
when it closes. Queries never read Python locals by name (DuckDB's replacement scan) and never create views that
outlive them, so two queries, or a derive between them, cannot collide on a name.

    with Scope(store.connection) as s:
        rows = s.rows(f"SELECT ... FROM {s.name(edges)} e JOIN {s.name(latent)} l ON ...", {"i": infection})

A relation built on a scope's names (`s.sql(...)`) is lazy: read or write it before the scope closes.
"""

import itertools

_counter = itertools.count()


class Scope:
    """Relations registered on `connection` under fresh names until the scope closes."""

    def __init__(self, connection):
        self.connection = connection
        self._names = []

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        while self._names:
            self.connection.unregister(self._names.pop())

    def name(self, relation):
        """Register `relation` (a DuckDB relation or an Arrow table) under a fresh name and return the name."""
        name = f"pt_{next(_counter)}"
        self.connection.register(name, relation)
        self._names.append(name)
        return name

    def sql(self, query):
        """`query` as a lazy relation: read it before the scope closes."""
        return self.connection.sql(query)

    def rows(self, query, params=None):
        """Every row of `query`, run with `params`."""
        return self.connection.execute(query, params).fetchall()

    def row(self, query, params=None):
        """The first row of `query`, run with `params`."""
        return self.connection.execute(query, params).fetchone()


def literal(text):
    """`text` as a SQL string literal, for a relation's `filter`, which takes no parameters."""
    return "'" + str(text).replace("'", "''") + "'"


def leftover_names(connection):
    """The temporary views on `connection`: every registered name still open. For tests."""
    return [name for (name,) in connection.execute(
        "SELECT view_name FROM duckdb_views() WHERE temporary AND NOT internal").fetchall()]
