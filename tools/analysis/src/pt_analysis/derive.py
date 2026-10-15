"""Derived tables: Parquet tables computed once from a Run's store and cached under `<run>/derived/`.

`derived/manifest.json` records the closed Windows per Phase the tables came from and the package version that wrote
them. `derive` recomputes every table when either differs from the store's closed Windows or this package's version,
and does nothing otherwise. Tables are computed from the sampling Phase, so every row comes from Draws;
`parent_set_members` is read over both Phases.

    set_classes   infection, set_id, n_observed, has_unobserved, class
    edges         iteration, child, parent, prob
    set_summary   infection, set_id, weight, rank, cumulative_mass
"""

# The SQL below names DuckDB relations held in local variables (`FROM parent_sets`); DuckDB finds them by name in
# the calling frame (a replacement scan), so those variables are used even though nothing else reads them.

import json
import os
from pathlib import Path

from pt_analysis import __version__
from pt_analysis.store import BOTH_PHASES

# Columns the derivations read, for a table with no rows in the closed Windows (a Run still in burn-in has no
# sampling `parent_sets`).
COLUMNS = {
    "parent_sets": {"iteration": "BIGINT", "infection": "VARCHAR", "set_id": "BIGINT", "prob": "DOUBLE"},
    "parent_set_members": {"infection": "VARCHAR", "set_id": "BIGINT", "member": "VARCHAR"},
}


def table(store, name):
    """`store.table(name)`, or an empty relation with the columns the derivations read."""
    relation = store.table(name)
    if relation is None:
        columns = ", ".join(f'NULL::{t} AS "{c}"' for c, t in COLUMNS[name].items())
        relation = store.connection.sql(f"SELECT {columns} WHERE false")
    return relation


# The composition class of a Parent set: `unobserved_only` when it has no observed parent (its only member is the
# Latent parent; there is never an empty set), else `observed_<n>`, suffixed `_unobserved` when the Latent parent is a
# member too.
CLASS_LABEL = """CASE WHEN n_observed = 0 THEN 'unobserved_only'
                     ELSE 'observed_' || n_observed || CASE WHEN has_unobserved THEN '_unobserved' ELSE '' END END"""


def set_classes(store):
    """Every Parent set of the Run with its observed-parent count, whether the unobserved source (the Infection's
    Latent parent) is a member, and its composition class."""
    members = table(store, "parent_set_members")
    infections = store.structure("infections")
    sets = store.connection.sql("""
        SELECT m.infection, m.set_id,
               count(*) FILTER (WHERE m.member <> i.latent_parent) AS n_observed,
               bool_or(m.member = i.latent_parent) AS has_unobserved
        FROM members m JOIN infections i USING (infection)
        GROUP BY ALL""")
    return store.connection.sql(f'SELECT *, {CLASS_LABEL} AS "class" FROM sets ORDER BY infection, set_id')


def edges(store):
    """Per Draw, each entry of the conditional transmission matrix: the summed weight of the child's Parent sets that
    contain the parent. The unobserved source is the child's Latent parent id. A parent with no row at a Draw has
    weight 0 there: the store drops sets below the Run's floor."""
    parent_sets = table(store, "parent_sets")
    members = table(store, "parent_set_members")
    return store.connection.sql("""
        SELECT p.iteration, p.infection AS child, m.member AS parent, sum(p.prob) AS prob
        FROM parent_sets p JOIN members m USING (infection, set_id)
        GROUP BY ALL ORDER BY iteration, child, parent""")


def set_summary(store):
    """Each Infection's Parent sets by pooled weight: the sum of the set's weight over Draws divided by the number of
    Draws (a set with no row at a Draw has weight 0 there), ranked from 1 by weight descending then set_id, with the
    cumulative weight down the ranking. Sets with no weight at any Draw are left out."""
    parent_sets = table(store, "parent_sets")
    draws = store.draws()
    return store.connection.sql(f"""
        WITH pooled AS (SELECT infection, set_id, sum(prob) / {draws} AS weight FROM parent_sets GROUP BY ALL),
        ranked AS (SELECT *, row_number() OVER (PARTITION BY infection ORDER BY weight DESC, set_id) AS rank
                   FROM pooled)
        SELECT *, sum(weight) OVER (PARTITION BY infection ORDER BY rank) AS cumulative_mass
        FROM ranked ORDER BY infection, rank""")


TABLES = {"set_classes": set_classes, "edges": edges, "set_summary": set_summary}


def source(store):
    """What the derived tables of `store` are computed from, as `derived/manifest.json` records it."""
    windows = {phase: [w.number for w in store.windows((phase,))] for phase in BOTH_PHASES}
    return {"windows": windows, "package_version": __version__}


def derived_dir(store):
    return Path(store.path) / "derived"


def is_current(store):
    """Whether `derived/` holds every table, computed from the store's closed Windows by this package version."""
    out = derived_dir(store)
    try:
        recorded = json.loads((out / "manifest.json").read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return False
    return recorded == source(store) and all((out / f"{name}.parquet").is_file() for name in TABLES)


def derive(store):
    """Write the derived tables of `store` under `<run>/derived/` unless they are current. True if they were written.
    Each table is written beside its final name and renamed into place, and the manifest is written last, so an
    interrupted derive leaves no manifest and is redone."""
    if is_current(store):
        return False
    out = derived_dir(store)
    out.mkdir(exist_ok=True)
    manifest = out / "manifest.json"
    manifest.unlink(missing_ok=True)
    for name, compute in TABLES.items():
        partial = out / f".{name}.parquet"
        compute(store).write_parquet(str(partial))
        os.replace(partial, out / f"{name}.parquet")
    partial = out / ".manifest.json"
    partial.write_text(json.dumps(source(store), indent=2) + "\n")
    os.replace(partial, manifest)
    return True


def derived(store, name):
    """Derived table `name` of `store` as a DuckDB relation, deriving the tables first if they are not current."""
    derive(store)
    return store.connection.read_parquet(str(derived_dir(store) / f"{name}.parquet"))
