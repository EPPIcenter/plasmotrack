"""Edges' weights per Draw per Chain, compact: the edge table's sparklines (the React app's Edges workspace).

An edge's weight at a Draw is its `edges` weight there, 0 at a Draw without a row, as `pt_analysis.edge_stats`
counts it. With `points`, each Chain's series is cut into that many runs of adjacent Draws (as even as the Draws
allow) and each run is replaced by its mean; a Chain with no more Draws than `points` keeps every Draw.
"""

import pyarrow as pa

from pt_analysis.derive import derive_analysis, derived, table
from pt_analysis.sql import Scope


def binned(values, points):
    """`values` as `points` means of adjacent runs, or as they are when there are no more than `points`."""
    n = len(values)
    if points is None or n <= points:
        return values
    edges = [k * n // points for k in range(points + 1)]
    return [sum(values[a:b]) / (b - a) for a, b in zip(edges, edges[1:])]


def edge_series(analysis, keys, points=None):
    """{"chains": [label], "series": [{child, parent, chains: [[weight per Draw] per Chain]}]} for `keys`, a list of
    (child, parent), in their order."""
    derive_analysis(analysis)
    wanted = pa.table({"child": pa.array([c for c, _ in keys], pa.string()),
                       "parent": pa.array([p for _, p in keys], pa.string())})
    per_chain = []
    for store in analysis.runs.values():
        iterations = sorted(i for (i,) in table(store, "run_scalars").aggregate("iteration").fetchall())
        position = {it: k for k, it in enumerate(iterations)}
        dense = {key: [0.0] * len(iterations) for key in keys}
        if keys and iterations:
            with Scope(store.connection) as s:
                rows = s.rows(f"""
                    SELECT e.child, e.parent, e.iteration, e.prob FROM {s.name(derived(store, "edges"))} e
                    JOIN {s.name(wanted)} w USING (child, parent)""")
            for child, parent, iteration, prob in rows:
                if iteration in position:
                    dense[(child, parent)][position[iteration]] = prob
        per_chain.append(dense)
    return {"chains": list(analysis.runs),
            "series": [{"child": child, "parent": parent,
                        "chains": [binned(dense[(child, parent)], points) for dense in per_chain]}
                       for child, parent in keys]}
