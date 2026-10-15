"""The network view's payload (U8): Infections as nodes and pooled Edge marginals as edges.

Nodes are the input's Infections, in input order, with their observation time and symptomatic flag from the input,
the pooled unobserved-source inclusion q, and the posterior infection time: the observation time minus the 97.5%,
50% and 2.5% quantiles of `infection_duration` over every Draw of every Chain. Edges are the pooled Edge marginals
from observed parents at or above a floor, with each Chain's mean and the lag (child's observation time minus the
parent's). Edges from the unobserved source are left out; q carries them.

Pooling is over Chain x Draw: a pooled mean is the sum over every Chain's Draws divided by the number of Draws in
all, and quantiles are over the pooled per-Draw values (R type 7).
"""

import gzip
import json

from pt_analysis.derive import derive_analysis, derived, table

DEFAULT_FLOOR = 0.01


def read_input(path):
    """The engine input JSON at `path`, gzipped or not."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return json.load(f)


def _views(analysis, name, relation_of):
    """Create one DuckDB view per Chain over `relation_of(store)` and return a UNION ALL query over them, with the
    Chain's position in the Analysis as column `chain`."""
    parts = []
    for k, store in enumerate(analysis.runs.values()):
        view = f"pt_network_{name}_{k}"
        relation_of(store).create_view(view, replace=True)
        parts.append(f"SELECT {k} AS chain, * FROM {view}")
    return " UNION ALL ".join(parts)


def network(analysis, floor=DEFAULT_FLOOR):
    """The network payload of `analysis`: {"chains", "floor", "nodes", "edges"}. See the module docstring."""
    derive_analysis(analysis)
    connection = next(iter(analysis.runs.values())).connection
    draws_query = _views(analysis, "infection_draws", lambda s: derived(s, "infection_draws"))
    q = dict(connection.sql(f"SELECT infection, avg(q) FROM ({draws_query}) GROUP BY infection").fetchall())
    durations_query = _views(analysis, "infections", lambda s: table(s, "infections"))
    durations = {infection: quantiles for infection, *quantiles in connection.sql(f"""
        SELECT infection, quantile_cont(infection_duration, 0.975), quantile_cont(infection_duration, 0.5),
               quantile_cont(infection_duration, 0.025)
        FROM ({durations_query}) GROUP BY infection""").fetchall()}

    nodes = []
    for node in read_input(analysis.input)["nodes"]:
        time = node["observation_time"]
        quantiles = durations.get(node["id"])
        interval = None if quantiles is None else dict(zip(("lower", "median", "upper"), (time - d for d in quantiles)))
        nodes.append({"id": node["id"], "observation_time": time, "symptomatic": node.get("symptomatic", True),
                      "q": q.get(node["id"]), "infection_time": interval})

    observation_time = {n["id"]: n["observation_time"] for n in nodes}
    return {"chains": list(analysis.runs), "floor": floor, "nodes": nodes,
            "edges": _edges(analysis, connection, floor, observation_time)}


def _edges(analysis, connection, floor, observation_time):
    """Edges from observed parents whose pooled marginal is at least `floor`, largest first."""
    draws = [store.draws() for store in analysis.runs.values()]
    if not sum(draws):
        return []
    edges_query = _views(analysis, "edges", lambda s: derived(s, "edges"))
    latent_query = _views(analysis, "latent", lambda s: s.structure("infections").select("latent_parent"))
    rows = connection.execute(f"""
        WITH totals AS (
            SELECT chain, child, parent, sum(prob) AS total FROM ({edges_query})
            WHERE parent NOT IN (SELECT latent_parent FROM ({latent_query}))
            GROUP BY ALL)
        SELECT child, parent, sum(total) / ? AS marginal, map(list(chain), list(total))
        FROM totals GROUP BY child, parent HAVING sum(total) / ? >= ?
        ORDER BY marginal DESC, child, parent""", [sum(draws), sum(draws), floor]).fetchall()
    edges = []
    for child, parent, marginal, totals in rows:
        chains = [totals.get(k, 0.0) / n if n else None for k, n in enumerate(draws)]
        lag = (observation_time[child] - observation_time[parent]
               if child in observation_time and parent in observation_time else None)
        edges.append({"child": child, "parent": parent, "marginal": marginal, "chains": chains, "lag": lag})
    return edges
