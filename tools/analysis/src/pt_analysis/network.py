"""The network view's payload (U8): Infections as nodes and pooled Edge marginals as edges.

Nodes are the input's Infections, in input order, with their observation time, symptomatic flag and candidate parents
(`allowed_parents`, the Infections the input lets be a parent of it) from the input, the pooled unobserved-source inclusion q, and the posterior infection time: the observation time minus the 97.5%,
50% and 2.5% quantiles of `infection_duration` over every Draw of every Chain. Edges are the pooled Edge marginals
from observed parents at or above a floor, with each Chain's mean and the lag (child's observation time minus the
parent's). Edges from the unobserved source are left out; q carries them.

With epi data (`pt_analysis.epi`), `epi_columns` lists its columns with their type and stratum values (as
`Epi.stratum_values` orders them) and each node carries its value of every column as text, None where it has none;
without epi data the list is empty and every node's values are `{}`.

Pooling is over Chain x Draw: a pooled mean is the sum over every Chain's Draws divided by the number of Draws in
all, and quantiles are over the pooled per-Draw values (R type 7).
"""

from functools import reduce

from pt_analysis import inputs
from pt_analysis.derive import derive_analysis, derived
from pt_analysis.epi import load_epi
from pt_analysis.sql import Scope
from pt_analysis.times import infection_times, lag, observation_times

DEFAULT_FLOOR = 0.01


def _union(analysis, relation_of):
    """`relation_of(store)` of every Chain, with the Chain's position in the Analysis as column `chain`, as one
    relation (UNION ALL)."""
    return reduce(lambda a, b: a.union(b), (relation_of(store).project(f"{k} AS chain, *")
                                            for k, store in enumerate(analysis.runs.values())))


def network(analysis, floor=DEFAULT_FLOOR):
    """The network payload of `analysis`: {"chains", "floor", "nodes", "edges", "epi_columns"}. See the module
    docstring."""
    derive_analysis(analysis)
    connection = next(iter(analysis.runs.values())).connection
    q = dict(_union(analysis, lambda s: derived(s, "infection_draws")).aggregate("infection, avg(q)", "infection")
             .fetchall())
    intervals = infection_times(analysis)
    epi = load_epi(analysis)
    columns = list(epi.columns) if epi is not None else []
    values = {column: epi.values(column) for column in columns}
    candidates = inputs.candidate_parents(analysis.input)
    nodes = []
    for id, node in inputs.nodes(analysis.input).items():
        nodes.append({"id": id, "observation_time": node["observation_time"],
                      "symptomatic": node.get("symptomatic", True), "candidate_parents": candidates[id], "q": q.get(id),
                      "infection_time": intervals.get(id),
                      "epi": {column: values[column].get(id) for column in columns}})
    epi_columns = [{"name": c, "type": epi.columns[c], "values": epi.stratum_values(c)} for c in columns]
    return {"chains": list(analysis.runs), "floor": floor, "nodes": nodes, "epi_columns": epi_columns,
            "edges": _edges(analysis, connection, floor, observation_times(analysis))}


def _edges(analysis, connection, floor, times):
    """Edges from observed parents whose pooled marginal is at least `floor`, largest first."""
    draws = [store.draws() for store in analysis.runs.values()]
    if not sum(draws):
        return []
    with Scope(connection) as s:
        edges = s.name(_union(analysis, lambda st: derived(st, "edges")))
        latent = s.name(_union(analysis, lambda st: st.structure("infections").select("latent_parent")))
        rows = s.rows(f"""
            WITH totals AS (
                SELECT chain, child, parent, sum(prob) AS total FROM {edges}
                WHERE parent NOT IN (SELECT latent_parent FROM {latent})
                GROUP BY ALL)
            SELECT child, parent, sum(total) / $draws AS marginal, map(list(chain), list(total))
            FROM totals GROUP BY child, parent HAVING sum(total) / $draws >= $floor
            ORDER BY marginal DESC, child, parent""", {"draws": sum(draws), "floor": floor})
    return [{"child": child, "parent": parent, "marginal": marginal,
             "chains": [totals.get(k, 0.0) / n if n else None for k, n in enumerate(draws)],
             "lag": lag(times, parent, child)}
            for child, parent, marginal, totals in rows]
