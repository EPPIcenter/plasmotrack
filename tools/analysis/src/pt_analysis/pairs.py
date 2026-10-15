"""The pairs worth looking at (the React app's Pairs workspace): pairs of Infections whose order in the Ordering is
uncertain or whose Chains disagree on it, with their Edge marginals both ways.

From each Chain's derived `order_pairs` (each pair once, ids in text order as A and B): the Draws with both in the
Ordering, P(A before B) and the switches per Chain, and pooled over Chains (Draws summed). A pair is listed when its
pooled P(A before B) lies in [lo, hi], or when its Chains' probabilities differ by at least `disagree` (largest minus
smallest over the Chains with Draws of the pair). Its Edge marginals are the pooled rows of the Analysis's
`edge_stats`: `a_to_b` with A the parent of B, `b_to_a` the reverse, 0 without a row. Pairs whose larger Edge
marginal is under `edge_min` are left out. Rows come with the larger Edge marginal first, then the order nearest an
even split, then by ids.
"""

from functools import reduce

from pt_analysis.derive import derive_analysis, derived, derived_analysis
from pt_analysis.sql import Scope

DEFAULT_BAND = (0.1, 0.9)
DEFAULT_LIMIT = 200


def uncertain_pairs(analysis, lo=DEFAULT_BAND[0], hi=DEFAULT_BAND[1], disagree=None, edge_min=0.0,
                    limit=DEFAULT_LIMIT):
    """{"chains", "total", "rows": [{a, b, draws, probability, switches, spread, chains: [{chain, draws, probability,
    switches}], a_to_b, b_to_a}]}: the first `limit` listed pairs and how many there are in all. `disagree` None
    lists no pair for disagreement alone."""
    derive_analysis(analysis)
    labels = list(analysis.runs)
    connection = next(iter(analysis.runs.values())).connection
    order = reduce(lambda x, y: x.union(y), (derived(store, "order_pairs").project(f"{k} AS chain, *")
                                              for k, store in enumerate(analysis.runs.values())))
    with Scope(connection) as s:
        pairs, edges = s.name(order), s.name(derived_analysis(analysis, "edge_stats"))
        rows = s.rows(f"""
            WITH per_pair AS (
                SELECT i, j, sum(draws) AS draws, sum(before) / sum(draws) AS p, sum(switches) AS switches,
                       max(before / draws) - min(before / draws) AS spread,
                       map(list(chain), list(struct_pack(draws, before, switches))) AS chains
                FROM {pairs} WHERE draws > 0 GROUP BY i, j),
            pooled AS (SELECT child, parent, mean FROM {edges} WHERE chain IS NULL AND NOT unobserved),
            listed AS (
                SELECT per_pair.*, coalesce(ab.mean, 0) AS a_to_b, coalesce(ba.mean, 0) AS b_to_a
                FROM per_pair
                LEFT JOIN pooled ab ON ab.parent = i AND ab.child = j
                LEFT JOIN pooled ba ON ba.parent = j AND ba.child = i
                WHERE (p BETWEEN $lo AND $hi OR ($disagree IS NOT NULL AND spread >= $disagree - 1e-12))
                  AND greatest(coalesce(ab.mean, 0), coalesce(ba.mean, 0)) >= $edge_min)
            SELECT *, count(*) OVER () AS total FROM listed
            ORDER BY greatest(a_to_b, b_to_a) DESC, abs(p - 0.5), i, j
            LIMIT $limit""", {"lo": lo, "hi": hi, "disagree": disagree, "edge_min": edge_min, "limit": limit})
    out = []
    total = 0
    for i, j, draws, p, switches, spread, chains, a_to_b, b_to_a, total in rows:
        per_chain = []
        for k, label in enumerate(labels):
            c = chains.get(k)
            per_chain.append({"chain": label, "draws": c["draws"] if c else 0,
                              "probability": c["before"] / c["draws"] if c and c["draws"] else None,
                              "switches": c["switches"] if c else 0})
        out.append({"a": i, "b": j, "draws": draws, "probability": p, "switches": switches, "spread": spread,
                    "chains": per_chain, "a_to_b": a_to_b, "b_to_a": b_to_a})
    if not rows and limit == 0:
        total = uncertain_pairs(analysis, lo, hi, disagree, edge_min, 1)["total"]
    return {"chains": labels, "total": total, "rows": out}
