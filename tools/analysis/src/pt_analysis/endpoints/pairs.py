"""The React app's Pairs workspace and the Edges workspace's sparklines (`web/DESIGN.md`).

    GET /api/pairs[?min=0.1&max=0.9&disagree=D&edge_min=0&limit=200]
        {chains, min, max, disagree, edge_min, limit, total,
         rows: [{a, b, draws, probability, switches, spread, chains: [{chain, draws, probability, switches}],
                 a_to_b, b_to_a, differing_loci}]}
    GET /api/edge-series?edges=child<-parent,...[&points=N]
        {chains, points, series: [{child, parent, chains: [[weight per Draw, or per run of Draws] per Chain]}]}

`/api/pairs` lists pairs whose pooled P(A before B) lies in [min, max], or whose Chains differ by at least `disagree`
(left out: no pair for disagreement alone), and whose larger Edge marginal is at least `edge_min`; `total` counts
them all, `rows` holds the first `limit`, each with `differing_loci`, the Loci where A's and B's pooled Latent genotype
modes differ (as `/api/pair`'s `genotypes.differing`); the list is remembered until a Chain closes another Window. See
`pt_analysis.pairs`. `/api/edge-series` gives at most 500 edges' weights per Draw, binned to `points` means per Chain
when given; see `pt_analysis.edge_series`. Both read the Analysis named by `analysis=` (the primary without it).
"""

import math

from starlette.routing import Route

from pt_analysis import pairs as compute
from pt_analysis.edge_series import edge_series
from pt_analysis.endpoints import error, json, memo, requested_analysis

MAX_EDGES = 500


def number(params, name, default, low=0.0, high=1.0):
    """The query parameter `name` as a number in [low, high]; `default` when absent; ValueError otherwise."""
    value = params.get(name, "")
    if value == "":
        return default
    x = float(value)
    if math.isnan(x) or not low <= x <= high:
        raise ValueError(name)
    return x


async def pairs(request):
    params = request.query_params
    try:
        lo = number(params, "min", compute.DEFAULT_BAND[0])
        hi = number(params, "max", compute.DEFAULT_BAND[1])
        disagree = number(params, "disagree", None)
        edge_min = number(params, "edge_min", 0.0)
        limit = int(number(params, "limit", compute.DEFAULT_LIMIT, 0, math.inf))
    except ValueError:
        return error("min, max, disagree and edge_min are probabilities in [0, 1] and limit a count of rows")
    if lo > hi:
        return error(f"min ({lo}) is above max ({hi})")
    analysis = requested_analysis(request)
    found = memo(request, ("pairs", lo, hi, disagree, edge_min, limit),
                 lambda: compute.uncertain_pairs(analysis, lo, hi, disagree, edge_min, limit))
    return json({"min": lo, "max": hi, "disagree": disagree, "edge_min": edge_min, "limit": limit, **found})


def edge_keys(text):
    """[(child, parent)] of `child<-parent,...`; ValueError for a key without `<-`."""
    keys = []
    for key in filter(None, text.split(",")):
        child, sep, parent = key.partition("<-")
        if not sep or not child or not parent:
            raise ValueError(key)
        keys.append((child, parent))
    return keys


async def series(request):
    params = request.query_params
    try:
        keys = edge_keys(params.get("edges", ""))
        points = int(params["points"]) if params.get("points") else None
    except ValueError:
        return error("edges are `child<-parent` keys, comma-separated, and points a whole number")
    if len(keys) > MAX_EDGES:
        return error(f"at most {MAX_EDGES} edges per request, not {len(keys)}")
    if points is not None and points < 1:
        return error("points is at least 1")
    return json({"points": points, **edge_series(requested_analysis(request), keys, points)})


ROUTES = [Route("/api/pairs", pairs), Route("/api/edge-series", series)]
