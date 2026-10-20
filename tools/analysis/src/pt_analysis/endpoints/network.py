"""The network view's nodes and edges (U8).

    GET /api/network[?floor=F]   {chains, floor, nodes: [{id, observation_time, symptomatic, candidate_parents, q,
                                                          infection_time: {lower, median, upper},
                                                          epi: {column: value}}],
                                  edges: [{child, parent, marginal, chains: [mean per Chain], lag},
                                  epi_columns: [{name, type, values}]}

Edges are those from observed parents with pooled marginal at least `floor` (default 0.01, in [0, 1]), largest first;
the view's threshold slider filters them further in the browser. `epi_columns` lists the epi data's columns with
their type (categorical, numeric or date) and stratum values, and each node's `epi` holds its value of each as text
(null where it has none); both are empty without epi data, and unreadable epi data is a 400. See
`pt_analysis.network`.

    GET /api/collection-edges   {edges: [{parent, child, entry}]}

`entry` is the pooled mass of the Parent sets ranked above the first one containing that parent, so the parent is
inside the child's credible collection at every level above it. The unobserved source is not an edge.
"""

import math

from starlette.routing import Route

from pt_analysis.endpoints import error, json, memo, requested_analysis
from pt_analysis.epi import EpiError
from pt_analysis.network import DEFAULT_FLOOR, network
from pt_analysis.parent_sets import collection_entries


async def network_payload(request):
    value = request.query_params.get("floor", "")
    try:
        floor = DEFAULT_FLOOR if value == "" else float(value)
    except ValueError:
        floor = math.nan
    if not 0 <= floor <= 1:
        return error(f"floor is an Edge marginal in [0, 1], not '{value}'")
    try:
        return json(network(requested_analysis(request), floor))
    except EpiError as e:
        return error(str(e))


async def collection(request):
    """Edges the credible collection retains, with the level at which each parent enters. The Network filters them."""
    try:
        edges = memo(request, "collection-edges", lambda: collection_entries(requested_analysis(request)))
    except (OSError, KeyError, ValueError) as e:
        return error(str(e))
    return json({"edges": [{"parent": p, "child": c, "entry": e} for p, c, e in edges]})


ROUTES = [Route("/api/network", network_payload), Route("/api/collection-edges", collection)]
