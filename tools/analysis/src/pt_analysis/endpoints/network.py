"""The network view's nodes and edges (U8).

    GET /api/network[?floor=F]   {chains, floor, nodes: [{id, observation_time, symptomatic, q,
                                                          infection_time: {lower, median, upper}}],
                                  edges: [{child, parent, marginal, chains: [mean per Chain], lag}]}

Edges are those from observed parents with pooled marginal at least `floor` (default 0.01, in [0, 1]), largest first;
the view's threshold slider filters them further in the browser. See `pt_analysis.network`.
"""

import math

from starlette.routing import Route

from pt_analysis.endpoints import error, json
from pt_analysis.network import DEFAULT_FLOOR, network


async def network_payload(request):
    value = request.query_params.get("floor", "")
    try:
        floor = DEFAULT_FLOOR if value == "" else float(value)
    except ValueError:
        floor = math.nan
    if not 0 <= floor <= 1:
        return error(f"floor is an Edge marginal in [0, 1], not '{value}'")
    return json(network(request.app.state.analysis, floor))


ROUTES = [Route("/api/network", network_payload)]
