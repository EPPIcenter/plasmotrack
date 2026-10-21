"""Several quantities' Draw values per Chain at once, compact: the convergence table's sparklines.

    GET /api/series?q=GROUP:NAME[&q=GROUP:NAME...][&points=N]
        {chains: [label], points, series: [{group, name, chains: [[value per Draw, or per run of Draws] per Chain]
                                           or null}]}

Each `q` is a quantity of the convergence table, split at its first ':' into group and name; at most 200, answered in
request order. Values are the sampling Draws of `/api/trace` (`traces.chain_traces`); with `points`, each Chain's
series is binned into that many means of adjacent Draws as `/api/edge-series` bins them (`edge_series.binned`).
`chains` is null for a quantity no Chain has or that is malformed for its group. Functionals of a stratum are not
served here. Results are remembered per quantity list and `points` until a Chain closes another Window.
"""

from starlette.routing import Route

from pt_analysis.edge_series import binned
from pt_analysis.endpoints import error, json, memo, requested_analysis
from pt_analysis.traces import TraceError, chain_traces

MAX_QUANTITIES = 200


def compact_series(analysis, quantities, points=None):
    """{"chains", "series"} of `quantities`, a list of (group, name); see the module docstring."""
    out = []
    for group, name in quantities:
        try:
            traces, _ = chain_traces(analysis, group, name)
            chains = [binned(s["value"], points) for s in traces.values()]
        except TraceError:  # UnknownQuantity too
            chains = None
        out.append({"group": group, "name": name, "chains": chains})
    return {"chains": list(analysis.runs), "series": out}


async def series(request):
    params = request.query_params
    quantities = []
    for q in params.getlist("q"):
        group, sep, name = q.partition(":")
        if not sep:
            return error(f"a quantity is GROUP:NAME, not '{q}'")
        quantities.append((group, name))
    if len(quantities) > MAX_QUANTITIES:
        return error(f"at most {MAX_QUANTITIES} quantities per request, not {len(quantities)}")
    try:
        points = int(params["points"]) if params.get("points") else None
    except ValueError:
        return error("points is a whole number")
    if points is not None and points < 1:
        return error("points is at least 1")
    analysis = requested_analysis(request)
    found = memo(request, ("series", tuple(quantities), points),
                 lambda: compact_series(analysis, quantities, points))
    return json({"points": points, **found})


ROUTES = [Route("/api/series", series)]
