"""Per-Chain traces of one quantity and ESS and R-hat over a Draw range (U1).

    GET /api/trace?group=G&name=N[&burn_in=true]   {group, name, burn_in, chains: [{chain, iteration, phase, value}]}
    GET /api/range?group=G&name=N[&first=I][&last=I]
                                                   {group, name, first, last, chains: [{chain, draws, ess}], ess, rhat}

`first` and `last` are Sampling iterations, inclusive; only Draws (sampling-Phase rows) count.
"""

from starlette.routing import Route

from pt_analysis.endpoints import error, json
from pt_analysis.traces import TraceError, UnknownQuantity, chain_traces, range_diagnostics


def quantity(request):
    group, name = request.query_params.get("group"), request.query_params.get("name")
    if group is None or name is None:
        raise TraceError("a quantity needs both group and name")
    return group, name


def iteration(request, key):
    value = request.query_params.get(key)
    if value in (None, ""):
        return None
    try:
        return int(float(value))
    except ValueError:
        raise TraceError(f"{key} is a Sampling iteration, not '{value}'") from None


async def trace(request):
    try:
        group, name = quantity(request)
        burn_in = request.query_params.get("burn_in", "false").lower() in ("1", "true", "yes")
        traces, with_burn_in = chain_traces(request.app.state.analysis, group, name, burn_in)
    except UnknownQuantity as e:
        return error(str(e), 404)
    except TraceError as e:
        return error(str(e))
    return json({"group": group, "name": name, "burn_in": with_burn_in,
                 "chains": [{"chain": label, **s} for label, s in traces.items()]})


async def draw_range(request):
    try:
        group, name = quantity(request)
        first, last = iteration(request, "first"), iteration(request, "last")
        result = range_diagnostics(request.app.state.analysis, group, name, first, last)
    except UnknownQuantity as e:
        return error(str(e), 404)
    except TraceError as e:
        return error(str(e))
    return json({"group": group, "name": name, "first": first, "last": last, **result})


ROUTES = [Route("/api/trace", trace), Route("/api/range", draw_range)]
