"""Per-Chain traces of one quantity and ESS and R-hat over a Draw range (U1).

    GET /api/trace?group=G&name=N[&burn_in=true]   {group, name, burn_in, chains: [{chain, iteration, phase, value}]}
    GET /api/range?group=G&name=N[&first=I][&last=I]
                                                   {group, name, first, last, chains: [{chain, draws, ess}], ess, rhat}

`first` and `last` are Sampling iterations, inclusive; only Draws (sampling-Phase rows) count.

A functional of one stratum (U9, `pt_analysis.strata`) is group `functionals` with `stratum_column=C` (empty for every
Infection) and `stratum=V` (left out for the Infections with no value of C); both are echoed in the response.

With several Analyses open (U5), `analysis=LABEL` reads that Analysis instead of the primary one, and is echoed.
"""

from starlette.routing import Route

from pt_analysis.endpoints import UnknownAnalysis, error, json, requested_analysis
from pt_analysis.traces import TraceError, UnknownQuantity, chain_traces, range_diagnostics


def quantity(request):
    group, name = request.query_params.get("group"), request.query_params.get("name")
    if group is None or name is None:
        raise TraceError("a quantity needs both group and name")
    return group, name


def stratum(request):
    """(column or None, value or None) when the query names a stratum column, else None."""
    params = request.query_params
    if "stratum_column" not in params:
        return None
    return params["stratum_column"] or None, params.get("stratum")


def echo(request):
    params = request.query_params
    out = {"analysis": params["analysis"]} if params.get("analysis") else {}
    if "stratum_column" not in params:
        return out
    return {**out, "stratum_column": params["stratum_column"], "stratum": params.get("stratum")}


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
        traces, with_burn_in = chain_traces(requested_analysis(request), group, name, burn_in, stratum(request))
    except (UnknownQuantity, UnknownAnalysis) as e:
        return error(str(e), 404)
    except TraceError as e:
        return error(str(e))
    return json({"group": group, "name": name, **echo(request), "burn_in": with_burn_in,
                 "chains": [{"chain": label, **s} for label, s in traces.items()]})


async def draw_range(request):
    try:
        group, name = quantity(request)
        first, last = iteration(request, "first"), iteration(request, "last")
        result = range_diagnostics(requested_analysis(request), group, name, first, last, stratum(request))
    except (UnknownQuantity, UnknownAnalysis) as e:
        return error(str(e), 404)
    except TraceError as e:
        return error(str(e))
    return json({"group": group, "name": name, **echo(request), "first": first, "last": last, **result})


ROUTES = [Route("/api/trace", trace), Route("/api/range", draw_range)]
