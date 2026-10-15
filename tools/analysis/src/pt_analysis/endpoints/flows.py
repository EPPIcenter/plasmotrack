"""Flows between strata, lags and spatial scale for the functionals panel (U9). See `pt_analysis.flows` for
the definitions.

    GET /api/flows[?column=C][&min_infections=N]
        {column, strata: [{value, n}],
         attribution: [{recipient, source, unobserved, mass, share}],
         enrichment: [{source, recipient, expected_edges, share, candidates, candidate_share, enrichment}],
         net_flow: {min_infections, hidden, strata: [{value, n, out_edges, in_edges, net}]}}
    GET /api/lags[?max_lag=L][&column=C][&width=W]
        {width, total, median_duration, bins: [{x0, x1, mass}],
         reassign: null without max_lag, else {max_lag, moved_mass, moved_share, overall, strata}}
            where overall and each stratum are {value, n, psi_edge, psi_edge_moved, mean_q, mean_q_moved}
    GET /api/spatial?unit=U&parent=P[&max_lag=L]
        {unit, parent, max_lag, total, shares: [{relation, mass, share}]}   relation: same_unit, same_parent, other

Stratum values are null for the Infections with no value of the column (and for the one stratum without a column);
in attribution the unobserved source is the row with `unobserved` true. Intervals are {mean, lower, upper}. A
column the epi data does not have, or a cutoff or width that is not a number, is a 400 naming it.
"""

from dataclasses import asdict

from starlette.routing import Route

from pt_analysis.endpoints import error, json
from pt_analysis.epi import EpiError
from pt_analysis.flows import (FlowsError, attribution, enrichment, lag_histogram, median_duration, nice_width,
                               net_flow, reassign, spatial_scale)
from pt_analysis.strata import strata


def _number(request, name, default=None):
    text = request.query_params.get(name)
    if not text:
        return default
    try:
        return float(text)
    except ValueError:
        raise FlowsError(f"{name} is a number, not '{text}'") from None


async def flows(request):
    analysis = request.app.state.analysis
    column = request.query_params.get("column") or None
    try:
        minimum = int(_number(request, "min_infections", 1))
        shown, hidden = net_flow(analysis, column, min_infections=minimum)
        return json({"column": column,
                     "strata": [{"value": v, "n": len(i)} for v, i in strata(analysis, column)],
                     "attribution": [asdict(c) for c in attribution(analysis, column)],
                     "enrichment": [asdict(c) for c in enrichment(analysis, column)],
                     "net_flow": {"min_infections": minimum, "hidden": hidden, "strata": [asdict(s) for s in shown]}})
    except (EpiError, FlowsError) as e:
        return error(str(e))


async def lags(request):
    analysis = request.app.state.analysis
    column = request.query_params.get("column") or None
    try:
        max_lag = _number(request, "max_lag")
        width = _number(request, "width")
        bins = lag_histogram(analysis, width)
        moved = None if max_lag is None else asdict(reassign(analysis, max_lag, column))
    except (EpiError, FlowsError) as e:
        return error(str(e))
    return json({"width": bins[0].x1 - bins[0].x0 if bins else width or nice_width(0),
                 "total": sum(b.mass for b in bins), "median_duration": median_duration(analysis),
                 "bins": [asdict(b) for b in bins], "reassign": moved})


async def spatial(request):
    analysis = request.app.state.analysis
    unit, parent = request.query_params.get("unit"), request.query_params.get("parent")
    if not unit or not parent:
        return error("spatial scale needs a unit column and its parent unit column: ?unit=...&parent=...")
    try:
        max_lag = _number(request, "max_lag")
        shares = spatial_scale(analysis, unit, parent, max_lag)
    except (EpiError, FlowsError) as e:
        return error(str(e))
    return json({"unit": unit, "parent": parent, "max_lag": max_lag, "total": sum(s.mass for s in shares),
                 "shares": [asdict(s) for s in shares]})


ROUTES = [Route("/api/flows", flows), Route("/api/lags", lags), Route("/api/spatial", spatial)]
