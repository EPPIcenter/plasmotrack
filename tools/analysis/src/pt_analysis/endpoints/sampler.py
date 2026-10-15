"""The sampler panel (U4): one Chain's acceptance and swap acceptance per Window, temperature ladder, burn-in-end
proposal variances and the Run events marked on its time axes.

    GET /api/sampler[?chain=LABEL]   {chain, windows, acceptance, swaps, ladder, burn_in_ended, events}
    GET /api/sampler?type=T[&chain=LABEL][&limit=N]
                                     {chain, type, targets, shown, acceptance, proposal_variances}

`chain` is a Chain label of the Analysis; the first Chain by default. Without `type`, acceptance and the burn-in-end
proposal variances are per sampler type (`pt_analysis.sampler.panel`); a Run can have hundreds of thousands of
samplers. With `type`, that type's targets, at most `limit` (default 500) with the lowest acceptance first
(`pt_analysis.sampler.type_panel`); a type the Chain has no sampler of is a 404. Many rows are columnar: one array
per field.
"""

from starlette.routing import Route

from pt_analysis import sampler
from pt_analysis.endpoints import error, json


async def sampler_panel(request):
    runs = request.app.state.analysis.runs
    chain = request.query_params.get("chain") or next(iter(runs))
    if chain not in runs:
        return error(f"no Chain '{chain}'; the Chains are {', '.join(runs)}", 404)
    kind = request.query_params.get("type")
    if not kind:
        return json({"chain": chain, **sampler.panel(runs[chain])})
    raw = request.query_params.get("limit", str(sampler.LIMIT))
    if not raw.isdigit() or int(raw) < 1:
        return error(f"limit is a positive number of targets, not '{raw}'")
    got = sampler.type_panel(runs[chain], kind, int(raw))
    if not got["targets"]:
        return error(f"Chain '{chain}' has no sampler of type '{kind}'", 404)
    return json({"chain": chain, **got})


ROUTES = [Route("/api/sampler", sampler_panel)]
