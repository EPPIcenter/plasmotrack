"""The sampler panel (U4): one Chain's acceptance and swap acceptance per Window, temperature ladder, burn-in-end
proposal variances and the Run events marked on its time axes.

    GET /api/sampler[?chain=LABEL]   {chain, windows, acceptance, swaps, ladder, burn_in_ended, events}

`chain` is a Chain label of the Analysis; the first Chain by default. The fields are those of
`pt_analysis.sampler.panel`.
"""

from starlette.routing import Route

from pt_analysis import sampler
from pt_analysis.endpoints import error, json


async def sampler_panel(request):
    runs = request.app.state.analysis.runs
    chain = request.query_params.get("chain") or next(iter(runs))
    if chain not in runs:
        return error(f"no Chain '{chain}'; the Chains are {', '.join(runs)}", 404)
    return json({"chain": chain, **sampler.panel(runs[chain])})


ROUTES = [Route("/api/sampler", sampler_panel)]
