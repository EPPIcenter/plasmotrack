"""The React app's Overview workspace and status bar (`web/DESIGN.md`); see `pt_analysis.overview`.

    GET /api/status     {derive: {current, chains: {label: {burn_in: [window], sampling: [window]}}},
                         capabilities: {epi, truth, cohort, pairwise_scores, analyses}, mixing: {quantities,
                         rhat_above_1_1, rhat_above_1_2, ess_below_100, ess_below_20}, events: [{chain, iteration,
                         phase, type}], attention: {workspace: count}}
    GET /api/overview[?band_lo=L&band_hi=H]
                        {groups: [{group, n, ess, rhat, ess_min, ess_median, rhat_max, rhat_above}], worst, edge_spread,
                         disagreement, band, samplers, events, headline: {psi_edge, mean_q, classes, q, effective_size,
                         top_edges}, thumbnail: {threshold, nodes, edges}}

`band_lo` and `band_hi` set the acceptance band sampler types are checked against (default 0.15 to 0.5). The overview
is remembered until a Chain closes another Window.
"""

from starlette.routing import Route

from pt_analysis import overview as compute
from pt_analysis.endpoints import error, json, memo, requested_analysis


async def status(request):
    return json(compute.status(requested_analysis(request), len(getattr(request.app.state, "analyses", {}) or [1])))


async def overview(request):
    params = request.query_params
    try:
        band = (float(params.get("band_lo", compute.BAND[0])), float(params.get("band_hi", compute.BAND[1])))
    except ValueError:
        return error("band_lo and band_hi are acceptance rates between 0 and 1")
    analysis = requested_analysis(request)
    return json(memo(request, ("overview", band), lambda: compute.overview(analysis, band)))


ROUTES = [Route("/api/status", status), Route("/api/overview", overview)]
