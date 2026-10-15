"""The React app's Overview workspace and status bar (`web/DESIGN.md`); see `pt_analysis.overview`.

    GET /api/status     {derive: {current, chains: {label: {burn_in: [window], sampling: [window]}}},
                         capabilities: {epi, truth, cohort, pairwise_scores, analyses}, mixing: {quantities,
                         rhat_above_1_1, rhat_above_1_2, ess_below_100, ess_below_20}, events: [{chain, iteration,
                         phase, type}], attention: {workspace: count}}
                        with a truth file, `attention.errors` counts the posterior's error rows at the default
                        threshold (as `/api/errors` lists them)
    GET /api/overview[?band_lo=L&band_hi=H]
                        {groups: [{group, n, ess, rhat, ess_min, ess_median, rhat_max, rhat_above}], worst, edge_spread,
                         disagreement, set_disagreement: [{infection, tvd, top, p_top}],
                         edge_counts: {between_infections, from_unobserved_source}, band, samplers, events,
                         headline: {psi_edge, mean_q, classes, q, effective_size, top_edges},
                         thumbnail: {threshold, nodes, edges}}

`set_disagreement` lists the ten Infections whose Parent-set distributions differ most between Chains: `tvd` is the
largest total variation distance between two Chains' set weights, `top` the pooled MAP set's observed parents and
`p_top` its pooled weight. `edge_counts` counts the edge table's rows (every edge with weight) by whether the parent is
an Infection or the unobserved source.

`band_lo` and `band_hi` set the acceptance band sampler types are checked against (default 0.15 to 0.5). The overview
is remembered until a Chain closes another Window.
"""

from starlette.routing import Route

from pt_analysis import overview as compute
from pt_analysis import truth_overlay
from pt_analysis.endpoints import error, json, memo, requested_analysis
from pt_analysis.endpoints.truth import FAILURES


def error_count(request):
    """The posterior's error rows against truth at the default threshold, remembered as `/api/errors` remembers
    them; None without truth or when they cannot be computed."""
    analysis = requested_analysis(request)
    if analysis.truth is None:
        return None
    threshold = truth_overlay.DEFAULT_THRESHOLD
    try:
        found = memo(request, ("errors", "posterior", threshold),
                     lambda: truth_overlay.errors(analysis, "posterior", threshold))
    except FAILURES:
        return None
    return len(found["rows"])


async def status(request):
    result = compute.status(requested_analysis(request), len(getattr(request.app.state, "analyses", {}) or [1]))
    errors = error_count(request)
    if errors is not None:
        result["attention"]["errors"] = errors
    return json(result)


async def overview(request):
    params = request.query_params
    try:
        band = (float(params.get("band_lo", compute.BAND[0])), float(params.get("band_hi", compute.BAND[1])))
    except ValueError:
        return error("band_lo and band_hi are acceptance rates between 0 and 1")
    analysis = requested_analysis(request)
    return json(memo(request, ("overview", band), lambda: compute.overview(analysis, band)))


ROUTES = [Route("/api/status", status), Route("/api/overview", overview)]
