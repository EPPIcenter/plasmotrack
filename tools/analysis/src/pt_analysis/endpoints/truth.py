"""The network view's truth overlay and the error table (U6, part 1), computed by `pt_analysis.truth_overlay`.

    GET /api/truth-network[?floor=0.01]
        {floor, truth: {n_true_edges, unmatched_infections, unobserved_truth_infections} or null,
         sources: [{source, label, directed, edges: [{parent, child, marginal, true, reverse_true, reverse_marginal,
                                                      lag}]}]}
    GET /api/errors[?source=posterior&threshold=0.5]
        {source, label, directed, threshold, truth, sources: [source names],
         rows: [{kind: missed | false | direction, parent, child, marginal, reverse_true, reverse_marginal, lag,
                 confidence, parent_true_infection_time, child_true_infection_time}]}
        404 without a truth file or when the manifest does not name that source

`sources` holds the posterior, then the null-model Analysis (`null`) and the pairwise scores (`ibd`, undirected) when
the manifest names them. Without truth the overlay's `true` and `reverse_true` are null. Results are remembered per
parameter (`endpoints.memo`) until a Chain closes another Window.
"""

import math

from starlette.routing import Route

from pt_analysis import metrics, truth_overlay
from pt_analysis.analysis import AnalysisError
from pt_analysis.endpoints import error, json, memo

FAILURES = (OSError, KeyError, ValueError, AnalysisError, metrics.MetricsError)


def probability(request, name, default):
    """The query parameter `name` as a number in [0, 1], or None."""
    value = request.query_params.get(name, "")
    try:
        x = default if value == "" else float(value)
    except ValueError:
        return None
    return x if 0 <= x <= 1 and not math.isnan(x) else None


async def truth_network(request):
    floor = probability(request, "floor", truth_overlay.DEFAULT_FLOOR)
    if floor is None:
        return error(f"floor is a score in [0, 1], not '{request.query_params['floor']}'")
    analysis = request.app.state.analysis
    try:
        return json(memo(request, ("truth-network", floor), lambda: truth_overlay.overlay(analysis, floor)))
    except FAILURES as e:
        return error(str(e))


async def errors(request):
    threshold = probability(request, "threshold", truth_overlay.DEFAULT_THRESHOLD)
    if threshold is None:
        return error(f"threshold is a score in [0, 1], not '{request.query_params['threshold']}'")
    source = request.query_params.get("source", "posterior")
    if source not in truth_overlay.SOURCES:
        return error(f"source is one of {', '.join(truth_overlay.SOURCES)}, not '{source}'")
    analysis = request.app.state.analysis
    if analysis.truth is None:
        return error("the Analysis manifest names no truth file, so there are no errors to list", 404)
    try:
        return json(memo(request, ("errors", source, threshold),
                           lambda: truth_overlay.errors(analysis, source, threshold)))
    except truth_overlay.UnknownSource as e:
        return error(str(e), 404)
    except FAILURES as e:
        return error(str(e))


ROUTES = [Route("/api/truth-network", truth_network), Route("/api/errors", errors)]
