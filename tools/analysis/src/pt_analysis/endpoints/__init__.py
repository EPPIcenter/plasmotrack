"""JSON endpoints of `pt-analysis serve`, one module per view. See `pt_analysis.server`."""

import math
from json import dumps as json_text

from pt_analysis.derive import analysis_source

from starlette.responses import JSONResponse


def finite(value):
    """`value` with non-finite floats as the strings "Infinity", "-Infinity" and "NaN", which JSON cannot hold; the
    frontend's API client turns them back into numbers. Recurses into lists and dicts."""
    if isinstance(value, float) and not math.isfinite(value):
        return "NaN" if math.isnan(value) else ("Infinity" if value > 0 else "-Infinity")
    if isinstance(value, dict):
        return {k: finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite(v) for v in value]
    return value


def json(value, status_code=200):
    return JSONResponse(finite(value), status_code=status_code)


def memo(request, key, compute):
    """`compute()`, remembered per app under `key` for the request's Analysis (`requested_analysis`) while its Chains'
    closed Windows stay the same: a Chain that closes another Window makes every result remembered for it stale."""
    analysis = requested_analysis(request)
    source = json_text(analysis_source(analysis))
    caches = getattr(request.app.state, "memo", None)
    if caches is None:
        caches = request.app.state.memo = {}
    cache = caches.get(analysis.label)
    if cache is None or cache["source"] != source:
        cache = caches[analysis.label] = {"source": source, "results": {}}
    if key not in cache["results"]:
        cache["results"][key] = compute()
    return cache["results"][key]


def error(message, status_code=400):
    return json({"error": message}, status_code)


class UnknownAnalysis(Exception):
    """The request names an Analysis the server was not started with: a 404 naming it (`unknown_analysis`)."""


async def unknown_analysis(request, exc):
    """The app's handler of `UnknownAnalysis`, raised from any endpoint."""
    return error(str(exc), 404)


def requested_analysis(request):
    """The Analysis named by the request's `analysis` parameter (a label, U5), or the primary Analysis without one.
    Every endpoint that reads one Analysis reads it here, so each takes `analysis=LABEL`; an unknown label raises
    `UnknownAnalysis`, which the app answers with a 404."""
    label = request.query_params.get("analysis")
    if label in (None, ""):
        return request.app.state.analysis
    analyses = getattr(request.app.state, "analyses", {})
    if label not in analyses:
        raise UnknownAnalysis(f"no Analysis '{label}' is open; open are: {', '.join(analyses)}")
    return analyses[label]
