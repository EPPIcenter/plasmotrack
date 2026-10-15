"""JSON endpoints of `pt-analysis serve`, one module per view. See `pt_analysis.server`."""

import math

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


def error(message, status_code=400):
    return json({"error": message}, status_code)


class UnknownAnalysis(Exception):
    """The request names an Analysis the server was not started with."""


def requested_analysis(request):
    """The Analysis named by the request's `analysis` parameter (a label, U5), or the primary Analysis without one."""
    label = request.query_params.get("analysis")
    if label in (None, ""):
        return request.app.state.analysis
    analyses = getattr(request.app.state, "analyses", {})
    if label not in analyses:
        raise UnknownAnalysis(f"no Analysis '{label}' is open; open are: {', '.join(analyses)}")
    return analyses[label]
