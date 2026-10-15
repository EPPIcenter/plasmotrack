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
