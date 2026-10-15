"""The epi data's columns and the functionals panel (U9).

    GET /api/epi                    {key, columns: [{name, type}], unmatched, missing, infections}
    GET /api/functionals[?column=C][&level=L]
                                    {column, level, classes, strata: [{value, n, psi_map, size_range, psi_edge,
                                     mean_q, classes: {label: interval}, multi_source_augmented,
                                     multi_source_observed}]}

An interval is {mean, lower, upper}. Without `column`, one stratum of every Infection; a stratum's `value` is null
for it and for the Infections with no value of the column. A column the epi data does not have is a 400 naming it.
The per-Draw series behind each cell come from `/api/trace` (see `pt_analysis.endpoints.traces`).
"""

from dataclasses import asdict

from starlette.routing import Route

from pt_analysis.endpoints import error, json
from pt_analysis.epi import EpiError, load_epi
from pt_analysis.strata import LEVEL, labels, stratified_functionals


async def epi(request):
    """The epi data's key and columns with their types, the keys that are no Infection (`unmatched`), the
    Infections without a row (`missing`) and the number of Infections with a row."""
    try:
        data = load_epi(request.app.state.analysis)
    except EpiError as e:
        return error(str(e))
    if data is None:
        return json({"key": None, "columns": [], "unmatched": [], "missing": [], "infections": 0})
    return json({"key": data.key, "columns": [{"name": n, "type": t} for n, t in data.columns.items()],
                 "unmatched": data.unmatched, "missing": data.missing,
                 "infections": len(data.infections) - len(data.missing)})


async def functionals(request):
    analysis = request.app.state.analysis
    column = request.query_params.get("column") or None
    try:
        level = float(request.query_params.get("level", LEVEL))
    except ValueError:
        return error(f"level is a probability, not '{request.query_params['level']}'")
    try:
        strata = stratified_functionals(analysis, column, level=level)
    except EpiError as e:
        return error(str(e))
    return json({"column": column, "level": level, "classes": labels(analysis),
                 "strata": [asdict(s) for s in strata]})


ROUTES = [Route("/api/epi", epi), Route("/api/functionals", functionals)]
