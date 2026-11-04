"""The pair inspector (U3).

    GET /api/pair?a=A&b=B   `pair.inspect_pair` for Infections A and B: {a, b, order, edges, co_inclusion, genotypes,
                            times}

400 without both ids or with the same id twice; 404 for an id that is no Infection of the Analysis.
"""

from starlette.routing import Route

from pt_analysis.endpoints import error, json
from pt_analysis.pair import PairError, UnknownInfection, inspect_pair


async def pair(request):
    a, b = request.query_params.get("a", "").strip(), request.query_params.get("b", "").strip()
    if not a or not b:
        return error("a pair needs two Infection ids, a and b")
    try:
        return json(inspect_pair(request.app.state.analysis, a, b))
    except UnknownInfection as e:
        return error(str(e), 404)
    except PairError as e:
        return error(str(e))


ROUTES = [Route("/api/pair", pair)]
