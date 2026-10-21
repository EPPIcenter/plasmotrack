"""Genotypes of a group of Infections, one column each.

    GET /api/genotypes?ids=A,B,C   `pair.group_genotypes`

400 without ids, with more than 12, or with an id that is no Infection of the Analysis.
"""

from starlette.routing import Route

from pt_analysis.endpoints import error, json, memo, requested_analysis
from pt_analysis.pair import PairError, UnknownInfection, group_genotypes


async def genotypes(request):
    ids = [i.strip() for i in request.query_params.get("ids", "").split(",") if i.strip()]
    if not ids:
        return error("name the Infections with ids")
    try:
        return json(memo(request, f"genotypes:{','.join(ids)}", lambda: group_genotypes(requested_analysis(request), ids)))
    except UnknownInfection as e:
        return error(str(e), 404)
    except PairError as e:
        return error(str(e))


ROUTES = [Route("/api/genotypes", genotypes)]
