"""The infection card (U10).

    GET /api/infection?id=ID   `infection.card` for Infection ID: {infection, observation_time, quantities, coi,
                               genotype, parent_sets, parents, children, summary}

400 without an id; 404 for an id that is no Infection of the Analysis.
"""

from starlette.routing import Route

from pt_analysis.endpoints import error, json
from pt_analysis.infection import UnknownInfection, card


async def infection(request):
    infection_id = request.query_params.get("id", "").strip()
    if not infection_id:
        return error("the infection card needs an Infection id")
    try:
        return json(card(request.app.state.analysis, infection_id))
    except UnknownInfection as e:
        return error(str(e), 404)


ROUTES = [Route("/api/infection", infection)]
