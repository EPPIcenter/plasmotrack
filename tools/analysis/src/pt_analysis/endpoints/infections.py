"""The Infections workspace's table (`web/DESIGN.md`, Infections); see `pt_analysis.infections`.

    GET /api/infections   {chains, columns: {column: [value per Infection]}, q_chains: {Chain: [mean q]},
                           epi_columns: [{name, type}]}

Columnar: every list in `columns` has one value per Infection, in input order. `analysis=LABEL` names an open
Analysis other than the primary. Remembered per Analysis until a Chain closes
another Window.
"""

from starlette.routing import Route

from pt_analysis.endpoints import json, memo, requested_analysis
from pt_analysis.infections import infections as compute


async def infections(request):
    analysis = requested_analysis(request)
    return json(memo(request, ("infections",), lambda: compute(analysis)))


ROUTES = [Route("/api/infections", infections)]
