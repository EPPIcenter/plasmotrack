"""The edge table (U2): per (child, parent), the pooled Edge marginal, each Chain's mean and switch count, the
between-Chain spread and the ESS, from the derived `edge_stats` table (`pt_analysis.edge_stats`).

    GET /api/edges   {chains: [label], rows: [{child, parent, unobserved, group, name, pooled, spread, switches, ess,
                                               chains: {label: {mean, switches}}}]}

`unobserved` marks the unobserved source (the child's Latent parent id) as parent; `group` and `name` are the edge's
quantity in the convergence table and the trace panel (`edge`, `child<-parent`; or `unobserved-source inclusion`,
child). `switches` sums the Chains' switch counts.
"""

from starlette.routing import Route

from pt_analysis.edge_stats import edge_table
from pt_analysis.endpoints import json, requested_analysis


async def edges(request):
    analysis = requested_analysis(request)
    rows = [{"child": e.child, "parent": e.parent, "unobserved": e.unobserved, "group": e.group, "name": e.name,
             "pooled": e.pooled, "spread": e.spread, "switches": e.switches, "ess": e.ess,
             "chains": {label: {"mean": c.mean, "switches": c.switches} for label, c in e.chains.items()}}
            for e in edge_table(analysis)]
    return json({"chains": list(analysis.runs), "rows": rows})


ROUTES = [Route("/api/edges", edges)]
