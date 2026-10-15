"""The Analysis summary and its convergence table (U1)."""

from starlette.routing import Route

from pt_analysis.derive import analysis_dir
from pt_analysis.endpoints import json


async def summary(request):
    """Label, Draws per Chain and in all, and the number of Infections and Loci (from the first Chain's structure)."""
    analysis = request.app.state.analysis
    chains = [{"chain": label, "draws": store.draws()} for label, store in analysis.runs.items()]
    store = next(iter(analysis.runs.values()))
    infections = store.structure("infections").aggregate("count(*)").fetchone()[0]
    loci = store.structure("loci").aggregate("count(distinct locus)").fetchone()[0]
    return json({"label": analysis.label, "chains": chains, "draws": sum(c["draws"] for c in chains),
                 "infections": infections, "loci": loci})


async def convergence(request):
    """The convergence table as rows {group, name, chain, ess, rhat}, in its stored order (groups in display order,
    names sorted, each quantity's Chain rows then its pooled row with chain null)."""
    analysis = request.app.state.analysis
    store = next(iter(analysis.runs.values()))
    table = store.connection.read_parquet(str(analysis_dir(analysis) / "convergence.parquet")).to_arrow_table()
    return json({"rows": table.to_pylist()})


ROUTES = [Route("/api/analysis", summary), Route("/api/convergence", convergence)]
