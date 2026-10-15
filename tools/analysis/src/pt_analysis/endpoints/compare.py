"""Comparison of Analyses of one input (U5), served by `pt-analysis serve a.json b.json ...`; see `pt_analysis.compare`.

    GET /api/analyses
        {input_sha256, analyses: [{label, primary, path, chains: [{chain, draws}], draws, sampling_seconds}]}
    GET /api/compare/convergence
        {analyses: [label], timed: [label], rows: [{group, name, analyses: {label: {ess, rhat[, ess_per_second]}}}]}

The primary Analysis comes first. `sampling_seconds` is null, and `ess_per_second` is absent, for an Analysis
without timings; `timed` lists the Analyses with them. The trace endpoints take `analysis=LABEL` to read any of them.
"""

from starlette.routing import Route

from pt_analysis.analysis import AnalysisError
from pt_analysis.compare import compare_convergence, input_hash, sampling_seconds
from pt_analysis.endpoints import json


def open_analyses(request):
    return request.app.state.analyses


async def analyses(request):
    opened = open_analyses(request)
    try:
        sha = input_hash(request.app.state.analysis)
    except AnalysisError:
        sha = None
    listed = []
    for label, analysis in opened.items():
        chains = [{"chain": chain, "draws": store.draws()} for chain, store in analysis.runs.items()]
        listed.append({"label": label, "primary": analysis is request.app.state.analysis, "path": str(analysis.path),
                       "chains": chains, "draws": sum(c["draws"] for c in chains),
                       "sampling_seconds": sampling_seconds(analysis)})
    return json({"input_sha256": sha, "analyses": listed})


async def convergence(request):
    opened = open_analyses(request)
    rows = compare_convergence(opened)
    timed = [label for label, a in opened.items() if sampling_seconds(a) is not None]
    return json({"analyses": list(opened), "timed": timed, "rows": rows})


ROUTES = [Route("/api/analyses", analyses), Route("/api/compare/convergence", convergence)]
