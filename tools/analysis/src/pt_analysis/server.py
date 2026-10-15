"""The local HTTP backend of `pt-analysis serve`: the built frontend as static files, and a JSON API over an Analysis.

Each module named in `ENDPOINTS` serves one view's data: it exports `ROUTES`, a list of Starlette routes under
`/api/`, whose handlers find the Analysis at `request.app.state.analysis`. A new view's backend is one module in
`pt_analysis/endpoints/` plus its line under its use case's comment in `ENDPOINTS`.
"""

import importlib
import threading
import webbrowser
from pathlib import Path

from starlette.applications import Starlette
from starlette.exceptions import HTTPException
from starlette.middleware import Middleware
from starlette.middleware.gzip import GZipMiddleware
from starlette.responses import PlainTextResponse
from starlette.routing import Mount, Route
from starlette.staticfiles import StaticFiles

from pt_analysis.derive import derive_analysis

# The frontend's build output, `tools/analysis/web/dist/` (`pnpm build` in `tools/analysis/web/`).
STATIC_DIR = Path(__file__).resolve().parents[2] / "web" / "dist"

ENDPOINTS = [
    "pt_analysis.endpoints.analysis",  # U1: Analysis summary, convergence table
    "pt_analysis.endpoints.traces",  # U1: trace panel, ESS and R-hat over a Draw range
    # U2 edge table
    "pt_analysis.endpoints.edges",  # U2: edge table

    # U8 network view
    "pt_analysis.endpoints.network",  # U8: network view's nodes and edges

    # U3 pair inspector
    "pt_analysis.endpoints.pair",  # U3: order, edges, co-inclusion, genotypes and infection times of a pair

    # U9 functionals panel
    "pt_analysis.endpoints.functionals",  # U9: epi columns, functionals overall and per stratum
    "pt_analysis.endpoints.flows",  # U9: attribution, enrichment, net flow, lags, spatial scale

    # U10 infection card
    "pt_analysis.endpoints.infection",  # U10: one Infection's posteriors, genotype, Parent sets and edges

    # U6 truth overlay, ranking
    "pt_analysis.endpoints.ranking",  # U6: ranking and parent-set metrics against truth
    "pt_analysis.endpoints.truth",  # U6: truth overlay of the network view, error table

    # U7 cohort view
    "pt_analysis.endpoints.cohort",  # U7: cohort of fit summaries, calibration counts, opening a fit's Analysis

    # U4 sampler panel
    "pt_analysis.endpoints.sampler",  # U4: acceptance, swaps, ladder, proposal variances, Run events

    # U5 comparison of Analyses
    "pt_analysis.endpoints.compare",  # U5: the open Analyses, their convergence side by side

    # React app (web/DESIGN.md): Overview workspace, status bar
    "pt_analysis.endpoints.overview",  # status bar and rail badges; Overview's health check and headline
    "pt_analysis.endpoints.pairs",  # Pairs workspace's list; Edges workspace's per-Draw series

]


class AppFiles(StaticFiles):
    """The built frontend, with `index.html` for any other path that is not under `/api/` and has no file extension:
    the app routes on the client (`/convergence?sel=...`), so a reload or a shared link of such a path gets the app."""

    async def get_response(self, path, scope):
        try:
            return await super().get_response(path, scope)
        except HTTPException as e:
            last = path.rsplit("/", 1)[-1]
            if e.status_code != 404 or path.startswith("api/") or path == "api" or "." in last:
                raise
            return await super().get_response("index.html", scope)


def not_built(request):
    return PlainTextResponse(f"The frontend is not built: run `pnpm install && pnpm build` in {STATIC_DIR.parent}.\n"
                             "The API is at /api/.", status_code=404)


def create_app(analysis, static_dir=STATIC_DIR, analyses=None):
    """A Starlette app over `analysis`, deriving its tables first if they are not current. `analyses`, label ->
    Analysis with `analysis` first, opens further Analyses of the same input for comparison (U5, `pt_analysis.compare`);
    every endpoint reads `analysis`, the primary, unless it takes an `analysis` parameter."""
    analyses = dict(analyses) if analyses else {analysis.label: analysis}
    for other in analyses.values():
        derive_analysis(other)
    routes = [route for module in ENDPOINTS for route in importlib.import_module(module).ROUTES]
    static_dir = Path(static_dir)
    if (static_dir / "index.html").is_file():
        routes.append(Mount("/", AppFiles(directory=static_dir, html=True)))
    else:
        routes.append(Route("/", not_built))
    app = Starlette(routes=routes, middleware=[Middleware(GZipMiddleware, minimum_size=1000)])
    app.state.analysis = analysis
    app.state.analyses = analyses
    return app


def serve(analysis, port=8765, open_browser=True, static_dir=STATIC_DIR, analyses=None):
    """Derive what is not current, then serve `analysis` (and `analyses`, as in `create_app`) on 127.0.0.1 until
    interrupted, opening the browser once the server listens unless `open_browser` is false."""
    import uvicorn

    app = create_app(analysis, static_dir, analyses)
    server = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=port, log_level="warning"))

    def announce():
        while not server.started and not server.should_exit:
            threading.Event().wait(0.05)
        if server.started:
            sockets = [s for srv in server.servers for s in srv.sockets]
            url = f"http://127.0.0.1:{sockets[0].getsockname()[1]}/"
            print(f"Serving {analysis.label} at {url} (Ctrl-C to stop)", flush=True)
            if open_browser:
                webbrowser.open(url)

    threading.Thread(target=announce, daemon=True).start()
    server.run()
