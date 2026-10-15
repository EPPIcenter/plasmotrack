"""The local HTTP backend of `pt-analysis serve`: the built frontend as static files, and a JSON API over an Analysis.

Each module named in `ENDPOINTS` serves one view's data: it exports `ROUTES`, a list of Starlette routes under
`/api/`, whose handlers find the Analysis at `request.app.state.analysis`. A new view's backend is one module in
`pt_analysis/endpoints/` plus its line in its ticket's slot of `ENDPOINTS`.
"""

import importlib
import threading
import webbrowser
from pathlib import Path

from starlette.applications import Starlette
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
    # #83 U2 edge table
    "pt_analysis.endpoints.edges",  # U2: edge table

    # #84 U8 network view

    # #85 U3 pair inspector
    "pt_analysis.endpoints.pair",  # U3: order, edges, co-inclusion, genotypes and infection times of a pair

    # #86 #87 U9 functionals panel

    # #88 U10 infection card
    "pt_analysis.endpoints.infection",  # U10: one Infection's posteriors, genotype, Parent sets and edges

    # #89 #90 U6 truth overlay, ranking

    # #91 U7 cohort view

    # #92 U4 sampler panel
    "pt_analysis.endpoints.sampler",  # U4: acceptance, swaps, ladder, proposal variances, Run events

    # #94 U5 comparison of Analyses

]


def not_built(request):
    return PlainTextResponse(f"The frontend is not built: run `pnpm install && pnpm build` in {STATIC_DIR.parent}.\n"
                             "The API is at /api/.", status_code=404)


def create_app(analysis, static_dir=STATIC_DIR):
    """A Starlette app over `analysis`, deriving its tables first if they are not current."""
    derive_analysis(analysis)
    routes = [route for module in ENDPOINTS for route in importlib.import_module(module).ROUTES]
    static_dir = Path(static_dir)
    if (static_dir / "index.html").is_file():
        routes.append(Mount("/", StaticFiles(directory=static_dir, html=True)))
    else:
        routes.append(Route("/", not_built))
    app = Starlette(routes=routes, middleware=[Middleware(GZipMiddleware, minimum_size=1000)])
    app.state.analysis = analysis
    return app


def serve(analysis, port=8765, open_browser=True, static_dir=STATIC_DIR):
    """Derive what is not current, then serve `analysis` on 127.0.0.1 until interrupted, opening the browser once the
    server listens unless `open_browser` is false."""
    import uvicorn

    app = create_app(analysis, static_dir)
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
