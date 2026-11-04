"""The cohort view (U7): a cohort of fit summaries (`pt_analysis.cohort`) and click-through to one fit's Analysis.

    GET  /api/cohort[?file=PATH][&bins=10]
                            {file, n, columns: [{name, kind}], data: {column: [values]},
                             calibration: {bins, row, bin, n, n_true, score_sum} | null}
                            the cohort the manifest's `cohort` names, or `file` (relative to the manifest's
                            directory); 404 without one
    POST /api/cohort/open?row=N[&file=PATH]
                            {url, manifest}: serves the Analysis of row N's `manifest` with `pt-analysis serve` on
                            a free port (once per manifest) and returns its address

The calibration counts are per cohort row and score bin, so the browser sums them over any selection of fits. A
fit's Analysis is served by a child process that lives as long as this server; only manifests the cohort names can
be opened.
"""

import atexit
import re
import subprocess
import sys
import threading

from starlette.concurrency import run_in_threadpool
from starlette.routing import Route

from pt_analysis import cohort
from pt_analysis.endpoints import error, json

ADDRESS = re.compile(r" at (http://127\.0\.0\.1:[0-9]+/)")


class Servers:
    """Child `pt-analysis serve` processes by manifest path."""

    def __init__(self):
        self.lock = threading.Lock()
        self.running = {}  # manifest -> (process, url)
        atexit.register(self.close)

    def url(self, manifest):
        with self.lock:
            found = self.running.get(manifest)
            if found and found[0].poll() is None:
                return found[1]
            process = subprocess.Popen(
                [sys.executable, "-m", "pt_analysis.cli", "serve", str(manifest), "--port", "0", "--no-browser"],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            said = []
            for line in process.stdout:
                match = ADDRESS.search(line)
                if match:
                    # Keep reading so that the child never blocks on a full pipe.
                    threading.Thread(target=process.stdout.read, daemon=True).start()
                    self.running[manifest] = (process, match.group(1))
                    return match.group(1)
                said.append(line)
            process.wait()
            raise RuntimeError(f"pt-analysis serve {manifest} stopped: {''.join(said).strip()}")

    def close(self):
        with self.lock:
            for process, _ in self.running.values():
                if process.poll() is None:
                    process.terminate()
            for process, _ in self.running.values():
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
            self.running.clear()


def cohort_file(request):
    """The cohort file of the request (`file`, relative to the manifest's directory) or the manifest's; None."""
    analysis = request.app.state.analysis
    given = request.query_params.get("file")
    if given:
        return (analysis.path.parent / given).resolve()
    return analysis.cohort


def servers(request):
    state = request.app.state
    if getattr(state, "cohort_servers", None) is None:
        state.cohort_servers = Servers()
    return state.cohort_servers


async def get_cohort(request):
    path = cohort_file(request)
    if path is None:
        return error("no cohort: the manifest names no 'cohort' file and the request gives no 'file'", 404)
    try:
        bins = int(request.query_params.get("bins", cohort.BINS))
    except ValueError:
        return error(f"bins is a whole number, not '{request.query_params['bins']}'")
    if not 1 <= bins <= 1000:
        return error(f"bins must be in 1..1000, not {bins}")
    if not path.is_file():
        return error(f"no cohort file at {path}", 404)
    try:
        return json(await run_in_threadpool(cohort.load, path, bins))
    except cohort.CohortError as e:
        return error(str(e))


async def open_fit(request):
    path = cohort_file(request)
    if path is None or not path.is_file():
        return error("no cohort file to open a fit from", 404)
    try:
        row = int(request.query_params.get("row", ""))
    except ValueError:
        return error(f"row is a cohort row number, not '{request.query_params.get('row', '')}'")
    manifests = await run_in_threadpool(cohort.manifests, path)
    if not 0 <= row < len(manifests) or manifests[row] is None:
        return error(f"cohort row {row} names no Analysis manifest")
    manifest = (path.parent / manifests[row]).resolve()
    if not manifest.is_file():
        return error(f"cohort row {row}'s manifest is not on disk: {manifest}", 404)
    try:
        url = await run_in_threadpool(servers(request).url, manifest)
    except RuntimeError as e:
        return error(str(e), 500)
    return json({"url": url, "manifest": str(manifest)})


ROUTES = [Route("/api/cohort", get_cohort), Route("/api/cohort/open", open_fit, methods=["POST"])]
