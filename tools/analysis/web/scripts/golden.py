"""The golden Analysis for the web app's tests: its endpoints' JSON as fixtures, or served for the Playwright smoke test.

The golden two-Chain Analysis (`tests/conftest.py`, `golden_analysis`) is copied into a temporary directory with an epi
CSV, a truth file and a cohort of three fits of it (full, partial and null; `cohort`), so every endpoint has data. Run from `tools/analysis/`:

    uv run python web/scripts/golden.py fixtures web/src/test/fixtures   # write <name>.json per endpoint
    uv run python web/scripts/golden.py serve --port 8899                 # serve it (and web/dist/) until stopped
"""

import argparse
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "tests"))

from conftest import copy_golden, open_golden, write_manifest  # noqa: E402

# Infections 1 to 3 in two districts of one region.
EPI = "infection,district,region,observed\n1,North,East,2020-01-05\n2,South,East,2020-02-01\n3,North,East,2020-03-10\n"
TRUTH = {"network": [{"src": "1", "dest": "3"}, {"src": "1", "dest": "2"}],
         "nodes": [{"id": i, "infection_time": 0.0, "flat_latent_genotype": [], "observation_time": 10.0}
                   for i in ("1", "2", "3")]}

# Fixture name -> request. One per endpoint the app reads, with the parameters its tests use.
ENDPOINTS = {
    "analysis": "/api/analysis",
    "analyses": "/api/analyses",
    "status": "/api/status",
    "overview": "/api/overview",
    "convergence": "/api/convergence",
    "trace-log-posterior": "/api/trace?group=log+posterior&name=log+posterior",
    "trace-log-posterior-burn-in": "/api/trace?group=log+posterior&name=log+posterior&burn_in=true",
    "trace-edge": "/api/trace?group=edge&name=3%3C-1",
    "range-log-posterior": "/api/range?group=log+posterior&name=log+posterior",
    "infection-3": "/api/infection?id=3",
    "edges": "/api/edges",
    "network": "/api/network",
    "pair-1-3": "/api/pair?a=1&b=3",
    "epi": "/api/epi",
    "functionals": "/api/functionals",
    "functionals-district": "/api/functionals?column=district",
    "flows-district": "/api/flows?column=district",
    "lags": "/api/lags?column=district&max_lag=30",
    "spatial": "/api/spatial?unit=district&parent=region",
    "sampler": "/api/sampler",
    "sampler-chain2": "/api/sampler?chain=chain2",
    "sampler-type": "/api/sampler?chain=chain1&type=ConstrainedContinuousRandomWalk&limit=5",
    "truth-network": "/api/truth-network",
    "errors": "/api/errors",
    "ranking": "/api/ranking",
    "parent-sets": "/api/parent-sets",
    "compare-convergence": "/api/compare/convergence",
    "pairs": "/api/pairs?min=0&max=1",
    "edge-series": "/api/edge-series?edges=1%3C-1_copy,2%3C-1,2%3C-2_copy,3%3C-1,3%3C-2,3%3C-3_copy&points=40",
    "cohort": "/api/cohort",
}

# The cohort's fits: (fit_type, missingness_condition), each a copy of the golden Analysis with these factors.
COHORT_FITS = [("posterior", "full"), ("posterior", "missing_50"), ("null_prior", "full")]
FACTORS = {"run_key": "biological__fp0.05__sp0__r00.5__coi1__high__sim0__s0", "simulation_model": "biological",
           "total_nodes": 3, "founder_probability": 0.05, "superinfection_probability": 0, "secondary_case_rate": 0.5,
           "founder_coi_rate": 1, "asymptomatic_probability": 0.5, "panel_name": "high",
           "panel_alleles_per_locus": 20, "simulation_index": 0, "seed": 0}


def cohort(directory):
    """A cohort of the golden Analysis's fit summaries, one per `COHORT_FITS`, with per-edge rows: `cohort.parquet`."""
    from pt_analysis import cohort as cohorts
    from pt_analysis import fit_summary

    rows = directory / "rows"
    rows.mkdir()
    for fit_type, missing in COHORT_FITS:
        fit = directory / f"{fit_type}-{missing}"
        fit.mkdir()
        (fit / "truth.json").write_text(json.dumps(TRUTH))
        analysis = open_golden(fit, label=fit.name, truth="truth.json",
                               factors=FACTORS | {"fit_type": fit_type, "missingness_condition": missing})
        fit_summary.write_row(fit_summary.summarize(analysis), rows / f"{fit.name}.parquet", manifest=analysis.path)
        fit_summary.write_edges(analysis, cohorts.edges_path(rows / f"{fit.name}.parquet"))
    cohorts.concatenate([rows], directory / "cohort.parquet")


# Fixture name -> request, against two Analyses open side by side (`pt-analysis serve single.json replica.json`): the
# golden single-chain store alone as "single" (primary, with timings) and the replica-exchange store alone as
# "replica", as in tests/test_compare_endpoint.py.
COMPARE_ENDPOINTS = {
    "compare-analyses": "/api/analyses",
    "compare-two": "/api/compare/convergence",
}
SINGLE_TIMINGS = {"chain1": {"sampling_seconds": 4.0, "instructions_retired": 2.0e9}}


def golden(directory):
    """The golden Analysis with epi data, truth and a cohort, its copies in `directory`."""
    directory = Path(directory).resolve()
    (directory / "epi.csv").write_text(EPI)
    (directory / "truth.json").write_text(json.dumps(TRUTH))
    cohort(directory)
    return open_golden(directory, label="golden", epi="epi.csv", truth="truth.json", cohort="cohort.parquet")


def comparison(directory):
    """The golden stores as two one-Chain Analyses of one input, "single" (primary) and "replica"."""
    from pt_analysis.compare import open_comparison

    directory = Path(directory).resolve()
    copy_golden(directory)
    manifests = [write_manifest(directory / "single.json", label="single", runs={"chain1": "single-chain"},
                                timings=SINGLE_TIMINGS),
                 write_manifest(directory / "replica.json", label="replica", runs={"chain1": "replica-exchange"})]
    return open_comparison(manifests)


def capture(client, endpoints, out):
    for name, path in endpoints.items():
        response = client.get(path)
        if response.status_code != 200:
            raise SystemExit(f"{path}: {response.status_code} {response.text}")
        body = response.json()
        if "analyses" in body and isinstance(body["analyses"][0], dict):  # temporary paths change on every capture
            for a in body["analyses"]:
                a["path"] = f"golden/{a['label'] if a['label'] != 'golden' else 'analysis'}.json"
        if name == "cohort":
            body["file"] = "golden/cohort.parquet"
        (out / f"{name}.json").write_text(json.dumps(body, indent=1, sort_keys=True) + "\n")
        print(f"{name}.json  {path}")


def fixtures(out):
    from starlette.testclient import TestClient

    from pt_analysis.server import create_app

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        analysis = golden(tmp)
        capture(TestClient(create_app(analysis, static_dir=Path(tmp) / "no-dist")), ENDPOINTS, out)
    with tempfile.TemporaryDirectory() as tmp:
        compared = comparison(tmp)
        app = create_app(compared.primary, static_dir=Path(tmp) / "no-dist", analyses=compared.analyses)
        capture(TestClient(app), COMPARE_ENDPOINTS, out)


def serve(port):
    from pt_analysis.server import serve as run

    with tempfile.TemporaryDirectory() as tmp:
        run(golden(tmp), port=port, open_browser=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("fixtures").add_argument("out")
    commands.add_parser("serve").add_argument("--port", type=int, default=8899)
    args = parser.parse_args()
    if args.command == "fixtures":
        fixtures(args.out)
    else:
        serve(args.port)


if __name__ == "__main__":
    main()
