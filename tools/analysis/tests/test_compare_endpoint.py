"""The comparison endpoints of `pt-analysis serve a.json b.json` (U5): the Analyses listed, their convergence side by
side, and the trace endpoint reading any of them. Fixture as in test_compare: the golden single-chain store alone as
Analysis "single" (primary), the replica-exchange store alone as "replica"."""

import pytest
from starlette.testclient import TestClient
from test_compare import one_run_analysis

from pt_analysis.compare import open_comparison
from pt_analysis.server import create_app

SINGLE_TIMINGS = {"chain1": {"exit": 0, "seed": 42, "sampling_seconds": 4.0, "samples": 50}}


def client_for(tmp_path, comparison):
    static = tmp_path / "dist"
    static.mkdir(exist_ok=True)
    (static / "index.html").write_text("<!doctype html><title>pt-analysis</title>")
    return TestClient(create_app(comparison.primary, static_dir=static, analyses=comparison.analyses))


@pytest.fixture
def client(tmp_path):
    comparison = open_comparison([one_run_analysis(tmp_path, "single-chain", "single", timings=SINGLE_TIMINGS),
                                  one_run_analysis(tmp_path, "replica-exchange", "replica")])
    return client_for(tmp_path, comparison)


def test_the_analyses_are_listed_primary_first_with_their_input_hash_and_timings(client):
    got = client.get("/api/analyses").json()
    assert got["input_sha256"] == "cd84e24b6713ef19bc3e27d4a45b9d579634262c808ee3530375e1ac409774ca"
    assert [(a["label"], a["primary"], a["draws"], a["sampling_seconds"], a["instructions_retired"])
            for a in got["analyses"]] == [("single", True, 10, 4.0, None), ("replica", False, 10, None, None)]
    assert got["analyses"][0]["chains"] == [{"chain": "chain1", "draws": 10}]


def test_existing_endpoints_read_the_primary(client):
    assert client.get("/api/analysis").json()["label"] == "single"
    rows = client.get("/api/convergence").json()["rows"]
    assert {r["chain"] for r in rows} == {"chain1", None}
    assert rows[0]["ess"] == pytest.approx(10.0)


def test_compared_convergence_puts_each_analysis_side_by_side(client):
    got = client.get("/api/compare/convergence").json()
    assert got["analyses"] == ["single", "replica"]
    assert got["timed"] == ["single"]
    assert got["instructions"] == []
    assert len(got["rows"]) == 42
    first = got["rows"][0]
    assert (first["group"], first["name"]) == ("log posterior", "log posterior")
    assert first["analyses"]["single"] == {"ess": pytest.approx(10.0), "rhat": None,
                                           "ess_per_second": pytest.approx(2.5)}
    assert first["analyses"]["replica"] == {"ess": pytest.approx(7.44641539328075), "rhat": None}


def test_a_trace_reads_the_named_analysis_and_defaults_to_the_primary(client):
    params = {"group": "log posterior", "name": "log posterior"}
    primary = client.get("/api/trace", params=params).json()
    single = client.get("/api/trace", params={**params, "analysis": "single"}).json()
    replica = client.get("/api/trace", params={**params, "analysis": "replica"}).json()
    assert primary["chains"] == single["chains"]
    assert "analysis" not in primary and single["analysis"] == "single"
    assert replica["analysis"] == "replica"
    assert replica["chains"][0]["value"] != single["chains"][0]["value"]
    ess = client.get("/api/range", params={**params, "analysis": "replica"}).json()["ess"]
    assert ess == pytest.approx(7.44641539328075)


def test_an_unknown_analysis_is_a_404_naming_it(client):
    response = client.get("/api/trace", params={"group": "log posterior", "name": "log posterior", "analysis": "nope"})
    assert response.status_code == 404
    assert "nope" in response.json()["error"]


def test_a_single_analysis_server_lists_itself(tmp_path):
    comparison = open_comparison([one_run_analysis(tmp_path, "single-chain", "single")])
    static = tmp_path / "dist"
    static.mkdir()
    (static / "index.html").write_text("")
    client = TestClient(create_app(comparison.primary, static_dir=static))
    got = client.get("/api/analyses").json()
    assert [a["label"] for a in got["analyses"]] == ["single"]
    assert client.get("/api/compare/convergence").json()["timed"] == []
