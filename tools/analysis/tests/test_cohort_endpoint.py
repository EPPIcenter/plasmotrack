"""The cohort view's endpoints (U7): the cohort the manifest names or `?file=` gives, and opening a fit's Analysis."""

import json
import urllib.request

import pytest
from test_cohort import rows  # noqa: F401  (fixture)
from test_metrics import open_analysis, stores  # noqa: F401  (fixture)
from test_ranking_endpoint import client_for

from pt_analysis import cohort


@pytest.fixture
def cohort_file(rows, tmp_path):  # noqa: F811
    out = tmp_path / "cohorts" / "cohort.parquet"
    cohort.concatenate([rows], out)
    return out


def test_the_manifests_cohort_is_served(tmp_path, stores, cohort_file):  # noqa: F811
    analysis = open_analysis(tmp_path, stores, name="viewer.json", cohort="cohorts/cohort.parquet")
    got = client_for(analysis, tmp_path).get("/api/cohort").json()
    assert got["n"] == 2
    assert got["data"]["manifest"] == ["../null_prior.json", "../posterior.json"]
    assert got["calibration"]["bins"] == 10


def test_a_file_query_overrides_the_manifest_and_sets_the_bins(tmp_path, stores, cohort_file):  # noqa: F811
    client = client_for(open_analysis(tmp_path, stores, name="viewer.json"), tmp_path)
    got = client.get("/api/cohort", params={"file": str(cohort_file), "bins": "4"}).json()
    assert got["n"] == 2 and got["calibration"]["bins"] == 4


def test_without_a_cohort_it_is_not_found(tmp_path, stores):  # noqa: F811
    response = client_for(open_analysis(tmp_path, stores, name="viewer.json"), tmp_path).get("/api/cohort")
    assert response.status_code == 404
    assert "cohort" in response.json()["error"]


def test_a_missing_file_is_an_error(tmp_path, stores):  # noqa: F811
    client = client_for(open_analysis(tmp_path, stores, name="viewer.json"), tmp_path)
    response = client.get("/api/cohort", params={"file": str(tmp_path / "none.parquet")})
    assert response.status_code == 404


def test_open_serves_a_fits_analysis_and_returns_its_address(tmp_path, stores, cohort_file):  # noqa: F811
    analysis = open_analysis(tmp_path, stores, name="viewer.json", cohort="cohorts/cohort.parquet")
    client = client_for(analysis, tmp_path)
    try:
        got = client.post("/api/cohort/open", params={"row": "1"}).json()
        assert got["url"].startswith("http://127.0.0.1:")
        with urllib.request.urlopen(got["url"] + "api/analysis") as response:
            assert json.load(response)["label"] == "golden"
        # Opening the same fit again reuses its server.
        assert client.post("/api/cohort/open", params={"row": "1"}).json()["url"] == got["url"]
    finally:
        client.app.state.cohort_servers.close()


def test_open_needs_a_row_with_a_manifest(tmp_path, stores, cohort_file):  # noqa: F811
    client = client_for(open_analysis(tmp_path, stores, name="viewer.json", cohort="cohorts/cohort.parquet"),
                        tmp_path)
    assert client.post("/api/cohort/open", params={"row": "7"}).status_code == 400
    assert client.post("/api/cohort/open", params={"row": "x"}).status_code == 400
