"""The functionals panel's endpoints (U9): the epi columns, functionals overall and per stratum, and the per-Draw
series of a stratum's functional through the trace panel's endpoints."""

import pytest
from starlette.testclient import TestClient
from test_strata import PSI_EDGE, open_analysis

from pt_analysis.server import create_app


def client_for(analysis, tmp_path):
    static = tmp_path / "dist"
    static.mkdir(exist_ok=True)
    (static / "index.html").write_text("<!doctype html><title>pt-analysis</title>")
    return TestClient(create_app(analysis, static_dir=static))


@pytest.fixture
def client(tmp_path):
    return client_for(open_analysis(tmp_path), tmp_path)


def test_epi_lists_its_columns_with_types_and_the_ids_that_do_not_match(client):
    got = client.get("/api/epi").json()
    assert got["key"] == "infection"
    assert got["columns"] == [{"name": "district", "type": "categorical"}, {"name": "observed", "type": "date"}]
    assert (got["unmatched"], got["missing"], got["infections"]) == (["9"], [], 3)


def test_without_epi_data_there_are_no_columns(tmp_path):
    client = client_for(open_analysis(tmp_path, epi=None), tmp_path)
    assert client.get("/api/epi").json() == {"key": None, "columns": [], "unmatched": [], "missing": [],
                                            "infections": 0}


def test_functionals_overall(client):
    got = client.get("/api/functionals").json()
    assert (got["column"], got["level"]) == (None, 0.8)
    assert got["classes"] == ["observed_1", "observed_1_unobserved", "observed_2", "observed_2_unobserved",
                              "unobserved_only"]
    (overall,) = got["strata"]
    assert (overall["value"], overall["n"]) == (None, 3)
    assert overall["psi_map"] == pytest.approx(1 / 3)
    assert overall["size_range"] == pytest.approx([1 / 3, 1 / 3])
    assert overall["psi_edge"]["mean"] == pytest.approx((0.299024 + 0.260587) / 2, rel=1e-5)
    assert set(overall["psi_edge"]) == {"mean", "lower", "upper"}
    assert set(overall) >= {"mean_q", "classes", "multi_source_augmented", "multi_source_observed"}
    assert list(overall["classes"]) == got["classes"]


def test_functionals_per_stratum_and_a_missing_column(client):
    got = client.get("/api/functionals", params={"column": "district"}).json()
    assert [(s["value"], s["n"]) for s in got["strata"]] == [("North", 2), ("South", 1)]
    response = client.get("/api/functionals", params={"column": "travel"})
    assert response.status_code == 400
    assert "no column 'travel'" in response.json()["error"]


def test_the_trace_of_a_stratums_functional(client):
    got = client.get("/api/trace", params={"group": "functionals", "name": "psi_edge",
                                           "stratum_column": "", "burn_in": "true"}).json()
    assert got["burn_in"] is False
    for chain, expected in zip(got["chains"], PSI_EDGE):
        assert chain["value"] == pytest.approx(expected, rel=2e-3)
        assert set(chain["phase"]) == {"sampling"}
    south = client.get("/api/trace", params={"group": "functionals", "name": "mean_q", "stratum_column": "district",
                                             "stratum": "South"}).json()
    assert [c["value"] for c in south["chains"]] == [pytest.approx([1.0] * 10, abs=1e-12)] * 2
    assert (south["stratum_column"], south["stratum"]) == ("district", "South")


def test_the_range_of_a_stratums_functional(client):
    whole = client.get("/api/range", params={"group": "functionals", "name": "mean_unobserved_only",
                                             "stratum_column": "district", "stratum": "North"}).json()
    assert [c["draws"] for c in whole["chains"]] == [10, 10]
    some = client.get("/api/range", params={"group": "functionals", "name": "mean_unobserved_only",
                                            "stratum_column": "district", "stratum": "North", "first": 50}).json()
    assert [c["draws"] for c in some["chains"]] == [5, 5]


def test_a_stratum_of_another_group_or_an_unknown_value_is_an_error(client):
    for params in ({"group": "edge", "name": "3<-1", "stratum_column": "district", "stratum": "North"},
                   {"group": "functionals", "name": "psi_edge", "stratum_column": "district", "stratum": "East"},
                   {"group": "functionals", "name": "psi_edge", "stratum_column": "travel", "stratum": "x"}):
        assert client.get("/api/trace", params=params).status_code == 400

