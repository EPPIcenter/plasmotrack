"""The flows, lags and spatial-scale endpoints of the functionals panel (U9, #87)."""

import pytest
from test_flows import E12, E13, E23, EPI, MEAN_Q, PSI_EDGE
from test_functionals_endpoint import client_for
from test_strata import open_analysis


@pytest.fixture
def client(tmp_path):
    return client_for(open_analysis(tmp_path, "both", EPI), tmp_path)


def test_flows_serve_attribution_enrichment_and_net_flow_per_stratum(client):
    got = client.get("/api/flows", params={"column": "district", "min_infections": "2"}).json()
    assert got["column"] == "district"
    assert got["strata"] == [{"value": "North", "n": 2}, {"value": "South", "n": 1}]
    north = [c for c in got["attribution"] if c["recipient"] == "North"]
    assert [(c["source"], c["unobserved"]) for c in north] == [("North", False), ("South", False), (None, True)]
    assert sum(c["share"] for c in north) == pytest.approx(1.0)
    cell = next(c for c in got["enrichment"] if (c["source"], c["recipient"]) == ("North", "North"))
    assert (cell["candidates"], cell["expected_edges"]) == (1, pytest.approx(E13, rel=1e-3))
    assert got["net_flow"] == {"min_infections": 2, "hidden": 1, "strata": [
        {"value": "North", "n": 2, "out_edges": pytest.approx(E12, rel=1e-2),
         "in_edges": pytest.approx(E23, rel=1e-2), "net": pytest.approx(E12 - E23, abs=1e-4)}]}


def test_flows_name_a_column_the_epi_data_lacks(client):
    response = client.get("/api/flows", params={"column": "parish"})
    assert response.status_code == 400 and "parish" in response.json()["error"]


def test_lags_serve_the_histogram_the_median_duration_and_no_reassignment_without_a_cutoff(client):
    got = client.get("/api/lags", params={"width": "10"}).json()
    assert got["width"] == 10
    assert sum(b["mass"] for b in got["bins"]) == pytest.approx(E12 + E13 + E23, rel=1e-3)
    assert got["total"] == pytest.approx(E12 + E13 + E23, rel=1e-3)
    assert got["median_duration"] == pytest.approx(14.02945, rel=1e-5)
    assert got["reassign"] is None


def test_lags_with_a_cutoff_report_the_moved_functionals_overall_and_per_stratum(client):
    got = client.get("/api/lags", params={"max_lag": "100", "column": "district"}).json()["reassign"]
    assert got["max_lag"] == 100
    assert got["moved_mass"] == pytest.approx(E13, rel=1e-3)
    assert got["overall"]["psi_edge"]["mean"] == pytest.approx(PSI_EDGE, rel=1e-5)
    assert got["overall"]["mean_q_moved"]["mean"] == pytest.approx(MEAN_Q + E13 / 3, rel=1e-3)
    assert [s["value"] for s in got["strata"]] == ["North", "South"]


@pytest.mark.parametrize("params", [{"max_lag": "soon"}, {"width": "-1"}, {"width": "wide"}])
def test_lags_reject_a_cutoff_or_width_that_is_not_a_number(client, params):
    assert client.get("/api/lags", params=params).status_code == 400


def test_spatial_serves_the_shares_by_relation(client):
    got = client.get("/api/spatial", params={"unit": "shehia", "parent": "district", "max_lag": "200"}).json()
    assert (got["unit"], got["parent"], got["max_lag"]) == ("shehia", "district", 200)
    assert [s["relation"] for s in got["shares"]] == ["same_unit", "same_parent", "other"]
    assert sum(s["share"] for s in got["shares"]) == pytest.approx(1.0)
    assert got["total"] == pytest.approx(E12 + E13 + E23, rel=1e-3)


@pytest.mark.parametrize("params", [{"unit": "shehia"}, {"unit": "shehia", "parent": "shehia"},
                                    {"unit": "parish", "parent": "district"}])
def test_spatial_needs_two_columns_the_epi_data_has(client, params):
    assert client.get("/api/spatial", params=params).status_code == 400
