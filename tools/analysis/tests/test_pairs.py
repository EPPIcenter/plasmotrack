"""The Pairs workspace's list (`/api/pairs`) and the edge table's per-Draw series (`/api/edge-series`), on the golden
stores as Chains chain1 and chain2 of one Analysis."""

import pytest
from starlette.testclient import TestClient
from test_pair import swap_1_and_2

from pt_analysis.edge_series import edge_series
from pt_analysis.pairs import uncertain_pairs
from pt_analysis.server import create_app


def client_of(analysis):
    return TestClient(create_app(analysis, static_dir=analysis.path.parent / "no-dist"))


# --- the pair list ----------------------------------------------------------------------------------------------------

def test_every_golden_order_is_certain_so_no_pair_is_listed(golden_analysis):
    # Infection 1 precedes 2 and 3, and 2 precedes 3, at every Draw of both Chains.
    got = uncertain_pairs(golden_analysis)
    assert (got["rows"], got["total"]) == ([], 0)


def test_a_pair_whose_order_flips_is_listed_with_its_order_per_chain_and_both_edge_marginals(golden_analysis):
    # chain1 has 2 before 1 at two of its ten Draws: P(1 before 2) is 0.8 there, 1 in chain2, 0.9 pooled.
    swap_1_and_2(golden_analysis)
    got = uncertain_pairs(golden_analysis, lo=0.1, hi=0.9)
    assert got["total"] == 1
    [row] = got["rows"]
    assert (row["a"], row["b"], row["draws"]) == ("1", "2", 20)
    assert row["probability"] == pytest.approx(0.9)
    assert row["chains"] == [
        {"chain": "chain1", "draws": 10, "probability": pytest.approx(0.8), "switches": 2},
        {"chain": "chain2", "draws": 10, "probability": 1.0, "switches": 0}]
    assert (row["switches"], row["spread"]) == (2, pytest.approx(0.2))
    # Edge 1 -> 2 (2's parent 1) has pooled marginal 0.0029 in the golden edge table; 2 -> 1 has no row.
    assert row["a_to_b"] == pytest.approx(0.002910808444662265)
    assert row["b_to_a"] == 0


def test_a_pair_outside_the_band_is_listed_when_its_chains_disagree(golden_analysis):
    swap_1_and_2(golden_analysis)
    assert uncertain_pairs(golden_analysis, lo=0.1, hi=0.85)["total"] == 0
    assert uncertain_pairs(golden_analysis, lo=0.1, hi=0.85, disagree=0.2)["total"] == 1


def test_pairs_without_an_edge_marginal_of_edge_min_either_way_are_left_out(golden_analysis):
    swap_1_and_2(golden_analysis)
    assert uncertain_pairs(golden_analysis, edge_min=0.002)["total"] == 1
    assert uncertain_pairs(golden_analysis, edge_min=0.01)["total"] == 0


def test_the_endpoint_serves_the_list_and_refuses_a_bad_band(golden_analysis):
    swap_1_and_2(golden_analysis)
    client = client_of(golden_analysis)
    body = client.get("/api/pairs?min=0.1&max=0.9&limit=5").json()
    assert body["chains"] == ["chain1", "chain2"]
    assert (body["min"], body["max"], body["limit"], body["total"]) == (0.1, 0.9, 5, 1)
    assert [(r["a"], r["b"]) for r in body["rows"]] == [("1", "2")]
    assert client.get("/api/pairs?limit=0").json()["rows"] == []
    assert client.get("/api/pairs?min=0.9&max=0.1").status_code == 400
    assert client.get("/api/pairs?min=x").status_code == 400
    assert client.get("/api/pairs?limit=-1").status_code == 400


# --- the per-Draw series of edges -------------------------------------------------------------------------------------

def test_an_edges_series_is_its_weight_per_draw_per_chain(golden_analysis):
    got = edge_series(golden_analysis, [("3", "1")])
    assert got["chains"] == ["chain1", "chain2"]
    [series] = got["series"]
    assert (series["child"], series["parent"]) == ("3", "1")
    chain1, chain2 = series["chains"]
    assert len(chain1) == len(chain2) == 10
    assert chain1[0] == pytest.approx(0.9638404138595045)
    assert chain2[:2] == [pytest.approx(0.003051631690787871), pytest.approx(0.003931889022720117)]


def test_a_draw_without_the_edge_counts_as_zero_and_an_unknown_edge_is_all_zero(golden_analysis):
    [series] = edge_series(golden_analysis, [("1", "3")])["series"]
    assert series["chains"] == [[0.0] * 10, [0.0] * 10]


def test_points_bins_each_series_into_means_of_adjacent_draws(golden_analysis):
    full = edge_series(golden_analysis, [("3", "1")])["series"][0]["chains"][0]
    binned = edge_series(golden_analysis, [("3", "1")], points=5)["series"][0]["chains"][0]
    assert binned == [pytest.approx((full[k] + full[k + 1]) / 2) for k in range(0, 10, 2)]


def test_the_series_endpoint_takes_edges_as_child_from_parent_keys(golden_analysis):
    client = client_of(golden_analysis)
    body = client.get("/api/edge-series", params={"edges": "3<-1,2<-1", "points": "4"}).json()
    assert [(s["child"], s["parent"]) for s in body["series"]] == [("3", "1"), ("2", "1")]
    assert all(len(c) == 4 for s in body["series"] for c in s["chains"])
    assert client.get("/api/edge-series", params={"edges": "3"}).status_code == 400
    assert client.get("/api/edge-series").json()["series"] == []


@pytest.mark.parametrize("path", ["/api/pairs", "/api/edge-series"])
def test_both_read_the_named_analysis_and_404_an_unknown_one(golden_analysis, path):
    client = client_of(golden_analysis)
    response = client.get(path, params={"analysis": "nope"})
    assert response.status_code == 404
    assert "nope" in response.json()["error"]
    assert client.get(path, params={"analysis": "golden"}).status_code == 200
