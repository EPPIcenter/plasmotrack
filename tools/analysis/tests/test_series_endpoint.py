"""`/api/series`: several quantities' per-Draw values per Chain at once, compact, for the convergence table's
sparklines, on the golden stores as Chains chain1 and chain2 of one Analysis."""

import pytest
from starlette.testclient import TestClient

from pt_analysis.server import create_app


@pytest.fixture
def client(golden_analysis):
    return TestClient(create_app(golden_analysis, static_dir=golden_analysis.path.parent / "no-dist"))


def sampling_values(client, group, name):
    """Each Chain's Draw values of a quantity as the trace panel gives them."""
    chains = client.get("/api/trace", params={"group": group, "name": name}).json()["chains"]
    return [c["value"] for c in chains]


def test_each_quantity_has_its_draw_values_per_chain_in_request_order(client):
    body = client.get("/api/series", params=[("q", "log posterior:log posterior"), ("q", "edge:3<-1"),
                                             ("q", "unobserved-source inclusion:3")]).json()
    assert body["chains"] == ["chain1", "chain2"]
    assert body["points"] is None
    assert [(s["group"], s["name"]) for s in body["series"]] == [
        ("log posterior", "log posterior"), ("edge", "3<-1"), ("unobserved-source inclusion", "3")]
    for s in body["series"]:
        assert s["chains"] == sampling_values(client, s["group"], s["name"])
    # The golden Chains have ten Draws each; Edge 3 <- 1 is 0.96384 at chain1's first.
    assert [len(c) for c in body["series"][1]["chains"]] == [10, 10]
    assert body["series"][1]["chains"][0][0] == pytest.approx(0.9638404138595045)


def test_points_bins_each_chain_into_means_of_adjacent_draws(client):
    full = sampling_values(client, "edge", "3<-1")[0]
    body = client.get("/api/series", params={"q": "edge:3<-1", "points": "5"}).json()
    assert body["points"] == 5
    assert body["series"][0]["chains"][0] == [pytest.approx((full[k] + full[k + 1]) / 2) for k in range(0, 10, 2)]


def test_a_name_is_split_at_the_first_colon_and_an_unknown_quantity_is_null(client):
    body = client.get("/api/series", params=[("q", "functionals:mean_q"), ("q", "edge:9<-1"),
                                             ("q", "no such group:x"), ("q", "infection duration:nobody")]).json()
    assert body["series"][0]["chains"] == sampling_values(client, "functionals", "mean_q")
    assert [s["chains"] for s in body["series"][1:]] == [None, None, None]
    assert [(s["group"], s["name"]) for s in body["series"][1:]] == [
        ("edge", "9<-1"), ("no such group", "x"), ("infection duration", "nobody")]


def test_bad_requests_are_refused(client):
    assert client.get("/api/series").json()["series"] == []
    assert client.get("/api/series", params={"q": "no colon"}).status_code == 400
    assert client.get("/api/series", params={"q": "edge:3<-1", "points": "0"}).status_code == 400
    assert client.get("/api/series", params={"q": "edge:3<-1", "points": "x"}).status_code == 400
    too_many = [("q", "edge:3<-1")] * 201
    assert client.get("/api/series", params=too_many).status_code == 400
    response = client.get("/api/series", params={"q": "edge:3<-1", "analysis": "nope"})
    assert response.status_code == 404
