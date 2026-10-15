"""The network view's payload (U8) on the golden Analysis: the two golden stores as Chains chain1 and chain2 of the
example input, whose three Infections are observed at 234.0737, 306.2297 and 350.5."""


import pytest
from starlette.testclient import TestClient

from pt_analysis.network import network
from pt_analysis.server import create_app


# Per-Draw Edge marginal 1 -> 3 and pooled q, computed with R arrow from the golden stores (exploration notes).
EDGE_3_FROM_1 = ([0.9638, 0.8808, 0.9799, 0.9337, 0.8457, 0.8624, 0.9414, 0.8481, 0.8313, 0.8267],
                 [0.0031, 0.0039, 0.9757, 0.9767, 0.9533, 0.9499, 0.9657, 0.9705, 0.9718, 0.9864])
Q = {"1": (1.0, 1.0), "2": (1.0, 1.0), "3": (0.1125, 0.2262)}

# Observation time minus the 97.5%, 50% and 2.5% quantiles (R type 7) of the 20 pooled `infection_duration` Draws,
# computed with Python's statistics.quantiles(method="inclusive") from the `infections` Trace.
INFECTION_TIME = {"1": (210.5691585473053, 223.2783343437396, 230.3022980395158),
                  "2": (281.18984163787144, 292.3976433463141, 301.20884910658725),
                  "3": (300.50056297780634, 328.20063946023583, 346.29913212204036)}


def mean(values):
    return sum(values) / len(values)


def test_nodes_are_the_input_infections_with_observation_time_symptomatic_flag_and_pooled_q(golden_analysis):
    nodes = network(golden_analysis)["nodes"]
    assert [n["id"] for n in nodes] == ["1", "2", "3"]
    assert [n["observation_time"] for n in nodes] == [234.0737, 306.2297, 350.5]
    assert [n["symptomatic"] for n in nodes] == [True, True, False]
    for n in nodes:
        assert n["q"] == pytest.approx(mean(Q[n["id"]]), abs=1e-4)


def test_infection_time_is_observation_time_minus_pooled_duration_quantiles(golden_analysis):
    for n in network(golden_analysis)["nodes"]:
        lower, median, upper = INFECTION_TIME[n["id"]]
        assert n["infection_time"] == {"lower": pytest.approx(lower), "median": pytest.approx(median),
                                       "upper": pytest.approx(upper)}


def test_at_the_default_floor_the_only_edge_is_1_to_3_with_its_pooled_and_per_chain_marginals_and_lag(golden_analysis):
    payload = network(golden_analysis)
    assert payload["chains"] == ["chain1", "chain2"]
    assert payload["floor"] == 0.01
    [edge] = payload["edges"]
    assert (edge["child"], edge["parent"]) == ("3", "1")
    assert edge["marginal"] == pytest.approx(mean(EDGE_3_FROM_1[0] + EDGE_3_FROM_1[1]), abs=1e-4)
    assert edge["chains"] == [pytest.approx(mean(EDGE_3_FROM_1[0]), abs=1e-4),
                              pytest.approx(mean(EDGE_3_FROM_1[1]), abs=1e-4)]
    assert edge["lag"] == pytest.approx(350.5 - 234.0737)


def test_a_lower_floor_admits_the_weak_edges_but_never_the_unobserved_source(golden_analysis):
    edges = network(golden_analysis, floor=0.001)["edges"]
    # Pooled 2 -> 3 is about 0.003 and 1 -> 2 about 0.0029 (exploration notes); both are below the default floor.
    assert {(e["child"], e["parent"]) for e in edges} == {("3", "1"), ("3", "2"), ("2", "1")}
    edge = next(e for e in edges if (e["child"], e["parent"]) == ("2", "1"))
    assert edge["marginal"] == pytest.approx((0.00267 + 0.00315) / 2, abs=1e-4)
    assert edge["lag"] == pytest.approx(306.2297 - 234.0737)
    assert edges == sorted(edges, key=lambda e: -e["marginal"])


@pytest.fixture
def client(golden_analysis, tmp_path):
    return TestClient(create_app(golden_analysis, static_dir=tmp_path / "no-dist"))


def test_the_endpoint_serves_the_payload_with_the_floor_as_a_query_parameter(client, golden_analysis):
    assert client.get("/api/network").json() == network(golden_analysis)
    got = client.get("/api/network", params={"floor": "0.001"}).json()
    assert got["floor"] == 0.001
    assert len(got["edges"]) == 3


@pytest.mark.parametrize("floor", ["x", "-0.1", "1.5", "nan"])
def test_a_floor_outside_0_to_1_is_refused(client, floor):
    response = client.get("/api/network", params={"floor": floor})
    assert response.status_code == 400
    assert "floor" in response.json()["error"]
