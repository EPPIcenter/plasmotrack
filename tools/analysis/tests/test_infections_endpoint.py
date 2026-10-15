"""The Infections workspace's table (`/api/infections`): one columnar payload with a row per Infection, on the two
golden stores as Chains chain1 and chain2 of one Analysis. Expected values come from the infection card, the
convergence table, the edge table and the network view, which compute them their own way."""

import math

import pytest
from conftest import open_golden
from starlette.testclient import TestClient

from pt_analysis.server import create_app


def client_of(analysis):
    return TestClient(create_app(analysis, static_dir=analysis.path.parent / "no-dist"))


@pytest.fixture
def client(golden_analysis):
    return client_of(golden_analysis)


def number(x):
    return float(x) if isinstance(x, str) else x


def rows(payload):
    """The columnar payload as {Infection: {column: value}}."""
    columns = payload["columns"]
    return {infection: {name: values[k] for name, values in columns.items()}
            for k, infection in enumerate(columns["infection"])}


@pytest.fixture
def table(client):
    return rows(client.get("/api/infections").json())


def test_one_row_per_infection_in_input_order_with_every_column_the_same_length(client):
    payload = client.get("/api/infections").json()
    assert payload["chains"] == ["chain1", "chain2"]
    assert payload["columns"]["infection"] == [n["id"] for n in client.get("/api/network").json()["nodes"]]
    assert {len(v) for v in payload["columns"].values()} == {3}
    assert all(len(v) == 3 for v in payload["q_chains"].values())
    assert payload["epi_columns"] == []


def test_observation_and_posterior_infection_times_match_the_network_view(client, table):
    for node in client.get("/api/network").json()["nodes"]:
        row = table[node["id"]]
        assert row["observation_time"] == node["observation_time"]
        assert row["symptomatic"] == node["symptomatic"]
        assert (row["infection_time_lower"], row["infection_time_median"], row["infection_time_upper"]) == \
            pytest.approx((node["infection_time"]["lower"], node["infection_time"]["median"],
                           node["infection_time"]["upper"]))
    assert table["3"]["observation_time"] == 350.5


def test_q_pooled_and_per_chain_with_its_ess_and_r_hat_from_the_convergence_table(client):
    payload = client.get("/api/infections").json()
    table = rows(payload)
    card = client.get("/api/infection?id=3").json()["summary"]
    assert table["3"]["q"] == pytest.approx(card["q"])
    k = payload["columns"]["infection"].index("3")
    assert {c: v[k] for c, v in payload["q_chains"].items()} == {
        c["chain"]: pytest.approx(c["q"]) for c in card["chains"]}
    convergence = {(r["group"], r["name"]): r for r in client.get("/api/convergence").json()["rows"]
                   if r["chain"] is None}
    q = convergence[("unobserved-source inclusion", "3")]
    assert number(table["3"]["q_ess"]) == pytest.approx(number(q["ess"]))
    assert number(table["3"]["q_rhat"]) == pytest.approx(number(q["rhat"]))


def test_expected_observed_parents_and_per_draw_entropy_and_effective_size_are_the_cards_pooled_means(client, table):
    for infection in ("1", "2", "3"):
        card = client.get(f"/api/infection?id={infection}").json()["summary"]
        assert table[infection]["expected_observed_parents"] == pytest.approx(card["expected_observed_parents"])
        assert table[infection]["entropy"] == pytest.approx(card["entropy"])
        assert table[infection]["effective_size"] == pytest.approx(card["effective_size"])


def test_the_effective_size_is_exp_of_the_entropy_of_the_pooled_parent_set_weights(client, table):
    for infection in ("1", "2", "3"):
        weights = [s["weight"] for s in client.get(f"/api/infection?id={infection}").json()["parent_sets"]["pooled"]]
        entropy = -sum(w * math.log(w) for w in weights if w > 0)
        assert table[infection]["posterior_entropy"] == pytest.approx(entropy)
        assert table[infection]["posterior_effective_size"] == pytest.approx(math.exp(entropy))
    # Infection 1 has one Parent set, the unobserved source alone.
    assert table["1"]["posterior_effective_size"] == pytest.approx(1.0)
    assert table["3"]["posterior_effective_size"] > 1.5


def test_the_top_parent_set_is_the_pooled_map_set_with_its_members_class_and_probability(table):
    # Infection 3's sets {1} weigh 0.885 and 0.771 in the two Chains (test_infection), so 0.828 pooled.
    assert table["3"]["top_set"] == ["1"]
    assert table["3"]["top_set_class"] == "observed_1"
    assert table["3"]["top_set_probability"] == pytest.approx((0.885 + 0.771) / 2, rel=1e-2)
    assert table["1"]["top_set"] == ["1_copy"]
    assert table["1"]["top_set_class"] == "unobserved_only"
    assert table["1"]["top_set_probability"] == pytest.approx(1.0)


def test_the_top_parent_is_the_observed_parent_with_the_largest_edge_marginal(client, table):
    edges = [e for e in client.get("/api/edges").json()["rows"] if e["child"] == "3" and not e["unobserved"]]
    best = max(edges, key=lambda e: e["pooled"])
    assert (table["3"]["top_parent"], table["3"]["top_parent_marginal"]) == ("1", pytest.approx(best["pooled"]))
    assert table["3"]["top_parent_marginal"] == pytest.approx((0.89138 + 0.77570) / 2, abs=1e-4)
    # Infection 1 has no observed parent.
    assert (table["1"]["top_parent"], table["1"]["top_parent_marginal"]) == (None, None)


def test_expected_out_degree_sums_the_infections_child_edge_marginals(client, table):
    children = client.get("/api/infection?id=1").json()["children"]
    assert table["1"]["out_degree"] == pytest.approx(sum(c["pooled"] for c in children))
    assert table["3"]["out_degree"] == pytest.approx(0.0)


def test_median_coi_and_multi_source_probabilities(table):
    assert [table[i]["coi_median"] for i in ("1", "2", "3")] == [1, 1, 1]
    # Infection 1's only set is the unobserved source alone; Infection 3 has sets with two members.
    assert table["1"]["multi_source"] == pytest.approx(0.0)
    assert 0 < table["3"]["multi_source"] < 0.3
    assert 0 <= table["3"]["multi_source_observed"] <= table["3"]["multi_source"]


def test_the_worst_r_hat_is_the_largest_over_the_infections_own_quantities_and_parent_edges(client, table):
    convergence = [r for r in client.get("/api/convergence").json()["rows"] if r["chain"] is None]
    for infection in ("1", "2", "3"):
        own = [r for r in convergence if r["rhat"] is not None and (
            (r["group"] in ("false positive rate", "false negative rate", "infection duration",
                            "unobserved-source inclusion") and r["name"] == infection)
            or (r["group"] == "edge" and r["name"].startswith(f"{infection}<-")))]
        worst = max(own, key=lambda r: number(r["rhat"]))
        row = table[infection]
        assert number(row["worst_rhat"]) == pytest.approx(number(worst["rhat"]))
        assert (row["worst_rhat_group"], row["worst_rhat_name"]) == (worst["group"], worst["name"])


def test_epi_columns_are_carried_as_text(tmp_path):
    (tmp_path / "e").mkdir()
    (tmp_path / "e" / "epi.csv").write_text("infection,district,age\n1,North,4\n3,South,12\n")
    payload = client_of(open_golden(tmp_path / "e", epi="epi.csv", label="golden")).get("/api/infections").json()
    assert payload["epi_columns"] == [{"name": "district", "type": "categorical"}, {"name": "age", "type": "numeric"}]
    table = rows(payload)
    assert (table["1"]["epi:district"], table["2"]["epi:district"], table["3"]["epi:age"]) == ("North", None, "12")
