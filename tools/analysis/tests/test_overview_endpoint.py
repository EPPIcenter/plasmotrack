"""The Overview workspace's data and the status bar's: `/api/status` and `/api/overview`, on the two golden stores as
Chains chain1 and chain2 of one Analysis."""

import pytest
from conftest import open_golden
from starlette.testclient import TestClient

from pt_analysis.server import create_app


def client_of(analysis):
    return TestClient(create_app(analysis, static_dir=analysis.path.parent / "no-dist"))


@pytest.fixture
def client(golden_analysis):
    return client_of(golden_analysis)


def pooled(client):
    return [r for r in client.get("/api/convergence").json()["rows"] if r["chain"] is None]


def number(x):
    return float(x) if isinstance(x, str) else x


# --- /api/status ------------------------------------------------------------------------------------------------------

def test_status_says_what_the_analysis_has(client, tmp_path):
    assert client.get("/api/status").json()["capabilities"] == {
        "epi": False, "truth": False, "cohort": False, "pairwise_scores": False, "analyses": 1}
    (tmp_path / "e").mkdir()
    (tmp_path / "e" / "epi.csv").write_text("infection,district\n1,North\n2,South\n3,North\n")
    with_epi = client_of(open_golden(tmp_path / "e", epi="epi.csv", label="golden"))
    assert with_epi.get("/api/status").json()["capabilities"]["epi"] is True


def test_status_gives_each_chains_closed_windows_as_the_sampler_panel_does(client):
    status = client.get("/api/status").json()
    for chain in ("chain1", "chain2"):
        windows = client.get(f"/api/sampler?chain={chain}").json()["windows"]
        assert status["derive"]["chains"][chain] == {
            "burn_in": [w["window"] for w in windows if w["phase"] == "burn_in"],
            "sampling": [w["window"] for w in windows if w["phase"] == "sampling"]}
    assert status["derive"]["current"] is True


def test_status_counts_quantities_by_the_designs_mixing_thresholds(client):
    rows = pooled(client)
    rhat = [number(r["rhat"]) for r in rows if r["rhat"] is not None]
    ess = [number(r["ess"]) for r in rows if r["ess"] is not None]
    assert client.get("/api/status").json()["mixing"] == {
        "quantities": len(rows),
        "rhat_above_1_1": sum(r > 1.1 for r in rhat), "rhat_above_1_2": sum(r > 1.2 for r in rhat),
        "ess_below_100": sum(e < 100 for e in ess), "ess_below_20": sum(e < 20 for e in ess)}


def test_status_lists_run_events_of_note_and_badges_workspaces_that_need_attention(client):
    status = client.get("/api/status").json()
    # The golden Runs end normally: burn-in ending and the Run ending are routine, not of note.
    assert status["events"] == []
    assert status["attention"]["convergence"] == status["mixing"]["rhat_above_1_1"]
    out_of_band = client.get("/api/overview").json()["samplers"]
    assert status["attention"]["samplers"] == len({s["type"] for s in out_of_band})


def test_status_badges_errors_with_the_error_rows_at_the_default_threshold_when_there_is_truth(client, tmp_path):
    assert "errors" not in client.get("/api/status").json()["attention"]
    (tmp_path / "t").mkdir()
    # 1 -> 2 is true but scores near 0 (missed); 2 -> 3 is not true.
    (tmp_path / "t" / "truth.json").write_text(
        '{"network": [{"src": "1", "dest": "3"}, {"src": "1", "dest": "2"}], "nodes": [{"id": "1"}, {"id": "2"}, {"id": "3"}]}')
    with_truth = client_of(open_golden(tmp_path / "t", truth="truth.json", label="golden"))
    errors = with_truth.get("/api/errors").json()["rows"]
    assert errors and with_truth.get("/api/status").json()["attention"]["errors"] == len(errors)


# --- /api/overview ----------------------------------------------------------------------------------------------------

def test_mixing_by_group_holds_every_pooled_value_in_convergence_group_order(client):
    rows = pooled(client)
    groups = client.get("/api/overview").json()["groups"]
    assert [g["group"] for g in groups] == list(dict.fromkeys(r["group"] for r in rows))
    for g in groups:
        mine = [r for r in rows if r["group"] == g["group"]]
        assert g["n"] == len(mine)
        assert sorted(g["ess"]) == pytest.approx(sorted(number(r["ess"]) for r in mine if r["ess"] is not None))
        assert len(g["rhat"]) == sum(r["rhat"] is not None for r in mine)


def test_the_worst_mixing_quantities_lead_with_the_largest_r_hat(client):
    rows = pooled(client)
    worst = client.get("/api/overview").json()["worst"]
    assert 0 < len(worst) <= 10
    largest = max(rows, key=lambda r: number(r["rhat"]) if r["rhat"] is not None else -1)
    assert (worst[0]["group"], worst[0]["name"]) == (largest["group"], largest["name"])
    assert [w["rhat"] for w in worst] == sorted((w["rhat"] for w in worst), reverse=True)


def test_edges_with_the_largest_between_chain_spread_come_from_the_edge_table(client):
    edges = client.get("/api/edges").json()["rows"]
    spread = client.get("/api/overview").json()["edge_spread"]
    expected = sorted((e for e in edges if not e["unobserved"]), key=lambda e: -e["spread"])[:len(spread)]
    assert [(s["parent"], s["child"]) for s in spread] == [(e["parent"], e["child"]) for e in expected]
    assert spread[0]["chains"] == {c: e["mean"] for c, e in expected[0]["chains"].items()}


def test_infections_whose_parents_disagree_most_between_chains_come_first(client):
    disagreement = client.get("/api/overview").json()["disagreement"]
    assert {d["infection"] for d in disagreement} == {"1", "2", "3"}
    assert [d["spread"] for d in disagreement] == sorted((d["spread"] for d in disagreement), reverse=True)
    # Infection 3's parent 1 has weight ~0.83 pooled, but chain2 starts without it (test_server's EDGE_3_FROM_1).
    assert disagreement[0]["infection"] == "3" and disagreement[0]["top"] == "1"


def test_sampler_types_outside_the_acceptance_band_are_listed_per_chain(client):
    assert client.get("/api/overview?band_lo=0&band_hi=1").json()["samplers"] == []
    out = client.get("/api/overview?band_lo=0.999&band_hi=1").json()["samplers"]
    types = set(client.get("/api/sampler?chain=chain1").json()["acceptance"]["type"])
    assert {s["type"] for s in out if s["chain"] == "chain1"} <= types
    assert all(s["rate"] < 0.999 for s in out)
    assert client.get("/api/overview?band_lo=x").status_code == 400


def test_the_headline_matches_the_functionals_panel_and_covers_every_infection(client):
    overview = client.get("/api/overview").json()
    headline = overview["headline"]
    everyone = client.get("/api/functionals").json()["strata"][0]
    assert headline["psi_edge"] == pytest.approx(everyone["psi_edge"])
    assert headline["mean_q"] == pytest.approx(everyone["mean_q"])
    assert headline["classes"].keys() == everyone["classes"].keys()
    for label, i in everyone["classes"].items():
        assert headline["classes"][label] == pytest.approx(i)
    assert sorted(headline["q"]) == pytest.approx(sorted(n["q"] for n in client.get("/api/network").json()["nodes"]))
    assert len(headline["effective_size"]) == 3
    top = headline["top_edges"]
    assert top[0] == {"parent": "1", "child": "3", "marginal": pytest.approx(0.83354, abs=1e-4)}


def test_the_network_thumbnail_has_every_infection_and_the_edges_at_or_above_one_half(client):
    thumbnail = client.get("/api/overview").json()["thumbnail"]
    assert sorted(n["id"] for n in thumbnail["nodes"]) == ["1", "2", "3"]
    assert [(e["parent"], e["child"]) for e in thumbnail["edges"]] == [("1", "3")]
