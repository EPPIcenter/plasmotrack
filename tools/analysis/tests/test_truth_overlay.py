"""The truth overlay and error table (U6, part 1) on the golden Analysis with a constructed truth.

Pooled Edge marginals of the golden Analysis (see test_metrics): 1->3 0.83354, 2->3 0.00296, 1->2 0.00291, and 0 for
2->1, 3->1 and 3->2. The example input observes Infections 1, 2 and 3 at 234.0737, 306.2297 and 350.5.
"""

import pytest
from starlette.testclient import TestClient
from conftest import write_manifest
from test_metrics import open_analysis, stores, write_truth  # noqa: F401 (fixture)

from pt_analysis import truth_overlay
from pt_analysis.server import create_app

P13, P23, P12 = 0.83354, 0.00296, 0.00291
LAG = {("1", "3"): 350.5 - 234.0737, ("2", "3"): 350.5 - 306.2297, ("3", "1"): 234.0737 - 350.5,
       ("1", "2"): 306.2297 - 234.0737}


def client_for(analysis, tmp_path):
    static = tmp_path / "dist"
    static.mkdir(exist_ok=True)
    (static / "index.html").write_text("<!doctype html>")
    return TestClient(create_app(analysis, static_dir=static))


@pytest.fixture
def with_baselines(tmp_path, stores):  # noqa: F811
    (tmp_path / "null").mkdir()
    write_manifest(tmp_path / "null" / "analysis.json", label="null", runs={"chain1": "../single-chain"})
    (tmp_path / "ibd.csv").write_text("a,b,score\n1,2,0.9\n3,1,0.8\n2,3,0.05\n")
    return open_analysis(tmp_path, stores, edges=[("1", "3"), ("2", "3")], null_model="null/analysis.json",
                         pairwise_scores="ibd.csv")


def key(rows):
    return [(r["kind"], r["parent"], r["child"]) for r in rows]


# --- matching truth to the Analysis's Infections ---------------------------------------------------------------------

def test_truth_ids_are_matched_to_infections_and_the_unmatched_reported(tmp_path, stores):  # noqa: F811
    write_truth(tmp_path / "sim.json", [("1", "3"), ("4", "3")], ids=("1", "3", "4", "5"))
    analysis = open_analysis(tmp_path, stores, edges=None, truth="sim.json")
    match = truth_overlay.match(analysis)
    assert match["n_true_edges"] == 1  # 4->3 has an end outside the Analysis
    assert match["unmatched_infections"] == ["2"]  # an Infection of the Analysis the truth has no node for
    assert match["unobserved_truth_infections"] == ["4", "5"]  # simulated Infections the input does not observe


# --- error table -----------------------------------------------------------------------------------------------------

def test_missed_true_edges_rank_by_how_low_their_marginal_is(tmp_path, stores):  # noqa: F811
    analysis = open_analysis(tmp_path, stores, edges=[("1", "3"), ("2", "3"), ("3", "2")])
    rows = truth_overlay.errors(analysis, "posterior", threshold=0.5)["rows"]
    # 1->3 is called; 3->2 (0) is missed more confidently than 2->3 (0.00296).
    assert key(rows) == [("missed", "3", "2"), ("missed", "2", "3")]
    assert rows[1]["marginal"] == pytest.approx(P23, abs=2e-5)
    assert rows[1]["lag"] == pytest.approx(LAG[("2", "3")])
    # The constructed truth infects every Infection at time 0.
    assert (rows[1]["parent_true_infection_time"], rows[1]["child_true_infection_time"]) == (0.0, 0.0)
    # Both directions of 2-3 are true, so each missed edge's reverse is true.
    assert [r["reverse_true"] for r in rows] == [True, True]


def test_a_false_edge_whose_reverse_is_true_is_a_direction_error(tmp_path, stores):  # noqa: F811
    analysis = open_analysis(tmp_path, stores, edges=[("3", "1"), ("2", "3")])
    got = truth_overlay.errors(analysis, "posterior", threshold=0.5)
    rows = got["rows"]
    # Confidence of being wrong: 3->1 missed at 0 (1.0), 2->3 missed at 0.00296 (0.997), 1->3 called at 0.834.
    assert key(rows) == [("missed", "3", "1"), ("missed", "2", "3"), ("direction", "1", "3")]
    direction = rows[2]
    assert direction["marginal"] == pytest.approx(P13, abs=1e-4)
    assert direction["reverse_true"] is True
    assert direction["reverse_marginal"] == 0.0
    assert direction["lag"] == pytest.approx(LAG[("1", "3")])
    assert rows[0]["reverse_marginal"] == pytest.approx(P13, abs=1e-4)
    assert (got["source"], got["directed"], got["threshold"]) == ("posterior", True, 0.5)


def test_false_edges_rank_by_marginal_descending_and_below_the_threshold_are_not_errors(tmp_path, stores):  # noqa: F811
    analysis = open_analysis(tmp_path, stores, edges=[("2", "1")])
    rows = truth_overlay.errors(analysis, "posterior", threshold=0.0029)["rows"]
    # 2->1 (0) missed; then the called edges by marginal: 1->3 0.834, 2->3 0.00296, and 1->2 0.00291, a direction
    # error as 2->1 is true. 3->1 and 3->2, at 0, are not called.
    assert key(rows) == [("missed", "2", "1"), ("false", "1", "3"), ("false", "2", "3"), ("direction", "1", "2")]
    assert [r["reverse_true"] for r in rows] == [False, False, False, True]


def test_pairwise_scores_are_undirected_so_no_direction_errors(tmp_path, with_baselines):
    rows = truth_overlay.errors(with_baselines, "ibd", threshold=0.5)["rows"]
    # True pairs {1,3} (0.8, called) and {2,3} (0.05, missed: 0.95 wrong); {1,2} 0.9 is false (0.9 wrong).
    assert key(rows) == [("missed", "2", "3"), ("false", "1", "2")]
    assert rows[0]["marginal"] == 0.05 and rows[1]["marginal"] == 0.9
    assert rows[0]["reverse_true"] is None and rows[0]["reverse_marginal"] is None


def test_the_null_model_ranks_its_own_marginals(tmp_path, with_baselines):
    got = truth_overlay.errors(with_baselines, "null", threshold=0.5)
    assert got["label"] == "null"
    assert all(r["kind"] in ("missed", "false", "direction") for r in got["rows"])


def test_a_baseline_the_manifest_does_not_name_is_absent(tmp_path, stores):  # noqa: F811
    analysis = open_analysis(tmp_path, stores, edges=[("1", "3")])
    assert [s.source for s in truth_overlay.edge_sets(analysis)] == ["posterior"]
    with pytest.raises(truth_overlay.UnknownSource):
        truth_overlay.errors(analysis, "ibd")


# --- overlay ---------------------------------------------------------------------------------------------------------

def test_overlay_has_every_true_edge_and_the_false_edges_at_or_above_the_floor(tmp_path, with_baselines):
    got = truth_overlay.overlay(with_baselines, floor=0.0029)
    assert [(s["source"], s["label"], s["directed"]) for s in got["sources"]] == [
        ("posterior", "golden", True), ("null", "null", True), ("ibd", "IBD", False)]
    posterior = {(e["parent"], e["child"]): e for e in got["sources"][0]["edges"]}
    assert set(posterior) == {("1", "3"), ("2", "3"), ("1", "2")}
    assert posterior[("1", "3")]["true"] is True and posterior[("2", "3")]["true"] is True
    assert posterior[("1", "2")]["true"] is False
    assert posterior[("1", "3")]["lag"] == pytest.approx(LAG[("1", "3")])
    # At a higher floor the false edge goes but a true edge stays, whatever its marginal.
    high = {(e["parent"], e["child"]) for e in truth_overlay.overlay(with_baselines, floor=0.5)["sources"][0]["edges"]}
    assert high == {("1", "3"), ("2", "3")}
    ibd = {(e["parent"], e["child"]): (e["marginal"], e["true"]) for e in got["sources"][2]["edges"]}
    # A true pair takes the truth's direction; a false pair its ids in order.
    assert ibd == {("1", "3"): (0.8, True), ("2", "3"): (0.05, True), ("1", "2"): (0.9, False)}
    assert got["truth"]["n_true_edges"] == 2


def test_overlay_without_truth_has_the_edge_sets_unmarked(tmp_path, stores):  # noqa: F811
    got = truth_overlay.overlay(open_analysis(tmp_path, stores, edges=None), floor=0.5)
    assert got["truth"] is None
    [posterior] = got["sources"]
    assert [(e["parent"], e["child"], e["true"]) for e in posterior["edges"]] == [("1", "3", None)]


# --- endpoints -------------------------------------------------------------------------------------------------------

def test_truth_network_and_errors_endpoints(tmp_path, with_baselines):
    client = client_for(with_baselines, tmp_path)
    got = client.get("/api/truth-network", params={"floor": "0.5"}).json()
    assert [s["source"] for s in got["sources"]] == ["posterior", "null", "ibd"]
    got = client.get("/api/errors", params={"source": "ibd", "threshold": "0.5"}).json()
    assert key(got["rows"]) == [("missed", "2", "3"), ("false", "1", "2")]
    assert got["sources"] == ["posterior", "null", "ibd"]
    assert client.get("/api/errors").json()["source"] == "posterior"
    assert client.get("/api/errors", params={"threshold": "2"}).status_code == 400
    assert client.get("/api/truth-network", params={"floor": "x"}).status_code == 400


def test_errors_without_truth_or_with_an_absent_source_are_not_found(tmp_path, stores):  # noqa: F811
    client = client_for(open_analysis(tmp_path, stores, edges=None), tmp_path)
    assert "truth" in client.get("/api/errors").json()["error"]
    assert client.get("/api/errors").status_code == 404
    client = client_for(open_analysis(tmp_path, stores, name="b.json"), tmp_path)
    assert client.get("/api/errors", params={"source": "null"}).status_code == 404
    assert client.get("/api/errors", params={"source": "nonsense"}).status_code == 400
