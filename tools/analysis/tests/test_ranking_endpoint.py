"""The ranking and parent-set endpoints (U6, part 2) on the golden Analysis with a constructed truth."""

import pytest
from starlette.testclient import TestClient
from conftest import write_manifest
from test_metrics import open_analysis, stores  # noqa: F401 (fixture)

from pt_analysis import metrics
from pt_analysis.server import create_app


def client_for(analysis, tmp_path):
    static = tmp_path / "dist"
    static.mkdir(exist_ok=True)
    (static / "index.html").write_text("<!doctype html>")
    return TestClient(create_app(analysis, static_dir=static))


@pytest.fixture
def with_baselines(tmp_path, stores):  # noqa: F811
    (tmp_path / "null").mkdir()
    write_manifest(tmp_path / "null" / "analysis.json", label="null", runs={"chain1": "../single-chain"})
    (tmp_path / "ibd.csv").write_text("a,b,score\n1,2,0.9\n1,3,0.8\n2,3,0.1\n")
    return open_analysis(tmp_path, stores, null_model="null/analysis.json", pairwise_scores="ibd.csv")


def test_ranking_returns_curves_and_metrics_for_the_posterior_and_each_baseline(tmp_path, with_baselines):
    got = client_for(with_baselines, tmp_path).get("/api/ranking").json()
    assert got["truth"] == {"n_true_cases": 3, "n_observed_cases": 3, "true_average_outdegree": pytest.approx(2 / 3),
                            "true_unobserved_source_proportion": pytest.approx(1 / 3)}
    assert [(f["fit"], f["label"]) for f in got["fits"]] == [("posterior", "golden"), ("null", "null"), ("ibd", "IBD")]
    posterior = got["fits"][0]
    assert posterior["directed"]["average_precision"] == pytest.approx(0.5 + 0.5 * 2 / 3)
    assert posterior["directed"]["n_edge_universe"] == 6
    assert [(p["recall"], p["precision"]) for p in posterior["directed"]["curve"]] == pytest.approx(
        [(0.5, 1), (0.5, 0.5), (1, 2 / 3), (1, 1 / 3)])
    assert posterior["undirected"]["n_edge_universe"] == 3
    ibd = got["fits"][2]
    assert ibd["directed"] is None
    assert ibd["undirected"]["average_precision"] == 1.0


def test_without_baselines_only_the_posterior_is_returned(tmp_path, stores):  # noqa: F811
    got = client_for(open_analysis(tmp_path, stores), tmp_path).get("/api/ranking").json()
    assert [f["fit"] for f in got["fits"]] == ["posterior"]
    got = client_for(open_analysis(tmp_path, stores), tmp_path).get("/api/parent-sets").json()
    assert [f["fit"] for f in got["fits"]] == ["posterior"]
    assert not any(key.startswith("null_") for row in got["rows"] for key in row)


def test_ranking_without_truth_is_not_found(tmp_path, stores):  # noqa: F811
    response = client_for(open_analysis(tmp_path, stores, edges=None), tmp_path).get("/api/ranking")
    assert response.status_code == 404
    assert "truth" in response.json()["error"]


def test_parent_sets_return_one_row_per_infection_beside_the_null_fit(tmp_path, with_baselines):
    got = client_for(with_baselines, tmp_path).get("/api/parent-sets", params={"level": "0.95"}).json()
    assert (got["level"], got["truth"]) == (0.95, True)
    assert [(f["fit"], f["label"]) for f in got["fits"]] == [("posterior", "golden"), ("null", "null")]
    assert got["fits"][0]["summary"]["credible_parent_set_coverage_a095"] == pytest.approx(2 / 3)
    row = got["rows"][2]
    assert row["infection"] == "3"
    assert (row["credible_size"], row["covered"]) == (2, True)
    # The null fit here is the single-chain store alone: Infection 3's {1} has 0.885, so 0.95 needs {3_copy} too.
    assert (row["null_credible_size"], row["null_covered"]) == (2, True)
    assert row["null_top_probability"] == pytest.approx(0.885, abs=1e-3)


def test_parent_sets_without_truth_leave_coverage_empty(tmp_path, stores):  # noqa: F811
    got = client_for(open_analysis(tmp_path, stores, edges=None), tmp_path).get("/api/parent-sets").json()
    assert got["truth"] is False
    assert [r["covered"] for r in got["rows"]] == [None, None, None]


def test_a_level_outside_the_unit_interval_is_refused(tmp_path, stores):  # noqa: F811
    client = client_for(open_analysis(tmp_path, stores), tmp_path)
    assert client.get("/api/parent-sets", params={"level": "1.5"}).status_code == 400
    assert client.get("/api/parent-sets", params={"level": "x"}).status_code == 400


def test_a_drawn_curve_keeps_only_the_points_where_true_positives_change():
    r = metrics.rank([0.9, 0.8, 0.7, 0.6, 0.5, 0.4], [True, False, False, False, True, False])
    # tp 1,1,1,1,2,2: the run of fp-only points at tp 1 keeps its ends.
    assert [(p.tp, p.fp) for p in metrics.drawn_curve(r)] == [(1, 0), (1, 3), (2, 3), (2, 4)]
