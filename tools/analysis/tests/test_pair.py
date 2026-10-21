"""The pair inspector's data (U3): for Infections A and B, their order per Chain, both Edge marginals, their
co-inclusion in third Infections' Parent sets, their genotypes per Locus and their infection times. On the golden
stores as Chains chain1 and chain2 of one Analysis, and on constructed data."""

import gzip
import json

import duckdb
import pytest
from conftest import open_golden
from starlette.testclient import TestClient
from test_derive import quantile_type_7

from pt_analysis.pair import MAX_GROUP, PairError, group_genotypes, inspect_pair, order_summary
from pt_analysis.server import create_app
from pt_analysis.summary import ordered_before, pair_order

DRAWS = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70]


def swap_1_and_2(analysis, at=(35, 40)):
    """Swap Infections 1 and 2 in chain1's Ordering at the Draws `at`."""
    window = analysis.runs["chain1"].path / "ordering" / "window-000001.parquet"
    swapped = window.with_name("swapped.parquet")
    duckdb.sql(f"""COPY (SELECT * REPLACE (CASE WHEN iteration IN {tuple(at)} AND infection IN ('1', '2')
                                               THEN 1 - rank ELSE rank END AS rank)
                        FROM '{window}') TO '{swapped}' (FORMAT parquet)""")
    swapped.replace(window)


def test_a_flipping_ordering_gives_the_fraction_of_draws_with_a_first_and_its_switches():
    # A before B at Draws 1, 2, 5 and 6; after it at 3 and 4: the order changes twice.
    iterations = [10, 20, 30, 40, 50, 60]
    before = [True, True, False, False, True, True]
    assert order_summary(iterations, before) == {
        "iteration": iterations, "before": [1, 1, 0, 0, 1, 1], "draws": 6, "before_count": 4,
        "probability": pytest.approx(4 / 6), "switches": 2}


def test_an_order_held_at_every_draw_has_probability_one_and_no_switches():
    assert order_summary([1, 2, 3], [False] * 3)["probability"] == 0
    got = order_summary([1, 2, 3], [True] * 3)
    assert (got["probability"], got["switches"]) == (1, 0)


def test_no_draws_give_no_probability():
    assert order_summary([], []) == {"iteration": [], "before": [], "draws": 0, "before_count": 0,
                                     "probability": None, "switches": 0}


def test_infection_1_precedes_3_at_every_draw_of_both_chains(golden_analysis):
    order = inspect_pair(golden_analysis, "1", "3")["order"]
    assert [c["chain"] for c in order["chains"]] == ["chain1", "chain2"]
    for chain in order["chains"]:
        assert chain["iteration"] == DRAWS
        assert chain["before"] == [1] * 10
        assert (chain["probability"], chain["switches"]) == (1.0, 0)
    assert order["probability"] == 1.0
    assert [c["probability"] for c in inspect_pair(golden_analysis, "3", "1")["order"]["chains"]] == [0.0, 0.0]


def test_the_order_and_switches_match_the_chain_summary_on_a_flipped_ordering(golden_analysis):
    swap_1_and_2(golden_analysis)
    order = inspect_pair(golden_analysis, "1", "2")["order"]
    chain1 = order["chains"][0]
    assert chain1["before"] == [1, 1, 0, 0, 1, 1, 1, 1, 1, 1]
    assert (chain1["before_count"], chain1["draws"], chain1["switches"]) == pair_order(
        ordered_before(golden_analysis.runs["chain1"], "1", "2")) == (8, 10, 2)
    assert order["probability"] == pytest.approx(18 / 20)


def test_an_unknown_infection_or_a_pair_of_one_infection_is_refused(golden_analysis):
    with pytest.raises(PairError, match="99"):
        inspect_pair(golden_analysis, "1", "99")
    with pytest.raises(PairError):
        inspect_pair(golden_analysis, "1", "1")


# Edge 1 -> 3 per Chain: the mean over 10 Draws of the per-Draw marginals in the exploration notes (R arrow).
EDGE_1_TO_3 = {"chain1": 0.89138, "chain2": 0.77570}


def test_both_edge_marginals_per_chain_with_an_absent_edge_at_zero(golden_analysis):
    edges = inspect_pair(golden_analysis, "1", "3")["edges"]
    assert [c["chain"] for c in edges["chains"]] == ["chain1", "chain2"]
    for chain in edges["chains"]:
        assert chain["a_to_b"] == pytest.approx(EDGE_1_TO_3[chain["chain"]], abs=1e-4)
        assert chain["b_to_a"] == 0
    assert edges["a_to_b"] == pytest.approx((0.89138 + 0.77570) / 2, abs=1e-4)
    assert edges["b_to_a"] == 0
    reversed_pair = inspect_pair(golden_analysis, "3", "1")["edges"]["chains"]
    assert [c["b_to_a"] for c in reversed_pair] == pytest.approx([0.89138, 0.77570], abs=1e-4)


def test_co_inclusion_sums_the_weight_of_a_third_infections_sets_holding_both(golden_analysis):
    """Infection 3's sets holding both 1 and 2 are {1, 2, 3_copy} and {1, 2}, the sets with two observed parents, so
    the co-inclusion is Infection 3's pooled multi-source (observed only) probability in the exploration notes."""
    (row,) = inspect_pair(golden_analysis, "1", "2")["co_inclusion"]
    assert row["infection"] == "3"
    assert [(c["chain"], c["weight"]) for c in row["chains"]] == [
        ("chain1", pytest.approx(0.00275, rel=5e-3)), ("chain2", pytest.approx(0.00273, rel=5e-3))]
    assert row["weight"] == pytest.approx((row["chains"][0]["weight"] + row["chains"][1]["weight"]) / 2)
    assert inspect_pair(golden_analysis, "1", "3")["co_inclusion"] == []


# examples/example_input.json: Infection 1 is observed 00001000 at AS1 and 000000010000 at AS11, Infection 2
# 00010000 and 000000010000, Infection 3 00001000 at AS1 only. The Latent genotypes are read from the stores (see
# test_derive): equal to the Observed ones in every Draw, except Infection 3 at AS11, which is 000000010000 in 18 of
# the 20 Draws and 000100000000 in 2.
def test_genotypes_per_locus_put_observed_and_latent_modes_side_by_side(golden_analysis):
    genotypes = inspect_pair(golden_analysis, "1", "2")["genotypes"]
    assert genotypes["truth"] is False
    assert genotypes["rows"] == [
        {"locus": "AS1", "observed_a": "00001000", "latent_a": "00001000", "agreement_a": 1.0,
         "latent_b": "00010000", "agreement_b": 1.0, "observed_b": "00010000", "differs": True},
        {"locus": "AS11", "observed_a": "000000010000", "latent_a": "000000010000", "agreement_a": 1.0,
         "latent_b": "000000010000", "agreement_b": 1.0, "observed_b": "000000010000", "differs": False}]
    assert genotypes["differing"] == 1


def test_the_latent_mode_pools_chains_and_an_unobserved_locus_has_no_observed_genotype(golden_analysis):
    rows = inspect_pair(golden_analysis, "1", "3")["genotypes"]["rows"]
    as11 = rows[1]
    assert (as11["locus"], as11["observed_b"], as11["latent_b"]) == ("AS11", None, "000000010000")
    assert as11["agreement_b"] == pytest.approx(18 / 20)
    assert [r["differs"] for r in rows] == [False, False]


def test_with_truth_the_true_latent_genotypes_are_shown_too(tmp_path):
    truth = tmp_path / "sim_0_s0.json.gz"
    nodes = [{"id": id_, "flat_latent_genotype": [{"locus": "AS1", "genotype": as1},
                                                  {"locus": "AS11", "genotype": "000000110000"}]}
             for id_, as1 in (("1", "00001100"), ("2", "00010000"), ("3", "00001000"))]
    with gzip.open(truth, "wt") as f:
        json.dump({"network": [{"src": "1", "dest": "3"}], "nodes": nodes, "full_nodes": nodes, "loci": []}, f)
    analysis = open_golden(tmp_path, truth=truth.name)
    genotypes = inspect_pair(analysis, "1", "3")["genotypes"]
    assert genotypes["truth"] is True
    assert [(r["true_a"], r["true_b"]) for r in genotypes["rows"]] == [
        ("00001100", "00001000"), ("000000110000", "000000110000")]


def test_a_truth_without_latent_genotypes_for_the_pair_shows_no_truth_columns(tmp_path):
    """A truth made for real data (edges and infection times only) has empty `flat_latent_genotype` lists: there is
    no true genotype to show, which is not the same as an unobserved one."""
    truth = tmp_path / "truth.json"
    nodes = [{"id": id_, "infection_time": 0.0, "flat_latent_genotype": []} for id_ in ("1", "2", "3")]
    truth.write_text(json.dumps({"network": [{"src": "1", "dest": "3"}], "nodes": nodes}))
    genotypes = inspect_pair(open_golden(tmp_path, truth=truth.name), "1", "3")["genotypes"]
    assert genotypes["truth"] is False
    assert all("true_a" not in r and "true_b" not in r for r in genotypes["rows"])


def test_infection_times_are_observation_time_minus_duration_quantiles_and_the_lag_is_observation_b_minus_a(
        golden_analysis):
    """examples/example_input.json observes Infection 1 at 234.0737 and 3 at 350.5. Durations are pooled over both
    Chains' Draws; the interval is the observation time minus the 97.5% and 2.5% duration quantiles (R type 7)."""
    times = inspect_pair(golden_analysis, "1", "3")["times"]
    assert times["lag"] == pytest.approx(350.5 - 234.0737)
    for infection, key, observed in (("1", "a", 234.0737), ("3", "b", 350.5)):
        durations = [r["infection_duration"] for store in golden_analysis.runs.values()
                     for r in store.rows("infections") if r["infection"] == infection]
        assert len(durations) == 20
        got = times[key]
        assert got["observation_time"] == observed
        assert got["median"] == pytest.approx(observed - quantile_type_7(durations, 0.5))
        assert (got["lower"], got["upper"]) == pytest.approx(
            (observed - quantile_type_7(durations, 0.975), observed - quantile_type_7(durations, 0.025)))
        assert got["lower"] < got["median"] < got["upper"] < observed


@pytest.fixture
def client(golden_analysis, tmp_path):
    static = tmp_path / "dist"
    static.mkdir()
    (static / "index.html").write_text("<!doctype html><title>pt-analysis</title>")
    return TestClient(create_app(golden_analysis, static_dir=static))


def test_the_endpoint_serves_every_panel_of_a_pair(client):
    response = client.get("/api/pair", params={"a": "1", "b": "3"})
    assert response.status_code == 200, response.text
    got = response.json()
    assert (got["a"], got["b"]) == ("1", "3")
    assert set(got) == {"a", "b", "order", "edges", "co_inclusion", "genotypes", "times"}
    assert [c["probability"] for c in got["order"]["chains"]] == [1.0, 1.0]
    assert [c["switches"] for c in got["order"]["chains"]] == [0, 0]
    assert [c["a_to_b"] for c in got["edges"]["chains"]] == pytest.approx([0.89138, 0.77570], abs=1e-4)
    assert got["genotypes"]["differing"] == 0
    assert got["times"]["lag"] == pytest.approx(116.4263)


def test_the_endpoint_refuses_a_missing_id_an_unknown_infection_and_a_pair_of_one(client):
    assert client.get("/api/pair", params={"a": "1"}).status_code == 400
    response = client.get("/api/pair", params={"a": "1", "b": "99"})
    assert response.status_code == 404
    assert "99" in response.json()["error"]
    assert client.get("/api/pair", params={"a": "2", "b": "2"}).status_code == 400


def test_a_group_lines_each_infection_up_at_every_locus(golden_analysis):
    got = group_genotypes(golden_analysis, ["3", "1", "3"])
    assert got["infections"] == ["3", "1"]
    assert got["differing"] == sum(locus["differs"] for locus in got["loci"])
    first = got["loci"][0]
    assert [member["id"] for member in first["members"]] == ["3", "1"]
    assert first["alleles"] >= len(first["members"][0]["mode"] or "")
    assert group_genotypes(golden_analysis, []) == {"infections": [], "loci": [], "differing": 0}
    with pytest.raises(PairError, match=str(MAX_GROUP)):
        group_genotypes(golden_analysis, [str(i) for i in range(MAX_GROUP + 1)])


def test_the_genotypes_endpoint_serves_the_group_and_refuses_a_bad_request(client, golden_analysis):
    response = client.get("/api/genotypes", params={"ids": "1,3"})
    assert response.status_code == 200, response.text
    assert response.json() == group_genotypes(golden_analysis, ["1", "3"])
    assert client.get("/api/genotypes").status_code == 400
    missing = client.get("/api/genotypes", params={"ids": "99"})
    assert missing.status_code == 404
    assert "99" in missing.json()["error"]
