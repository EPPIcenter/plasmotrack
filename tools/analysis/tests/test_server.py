"""The HTTP backend of `pt-analysis serve`: the Analysis summary, the convergence table, per-Chain traces and ESS and
R-hat over a Draw range, on the two golden stores as Chains chain1 and chain2 of one Analysis."""

import pytest
from starlette.testclient import TestClient

from pt_analysis.diagnostics import multichain_ess, split_rhat
from pt_analysis.server import create_app


@pytest.fixture
def client(golden_analysis, tmp_path):
    static = tmp_path / "dist"
    static.mkdir()
    (static / "index.html").write_text("<!doctype html><title>pt-analysis</title>")
    return TestClient(create_app(golden_analysis, static_dir=static))


def test_the_summary_names_the_analysis_and_counts_chains_draws_infections_and_loci(client):
    assert client.get("/api/analysis").json() == {
        "label": "golden",
        "chains": [{"chain": "chain1", "draws": 10}, {"chain": "chain2", "draws": 10}],
        "draws": 20, "infections": 3, "loci": 2}


def test_serving_derives_the_convergence_table_and_returns_its_rows(client, golden_analysis):
    assert (golden_analysis.path.parent / "analysis.derived" / "convergence.parquet").is_file()
    rows = client.get("/api/convergence").json()["rows"]
    # 1 + 3 + 3 + 3 + 3 + 20 + 3 + 3 + 3 quantities (test_convergence), each with two Chain rows and a pooled row.
    assert len(rows) == 42 * 3
    assert rows[:3] == [
        {"group": "log posterior", "name": "log posterior", "chain": "chain1", "ess": pytest.approx(10.0),
         "rhat": None},
        {"group": "log posterior", "name": "log posterior", "chain": "chain2", "ess": pytest.approx(7.44641539328075),
         "rhat": None},
        {"group": "log posterior", "name": "log posterior", "chain": None, "ess": pytest.approx(3.58385705875339),
         "rhat": pytest.approx(1.2643097719105674)}]


def test_the_convergence_table_follows_a_chain_that_closed_another_window_while_serving(client, golden_analysis):
    # Hide chain2's last sampling Window, so the served table is derived without it, then close it again.
    marker = golden_analysis.runs["chain2"].path / "_closed" / "window-000002"
    aside = marker.with_name("aside")
    marker.rename(aside)
    before = {r["chain"]: r["ess"] for r in client.get("/api/convergence").json()["rows"][:3]}
    aside.rename(marker)
    after = {r["chain"]: r["ess"] for r in client.get("/api/convergence").json()["rows"][:3]}
    assert after["chain2"] == pytest.approx(7.44641539328075)
    assert before["chain2"] != after["chain2"]


DRAWS = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70]

# Per-Draw values computed with R arrow from the golden stores (exploration notes), chain1 then chain2.
EDGE_3_FROM_1 = ([0.9638, 0.8808, 0.9799, 0.9337, 0.8457, 0.8624, 0.9414, 0.8481, 0.8313, 0.8267],
                 [0.0031, 0.0039, 0.9757, 0.9767, 0.9533, 0.9499, 0.9657, 0.9705, 0.9718, 0.9864])
PSI_EDGE = ([0.3217, 0.2956, 0.3332, 0.3125, 0.2828, 0.2888, 0.3155, 0.2844, 0.2788, 0.2769],
            [0.00216, 0.00289, 0.3260, 0.3259, 0.3223, 0.3167, 0.3241, 0.3262, 0.3283, 0.3311])


def trace(client, group, name, **params):
    response = client.get("/api/trace", params={"group": group, "name": name, **params})
    assert response.status_code == 200, response.text
    return response.json()


def test_a_trace_has_one_series_per_chain_over_its_draws(client):
    got = trace(client, "edge", "3<-1")
    assert (got["group"], got["name"], got["burn_in"]) == ("edge", "3<-1", False)
    assert [c["chain"] for c in got["chains"]] == ["chain1", "chain2"]
    for chain, expected in zip(got["chains"], EDGE_3_FROM_1):
        assert chain["iteration"] == DRAWS
        assert chain["phase"] == ["sampling"] * 10
        assert chain["value"] == pytest.approx(expected, abs=1e-4)


def test_functional_traces_come_from_functional_draws(client):
    got = trace(client, "functionals", "psi_edge")
    for chain, expected in zip(got["chains"], PSI_EDGE):
        assert chain["value"] == pytest.approx(expected, abs=1e-4)


def test_unobserved_source_inclusion_is_q_per_draw(client):
    chain1 = trace(client, "unobserved-source inclusion", "3")["chains"][0]
    assert sum(chain1["value"]) / 10 == pytest.approx(0.1125, abs=1e-4)


def test_burn_in_rows_come_first_on_request_for_traces_of_the_store(client):
    got = trace(client, "log posterior", "log posterior", burn_in="true")
    assert got["burn_in"] is True
    chain1 = got["chains"][0]
    assert chain1["iteration"] == [0, 5, 10, 15, 20] + DRAWS
    assert chain1["phase"] == ["burn_in"] * 5 + ["sampling"] * 10
    assert chain1["value"][0] == pytest.approx(17.6335 - 10.5023, abs=1e-4)
    rates = trace(client, "false positive rate", "1", burn_in="true")["chains"][0]
    assert rates["value"][1] == pytest.approx(0.01145, abs=1e-5)
    duration = trace(client, "infection duration", "1", burn_in="true")["chains"][0]
    assert duration["value"][1] == pytest.approx(17.910, abs=1e-3)
    frequency = trace(client, "allele frequency", "AS1[3]", burn_in="true")["chains"][0]
    assert frequency["value"][0] == pytest.approx(0.7533, abs=1e-4)
    scalar = trace(client, "global scalars", "mean_coi", burn_in="true")["chains"][0]
    assert scalar["value"][0] == pytest.approx(1.01, abs=1e-6)


def test_derived_traces_have_no_burn_in_rows(client):
    got = trace(client, "edge", "3<-1", burn_in="true")
    assert got["burn_in"] is False
    assert got["chains"][0]["iteration"] == DRAWS


def test_an_unknown_group_or_quantity_is_an_error(client):
    assert client.get("/api/trace", params={"group": "nope", "name": "x"}).status_code == 400
    assert client.get("/api/trace", params={"group": "edge", "name": "1<-3"}).status_code == 404
    assert client.get("/api/trace", params={"group": "edge"}).status_code == 400


def draw_range(client, group, name, **params):
    response = client.get("/api/range", params={"group": group, "name": name, **params})
    assert response.status_code == 200, response.text
    return response.json()


def test_the_whole_range_gives_the_convergence_tables_ess_and_rhat(client):
    got = draw_range(client, "log posterior", "log posterior")
    # BASELINE in test_convergence: the benchmark scripts' numbers on these stores.
    assert got["chains"] == [{"chain": "chain1", "draws": 10, "ess": pytest.approx(10.0)},
                             {"chain": "chain2", "draws": 10, "ess": pytest.approx(7.44641539328075)}]
    assert (got["ess"], got["rhat"]) == pytest.approx((3.58385705875339, 1.2643097719105674))
    edge = draw_range(client, "edge", "3<-1")
    assert (edge["ess"], edge["rhat"]) == pytest.approx((12.484076791343613, 1.1044701085052935))


def test_a_brushed_range_recomputes_ess_and_rhat_over_its_draws_only(client):
    # Iterations 45-70 are the last six Draws: replica-exchange's two Draws near 0 are left out.
    got = draw_range(client, "edge", "3<-1", first=45, last=70)
    assert (got["first"], got["last"]) == (45, 70)
    assert [c["draws"] for c in got["chains"]] == [6, 6]
    last_six = [values[4:] for values in EDGE_3_FROM_1]
    assert got["ess"] == pytest.approx(multichain_ess(last_six), rel=1e-2)
    assert got["rhat"] == pytest.approx(split_rhat(last_six), rel=1e-2)
    assert got["chains"][0]["ess"] == pytest.approx(multichain_ess(last_six[:1]), rel=1e-2)
    # Brush ends between Draws select the Draws inside.
    assert [c["draws"] for c in draw_range(client, "edge", "3<-1", first=43, last=71)["chains"]] == [6, 6]


def test_a_range_of_fewer_than_four_draws_has_no_statistics(client):
    got = draw_range(client, "functionals", "psi_edge", first=60, last=70)
    assert [c["draws"] for c in got["chains"]] == [3, 3]
    assert (got["ess"], got["rhat"]) == (None, None)
    assert [c["ess"] for c in got["chains"]] == [None, None]


def test_burn_in_rows_never_count_in_a_range(client):
    got = draw_range(client, "log posterior", "log posterior", first=0, last=70)
    assert [c["draws"] for c in got["chains"]] == [10, 10]


def test_the_built_frontend_is_served_at_the_root(client):
    response = client.get("/")
    assert response.status_code == 200 and "<title>pt-analysis</title>" in response.text


def test_without_a_built_frontend_the_root_says_how_to_build_it(golden_analysis, tmp_path):
    response = TestClient(create_app(golden_analysis, static_dir=tmp_path / "missing")).get("/")
    assert response.status_code == 404 and "pnpm build" in response.text
