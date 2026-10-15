"""The infection card's data (U10): for one Infection of an Analysis, its observation-process and duration posteriors,
COI per Draw (the derived `coi` table), its Latent genotype against the Observed one, its Parent-set distribution
with the credible collection, its parent and child Edge marginals, and its entropy, effective size and q. On the
golden stores as Chains chain1 and chain2 of one Analysis, and on constructed Latent genotypes."""

import shutil

import duckdb
import pytest
from conftest import GOLDEN

from pt_analysis.derive import derive
from pt_analysis.infection import UnknownInfection, card
from pt_analysis.store import Store

DRAWS = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70]


def read_coi(store):
    derive(store)
    return duckdb.sql(f"SELECT * FROM '{store.path / 'derived' / 'coi.parquet'}' ORDER BY iteration, infection"
                      ).fetchall()


def test_coi_is_one_at_every_draw_of_the_golden_stores(store):
    """Every golden Latent genotype has one Allele at each Locus (see test_derive)."""
    assert read_coi(store) == [(it, inf, 1) for it in DRAWS for inf in ("1", "2", "3")]


def rewrite_latent_genotypes(path, bitstrings):
    """Replace the Latent genotypes of Window 1 in the store at `path`: {(iteration, infection, locus): bitstring}."""
    window = path / "latent_genotypes" / "window-000001.parquet"
    rewritten = window.with_name("rewritten.parquet")
    cases = " ".join(f"WHEN iteration = {it} AND infection = '{inf}' AND locus = '{locus}' THEN '{bits}'"
                     for (it, inf, locus), bits in bitstrings.items())
    duckdb.sql(f"""COPY (SELECT * REPLACE (CASE {cases} ELSE alleles END AS alleles) FROM '{window}')
                   TO '{rewritten}' (FORMAT parquet)""")
    rewritten.replace(window)


def test_coi_is_the_largest_allele_count_over_loci_at_each_draw(tmp_path):
    shutil.copytree(GOLDEN / "single-chain", tmp_path / "single-chain")
    rewrite_latent_genotypes(tmp_path / "single-chain", {
        (25, "1", "AS1"): "01101000",       # 3 Alleles at AS1, 1 at AS11: COI 3
        (25, "2", "AS1"): "11000000",       # 2 at AS1, 4 at AS11: COI 4
        (25, "2", "AS11"): "100100010001",
        (30, "3", "AS11"): "000000000000",  # no Allele at AS11, 1 at AS1: COI 1
        (35, "3", "AS1"): "11111111",       # every Allele: COI 8
    })
    coi = {(it, inf): c for it, inf, c in read_coi(Store(tmp_path / "single-chain"))}
    assert (coi[(25, "1")], coi[(25, "2")], coi[(25, "3")]) == (3, 4, 1)
    assert (coi[(30, "3")], coi[(35, "3")], coi[(40, "1")]) == (1, 8, 1)
    assert len(coi) == 30


def by_chain(rows, key):
    return {r["chain"]: r[key] for r in rows}


# Infection 3 per Chain (exploration notes, R arrow): q and the expected observed-parent count, and the entropy, each
# the mean over its 10 Draws of `infection_draws`.
def test_q_entropy_and_effective_size_per_chain_and_pooled(golden_analysis):
    summary = card(golden_analysis, "3")["summary"]
    assert by_chain(summary["chains"], "q") == {"chain1": pytest.approx(0.1125, rel=5e-3),
                                                "chain2": pytest.approx(0.2262, rel=5e-3)}
    assert by_chain(summary["chains"], "entropy") == {"chain1": pytest.approx(0.3697, rel=5e-3),
                                                      "chain2": pytest.approx(0.1412, rel=5e-3)}
    assert by_chain(summary["chains"], "expected_observed_parents") == {
        "chain1": pytest.approx(0.8944, rel=5e-3), "chain2": pytest.approx(0.7786, rel=5e-3)}
    assert summary["q"] == pytest.approx((0.1125 + 0.2262) / 2, rel=5e-3)
    assert summary["entropy"] == pytest.approx((0.3697 + 0.1412) / 2, rel=5e-3)
    # exp is convex, so the mean effective size is at least exp of the mean entropy; it is at most the 6 sets.
    for chain in summary["chains"]:
        assert 1 <= chain["effective_size"] <= 6
    assert summary["effective_size"] >= 2.718281828 ** summary["entropy"]
    assert card(golden_analysis, "1")["summary"]["chains"] == [
        {"chain": c, "q": 1.0, "expected_observed_parents": 0.0, "entropy": 0.0, "effective_size": 1.0}
        for c in ("chain1", "chain2")]


# Infection 3's Parent sets by member, and their composition classes (see test_derive), with each Chain's ranking by
# pooled weight (exploration notes).
MEMBERS = {4: ("1",), 0: ("3_copy",), 1: ("1", "3_copy"), 5: ("1", "2"), 2: ("2", "3_copy"), 3: ("1", "2", "3_copy")}
CLASS = {4: "observed_1", 0: "unobserved_only", 1: "observed_1_unobserved", 5: "observed_2",
         2: "observed_1_unobserved", 3: "observed_2_unobserved"}
RANKING = {"chain1": [4, 0, 1, 5, 2, 3], "chain2": [4, 0, 5, 1, 2, 3]}


def test_each_chains_parent_sets_are_ranked_with_members_class_and_the_credible_collection(golden_analysis):
    sets = card(golden_analysis, "3")["parent_sets"]
    assert sets["level"] == 0.8
    assert [c["chain"] for c in sets["chains"]] == ["chain1", "chain2"]
    for chain in sets["chains"]:
        rows = chain["sets"]
        assert [r["set_id"] for r in rows] == RANKING[chain["chain"]]
        assert [r["rank"] for r in rows] == [1, 2, 3, 4, 5, 6]
        assert all((r["members"], r["class"]) == (list(MEMBERS[r["set_id"]]), CLASS[r["set_id"]]) for r in rows)
        assert [r["cumulative_mass"] for r in rows] == pytest.approx(
            [sum(r["weight"] for r in rows[:k + 1]) for k in range(6)])
    # The 80% credible collection: set 4 alone in chain1 (0.885), sets 4 and 0 in chain2 (0.771, then 0.995).
    assert [[r["set_id"] for r in c["sets"] if r["credible"]] for c in sets["chains"]] == [[4], [4, 0]]


def test_the_credible_collection_matches_the_derive_commands(golden_analysis):
    from pt_analysis import functionals

    for chain in card(golden_analysis, "3")["parent_sets"]["chains"]:
        expected = [s.set_id for s in functionals.credible_collection(golden_analysis.runs[chain["chain"]], "3")]
        assert [r["set_id"] for r in chain["sets"] if r["credible"]] == expected


def test_pooled_parent_sets_weight_each_chain_by_its_draws(golden_analysis):
    pooled = card(golden_analysis, "3")["parent_sets"]["pooled"]
    assert [r["set_id"] for r in pooled] == [4, 0, 1, 5, 2, 3]
    weights = {r["set_id"]: r["weight"] for r in pooled}
    assert weights[4] == pytest.approx((0.885 + 0.771) / 2, rel=1e-2)
    assert weights[0] == pytest.approx((0.108 + 0.224) / 2, rel=1e-2)
    assert [r["chains"] for r in pooled if r["set_id"] == 4] == [
        [{"chain": "chain1", "weight": pytest.approx(0.885, rel=1e-2)},
         {"chain": "chain2", "weight": pytest.approx(0.771, rel=1e-2)}]]
    assert [r["credible"] for r in pooled] == [True] + [False] * 5  # 0.828 reaches 0.8 alone
    assert all((r["members"], r["class"]) == (list(MEMBERS[r["set_id"]]), CLASS[r["set_id"]]) for r in pooled)


# Edge 1 -> 3 per Chain (see test_pair) and q of Infection 3 as the unobserved source's edge.
def test_parent_marginals_hold_observed_parents_and_the_unobserved_source(golden_analysis):
    parents = {r["parent"]: r for r in card(golden_analysis, "3")["parents"]}
    assert set(parents) == {"1", "2", "3_copy"}
    assert by_chain(parents["1"]["chains"], "mean") == {"chain1": pytest.approx(0.89138, abs=1e-4),
                                                        "chain2": pytest.approx(0.77570, abs=1e-4)}
    assert parents["1"]["pooled"] == pytest.approx((0.89138 + 0.77570) / 2, abs=1e-4)
    assert parents["2"]["pooled"] == pytest.approx(0.003, abs=1e-3)
    assert (parents["1"]["group"], parents["1"]["name"], parents["1"]["unobserved"]) == ("edge", "3<-1", False)
    q = parents["3_copy"]
    assert (q["group"], q["name"], q["unobserved"]) == ("unobserved-source inclusion", "3", True)
    assert q["pooled"] == pytest.approx((0.1125 + 0.2262) / 2, rel=5e-3)
    assert [r["parent"] for r in card(golden_analysis, "3")["parents"]] == ["1", "3_copy", "2"]  # largest first


def test_child_marginals_hold_the_infections_it_may_have_infected(golden_analysis):
    children = card(golden_analysis, "1")["children"]
    assert [r["child"] for r in children] == ["3", "2"]
    assert children[0]["pooled"] == pytest.approx((0.89138 + 0.77570) / 2, abs=1e-4)
    assert children[1]["pooled"] == pytest.approx(0.0029, abs=5e-4)
    assert (children[0]["group"], children[0]["name"]) == ("edge", "3<-1")
    assert card(golden_analysis, "3")["children"] == []


QUANTITIES = {"false_positives": "expected_false_positives", "false_negatives": "expected_false_negatives",
              "duration": "infection_duration"}


def test_observation_process_and_duration_per_chain_and_pooled_quantiles(golden_analysis):
    from test_derive import quantile_type_7

    quantities = card(golden_analysis, "3")["quantities"]
    rows = {label: [r for r in store.rows("infections") if r["infection"] == "3"]
            for label, store in golden_analysis.runs.items()}
    for key, column in QUANTITIES.items():
        q = quantities[key]
        assert [c["chain"] for c in q["chains"]] == ["chain1", "chain2"]
        for chain in q["chains"]:
            assert chain["iteration"] == DRAWS
            assert chain["value"] == [r[column] for r in rows[chain["chain"]]]
        pooled = [r[column] for label in rows for r in rows[label]]
        assert (q["lower"], q["median"], q["upper"]) == pytest.approx(
            tuple(quantile_type_7(pooled, p) for p in (0.025, 0.5, 0.975)))
    assert (quantities["false_positives"]["group"], quantities["false_positives"]["name"]) == ("false positive rate", "3")
    assert quantities["infection_time"]["group"] is None
    # Infection 3 is observed at 350.5; infection time is that minus the duration, so its quantiles swap ends.
    duration, time = quantities["duration"], quantities["infection_time"]
    assert time["chains"][0]["value"] == pytest.approx([350.5 - d for d in duration["chains"][0]["value"]])
    assert (time["lower"], time["median"], time["upper"]) == pytest.approx(
        (350.5 - duration["upper"], 350.5 - duration["median"], 350.5 - duration["lower"]))
    assert card(golden_analysis, "3")["observation_time"] == 350.5


def test_coi_per_draw_per_chain_and_its_pooled_distribution(golden_analysis):
    coi = card(golden_analysis, "3")["coi"]
    assert [(c["chain"], c["iteration"], c["coi"]) for c in coi["chains"]] == [
        ("chain1", DRAWS, [1] * 10), ("chain2", DRAWS, [1] * 10)]
    assert coi["distribution"] == [{"coi": 1, "draws": 20, "fraction": 1.0}]


def test_the_latent_genotype_has_posterior_presence_per_allele_and_the_observed_bits(golden_analysis):
    """Infection 3 is observed 00001000 at AS1 and not at AS11; its Latent genotype is 00001000 at AS1 in all 20 Draws
    and at AS11 000000010000 in 18, 000100000000 in 2 (see test_derive)."""
    loci = card(golden_analysis, "3")["genotype"]["loci"]
    assert [(locus["locus"], locus["observed"]) for locus in loci] == [("AS1", "00001000"), ("AS11", None)]
    assert loci[0]["presence"] == [0, 0, 0, 0, 1, 0, 0, 0]
    assert loci[1]["presence"] == pytest.approx([0, 0, 0, 0.1, 0, 0, 0, 0.9, 0, 0, 0, 0])


def test_an_unknown_infection_is_refused(golden_analysis):
    with pytest.raises(UnknownInfection, match="99"):
        card(golden_analysis, "99")


@pytest.fixture
def client(golden_analysis, tmp_path):
    from starlette.testclient import TestClient

    from pt_analysis.server import create_app

    static = tmp_path / "dist"
    static.mkdir()
    (static / "index.html").write_text("<!doctype html><title>pt-analysis</title>")
    return TestClient(create_app(golden_analysis, static_dir=static))


def test_the_endpoint_serves_every_panel_of_an_infection(client):
    response = client.get("/api/infection", params={"id": "3"})
    assert response.status_code == 200, response.text
    got = response.json()
    assert set(got) == {"infection", "observation_time", "quantities", "coi", "genotype", "parent_sets", "parents",
                        "children", "summary"}
    assert got["infection"] == "3"
    assert [len([r for r in c["sets"] if r["credible"]]) for c in got["parent_sets"]["chains"]] == [1, 2]
    assert got["summary"]["q"] == pytest.approx(0.16935, rel=5e-3)


def test_the_endpoint_refuses_a_missing_or_unknown_id(client):
    assert client.get("/api/infection").status_code == 400
    response = client.get("/api/infection", params={"id": "99"})
    assert response.status_code == 404
    assert "99" in response.json()["error"]
