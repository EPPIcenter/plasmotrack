"""Edge statistics (U2): per (child, parent), the pooled Edge marginal, each Chain's mean, the between-Chain spread,
each Chain's switch count (adjacent Draws on either side of one half) and the multi-chain ESS."""

import duckdb
import pytest
from starlette.testclient import TestClient

from pt_analysis.derive import derived_analysis
from pt_analysis.edge_stats import ChainEdge, EdgeStats, chain_edge_stats, edge_stats, edge_table
from pt_analysis.server import create_app

DRAWS = [10, 20, 30, 40, 50, 60]


def edges(*rows):
    """An in-memory `edges` table of (iteration, child, parent, prob) rows."""
    connection = duckdb.connect()
    return connection.sql("SELECT * FROM (VALUES " + ", ".join(
        f"({i}, '{c}', '{p}', {w}::DOUBLE)" for i, c, p, w in rows) + ") t(iteration, child, parent, prob)")


def test_a_chain_gives_each_edge_its_mean_over_draws_and_its_crossings_of_one_half():
    got = chain_edge_stats(edges(
        # 9 <- 1 locked near 1: no crossing.
        *[(d, "9", "1", 0.99) for d in DRAWS],
        # 9 <- 2 above, below, above, below, above, above: four crossings.
        *[(d, "9", "2", w) for d, w in zip(DRAWS, [0.8, 0.2, 0.7, 0.4, 0.6, 0.9])],
    ), DRAWS)
    assert got == {("9", "1"): ChainEdge(mean=pytest.approx(0.99), switches=0),
                   ("9", "2"): ChainEdge(mean=pytest.approx(3.6 / 6), switches=4)}


def test_a_draw_without_a_row_is_zero_for_the_mean_and_the_crossings():
    # Rows only at 10 and 30: 0.9, (0), 0.9, (0), (0), (0) has three crossings and mean 1.8 / 6.
    got = chain_edge_stats(edges((10, "9", "1", 0.9), (30, "9", "1", 0.9)), DRAWS)
    assert got == {("9", "1"): ChainEdge(mean=pytest.approx(0.3), switches=3)}


def test_exactly_one_half_is_not_above_it():
    got = chain_edge_stats(edges(*[(d, "9", "1", w) for d, w in zip(DRAWS, [0.4, 0.5, 0.4, 0.5, 0.6, 0.5])]), DRAWS)
    assert got[("9", "1")].switches == 2


def test_chains_that_disagree_have_a_spread_and_a_chain_without_the_edge_counts_as_zero():
    # He01-style: 9 <- 1 locked at 1 in chain A and absent (0) in chain B; 9 <- 2 agreed on; 9 <- 3 only in A.
    per_chain = {"A": {("9", "1"): ChainEdge(0.98, 0), ("9", "2"): ChainEdge(0.4, 3), ("9", "3"): ChainEdge(0.1, 0),
                       ("9", "9_copy"): ChainEdge(0.6, 1)},
                 "B": {("9", "2"): ChainEdge(0.5, 2), ("9", "9_copy"): ChainEdge(1.0, 0)}}
    ess = {("edge", "9<-1"): 2.5, ("edge", "9<-2"): 11.0, ("unobserved-source inclusion", "9"): 4.0}
    got = {(e.child, e.parent): e for e in edge_stats(per_chain, {"A": 6, "B": 4}, {"9_copy"}, ess)}

    locked = got[("9", "1")]
    assert locked == EdgeStats(child="9", parent="1", unobserved=False, pooled=pytest.approx(0.98 * 6 / 10),
                               spread=pytest.approx(0.98), switches=0, ess=2.5,
                               chains={"A": ChainEdge(0.98, 0), "B": ChainEdge(0.0, 0)})
    assert (locked.group, locked.name) == ("edge", "9<-1")

    agreed = got[("9", "2")]
    assert (agreed.pooled, agreed.spread, agreed.switches) == (pytest.approx(0.44), pytest.approx(0.1), 5)
    assert got[("9", "3")].ess is None

    source = got[("9", "9_copy")]
    assert (source.unobserved, source.group, source.name, source.ess) == (True, "unobserved-source inclusion", "9", 4.0)
    assert source.pooled == pytest.approx((0.6 * 6 + 4) / 10)


def test_edges_are_in_child_then_parent_order_and_chains_without_draws_are_left_out():
    per_chain = {"A": {("9", "2"): ChainEdge(0.5, 1), ("10", "9"): ChainEdge(0.2, 0)}, "B": {}}
    got = edge_stats(per_chain, {"A": 6, "B": 0}, set(), {})
    assert [(e.child, e.parent) for e in got] == [("10", "9"), ("9", "2")]
    assert got[1].chains == {"A": ChainEdge(0.5, 1)}
    assert (got[1].pooled, got[1].spread) == (pytest.approx(0.5), 0.0)


def test_golden_edge_3_from_1_sits_high_in_chain1_and_crosses_one_half_once_in_chain2(golden_analysis):
    got = {(e.child, e.parent): e for e in edge_table(golden_analysis)}
    edge = got[("3", "1")]
    # Per-Draw values from R arrow (exploration notes): chain1 means 0.89138, chain2 0.7757; chain2 starts at 0.0031,
    # 0.0039 then stays above 0.97.
    assert edge.chains["chain1"].mean == pytest.approx(0.89138, abs=1e-4)
    assert edge.chains["chain2"].mean == pytest.approx(0.7757, abs=1e-4)
    assert (edge.chains["chain1"].switches, edge.chains["chain2"].switches, edge.switches) == (0, 1, 1)
    assert edge.pooled == pytest.approx(0.83354, abs=1e-4)
    assert edge.spread == pytest.approx(0.11568, abs=1e-4)
    assert not edge.unobserved
    pooled_ess = derived_analysis(golden_analysis, "convergence").filter(
        "\"group\" = 'edge' AND name = '3<-1' AND chain IS NULL").project("ess").fetchone()[0]
    assert edge.ess == pytest.approx(pooled_ess)


def test_golden_unobserved_source_edges_are_flagged_and_carry_q(golden_analysis):
    got = {(e.child, e.parent): e for e in edge_table(golden_analysis)}
    assert {key for key, e in got.items() if e.unobserved} == {("1", "1_copy"), ("2", "2_copy"), ("3", "3_copy")}
    # Pooled q of Infection 3: 0.1125 in chain1 and 0.2262 in chain2 (exploration notes).
    assert got[("3", "3_copy")].pooled == pytest.approx((0.1125 + 0.2262) / 2, abs=1e-3)
    assert got[("1", "1_copy")].pooled == pytest.approx(1.0)
    assert got[("1", "1_copy")].ess is None  # constant: no ESS
    assert (got[("3", "3_copy")].group, got[("3", "3_copy")].name) == ("unobserved-source inclusion", "3")


def test_the_edge_stats_table_is_derived_beside_the_convergence_table(golden_analysis):
    edge_table(golden_analysis)
    assert (golden_analysis.path.parent / "analysis.derived" / "edge_stats.parquet").is_file()


def test_the_edges_endpoint_gives_a_row_per_edge_with_its_trace_quantity(golden_analysis, tmp_path):
    client = TestClient(create_app(golden_analysis, static_dir=tmp_path / "no-dist"))
    body = client.get("/api/edges").json()
    assert body["chains"] == ["chain1", "chain2"]
    rows = {(r["child"], r["parent"]): r for r in body["rows"]}
    edge = rows[("3", "1")]
    assert (edge["group"], edge["name"], edge["unobserved"], edge["switches"]) == ("edge", "3<-1", False, 1)
    assert edge["pooled"] == pytest.approx(0.83354, abs=1e-4)
    assert edge["spread"] == pytest.approx(0.11568, abs=1e-4)
    assert edge["chains"]["chain1"] == {"mean": pytest.approx(0.89138, abs=1e-4), "switches": 0}
    assert edge["chains"]["chain2"] == {"mean": pytest.approx(0.7757, abs=1e-4), "switches": 1}
    assert isinstance(edge["ess"], float)
    source = rows[("3", "3_copy")]
    assert (source["group"], source["name"], source["unobserved"]) == ("unobserved-source inclusion", "3", True)
    assert rows[("1", "1_copy")]["ess"] is None
