"""Transmission-network summaries of one Run, from its derived tables, against values computed with R arrow."""

import pytest

from pt_analysis.functionals import credible_collection, map_graph, psi_map, size_range, top_edges

CREDIBLE_SIZES_80 = {"single-chain": {"1": 1, "2": 1, "3": 1}, "replica-exchange": {"1": 1, "2": 1, "3": 2}}
SIZE_RANGE_80 = {"single-chain": (1 / 3, 1 / 3), "replica-exchange": (0, 1 / 3)}
# Pooled edge 1->3: the mean of the per-Draw values in the exploration notes.
EDGE_1_TO_3 = {"single-chain": 0.89138, "replica-exchange": 0.77570}
Q_3 = {"single-chain": 0.1125, "replica-exchange": 0.2262}  # unobserved-source inclusion of inf 3
EDGE_1_TO_2 = {"single-chain": 0.00267, "replica-exchange": 0.00315}  # inf 2's set {1, 2_copy}


def test_the_map_graph_gives_every_infection_its_modal_set(store):
    graph = map_graph(store)
    assert {inf: s.set_id for inf, s in graph.items()} == {"1": 0, "2": 0, "3": 4}
    assert graph["3"].members == ("1",)
    assert graph["1"].members == ("1_copy",)
    assert psi_map(store) == pytest.approx(1 / 3)


def test_the_parent_set_size_range_at_0_8(store):
    assert size_range(store) == pytest.approx(SIZE_RANGE_80[store.path.name])


def test_credible_collection_sizes_at_0_8(store):
    sizes = {inf: len(credible_collection(store, inf)) for inf in ("1", "2", "3")}
    assert sizes == CREDIBLE_SIZES_80[store.path.name]


def test_a_credible_collection_lists_its_sets_with_members_in_rank_order(store):
    collection = credible_collection(store, "3", level=0.999)
    assert [s.set_id for s in collection] == {"single-chain": [4, 0, 1, 5], "replica-exchange": [4, 0, 5, 1]}[
        store.path.name]
    assert [s.members for s in collection][:3] == [("1",), ("3_copy",)] + {
        "single-chain": [("1", "3_copy")], "replica-exchange": [("1", "2")]}[store.path.name]
    assert [s.composition for s in collection[:2]] == ["observed_1", "unobserved_only"]
    assert [s.rank for s in collection] == [1, 2, 3, 4]
    assert collection[-1].cumulative_mass >= 0.999


def test_top_edges_of_the_posterior_transmission_matrix(store):
    (top,) = top_edges(store, n=1)
    assert (top.parent, top.child) == ("1", "3")
    assert top.weight == pytest.approx(EDGE_1_TO_3[store.path.name], abs=1e-4)
    edges = {(e.parent, e.child): e.weight for e in top_edges(store)}
    assert set(edges) == {("1", "2"), ("1", "3"), ("2", "3")}
    assert edges[("1", "2")] == pytest.approx(EDGE_1_TO_2[store.path.name], rel=1e-2)


def test_top_edges_can_include_the_unobserved_source(store):
    edges = {(e.parent, e.child): e.weight for e in top_edges(store, unobserved=True)}
    assert edges[("1_copy", "1")] == pytest.approx(1.0)
    assert edges[("3_copy", "3")] == pytest.approx(Q_3[store.path.name], abs=1e-4)
