"""One MAP rule and one credible-collection rule for every summary of Parent sets (`pt_analysis.parent_sets`)."""

from conftest import open_golden

from pt_analysis import functionals, metrics, strata
from pt_analysis.parent_sets import credible_flags, map_key, pooled_sets


def test_the_map_rule_ranks_by_weight_then_members():
    sets = [(("2", "3_copy"), 0.4), (("1",), 0.4), (("3_copy",), 0.2), (("1", "2"), 0.4)]
    assert [m for m, w in sorted(sets, key=lambda s: map_key(s[1], s[0]))] == [
        ("1",), ("1", "2"), ("2", "3_copy"), ("3_copy",)]


def test_the_credible_collection_is_the_prefix_whose_mass_reaches_the_level():
    assert credible_flags([0.5, 0.3, 0.2], 0.8) == [True, True, False]
    assert credible_flags([0.5, 0.29, 0.2], 0.8) == [True, True, True]
    assert credible_flags([0.8 - 1e-13, 0.2], 0.8) == [True, False]  # reached within the tolerance
    assert credible_flags([0.3, 0.3], 0.8) == [True, True]  # never reached: every set


def test_one_chains_modal_sets_and_credible_collections_agree_across_summaries(tmp_path):
    analysis = open_golden(tmp_path, runs={"chain1": "replica-exchange"})
    store = analysis.runs["chain1"]
    pooled = pooled_sets(analysis)
    graph = functionals.map_graph(store)
    assert {i: sets[0].members for i, sets in pooled.items()} == {i: s.members for i, s in graph.items()}
    for infection, sets in pooled.items():
        assert [s.members for s in sets if s.credible] == [
            s.members for s in functionals.credible_collection(store, infection)]
    assert strata.stratified_functionals(analysis)[0].psi_map == functionals.psi_map(store)
    rows = {r.infection: r for r in metrics.parent_set_rows(analysis)}
    assert {i: r.top_probability for i, r in rows.items()} == {i: s.weight for i, s in graph.items()}
