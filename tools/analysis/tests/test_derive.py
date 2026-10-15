"""`derive` writes a Run's derived tables under `<run>/derived/`, keyed by the closed Windows they came from."""

import json
import sys
from collections import defaultdict

import duckdb
import pytest
from conftest import GOLDEN

from pt_analysis import __version__
from pt_analysis.derive import derive

sys.path.insert(0, str(GOLDEN.parents[2] / "tools" / "benchmark"))
import mixing_benchmark as mb  # noqa: E402


def read(store, table, order):
    return duckdb.sql(f"SELECT * FROM '{store.path / 'derived' / f'{table}.parquet'}' ORDER BY {order}").fetchall()


def test_set_classes_label_every_set_by_observed_count_and_unobserved_source(store):
    derive(store)
    assert read(store, "set_classes", "infection, set_id") == [
        ("1", 0, 0, True, "unobserved_only"),
        ("2", 0, 0, True, "unobserved_only"),
        ("2", 1, 1, True, "observed_1_unobserved"),
        ("3", 0, 0, True, "unobserved_only"),
        ("3", 1, 1, True, "observed_1_unobserved"),
        ("3", 2, 1, True, "observed_1_unobserved"),
        ("3", 3, 2, True, "observed_2_unobserved"),
        ("3", 4, 1, False, "observed_1"),
        ("3", 5, 2, False, "observed_2"),
    ]


# From the exploration notes, computed independently with R arrow.
EDGE_1_TO_3 = {
    "single-chain": [0.9638, 0.8808, 0.9799, 0.9337, 0.8457, 0.8624, 0.9414, 0.8481, 0.8313, 0.8267],
    "replica-exchange": [0.0031, 0.0039, 0.9757, 0.9767, 0.9533, 0.9499, 0.9657, 0.9705, 0.9718, 0.9864],
}


def test_edges_hold_each_draws_edge_marginal(store):
    derive(store)
    edges = read(store, "edges", "iteration")
    edge_1_to_3 = [(it, round(prob, 4)) for it, child, parent, prob in edges if (child, parent) == ("3", "1")]
    assert edge_1_to_3 == list(zip(range(25, 75, 5), EDGE_1_TO_3[store.path.name]))


def test_edges_name_the_latent_parent_as_the_unobserved_source(store):
    derive(store)
    edges = read(store, "edges", "iteration, child, parent")
    assert {(child, parent) for _, child, parent, _ in edges} == {
        ("1", "1_copy"), ("2", "1"), ("2", "2_copy"), ("3", "1"), ("3", "2"), ("3", "3_copy")}
    assert all(prob == 1.0 for _, child, _, prob in edges if child == "1")
    assert {it for it, *_ in edges} == set(range(25, 75, 5))  # Draws only: no burn-in iteration


def test_edges_match_the_mixing_benchmarks_edge_traces(store):
    derive(store)
    groups = defaultdict(dict)
    mb.parent_set_traces(store, groups)
    draws = list(range(25, 75, 5))
    series = defaultdict(lambda: dict.fromkeys(draws, 0.0))
    for it, child, parent, prob in read(store, "edges", "iteration"):
        if not parent.endswith("_copy"):
            series[f"{child}<-{parent}"][it] = prob
    assert sorted(series) == sorted(groups["_edges"])
    for name, benchmark in groups["_edges"].items():
        assert list(series[name].values()) == pytest.approx(benchmark, abs=1e-12), name


# Pooled set weights (sum over Draws / 10 Draws) by Infection, in rank order, from the exploration notes.
POOLED_WEIGHTS = {
    "single-chain": {"1": [(0, 1.0)], "2": [(0, 0.997), (1, 0.00267)],
                     "3": [(4, 0.885), (0, 0.108), (1, 0.00387), (5, 0.00273), (2, 0.000269), (3, 1.5e-5)]},
    "replica-exchange": {"1": [(0, 1.0)], "2": [(0, 0.997), (1, 0.00315)],
                         "3": [(4, 0.771), (0, 0.224), (5, 0.00272), (1, 0.00186), (2, 0.000187), (3, 7.8e-6)]},
}


def test_set_summary_ranks_sets_by_pooled_weight_over_all_draws(store):
    derive(store)
    rows = read(store, "set_summary", "infection, rank")
    for infection, expected in POOLED_WEIGHTS[store.path.name].items():
        got = [(set_id, weight, rank, cumulative) for inf, set_id, weight, rank, cumulative in rows if inf == infection]
        assert [g[0] for g in got] == [set_id for set_id, _ in expected]
        assert [g[1] for g in got] == pytest.approx([w for _, w in expected], rel=1e-2)
        assert [g[2] for g in got] == list(range(1, len(expected) + 1))
        assert [g[3] for g in got] == pytest.approx([sum(g[1] for g in got[:k + 1]) for k in range(len(got))])


TABLES = ["set_classes", "edges", "set_summary"]


def test_the_manifest_records_the_closed_windows_per_phase_and_the_package_version(store):
    assert derive(store) is True
    assert json.loads((store.path / "derived" / "manifest.json").read_text()) == {
        "windows": {"burn_in": [0], "sampling": [1, 2]}, "package_version": __version__}
    assert sorted(p.name for p in (store.path / "derived").iterdir()) == sorted(
        [f"{t}.parquet" for t in TABLES] + ["manifest.json"])


def test_a_second_derive_writes_nothing(store):
    derive(store)
    derived = store.path / "derived"
    before = {p.name: p.stat().st_mtime_ns for p in derived.iterdir()}
    assert derive(store) is False
    assert {p.name: p.stat().st_mtime_ns for p in derived.iterdir()} == before


def test_a_store_that_gained_a_window_is_derived_again(store):
    """Derive with Window 2's marker set aside, as if the Run were still writing it, then put the marker back."""
    marker = store.path / "_closed" / "window-000002"
    aside = store.path / "window-000002.aside"
    marker.rename(aside)
    assert derive(store) is True
    assert {it for it, *_ in read(store, "edges", "iteration")} == set(range(25, 50, 5))

    aside.rename(marker)
    assert derive(store) is True
    manifest = json.loads((store.path / "derived" / "manifest.json").read_text())
    assert manifest["windows"] == {"burn_in": [0], "sampling": [1, 2]}
    assert {it for it, *_ in read(store, "edges", "iteration")} == set(range(25, 75, 5))
    weights = {(inf, set_id): w for inf, set_id, w, *_ in read(store, "set_summary", "infection, rank")}
    assert weights[("3", 4)] == pytest.approx(POOLED_WEIGHTS[store.path.name]["3"][0][1], rel=1e-2)
    assert sorted(p.name for p in (store.path / "derived").iterdir()) == sorted(
        [f"{t}.parquet" for t in TABLES] + ["manifest.json"])


def test_a_run_still_in_burn_in_has_no_edges_or_pooled_weights(store):
    for w in (1, 2):
        (store.path / "_closed" / f"window-{w:06d}").unlink()
    assert derive(store) is True
    assert read(store, "edges", "iteration") == []
    assert read(store, "set_summary", "infection") == []
    assert len(read(store, "set_classes", "infection")) == 9
