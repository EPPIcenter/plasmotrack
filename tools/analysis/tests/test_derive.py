"""`derive` writes a Run's derived tables under `<run>/derived/`, keyed by the closed Windows they came from."""

import json
import shutil
import sys
from collections import defaultdict

import duckdb
import pytest
from conftest import GOLDEN

from pt_analysis import __version__
from pt_analysis.derive import derive
from pt_analysis.store import Store

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


# Inf 3 per-Draw quantities pooled over the 10 Draws, from the exploration notes (R arrow).
INF_3 = {
    "single-chain": {"q": 0.1125, "expected_observed_parents": 0.8944, "entropy": 0.3697, "observed_only": 0.8875,
                     "mixed": 0.00415, "unobserved_only": 0.1084, "multi_source_augmented": 0.00688,
                     "multi_source_observed": 0.00275},
    "replica-exchange": {"q": 0.2262, "expected_observed_parents": 0.7786, "entropy": 0.1412, "observed_only": 0.7738,
                         "mixed": 0.00205, "unobserved_only": 0.2241, "multi_source_augmented": 0.00477,
                         "multi_source_observed": 0.00273},
}
CLASSES = ["observed_1", "observed_1_unobserved", "observed_2", "observed_2_unobserved", "unobserved_only"]


def infection_draws(store):
    derive(store)
    relation = duckdb.sql(f"SELECT * FROM '{store.path / 'derived' / 'infection_draws.parquet'}' "
                          "ORDER BY iteration, infection")
    return [dict(zip(relation.columns, row)) for row in relation.fetchall()]


def test_infection_draws_have_a_row_per_draw_and_infection_with_a_column_per_class(store):
    rows = infection_draws(store)
    assert list(rows[0]) == ["iteration", "infection", "q", "expected_observed_parents", "entropy", "effective_size",
                             *CLASSES, "multi_source_augmented", "multi_source_observed"]
    assert [(r["iteration"], r["infection"]) for r in rows] == [
        (it, inf) for it in range(25, 75, 5) for inf in ("1", "2", "3")]


def test_infection_draws_pool_to_the_oracle_values_for_infection_3(store):
    rows = [r for r in infection_draws(store) if r["infection"] == "3"]
    mean = lambda column: sum(r[column] for r in rows) / len(rows)  # noqa: E731
    expected = INF_3[store.path.name]
    for column in ("q", "expected_observed_parents", "entropy", "unobserved_only", "multi_source_augmented",
                   "multi_source_observed"):
        assert mean(column) == pytest.approx(expected[column], rel=5e-3), column
    assert mean("observed_1") + mean("observed_2") == pytest.approx(expected["observed_only"], rel=5e-3)
    assert mean("observed_1_unobserved") + mean("observed_2_unobserved") == pytest.approx(expected["mixed"], rel=5e-3)


def test_infection_draws_partition_each_draws_weight_into_classes(store):
    """The class columns sum to the Infection's total weight, which falls short of 1 by at most the floor for each of
    its sets the store dropped; q is the weight of the classes that include the unobserved source."""
    floor = store.manifest["parent_set_floor"]
    sets = {"1": 1, "2": 2, "3": 6}
    for r in infection_draws(store):
        assert sum(r[c] for c in CLASSES) == pytest.approx(1, abs=sets[r["infection"]] * floor)
        assert r["q"] == pytest.approx(r["observed_1_unobserved"] + r["observed_2_unobserved"] + r["unobserved_only"])
    infection_1 = [r for r in infection_draws(store) if r["infection"] == "1"]
    assert all((r["q"], r["unobserved_only"], r["entropy"], r["effective_size"]) == (1, 1, 0, 1) for r in infection_1)


def test_entropy_and_effective_size_of_a_hand_computed_infection(tmp_path):
    """Single-chain inf 2 at Draw 25 has two sets: {2_copy} with weight p = 0.9993735522251003 and {1, 2_copy} with
    1 - p = 6.264477748996795e-4. By hand: (1-p) ln(1-p) = -4.620332e-3 and p ln p = -6.262515e-4, so
    H = 5.246584e-3 and exp(H) = 1.0052604."""
    shutil.copytree(GOLDEN / "single-chain", tmp_path / "single-chain")
    rows = infection_draws(Store(tmp_path / "single-chain"))
    (row,) = [r for r in rows if (r["iteration"], r["infection"]) == (25, "2")]
    assert row["entropy"] == pytest.approx(5.246584e-3, abs=1e-8)
    assert row["effective_size"] == pytest.approx(1.0052604, abs=1e-7)
    assert row["q"] == pytest.approx(1.0)
    assert row["expected_observed_parents"] == pytest.approx(6.264477748996795e-4)
    assert row["observed_1_unobserved"] == pytest.approx(6.264477748996795e-4)
    assert row["multi_source_augmented"] == pytest.approx(6.264477748996795e-4)
    assert row["multi_source_observed"] == 0


# psi_edge per Draw, its mean and 95% interval (R type-7 quantiles), and mean q, from the exploration notes.
PSI_EDGE = {
    "single-chain": [0.3217, 0.2956, 0.3332, 0.3125, 0.2828, 0.2888, 0.3155, 0.2844, 0.2788, 0.2769],
    "replica-exchange": [0.00216, 0.00289, 0.3260, 0.3259, 0.3223, 0.3167, 0.3241, 0.3262, 0.3283, 0.3311],
}
PSI_EDGE_MEAN = {"single-chain": 0.299024, "replica-exchange": 0.260587}
PSI_EDGE_INTERVAL = {"single-chain": (0.27736, 0.33064), "replica-exchange": (0.002327, 0.3305)}
MEAN_Q = {"single-chain": 0.704166, "replica-exchange": 0.742057}


def functional_draws(store):
    derive(store)
    relation = duckdb.sql(f"SELECT * FROM '{store.path / 'derived' / 'functional_draws.parquet'}' ORDER BY iteration")
    return [dict(zip(relation.columns, row)) for row in relation.fetchall()]


def quantile_type_7(values, p):
    """R's default quantile: linear interpolation between order statistics at (n - 1) p."""
    values = sorted(values)
    h = (len(values) - 1) * p
    lo = int(h)
    return values[lo] + (h - lo) * (values[min(lo + 1, len(values) - 1)] - values[lo])


def test_functional_draws_hold_psi_edge_mean_q_mean_entropy_and_mean_classes_per_draw(store):
    rows = functional_draws(store)
    assert list(rows[0]) == ["iteration", "psi_edge", "mean_q", "mean_entropy", *(f"mean_{c}" for c in CLASSES)]
    assert [r["iteration"] for r in rows] == list(range(25, 75, 5))
    name = store.path.name
    assert [r["psi_edge"] for r in rows] == pytest.approx(PSI_EDGE[name], rel=2e-3)
    psi = [r["psi_edge"] for r in rows]
    assert sum(psi) / len(psi) == pytest.approx(PSI_EDGE_MEAN[name], rel=1e-5)
    interval = (quantile_type_7(psi, 0.025), quantile_type_7(psi, 0.975))
    assert interval == pytest.approx(PSI_EDGE_INTERVAL[name], rel=1e-3)
    assert sum(r["mean_q"] for r in rows) / len(rows) == pytest.approx(MEAN_Q[name], rel=1e-5)


def test_psi_edge_is_the_mean_over_infections_of_expected_observed_parents(store):
    by_draw = defaultdict(list)
    for r in infection_draws(store):
        by_draw[r["iteration"]].append(r)
    for f in functional_draws(store):
        infections = by_draw[f["iteration"]]
        assert len(infections) == 3
        assert f["psi_edge"] == pytest.approx(sum(r["expected_observed_parents"] for r in infections) / 3, abs=1e-15)
        assert f["mean_entropy"] == pytest.approx(sum(r["entropy"] for r in infections) / 3, abs=1e-15)
        assert f["mean_unobserved_only"] == pytest.approx(sum(r["unobserved_only"] for r in infections) / 3)


def test_out_degree_draws_give_each_infection_its_expected_observed_children_per_draw(store):
    derive(store)
    rows = read(store, "out_degree_draws", "iteration, infection")
    assert [(it, inf) for it, inf, _ in rows] == [(it, inf) for it in range(25, 75, 5) for inf in ("1", "2", "3")]
    by_draw = defaultdict(dict)
    for it, inf, out_degree in rows:
        by_draw[it][inf] = out_degree
    name = store.path.name
    # Infection 3 is no Infection's parent; Infection 1 is 3's parent with the oracle's edge 1->3, plus 2's at < 0.01.
    assert all(d["3"] == 0 for d in by_draw.values())
    edge_1_to_3 = [by_draw[it]["1"] for it in sorted(by_draw)]
    assert edge_1_to_3 == pytest.approx(EDGE_1_TO_3[name], abs=0.01)
    # Out-degrees and in-degrees count the same edges, so their mean over Infections is psi_edge.
    psi = [sum(by_draw[it].values()) / 3 for it in sorted(by_draw)]
    assert psi == pytest.approx(PSI_EDGE[name], rel=2e-3)


TABLES = ["set_classes", "edges", "set_summary", "infection_draws", "functional_draws", "genotype_summary",
          "genotype_mode", "out_degree_draws", "coi"]


# Latent genotype bitstrings per Draw, read from the golden stores: Infection 1 is 00001000 at AS1 in every Draw of
# both; Infection 3 at AS11 is 000000010000 in every single-chain Draw and in 8 of 10 replica-exchange Draws, and
# 000100000000 in the other 2.
def test_genotype_summary_holds_posterior_presence_per_allele(store):
    derive(store)
    rows = read(store, "genotype_summary", "infection, locus, allele")
    assert len(rows) == 3 * 20  # every Infection at every Allele of structure/loci (8 at AS1, 12 at AS11)
    presence = {(inf, locus, allele): p for inf, locus, allele, p in rows}
    assert [presence[("1", "AS1", k)] for k in range(8)] == [0, 0, 0, 0, 1, 0, 0, 0]
    as11 = [presence[("3", "AS11", k)] for k in range(12)]
    if store.path.name == "single-chain":
        assert as11 == [0] * 7 + [1] + [0] * 4
    else:
        assert as11 == pytest.approx([0, 0, 0, 0.2, 0, 0, 0, 0.8, 0, 0, 0, 0])


def test_genotype_mode_ranks_each_locus_bitstrings_with_the_mode_first(store):
    derive(store)
    rows = read(store, "genotype_mode", "infection, locus, rank")
    assert [r for r in rows if r[:2] == ("1", "AS1")] == [("1", "AS1", "00001000", 10, 1.0, 1)]
    got = [r for r in rows if r[:2] == ("3", "AS11")]
    if store.path.name == "single-chain":
        assert got == [("3", "AS11", "000000010000", 10, 1.0, 1)]
    else:
        assert got == [("3", "AS11", "000000010000", 8, 0.8, 1), ("3", "AS11", "000100000000", 2, 0.2, 2)]
    assert len({r[:2] for r in rows}) == 6


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
    assert read(store, "infection_draws", "iteration") == []
    assert read(store, "functional_draws", "iteration") == []
    assert read(store, "genotype_summary", "infection") == []
    assert read(store, "genotype_mode", "infection") == []
    assert read(store, "out_degree_draws", "iteration") == []
    assert len(read(store, "set_classes", "infection")) == 9
