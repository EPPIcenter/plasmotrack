"""Ranking and parent-set metrics against truth (U6, part 2), on the two golden stores as Chains of one Analysis.

Hand computation from the pooled weights in the exploration notes (sum over 20 Draws / 20):
- Edge marginals: 1->3 0.83354, 2->3 0.00296 (inf 3's sets {2,3_copy}, {1,2,3_copy}, {1,2}), 1->2 0.00291 (inf 2's
  set {1,2_copy}); 2->1, 3->1 and 3->2 have no rows, so 0.
- Infection 3's pooled set weights: {1} 0.828, {3_copy} 0.166, {1,3_copy} 0.00287, {1,2} 0.00273, {2,3_copy}
  0.000228, {1,2,3_copy} 1.2e-5. Infection 2's: {2_copy} 0.997, {1,2_copy} 0.00291. Infection 1's: {1_copy} 1.
"""

import gzip
import json
import shutil

import pytest
from conftest import GOLDEN
from test_analysis import write_manifest

from pt_analysis import metrics
from pt_analysis.analysis import Analysis


def node(id, infection_time=0.0):
    return {"id": id, "infection_time": infection_time, "flat_latent_genotype": [], "observation_time": 10.0}


def write_truth(path, edges, ids=("1", "2", "3"), gzipped=False):
    truth = {"network": [{"src": s, "dest": d} for s, d in edges], "nodes": [node(i) for i in ids]}
    data = json.dumps(truth).encode()
    path.write_bytes(gzip.compress(data) if gzipped else data)
    return path.name


@pytest.fixture
def stores(tmp_path):
    for name in ("single-chain", "replica-exchange"):
        shutil.copytree(GOLDEN / name, tmp_path / name)
    return {"chain1": "single-chain", "chain2": "replica-exchange"}


def open_analysis(tmp_path, stores, edges=(("1", "3"), ("1", "2")), name="analysis.json", **sections):
    truth = write_truth(tmp_path / "truth.json", edges) if edges is not None else None
    if truth:
        sections["truth"] = truth
    return Analysis.open(write_manifest(tmp_path / name, label="golden", runs=stores, **sections))


# --- truth -----------------------------------------------------------------------------------------------------------

def test_truth_reads_edges_parent_sets_and_realised_quantities(tmp_path):
    path = tmp_path / "truth.json"
    write_truth(path, [("1", "3"), ("1", "2")])
    truth = metrics.load_truth(path, ["1", "2", "3"])
    assert truth.edges == {("1", "3"), ("1", "2")}
    # Infection 1 is a root: its true set is the unobserved source, whose token is its own id.
    assert truth.parent_sets == {"1": frozenset({"1"}), "2": frozenset({"1"}), "3": frozenset({"1"})}
    assert truth.n_true_cases == 3
    assert truth.n_observed_cases == 3
    assert truth.true_average_outdegree == pytest.approx(2 / 3)
    assert truth.true_unobserved_source_proportion == pytest.approx(1 / 3)


def test_a_true_parent_outside_the_analysis_stands_for_the_unobserved_source(tmp_path):
    path = tmp_path / "truth.json.gz"
    write_truth(path, [("1", "3"), ("4", "2"), ("1", "2"), ("4", "5")], ids=("1", "2", "3", "4", "5"), gzipped=True)
    truth = metrics.load_truth(path, ["1", "2", "3"])
    assert truth.edges == {("1", "3"), ("1", "2")}  # 4->2 and 4->5 have an end outside the Analysis
    assert truth.parent_sets == {"1": frozenset({"1"}), "2": frozenset({"1", "2"}), "3": frozenset({"1"})}
    assert truth.n_true_cases == 5
    assert truth.n_observed_cases == 3
    assert truth.true_average_outdegree == pytest.approx(2 / 3)
    assert truth.true_unobserved_source_proportion == pytest.approx(2 / 3)  # 1 is a root, 2 has a parent outside


# --- ranking on constructed scores -----------------------------------------------------------------------------------

def test_ties_are_one_breakpoint_of_the_curve_and_all_enter_the_top_k():
    # Distinct scores 0.9 (true), 0.5 (true, false, false), 0.2 (true), 0 (false); 3 true edges, so k = 3 and the
    # 3rd highest score is 0.5, which selects all four pairs at 0.5 or above, two of them true.
    r = metrics.rank([0.9, 0.5, 0.5, 0.5, 0.2, 0.0], [True, False, True, False, True, False])
    assert [(p.score, p.tp, p.fp) for p in r.curve] == [(0.9, 1, 0), (0.5, 2, 2), (0.2, 3, 2), (0.0, 3, 3)]
    assert [p.recall for p in r.curve] == pytest.approx([1 / 3, 2 / 3, 1, 1])
    assert [p.precision for p in r.curve] == pytest.approx([1, 0.5, 0.6, 0.5])
    assert r.average_precision == pytest.approx(1 / 3 * 1 + 1 / 3 * 0.5 + 1 / 3 * 0.6)  # 0.7
    # Trapezoids from (0, 1): 1/3 * 1, 1/3 * 0.75, 1/3 * 0.55, 0.
    assert r.pr_auc == pytest.approx(0.766666666667)
    assert r.top_k_precision == pytest.approx(0.5)
    assert (r.n_edge_universe, r.n_true_edges, r.n_distinct_scores) == (6, 3, 4)


def test_scores_are_clipped_to_the_unit_interval_and_rounded_before_ties_are_found():
    r = metrics.rank([0.1 + 0.2, 0.3, 1.2, -0.1], [True, False, True, False])
    assert [(p.score, p.tp, p.fp) for p in r.curve] == [(1.0, 1, 0), (0.3, 2, 1), (0.0, 2, 2)]
    assert r.n_distinct_scores == 3


def test_without_a_true_edge_the_metrics_are_none():
    r = metrics.rank([0.5, 0.1], [False, False])
    assert (r.average_precision, r.pr_auc, r.top_k_precision, r.curve) == (None, None, None, [])
    assert (r.n_edge_universe, r.n_true_edges, r.n_distinct_scores) == (2, 0, 2)


# --- edge ranking on the golden Analysis ----------------------------------------------------------------------------

def test_edge_scores_pool_every_draw_of_every_chain(tmp_path, stores):
    scores = metrics.edge_scores(open_analysis(tmp_path, stores))
    # The per-Draw values of edge 1->3 in the notes, summed over both Chains, over 20 Draws.
    assert scores[("1", "3")] == pytest.approx((8.9138 + 7.757) / 20, abs=1e-4)
    assert scores[("2", "3")] == pytest.approx(0.00296, abs=2e-5)
    assert scores[("1", "2")] == pytest.approx(0.00291, abs=2e-5)
    assert set(scores) == {("1", "3"), ("2", "3"), ("1", "2")}  # no edge from an unobserved source


def test_directed_ranking_matches_the_hand_computed_metrics(tmp_path, stores):
    analysis = open_analysis(tmp_path, stores)
    r = metrics.edge_rankings(analysis, metrics.analysis_truth(analysis))["directed"]
    # Order: 1->3 (true), 2->3 (false), 1->2 (true), then 2->1, 3->1, 3->2 at 0 (false).
    assert [(p.tp, p.fp) for p in r.curve] == [(1, 0), (1, 1), (2, 1), (2, 4)]
    assert [p.precision for p in r.curve] == pytest.approx([1, 1 / 2, 2 / 3, 1 / 3])
    assert r.average_precision == pytest.approx(0.5 * 1 + 0.5 * 2 / 3)
    assert r.pr_auc == pytest.approx(0.5 * 1 + 0.5 * (1 / 2 + 2 / 3) / 2)
    assert r.top_k_precision == pytest.approx(0.5)  # k = 2: 1->3 and 2->3
    assert (r.n_edge_universe, r.n_true_edges, r.n_distinct_scores) == (6, 2, 4)


def test_undirected_ranking_sums_the_two_directions_and_a_pair_is_true_either_way(tmp_path, stores):
    analysis = open_analysis(tmp_path, stores, edges=[("1", "3"), ("2", "1")])
    r = metrics.edge_rankings(analysis, metrics.analysis_truth(analysis))["undirected"]
    # {1,3} 0.83354 (true), {2,3} 0.00296 (false), {1,2} 0.00291 + 0 (true: 2->1 is a true edge).
    assert [(p.score, p.tp, p.fp) for p in r.curve] == [
        (pytest.approx(0.83354, abs=1e-4), 1, 0), (pytest.approx(0.00296, abs=2e-5), 1, 1),
        (pytest.approx(0.00291, abs=2e-5), 2, 1)]
    assert r.average_precision == pytest.approx(0.5 + 0.5 * 2 / 3)
    assert r.pr_auc == pytest.approx(0.5 + 0.5 * (1 / 2 + 2 / 3) / 2)
    assert r.top_k_precision == pytest.approx(0.5)
    assert (r.n_edge_universe, r.n_true_edges, r.n_distinct_scores) == (3, 2, 3)


# --- baselines -------------------------------------------------------------------------------------------------------

def test_pairwise_scores_are_read_per_unordered_pair_and_ranked_undirected(tmp_path, stores):
    (tmp_path / "ibd.csv").write_text("a,b,score\n1,2,0.9\n2,1,0.95\n3,1,0.8\n2,3,0.1\n2,7,1.0\n")
    analysis = open_analysis(tmp_path, stores, pairwise_scores="ibd.csv")
    scores = metrics.read_pairwise_scores(analysis.pairwise_scores)
    assert scores[("1", "2")] == 0.9 and scores[("2", "1")] == 0.95
    r = metrics.pairwise_ranking(analysis, metrics.analysis_truth(analysis))
    # {1,2} 0.95 (the larger of its two rows, true), {1,3} 0.8 (true), {2,3} 0.1; 2-7 is outside the Analysis.
    assert [(p.score, p.tp, p.fp) for p in r.curve] == [(0.95, 1, 0), (0.8, 2, 0), (0.1, 2, 1)]
    assert (r.average_precision, r.pr_auc, r.top_k_precision) == (1.0, 1.0, 1.0)
    assert r.n_edge_universe == 3


def test_baselines_are_scored_like_the_posterior_and_left_out_when_not_named(tmp_path, stores):
    # The null-model Analysis here is the single-chain store alone: 1->3 0.89138, 2->3 0.00301, 1->2 0.00267.
    (tmp_path / "null").mkdir()
    write_manifest(tmp_path / "null" / "analysis.json", label="null", runs={"chain1": "../single-chain"})
    analysis = open_analysis(tmp_path, stores, null_model="null/analysis.json")
    fits = metrics.ranking_fits(analysis, metrics.analysis_truth(analysis))
    assert [(f.fit, f.label) for f in fits] == [("posterior", "golden"), ("null", "null")]
    null = fits[1]
    assert null.directed.curve[0].score == pytest.approx(8.9138 / 10, abs=1e-4)
    assert null.directed.average_precision == pytest.approx(0.5 + 0.5 * 2 / 3)

    plain = open_analysis(tmp_path, stores, name="plain.json")
    assert [f.fit for f in metrics.ranking_fits(plain, metrics.analysis_truth(plain))] == ["posterior"]


def test_the_pairwise_baseline_has_no_directed_ranking(tmp_path, stores):
    (tmp_path / "ibd.csv").write_text("a,b,score\n1,2,0.9\n")
    analysis = open_analysis(tmp_path, stores, pairwise_scores="ibd.csv")
    fits = metrics.ranking_fits(analysis, metrics.analysis_truth(analysis))
    assert [(f.fit, f.directed is None, f.undirected is None) for f in fits] == [
        ("posterior", False, False), ("ibd", True, False)]


# --- parent-set metrics ----------------------------------------------------------------------------------------------

def test_each_infection_has_effective_size_credible_size_coverage_and_top_probability(tmp_path, stores):
    analysis = open_analysis(tmp_path, stores)  # true edges 1->3 and 1->2
    rows = {r.infection: r for r in metrics.parent_set_rows(analysis, metrics.analysis_truth(analysis))}
    # Infection 1 has one set, {1_copy}, the unobserved source: its true set, as Infection 1 is a root.
    assert rows["1"] == metrics.InfectionParentSets("1", pytest.approx(1.0), 1, True, pytest.approx(1.0))
    # Infection 2's true set is {1}; its sets are {2_copy} 0.997 and {1,2_copy}, so it is never covered.
    assert (rows["2"].credible_size, rows["2"].covered) == (1, False)
    assert rows["2"].top_probability == pytest.approx(0.997, abs=1e-3)
    # Infection 3's true set {1} is its modal set, 0.828.
    assert (rows["3"].credible_size, rows["3"].covered) == (1, True)
    assert rows["3"].top_probability == pytest.approx(0.828, abs=1e-3)
    # exp of minus the sum of w log w over 0.828, 0.166, 0.00287, 0.00273, 0.000228, 1.2e-5.
    assert rows["3"].effective_size == pytest.approx(1.631, abs=2e-3)


def test_coverage_at_two_levels(tmp_path, stores):
    analysis = open_analysis(tmp_path, stores, edges=[("1", "2")])  # Infection 3 is a root: its true set is {3_copy}
    truth = metrics.analysis_truth(analysis)
    at_80 = metrics.parent_set_summary(metrics.parent_set_rows(analysis, truth, 0.8), 0.8)
    at_95 = metrics.parent_set_summary(metrics.parent_set_rows(analysis, truth, 0.95), 0.95)
    # At 0.8 Infection 3's collection is {1} alone; at 0.95 it adds {3_copy} (0.828 + 0.166 = 0.994).
    assert at_80 == {"mean_credible_parent_set_size_80": pytest.approx(1.0),
                     "median_credible_parent_set_size_80": 1,
                     "credible_parent_set_coverage_80": pytest.approx(1 / 3),
                     "mean_top_parent_set_probability": pytest.approx((1 + 0.997 + 0.828) / 3, abs=1e-3)}
    assert at_95["mean_credible_parent_set_size_a095"] == pytest.approx(4 / 3)
    assert at_95["median_credible_parent_set_size_a095"] == 1
    assert at_95["credible_parent_set_coverage_a095"] == pytest.approx(2 / 3)


def test_one_chains_credible_collection_is_the_functionals_one(tmp_path, stores):
    from pt_analysis import functionals

    analysis = open_analysis(tmp_path, stores)
    store = analysis.runs["chain2"]
    one_chain = Analysis.open(write_manifest(tmp_path / "one.json", runs={"chain2": "replica-exchange"}))
    rows = {r.infection: r for r in metrics.parent_set_rows(one_chain, level=0.8)}
    for infection in ("1", "2", "3"):
        assert rows[infection].credible_size == len(functionals.credible_collection(store, infection, 0.8))
    assert rows["3"].credible_size == 2  # the notes' 80% credible collection sizes for replica-exchange: 1, 1, 2


def test_without_truth_coverage_is_none(tmp_path, stores):
    analysis = open_analysis(tmp_path, stores, edges=None)
    rows = metrics.parent_set_rows(analysis)
    assert [r.covered for r in rows] == [None, None, None]
    assert metrics.parent_set_summary(rows)["credible_parent_set_coverage_80"] is None


def test_parent_set_fits_compare_the_posterior_with_the_null_fit(tmp_path, stores):
    (tmp_path / "null").mkdir()
    write_manifest(tmp_path / "null" / "analysis.json", label="null", runs={"chain1": "../replica-exchange"})
    analysis = open_analysis(tmp_path, stores, null_model="null/analysis.json")
    fits = metrics.parent_set_fits(analysis, metrics.analysis_truth(analysis), 0.8)
    assert [(f.fit, f.label) for f in fits] == [("posterior", "golden"), ("null", "null")]
    # replica-exchange alone: Infection 3's {1} has 0.771 < 0.8, so its collection has 2 sets.
    assert fits[1].summary["mean_credible_parent_set_size_80"] == pytest.approx(4 / 3)
    assert fits[1].summary["credible_parent_set_coverage_80"] == pytest.approx(2 / 3)


# --- fit-summary columns ---------------------------------------------------------------------------------------------

# The columns of `simulation_fit_summary.rds` these metrics fill, by name and in its order (exploration notes).
FIT_SUMMARY_COLUMNS = [
    "n_true_cases", "n_observed_cases", "true_average_outdegree", "true_unobserved_source_proportion",
    "directed_n_edge_universe", "directed_n_true_edges", "directed_n_distinct_scores", "directed_average_precision",
    "directed_pr_auc", "directed_top_k_precision", "undirected_n_edge_universe", "undirected_n_true_edges",
    "undirected_n_distinct_scores", "undirected_average_precision", "undirected_pr_auc", "undirected_top_k_precision",
    "mean_credible_parent_set_size_80", "median_credible_parent_set_size_80", "credible_parent_set_coverage_80",
    "mean_top_parent_set_probability"]


def test_fit_summary_metrics_carry_the_fit_summary_column_names(tmp_path, stores):
    analysis = open_analysis(tmp_path, stores)
    got = metrics.fit_summary_metrics(analysis, metrics.analysis_truth(analysis))
    assert list(got) == FIT_SUMMARY_COLUMNS
    assert got["n_true_cases"] == 3
    assert got["directed_n_edge_universe"] == 6
    assert got["directed_average_precision"] == pytest.approx(0.5 + 0.5 * 2 / 3)
    assert got["undirected_n_edge_universe"] == 3
    assert got["credible_parent_set_coverage_80"] == pytest.approx(2 / 3)
