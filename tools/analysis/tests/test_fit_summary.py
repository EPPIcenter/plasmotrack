"""The per-fit summary row (U7): the 70 columns of the thesis pipeline's `simulation_fit_summary.rds`, in order and
with its types, computed for the golden stores as two Chains of one Analysis against a constructed truth file."""


import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from test_metrics import open_analysis, stores  # noqa: F401  (fixture)

from pt_analysis import fit_summary, metrics, strata

# The RDS's columns in order with their R classes (exploration notes, "simulation_fit_summary.rds").
RDS_COLUMNS = """
run_key:character simulation_model:character total_nodes:numeric founder_probability:numeric
superinfection_probability:numeric secondary_case_rate:numeric founder_coi_rate:numeric
asymptomatic_probability:numeric panel_name:character panel_alleles_per_locus:numeric simulation_index:integer
seed:integer missingness_condition:character fit_type:character n_true_cases:integer n_observed_cases:integer
true_average_outdegree:numeric true_unobserved_source_proportion:numeric estimated_average_outdegree:numeric
average_outdegree_lower:numeric average_outdegree_upper:numeric estimated_unobserved_source_proportion:numeric
unobserved_source_lower:numeric unobserved_source_upper:numeric mean_observed_only_probability:numeric
mean_observed_only_probability_lower:numeric mean_observed_only_probability_upper:numeric
mean_mixed_observed_unobserved_probability:numeric mean_mixed_observed_unobserved_probability_lower:numeric
mean_mixed_observed_unobserved_probability_upper:numeric mean_unobserved_only_probability:numeric
mean_unobserved_only_probability_lower:numeric mean_unobserved_only_probability_upper:numeric
mean_no_parent_probability:numeric mean_no_parent_probability_lower:numeric mean_no_parent_probability_upper:numeric
mean_multi_source_probability_augmented:numeric mean_multi_source_probability_augmented_lower:numeric
mean_multi_source_probability_augmented_upper:numeric mean_multi_source_probability_observed_only:numeric
mean_multi_source_probability_observed_only_lower:numeric mean_multi_source_probability_observed_only_upper:numeric
mean_parent_set_entropy:numeric mean_parent_set_entropy_lower:numeric mean_parent_set_entropy_upper:numeric
mean_effective_parent_sets:numeric mean_effective_parent_sets_lower:numeric mean_effective_parent_sets_upper:numeric
directed_n_edge_universe:integer directed_n_true_edges:integer directed_n_distinct_scores:integer
directed_average_precision:numeric directed_pr_auc:numeric directed_top_k_precision:numeric
undirected_n_edge_universe:integer undirected_n_true_edges:integer undirected_n_distinct_scores:integer
undirected_average_precision:numeric undirected_pr_auc:numeric undirected_top_k_precision:numeric
mean_credible_parent_set_size_80:numeric median_credible_parent_set_size_80:numeric
credible_parent_set_coverage_80:numeric mean_top_parent_set_probability:numeric average_outdegree_error:numeric
average_outdegree_abs_error:numeric average_outdegree_sq_error:numeric unobserved_source_error:numeric
unobserved_source_abs_error:numeric unobserved_source_sq_error:numeric
""".split()
ARROW_TYPE = {"character": pa.string(), "integer": pa.int64(), "numeric": pa.float64()}

FACTORS = {"run_key": "biological__fp0.05__sp0__r00.5__coi1__high__sim0__s0", "simulation_model": "biological",
           "total_nodes": 200, "founder_probability": 0.05, "superinfection_probability": 0,
           "secondary_case_rate": 0.5, "founder_coi_rate": 1, "asymptomatic_probability": 0.5,
           "panel_name": "high", "panel_alleles_per_locus": 20, "simulation_index": 0, "seed": 0,
           "missingness_condition": "full", "fit_type": "posterior"}

# Pooled per-Draw psi_edge of the two golden stores (exploration notes): single-chain then replica-exchange.
PSI_EDGE = [0.3217, 0.2956, 0.3332, 0.3125, 0.2828, 0.2888, 0.3155, 0.2844, 0.2788, 0.2769,
            0.00216, 0.00289, 0.3260, 0.3259, 0.3223, 0.3167, 0.3241, 0.3262, 0.3283, 0.3311]


@pytest.fixture
def analysis(tmp_path, stores):  # noqa: F811
    return open_analysis(tmp_path, stores, factors=FACTORS)


def test_the_schema_is_the_rds_columns_in_order_with_their_types():
    assert fit_summary.SCHEMA.names == [c.split(":")[0] for c in RDS_COLUMNS]
    assert [str(t) for t in fit_summary.SCHEMA.types] == [str(ARROW_TYPE[c.split(":")[1]]) for c in RDS_COLUMNS]


def test_the_row_written_has_the_schema(analysis, tmp_path):
    out = tmp_path / "row.parquet"
    fit_summary.write_row(fit_summary.summarize(analysis), out)
    table = pq.read_table(out)
    assert table.schema.remove_metadata() == fit_summary.SCHEMA
    assert table.num_rows == 1


def quantile7(values, p):
    ordered = sorted(values)
    h = (len(ordered) - 1) * p
    lo = int(h)
    return ordered[lo] + (h - lo) * (ordered[min(lo + 1, len(ordered) - 1)] - ordered[lo])


def test_factors_come_from_the_manifest(analysis):
    row = fit_summary.summarize(analysis)
    assert {k: row[k] for k in FACTORS} == FACTORS
    assert isinstance(row["total_nodes"], float) and isinstance(row["seed"], int)


def test_psi_edge_and_mean_q_pool_the_draws_of_both_chains(analysis):
    row = fit_summary.summarize(analysis)
    # Means from the notes: psi_edge 0.299024 and 0.260587, mean q 0.704166 and 0.742057, 10 Draws each.
    assert row["estimated_average_outdegree"] == pytest.approx((0.299024 + 0.260587) / 2, abs=1e-5)
    assert row["average_outdegree_lower"] == pytest.approx(quantile7(PSI_EDGE, 0.025), abs=1e-3)
    assert row["average_outdegree_upper"] == pytest.approx(quantile7(PSI_EDGE, 0.975), abs=1e-3)
    assert row["estimated_unobserved_source_proportion"] == pytest.approx((0.704166 + 0.742057) / 2, abs=1e-5)


def test_functionals_agree_with_strata(analysis):
    row = fit_summary.summarize(analysis)
    (overall,) = strata.stratified_functionals(analysis)
    for column, prefix, interval in (
            ("estimated_average_outdegree", "average_outdegree", overall.psi_edge),
            ("estimated_unobserved_source_proportion", "unobserved_source", overall.mean_q),
            ("mean_unobserved_only_probability", "mean_unobserved_only_probability",
             overall.classes["unobserved_only"]),
            ("mean_multi_source_probability_augmented", "mean_multi_source_probability_augmented",
             overall.multi_source_augmented),
            ("mean_multi_source_probability_observed_only", "mean_multi_source_probability_observed_only",
             overall.multi_source_observed)):
        assert (row[column], row[f"{prefix}_lower"], row[f"{prefix}_upper"]) == pytest.approx(
            (interval.mean, interval.lower, interval.upper))
    observed_only = overall.classes["observed_1"].mean + overall.classes["observed_2"].mean
    mixed = overall.classes["observed_1_unobserved"].mean + overall.classes["observed_2_unobserved"].mean
    assert row["mean_observed_only_probability"] == pytest.approx(observed_only)
    assert row["mean_mixed_observed_unobserved_probability"] == pytest.approx(mixed)


def test_composition_classes_match_the_notes(analysis):
    row = fit_summary.summarize(analysis)
    # Infections 1 and 2 are unobserved-only at about 1; infection 3's classes (notes), averaged over 3 Infections.
    assert row["mean_observed_only_probability"] == pytest.approx((0.8875 + 0.7738) / 2 / 3 + 0.00267 / 3, abs=2e-3)
    assert row["mean_no_parent_probability"] == row["mean_no_parent_probability_upper"] == 0.0


def test_entropy_is_the_mean_over_infections_and_draws(analysis):
    row = fit_summary.summarize(analysis)
    values = []
    for store in analysis.runs.values():
        values += [h for (h,) in store.connection.sql(
            "SELECT avg(entropy) FROM read_parquet(?) GROUP BY iteration ORDER BY iteration",
            params=[str(store.path / "derived" / "infection_draws.parquet")]).fetchall()]
    assert row["mean_parent_set_entropy"] == pytest.approx(sum(values) / len(values))
    assert row["mean_parent_set_entropy_upper"] == pytest.approx(quantile7(values, 0.975))
    assert row["mean_effective_parent_sets"] >= 1.0


def test_truth_ranking_and_parent_set_columns_are_the_metrics(analysis):
    row = fit_summary.summarize(analysis)
    expected = metrics.fit_summary_metrics(analysis, metrics.analysis_truth(analysis))
    assert len(expected) == 20
    assert {k: row[k] for k in expected} == expected
    assert row["true_average_outdegree"] == pytest.approx(2 / 3)


def test_errors_are_estimate_minus_truth(analysis):
    row = fit_summary.summarize(analysis)
    error = row["estimated_average_outdegree"] - 2 / 3
    assert row["average_outdegree_error"] == pytest.approx(error)
    assert row["average_outdegree_abs_error"] == pytest.approx(abs(error))
    assert row["average_outdegree_sq_error"] == pytest.approx(error ** 2)
    assert row["unobserved_source_error"] == pytest.approx(row["estimated_unobserved_source_proportion"] - 1 / 3)


def test_without_truth_or_factors_those_columns_are_null(tmp_path, stores):  # noqa: F811
    row = fit_summary.summarize(open_analysis(tmp_path, stores, edges=None))
    assert row["run_key"] is None and row["directed_pr_auc"] is None and row["average_outdegree_error"] is None
    assert row["estimated_average_outdegree"] is not None
    assert row["credible_parent_set_coverage_80"] is None and row["mean_credible_parent_set_size_80"] is not None


def test_edge_rows_hold_scored_and_true_pairs(analysis):
    table = fit_summary.edge_rows(analysis)
    assert table.schema == fit_summary.EDGE_SCHEMA
    rows = {(r["parent"], r["child"]): (r["score"], r["truth"]) for r in table.to_pylist()}
    # Scored: 1->3 0.83354, 2->3 0.00296, 1->2 0.00291; true: 1->3 and 1->2 (test_metrics' hand computation).
    assert set(rows) == {("1", "3"), ("2", "3"), ("1", "2")}
    assert rows[("1", "3")][0] == pytest.approx(0.83354, abs=1e-4) and rows[("1", "3")][1]
    assert not rows[("2", "3")][1]
    assert set(table.column("run_key").to_pylist()) == {FACTORS["run_key"]}


def test_the_row_file_names_its_manifest(analysis, tmp_path):
    out = tmp_path / "row.parquet"
    fit_summary.write_row(fit_summary.summarize(analysis), out, analysis.path)
    assert pq.read_schema(out).metadata[fit_summary.MANIFEST_KEY].decode() == str(analysis.path)
