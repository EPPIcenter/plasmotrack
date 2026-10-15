"""Flows, lags and spatial scale per stratum of an epi column (U9), on the golden stores as a two-Chain
Analysis with a constructed epi CSV. Expected values are hand-computed from the exploration notes' R arrow oracles:
pooled over both Chains' 20 Draws, edge 1->3 is (0.89138 + 0.7757) / 2, edge 1->2 is (0.00267 + 0.00315) / 2, edge
2->3 is the weight of Infection 3's sets holding 2, and q is 1, 1 and (0.1125 + 0.2262) / 2.

The example input's observation times are 234.0737, 306.2297 and 350.5, so the lags are 1->2 72.156, 1->3 116.4263
and 2->3 44.2703; its Candidate parents are 2: [1] and 3: [1, 2].
"""

import pytest
from test_strata import open_analysis

from pt_analysis.epi import EpiError
from pt_analysis.flows import (FlowsError, attribution, enrichment, lag_histogram, median_duration, net_flow,
                               reassign, spatial_scale)

E13 = (0.89138 + 0.7757) / 2
E12 = (0.00267 + 0.00315) / 2
E23 = ((0.000269 + 1.5e-5 + 0.00273) + (0.000187 + 7.8e-6 + 0.00272)) / 2
Q = {"1": 1.0, "2": 1.0, "3": (0.1125 + 0.2262) / 2}
PSI_EDGE = (0.299024 + 0.260587) / 2
MEAN_Q = (0.704166 + 0.742057) / 2

# Two strata of district, North {1, 3} and South {2}; shehia nested in district, with shehia A in both districts.
EPI = "infection,district,shehia\n1,North,A\n2,South,B\n3,North,A\n"
NESTED = "infection,district,shehia\n1,North,A\n2,North,B\n3,South,A\n"


@pytest.fixture
def analysis(tmp_path):
    return open_analysis(tmp_path, "both", EPI)


def cells(rows, key):
    return {key(r): r for r in rows}


def test_attribution_rows_sum_to_one_per_recipient_stratum(analysis):
    rows = attribution(analysis, "district")
    by = cells(rows, lambda r: (r.recipient, r.source, r.unobserved))
    for recipient in ("North", "South"):
        assert sum(r.share for r in rows if r.recipient == recipient) == pytest.approx(1.0)
    # North's parent mass: edge 1->3 from North, 2->3 from South, q1 + q3 unobserved.
    north = E13 + E23 + Q["1"] + Q["3"]
    assert by[("North", "North", False)].mass == pytest.approx(E13, rel=1e-3)
    assert by[("North", "South", False)].mass == pytest.approx(E23, rel=1e-2)
    assert by[("North", None, True)].mass == pytest.approx(Q["1"] + Q["3"], rel=1e-3)
    assert by[("North", None, True)].share == pytest.approx((Q["1"] + Q["3"]) / north, rel=1e-3)
    assert by[("South", "North", False)].mass == pytest.approx(E12, rel=1e-2)
    assert by[("South", "South", False)].mass == 0


def test_enrichment_divides_by_the_share_of_candidate_pairs_in_each_cell(analysis):
    by = cells(enrichment(analysis, "district"), lambda r: (r.source, r.recipient))
    # Candidate pairs: 1->2 North->South, 1->3 North->North, 2->3 South->North, a third of them each.
    assert {k: r.candidates for k, r in by.items()} == {("North", "North"): 1, ("North", "South"): 1,
                                                      ("South", "North"): 1, ("South", "South"): 0}
    total = E13 + E12 + E23
    assert by[("North", "North")].expected_edges == pytest.approx(E13, rel=1e-3)
    assert by[("North", "North")].enrichment == pytest.approx((E13 / total) / (1 / 3), rel=1e-3)
    assert by[("North", "South")].enrichment == pytest.approx((E12 / total) / (1 / 3), rel=1e-2)
    assert by[("South", "South")].enrichment is None
    assert sum(r.share for r in by.values()) == pytest.approx(1.0)
    assert sum(r.candidate_share for r in by.values()) == pytest.approx(1.0)


def test_net_flow_counts_only_cross_stratum_edges(analysis):
    rows, hidden = net_flow(analysis, "district")
    by = cells(rows, lambda r: r.value)
    assert hidden == 0
    # 1->3 stays within North and counts in neither direction.
    assert (by["North"].out_edges, by["North"].in_edges) == pytest.approx((E12, E23), rel=1e-2)
    assert by["North"].net == pytest.approx(E12 - E23, abs=1e-4)
    assert by["South"].net == pytest.approx(E23 - E12, abs=1e-4)
    assert (by["North"].n, by["South"].n) == (2, 1)


def test_net_flow_hides_strata_under_the_minimum_infection_count(analysis):
    rows, hidden = net_flow(analysis, "district", min_infections=2)
    assert [r.value for r in rows] == ["North"]
    assert hidden == 1


def test_the_lag_histogram_bins_edge_mass_by_child_minus_parent_observation_time(analysis):
    bins = lag_histogram(analysis, width=10)
    assert sum(b.mass for b in bins) == pytest.approx(E13 + E12 + E23, rel=1e-3)
    holding = {lag: next(b for b in bins if b.x0 <= lag < b.x1) for lag in (44.2703, 72.156, 116.4263)}
    assert holding[116.4263].mass == pytest.approx(E13, rel=1e-3)
    assert holding[72.156].mass == pytest.approx(E12, rel=1e-2)
    assert all(b.x1 - b.x0 == pytest.approx(10) for b in bins)


def test_the_median_fitted_duration_is_the_median_of_each_infections_posterior_mean(analysis):
    # R arrow: per-Infection mean infection_duration over both Chains' Draws is 11.2, 14.03, 24.5.
    assert median_duration(analysis) == pytest.approx(14.02945, rel=1e-5)


def test_reassigning_mass_beyond_a_lag_moves_edge_1_to_3_to_the_unobserved_source(analysis):
    result = reassign(analysis, max_lag=100)
    assert result.moved_mass == pytest.approx(E13, rel=1e-3)
    assert result.moved_share == pytest.approx(E13 / (E13 + E12 + E23), rel=1e-3)
    assert result.overall.psi_edge.mean == pytest.approx(PSI_EDGE, rel=1e-5)
    assert result.overall.mean_q.mean == pytest.approx(MEAN_Q, rel=1e-5)
    assert result.overall.psi_edge_moved.mean == pytest.approx(PSI_EDGE - E13 / 3, rel=1e-2)
    assert result.overall.mean_q_moved.mean == pytest.approx(MEAN_Q + E13 / 3, rel=1e-3)


def test_reassigned_q_is_clamped_to_one_per_infection(analysis):
    # Beyond lag 0 every edge moves. Unclamped, Infection 3's q plus its moved in-edges exceeds 1 (a set holding 1
    # and 3_copy counts in both), as does Infection 2's (q 1 plus edge 1->2); clamped, all three reach exactly 1.
    moved = reassign(analysis, max_lag=0).overall.mean_q_moved
    assert (moved.mean, moved.lower, moved.upper) == (pytest.approx(1.0), pytest.approx(1.0), pytest.approx(1.0))
    assert moved.upper <= 1.0


def test_reassigning_beyond_the_largest_lag_moves_nothing(analysis):
    result = reassign(analysis, max_lag=116.4263)
    assert result.moved_mass == 0
    assert result.overall.psi_edge_moved == result.overall.psi_edge


def test_reassigned_functionals_per_stratum_count_sources_for_psi_edge_and_recipients_for_q(analysis):
    by = cells(reassign(analysis, max_lag=100, column="district").strata, lambda s: s.value)
    # Edge 1->3: its source 1 and its recipient 3 are both North, n = 2.
    assert by["North"].psi_edge.mean - by["North"].psi_edge_moved.mean == pytest.approx(E13 / 2, rel=1e-3)
    # q rises by the moved mass less the weight that would push Infection 3's q past 1: its sets holding both 1 and
    # 3_copy, {1, 3_copy} 0.00287 and {1, 2, 3_copy} 1.2e-5 pooled (test_metrics), counted in q and in edge 1->3.
    assert by["North"].mean_q_moved.mean - by["North"].mean_q.mean == pytest.approx(
        (E13 - 0.00287 - 1.2e-5) / 2, rel=1e-3)
    assert by["South"].psi_edge_moved == by["South"].psi_edge
    assert by["South"].mean_q_moved == by["South"].mean_q


def test_spatial_shares_sum_to_one_and_name_the_relation_of_source_to_recipient(tmp_path):
    analysis = open_analysis(tmp_path, "both", EPI)
    rows = {r.relation: r for r in spatial_scale(analysis, unit="shehia", parent="district")}
    assert list(rows) == ["same_unit", "same_parent", "other"]
    assert sum(r.share for r in rows.values()) == pytest.approx(1.0)
    # 1->3 within shehia A; 1->2 and 2->3 cross districts.
    assert rows["same_unit"].mass == pytest.approx(E13, rel=1e-3)
    assert rows["same_parent"].mass == 0
    assert rows["other"].mass == pytest.approx(E12 + E23, rel=1e-2)


def test_a_unit_is_known_by_its_parent_unit_too(tmp_path):
    # Shehia A is in North for 1 and in South for 3: 1->3 is in another district, 1->2 in the same one.
    analysis = open_analysis(tmp_path, "both", NESTED)
    rows = {r.relation: r for r in spatial_scale(analysis, unit="shehia", parent="district")}
    assert rows["same_unit"].mass == 0
    assert rows["same_parent"].mass == pytest.approx(E12, rel=1e-2)
    assert rows["other"].mass == pytest.approx(E13 + E23, rel=1e-2)


def test_spatial_shares_keep_only_edges_up_to_the_lag(analysis):
    rows = {r.relation: r for r in spatial_scale(analysis, unit="shehia", parent="district", max_lag=100)}
    assert rows["same_unit"].mass == 0
    assert rows["other"].share == pytest.approx(1.0)


def test_an_unknown_column_names_itself(analysis):
    with pytest.raises(EpiError, match="parish"):
        attribution(analysis, "parish")
    with pytest.raises(EpiError, match="parish"):
        spatial_scale(analysis, unit="parish", parent="district")


def test_a_unit_cannot_be_its_own_parent(analysis):
    with pytest.raises(FlowsError):
        spatial_scale(analysis, unit="district", parent="district")
