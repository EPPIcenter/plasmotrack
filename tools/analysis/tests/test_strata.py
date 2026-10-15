"""Functionals overall and per stratum of an epi column, on the golden stores as one-Chain and two-Chain Analyses,
against values computed with R arrow (exploration notes)."""

import pytest
from conftest import open_golden

from pt_analysis.epi import EpiError
from pt_analysis.strata import Interval, StrataError, interval, stratified_functionals, stratum_series

RUNS = {"single": {"chain1": "single-chain"}, "replica": {"chain1": "replica-exchange"},
        "both": {"chain1": "single-chain", "chain2": "replica-exchange"}}

PSI_EDGE_MEAN = {"single": 0.299024, "replica": 0.260587, "both": (0.299024 + 0.260587) / 2}
PSI_EDGE_INTERVAL = {"single": (0.27736, 0.33064), "replica": (0.002327, 0.3305)}
MEAN_Q = {"single": 0.704166, "replica": 0.742057, "both": (0.704166 + 0.742057) / 2}
PSI_EDGE = ([0.3217, 0.2956, 0.3332, 0.3125, 0.2828, 0.2888, 0.3155, 0.2844, 0.2788, 0.2769],
            [0.00216, 0.00289, 0.3260, 0.3259, 0.3223, 0.3167, 0.3241, 0.3262, 0.3283, 0.3311])

EPI = "infection,district,observed\n1,North,2020-01-05\n2,South,2020-02-01\n3,North,2020-03-10\n9,South,2020-02-02\n"


def open_analysis(tmp_path, runs="both", epi=EPI):
    sections = {}
    if epi is not None:
        (tmp_path / "epi.csv").write_text(epi)
        sections["epi"] = "epi.csv"
    return open_golden(tmp_path, f"{runs}.json", runs=RUNS[runs], **sections)


def test_an_interval_is_the_mean_and_type_7_quantiles():
    # Type 7: the quantile at p lies (n - 1) p of the way along the sorted values, interpolated linearly.
    got = interval([4.0, 1.0, 3.0, 2.0, 10.0])
    assert got == Interval(mean=4.0, lower=pytest.approx(1.1), upper=pytest.approx(9.4))
    assert interval([2.5]) == Interval(2.5, 2.5, 2.5)
    assert interval([]) is None


@pytest.mark.parametrize("runs", ["single", "replica", "both"])
def test_overall_psi_edge_and_mean_q_pool_every_draw_of_every_chain(tmp_path, runs):
    (overall,) = stratified_functionals(open_analysis(tmp_path, runs))
    assert (overall.value, overall.n) == (None, 3)
    assert overall.psi_edge.mean == pytest.approx(PSI_EDGE_MEAN[runs], rel=1e-5)
    assert overall.mean_q.mean == pytest.approx(MEAN_Q[runs], rel=1e-5)
    if runs in PSI_EDGE_INTERVAL:
        assert (overall.psi_edge.lower, overall.psi_edge.upper) == pytest.approx(PSI_EDGE_INTERVAL[runs], rel=1e-3)


def test_the_pooled_interval_is_over_the_draws_of_both_chains(tmp_path):
    (overall,) = stratified_functionals(open_analysis(tmp_path, "both"))
    pooled = sorted(PSI_EDGE[0] + PSI_EDGE[1])  # 20 Draws: (n - 1) p = 0.475 and 18.525
    lower = pooled[0] + 0.475 * (pooled[1] - pooled[0])
    upper = pooled[18] + 0.525 * (pooled[19] - pooled[18])
    assert (overall.psi_edge.lower, overall.psi_edge.upper) == pytest.approx((lower, upper), rel=2e-3)


def weighted(strata, get):
    return sum(s.n * get(s) for s in strata) / sum(s.n for s in strata)


FUNCTIONALS = {"psi_map": lambda s: s.psi_map, "psi_lo": lambda s: s.size_range[0],
               "psi_hi": lambda s: s.size_range[1], "psi_edge": lambda s: s.psi_edge.mean, "mean_q": lambda s: s.mean_q.mean,
               "multi_source_augmented": lambda s: s.multi_source_augmented.mean,
               "multi_source_observed": lambda s: s.multi_source_observed.mean,
               **{c: (lambda c: lambda s: s.classes[c].mean)(c) for c in
                  ("observed_1", "observed_1_unobserved", "observed_2", "observed_2_unobserved", "unobserved_only")}}


@pytest.mark.parametrize("column", ["district", "observed"])
def test_strata_partition_the_infections_and_average_back_to_overall(tmp_path, column):
    analysis = open_analysis(tmp_path)
    (overall,) = stratified_functionals(analysis)
    strata = stratified_functionals(analysis, column)
    assert sum(s.n for s in strata) == overall.n == 3
    for name, get in FUNCTIONALS.items():
        assert weighted(strata, get) == pytest.approx(get(overall), abs=1e-12), name


def test_out_degree_counts_an_infection_as_a_source_and_inclusion_as_a_recipient(tmp_path):
    strata = {s.value: s for s in stratified_functionals(open_analysis(tmp_path, "single"), "district")}
    assert [(v, s.n) for v, s in strata.items()] == [("North", 2), ("South", 1)]
    north, south = strata["North"], strata["South"]
    # Infection 1 is the parent of 3 (pooled 0.89138) and of 2 (0.00267); 3 is no Infection's parent.
    assert north.psi_edge.mean == pytest.approx((0.89138 + 0.00267) / 2, abs=1e-4)
    assert south.psi_edge.mean == pytest.approx(0.003, abs=5e-4)  # 2 -> 3
    # Infection 1 has only the unobserved source; 3's inclusion is 0.1125; 2's is 1.
    assert north.mean_q.mean == pytest.approx((1 + 0.1125) / 2, abs=1e-4)
    assert south.mean_q.mean == pytest.approx(1, abs=1e-12)
    assert north.classes["unobserved_only"].mean == pytest.approx((1 + 0.1084) / 2, abs=1e-4)
    assert north.multi_source_augmented.mean == pytest.approx(0.00688 / 2, abs=1e-5)
    assert north.multi_source_observed.mean == pytest.approx(0.00275 / 2, abs=1e-5)


def test_infections_with_no_value_are_a_stratum_of_their_own(tmp_path):
    analysis = open_analysis(tmp_path, epi="infection,district\n1,North\n2,\n")
    assert [(s.value, s.n) for s in stratified_functionals(analysis, "district")] == [("North", 1), (None, 2)]


def test_composition_classes_partition_each_stratums_weight(tmp_path):
    for s in stratified_functionals(open_analysis(tmp_path), "district"):
        assert list(s.classes) == ["observed_1", "observed_1_unobserved", "observed_2", "observed_2_unobserved",
                                   "unobserved_only"]
        # Each class column sums to the Infection's weight, short of 1 by at most the floor (1e-9) per dropped set.
        assert sum(c.mean for c in s.classes.values()) == pytest.approx(1, abs=1e-6)


def test_a_missing_column_is_a_named_error(tmp_path):
    with pytest.raises(EpiError, match="no column 'travel'; its columns are district, observed"):
        stratified_functionals(open_analysis(tmp_path), "travel")
    with pytest.raises(EpiError, match="no epi data"):
        stratified_functionals(open_analysis(tmp_path / "none", epi=None), "district")


SIZE_RANGE_80 = {"single": (1 / 3, 1 / 3), "replica": (0, 1 / 3),
                 # Pooled over both Chains, Infection 3's modal set {1} has 0.828 alone.
                 "both": (1 / 3, 1 / 3)}


@pytest.mark.parametrize("runs", ["single", "replica", "both"])
def test_overall_psi_map_and_size_range_come_from_weights_pooled_over_chains(tmp_path, runs):
    (overall,) = stratified_functionals(open_analysis(tmp_path, runs))
    assert overall.psi_map == pytest.approx(1 / 3)
    assert overall.size_range == pytest.approx(SIZE_RANGE_80[runs])


def test_psi_map_counts_an_infection_as_a_source_and_the_size_range_as_a_recipient(tmp_path):
    strata = {s.value: s for s in stratified_functionals(open_analysis(tmp_path, "replica"), "district")}
    # The MAP graph has one edge, 1 -> 3. North holds 1 and 3; South holds 2.
    assert (strata["North"].psi_map, strata["South"].psi_map) == (0.5, 0)
    # At 0.8 Infection 3's collection is {1} and {3_copy} (counts 1 and 0); 1 and 2 have one set with no parent.
    assert strata["North"].size_range == pytest.approx((0, 0.5))
    assert strata["South"].size_range == (0, 0)
    assert stratified_functionals(open_analysis(tmp_path, "replica"), "district", level=0.7)[0].size_range == \
        pytest.approx((0.5, 0.5))


DRAWS = list(range(25, 75, 5))


def test_the_overall_series_is_the_per_draw_functional_of_each_chain(tmp_path):
    series = stratum_series(open_analysis(tmp_path), "psi_edge")
    assert list(series) == ["chain1", "chain2"]
    for chain, expected in zip(series.values(), PSI_EDGE):
        assert chain["iteration"] == DRAWS
        assert chain["value"] == pytest.approx(expected, rel=2e-3)


@pytest.mark.parametrize("name", ["psi_edge", "mean_q", "mean_observed_1", "mean_unobserved_only",
                                  "mean_multi_source_augmented", "mean_multi_source_observed"])
def test_stratum_series_average_back_to_the_overall_series_draw_by_draw(tmp_path, name):
    analysis = open_analysis(tmp_path, epi="infection,district\n1,North\n2,South\n")  # 3 has no value
    overall = stratum_series(analysis, name)
    parts = {(value, s.n): stratum_series(analysis, name, "district", value)
             for value, s in ((s.value, s) for s in stratified_functionals(analysis, "district"))}
    assert [value for value, _ in parts] == ["North", "South", None]
    for chain in ("chain1", "chain2"):
        recombined = [sum(n * part[chain]["value"][i] for (_, n), part in parts.items()) / 3 for i in range(10)]
        assert recombined == pytest.approx(overall[chain]["value"], abs=1e-12)


def test_a_stratum_series_mean_over_draws_is_the_stratums_estimate(tmp_path):
    analysis = open_analysis(tmp_path)
    south = stratified_functionals(analysis, "district")[1]
    values = [v for chain in stratum_series(analysis, "mean_q", "district", "South").values() for v in chain["value"]]
    assert interval(values) == south.mean_q


def test_an_unknown_functional_or_stratum_value_is_an_error(tmp_path):
    analysis = open_analysis(tmp_path)
    with pytest.raises(StrataError, match="mean_entropy.*psi_edge, mean_q"):
        stratum_series(analysis, "mean_entropy")
    with pytest.raises(StrataError, match="no Infection has district 'East'"):
        stratum_series(analysis, "psi_edge", "district", "East")
