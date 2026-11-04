"""The Analysis-level `convergence` table: per quantity, each Chain's ESS and the pooled multi-chain ESS and split
R-hat, written under `<manifest stem>.derived/` beside the Analysis manifest."""

import json

import duckdb
import pytest

from pt_analysis import __version__
from pt_analysis.derive import derive_analysis

GROUPS = ["log posterior", "global scalars", "false positive rate", "false negative rate", "infection duration",
          "allele frequency", "edge", "unobserved-source inclusion", "functionals"]


def convergence(analysis):
    derive_analysis(analysis)
    path = analysis.path.parent / "analysis.derived" / "convergence.parquet"
    relation = duckdb.sql(f"SELECT * FROM '{path}'")
    return [dict(zip(relation.columns, row)) for row in relation.fetchall()]


def test_every_quantity_has_a_row_per_chain_and_a_pooled_row_in_group_order(golden_analysis):
    rows = convergence(golden_analysis)
    assert list(rows[0]) == ["group", "name", "chain", "ess", "rhat"]
    groups = list(dict.fromkeys(r["group"] for r in rows))
    assert groups == GROUPS
    names = {g: list(dict.fromkeys(r["name"] for r in rows if r["group"] == g)) for g in GROUPS}
    assert names["log posterior"] == ["log posterior"]
    assert names["global scalars"] == ["mean_coi", "mean_strains_tx", "parent_set_size_prob"]
    assert names["false positive rate"] == names["infection duration"] == ["1", "2", "3"]
    assert len(names["allele frequency"]) == 20 and names["allele frequency"][0] == "AS11[0]"
    assert names["edge"] == ["2<-1", "3<-1", "3<-2"]
    assert names["unobserved-source inclusion"] == ["1", "2", "3"]
    assert names["functionals"] == ["mean_entropy", "mean_q", "psi_edge"]
    for group, group_names in names.items():
        for name in group_names:
            assert [r["chain"] for r in rows if (r["group"], r["name"]) == (group, name)] == ["chain1", "chain2", None]


def by_name(rows, group, name):
    return {r["chain"]: (r["ess"], r["rhat"]) for r in rows if (r["group"], r["name"]) == (group, name)}


# Printed by mixing_benchmark's multichain_ess and split_rhat on its run_traces series, before they moved into the
# package: (ESS, R-hat) per Chain (None: pooled).
BASELINE = {
    ("log posterior", "log posterior"): {"chain1": 10.0, "chain2": 7.44641539328075,
                                         None: (3.58385705875339, 1.2643097719105674)},
    ("allele frequency", "AS1[3]"): {"chain1": 6.493931772122144, "chain2": 10.0,
                                     None: (17.17541738618183, 1.0208460468250458)},
    ("false negative rate", "3"): {"chain1": 10.0, "chain2": 10.0, None: (26.020599913279625, 0.9607583510545856)},
    ("edge", "3<-1"): {"chain1": 10.0, "chain2": 6.702635780124615, None: (12.484076791343613, 1.1044701085052935)},
}


def test_ess_and_rhat_reproduce_the_benchmark_series(golden_analysis):
    rows = convergence(golden_analysis)
    for (group, name), expected in BASELINE.items():
        got = by_name(rows, group, name)
        assert got["chain1"] == pytest.approx((expected["chain1"], None), rel=1e-9), (group, name)
        assert got["chain2"] == pytest.approx((expected["chain2"], None), rel=1e-9), (group, name)
        assert got[None] == pytest.approx(expected[None], rel=1e-9), (group, name)


def test_quantities_constant_in_every_chain_have_no_ess(golden_analysis):
    rows = convergence(golden_analysis)
    # Infection 1 has one set, {1_copy}; every set of Infection 2 contains 2_copy, so its q is 1 up to rounding.
    for infection in ("1", "2"):
        assert by_name(rows, "unobserved-source inclusion", infection) == {
            "chain1": (None, None), "chain2": (None, None), None: (None, None)}
    q_3 = by_name(rows, "unobserved-source inclusion", "3")
    assert q_3[None][0] is not None and q_3[None][1] is not None


def test_the_convergence_manifest_records_each_runs_windows_and_a_second_derive_writes_nothing(golden_analysis):
    assert derive_analysis(golden_analysis) is True
    out = golden_analysis.path.parent / "analysis.derived"
    assert json.loads((out / "manifest.json").read_text()) == {
        "runs": {label: {"path": str(store.path), "windows": {"burn_in": [0], "sampling": [1, 2]}}
                 for label, store in golden_analysis.runs.items()},
        "package_version": __version__}
    before = {p.name: p.stat().st_mtime_ns for p in out.iterdir()}
    assert derive_analysis(golden_analysis) is False
    assert {p.name: p.stat().st_mtime_ns for p in out.iterdir()} == before


def test_a_chain_that_gained_a_window_rederives_the_convergence_table(golden_analysis):
    store = golden_analysis.runs["chain2"]
    marker = store.path / "_closed" / "window-000002"
    aside = store.path / "window-000002.aside"
    marker.rename(aside)
    rows = convergence(golden_analysis)
    assert by_name(rows, "log posterior", "log posterior")["chain2"][0] is not None
    aside.rename(marker)
    assert by_name(convergence(golden_analysis), "log posterior", "log posterior")[None] == pytest.approx(
        BASELINE[("log posterior", "log posterior")][None], rel=1e-9)


def test_a_chain_still_in_burn_in_has_no_ess_and_is_left_out_of_the_pooled_row(golden_analysis):
    for w in (1, 2):
        (golden_analysis.runs["chain2"].path / "_closed" / f"window-{w:06d}").unlink()
    rows = convergence(golden_analysis)
    log_posterior = by_name(rows, "log posterior", "log posterior")
    assert log_posterior["chain2"] == (None, None)
    assert log_posterior[None] == pytest.approx((10.0, None))
    assert by_name(rows, "edge", "3<-1")["chain2"] == (None, None)
