"""An Analysis manifest names the input, the IDPs, the labelled Runs and optional truth, epi and score files."""

import json
import os

import pytest
from conftest import GOLDEN

from pt_analysis.analysis import Analysis, AnalysisError
from pt_analysis.store import Store


def write_manifest(path, **sections):
    """Write a manifest whose paths are relative to it, so that resolution against the manifest is exercised."""
    rel = lambda p: os.path.relpath(p, path.parent)  # noqa: E731
    manifest = {
        "input": rel(GOLDEN.parents[2] / "examples" / "example_input.json"),
        "symptomatic_idp": rel(GOLDEN / "symptomatic_idp.txt"),
        "asymptomatic_idp": rel(GOLDEN / "asymptomatic_idp.txt"),
        "runs": {"single": rel(GOLDEN / "single-chain"), "replica": rel(GOLDEN / "replica-exchange")},
        **sections,
    }
    path.write_text(json.dumps(manifest))
    return path


def test_runs_are_labelled_stores_in_manifest_order(tmp_path):
    analysis = Analysis.open(write_manifest(tmp_path / "analysis.json"))
    assert list(analysis.runs) == ["single", "replica"]
    assert all(isinstance(s, Store) for s in analysis.runs.values())
    assert analysis.runs["single"].path == GOLDEN / "single-chain"
    assert analysis.runs["replica"].draws() == 10


def test_paths_resolve_against_the_manifest(tmp_path):
    analysis = Analysis.open(write_manifest(tmp_path / "analysis.json"))
    assert analysis.input == GOLDEN.parents[2] / "examples" / "example_input.json"
    assert analysis.symptomatic_idp == GOLDEN / "symptomatic_idp.txt"
    assert analysis.asymptomatic_idp == GOLDEN / "asymptomatic_idp.txt"
    assert analysis.path == tmp_path / "analysis.json"


def test_optional_sections_default_to_none(tmp_path):
    analysis = Analysis.open(write_manifest(tmp_path / "analysis.json"))
    assert analysis.truth is None
    assert analysis.epi is None
    assert analysis.pairwise_scores is None
    assert analysis.null_model is None
    assert analysis.timings is None


def test_optional_sections_are_read_and_resolved(tmp_path):
    (tmp_path / "timings.json").write_text(json.dumps({"single": {"exit": 0, "seed": 42}}))
    analysis = Analysis.open(write_manifest(
        tmp_path / "analysis.json", truth="sim_0_s0.json.gz", epi="epi.csv", pairwise_scores="ibd.csv",
        null_model="null/analysis.json", timings="timings.json"))
    assert analysis.truth == tmp_path / "sim_0_s0.json.gz"
    assert analysis.epi == tmp_path / "epi.csv"
    assert analysis.pairwise_scores == tmp_path / "ibd.csv"
    assert analysis.null_model == tmp_path / "null" / "analysis.json"
    assert analysis.timings == {"single": {"exit": 0, "seed": 42}}


def test_inline_timings_are_kept(tmp_path):
    analysis = Analysis.open(write_manifest(tmp_path / "analysis.json", timings={"single": {"exit": 0}}))
    assert analysis.timings == {"single": {"exit": 0}}


def test_a_missing_run_directory_fails_naming_the_path(tmp_path):
    path = write_manifest(tmp_path / "analysis.json", runs={"single": "no-such-run"})
    with pytest.raises(AnalysisError, match=str(tmp_path / "no-such-run")):
        Analysis.open(path)


def test_a_manifest_without_runs_fails(tmp_path):
    path = tmp_path / "analysis.json"
    path.write_text(json.dumps({"input": "x.json", "symptomatic_idp": "s", "asymptomatic_idp": "a"}))
    with pytest.raises(AnalysisError, match="runs"):
        Analysis.open(path)
