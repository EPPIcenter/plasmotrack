"""The benchmark scripts read stores through the package and give the numbers they gave before the port."""

import json
import shutil
import subprocess
import sys

import pytest
from conftest import GOLDEN

from pt_analysis.cli import main

BENCHMARK = GOLDEN.parents[2] / "tools" / "benchmark"
sys.path.insert(0, str(BENCHMARK))
import mixing_benchmark as mb  # noqa: E402

# From the exploration notes, computed independently with R arrow.
EDGE_1_TO_3 = {
    "single-chain": [0.9638, 0.8808, 0.9799, 0.9337, 0.8457, 0.8624, 0.9414, 0.8481, 0.8313, 0.8267],
    "replica-exchange": [0.0031, 0.0039, 0.9757, 0.9767, 0.9533, 0.9499, 0.9657, 0.9705, 0.9718, 0.9864],
}
UNOBSERVED_ONLY_3 = {"single-chain": 0.1084, "replica-exchange": 0.2241}


@pytest.mark.parametrize("name", ["single-chain", "replica-exchange"])
def test_run_traces_edge_and_no_parent_series(name):
    traces = mb.run_traces(GOLDEN / name)
    assert [round(v, 4) for v in traces["_edges"]["3<-1"]] == EDGE_1_TO_3[name]
    assert sum(traces["P(no sampled parent)"]["3"]) / 10 == pytest.approx(UNOBSERVED_ONLY_3[name], abs=1e-4)
    assert traces["P(no sampled parent)"]["1"] == [1.0] * 10


def test_run_traces_groups_and_series_counts():
    traces = mb.run_traces(GOLDEN / "single-chain")
    assert list(traces) == ["log posterior", "global scalars", "false positive rate", "false negative rate",
                            "infection duration", "allele frequency", "P(no sampled parent)", "_edges"]
    assert len(traces["log posterior"]["log posterior"]) == 10
    assert sorted(traces["global scalars"]) == ["mean_coi", "mean_strains_tx", "parent_set_size_prob"]
    assert sorted(traces["infection duration"]) == ["1", "2", "3"]
    assert len(traces["allele frequency"]) == 20
    assert sorted(traces["_edges"]) == ["2<-1", "3<-1", "3<-2"]
    assert round(traces["log posterior"]["log posterior"][0], 4) == round(
        sum(mb.Store(GOLDEN / "single-chain").rows("run_scalars")[0][k] for k in ("prior", "likelihood")), 4)


def chain_summary(*args):
    return subprocess.run([sys.executable, str(BENCHMARK / "chain_summary.py"), *map(str, args)],
                          capture_output=True, text=True)


def test_chain_summary_prints_the_package_summary_of_run_dir_chains(tmp_path, capsys):
    """RUN_DIR/chain1 and chain2 become the Chains of a generated manifest, RUN_DIR/chain_summary.json, whose
    summary the script prints: the same as `pt-analysis summary` on it, starting with the diagnostics table."""
    shutil.copytree(GOLDEN / "single-chain", tmp_path / "chain1")
    shutil.copytree(GOLDEN / "replica-exchange", tmp_path / "chain2")
    (tmp_path / "chain2-notes.txt").write_text("not a Chain")
    result = chain_summary(tmp_path, "--pairs", "1:2", "--show", "2")
    assert result.returncode == 0, result.stderr
    manifest = json.loads((tmp_path / "chain_summary.json").read_text())
    assert manifest["runs"] == {"chain1": "chain1", "chain2": "chain2"}
    assert (manifest["input"], manifest["symptomatic_idp"], manifest["asymptomatic_idp"]) == (
        "<input>", "<golden>/symptomatic_idp.txt", "<golden>/asymptomatic_idp.txt")
    assert main(["summary", str(tmp_path / "chain_summary.json"), "--pairs", "1:2", "--show", "2"]) == 0
    assert result.stdout == capsys.readouterr().out
    assert result.stdout.startswith("2 Chains: chain1, chain2\n\ngroup")
    assert "log posterior                        1         3.6       3.6       100%       1.26\n" in result.stdout
    assert "chain1: 1 before 2: 10/10 (changes 0)\n" in result.stdout
    assert chain_summary(tmp_path, "--pairs", "1:2", "--show", "2").stdout == result.stdout  # beside its own output


def test_chain_summary_selects_chains_by_number(tmp_path):
    for k, name in ((1, "single-chain"), (2, "replica-exchange"), (10, "single-chain")):
        shutil.copytree(GOLDEN / name, tmp_path / f"chain{k}")
    result = chain_summary(tmp_path, "--chains", "10,2")
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("2 Chains: chain2, chain10\n")


def test_mixing_benchmark_sampling_time_reads_the_store():
    assert mb.Store(GOLDEN / "replica-exchange").sampling_time() == {"sampling_seconds": 0.260588, "samples": 50}
