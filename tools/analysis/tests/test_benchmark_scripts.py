"""The benchmark scripts read stores through the package and give the numbers they gave before the port."""

import subprocess
import sys
from pathlib import Path

import pytest
from conftest import GOLDEN

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


def test_chain_summary_diagnostics_table_is_unchanged(tmp_path):
    """The ESS and R-hat table on the two golden stores as two Chains, as the script printed it before the port.
    The pairs and agreement sections still read the old gzip CSV layout; the agreement section fails after the
    table, so the exit status is not checked."""
    (tmp_path / "chain1").symlink_to(GOLDEN / "single-chain")
    (tmp_path / "chain2").symlink_to(GOLDEN / "replica-exchange")
    result = subprocess.run([sys.executable, str(BENCHMARK / "chain_summary.py"), str(tmp_path)],
                            capture_output=True, text=True)
    assert result.stdout.startswith("""2 chains: chain1, chain2

group                     series  ESS median   ESS min  R-hat>1.1  R-hat max
log posterior                  1         3.6       3.6       100%       1.26
global scalars                 3        26.0      19.5         0%      0.938
false positive rate            3        26.0      21.3         0%      0.999
false negative rate            3        26.0       3.2        33%       1.11
infection duration             3        26.0      24.1         0%       1.01
allele frequency              20        10.6       1.4        65%       2.43
P(no sampled parent)           2        19.3      12.5        50%        1.1
top parent edge                2        19.3      12.5        50%        1.1
""")


def test_mixing_benchmark_sampling_time_reads_the_store():
    assert mb.Store(GOLDEN / "replica-exchange").sampling_time() == {"sampling_seconds": 0.260588, "samples": 50}
