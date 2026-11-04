"""The thesis-rerun loader: one Analysis manifest per (model, simulation) from `<root>/<model>/<sim>/chainN`, with
the simulation's inputs and truth from the inputs directory and its factors from the path. Fixtures are built at
test time (`.gitignore` ignores `*.gz`): golden stores copied as chain1 and chain2, a fake input and truth file."""

import gzip
import json
import shutil

import pytest
from conftest import GOLDEN

from pt_analysis import rerun
from pt_analysis.analysis import Analysis
from pt_analysis.cli import main

SCENARIO = ("total_nodes_200_founder_prop_0.05_super_infected_prop_0_r0_1.5_mean_moi_2_asymptomatic_prop_0.5")
SIM = f"{SCENARIO}/simple_tx_network_med_div/sim_0_s7"
OTHER = f"{SCENARIO}/tx_network_high_div/sim_0_s3"


def gz(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(json.dumps(value).encode()))


@pytest.fixture
def layout(tmp_path):
    root, inputs = tmp_path / "results", tmp_path / "inputs"
    example = json.loads((GOLDEN.parents[2] / "examples" / "example_input.json").read_text())
    truth = {"network": [{"src": "1", "dest": "3"}], "nodes": [{"id": i} for i in ("1", "2", "3", "4")]}
    for sim in (SIM, OTHER):
        gz(inputs / f"{sim}_to_run.json.gz", example)
        gz(inputs / f"{sim}_to_run_partially_obs.json.gz", example)
        gz(inputs / f"{sim}.json.gz", truth)
    (tmp_path / "symptomatic_idp.csv").write_text("0.5\n0.5\n")
    (tmp_path / "asymptomatic_idp.csv").write_text("0.5\n0.5\n")
    for model in ("full", "partial", "null"):
        for chain, golden in (("chain1", "single-chain"), ("chain2", "replica-exchange")):
            shutil.copytree(GOLDEN / golden, root / model / SIM / chain)
    shutil.copytree(GOLDEN / "single-chain", root / "full" / OTHER / "chain1")  # no partial or null fit yet
    (root / "full" / OTHER / "chain2").mkdir()  # a Chain that has not started: no store yet
    return root, inputs, tmp_path


def test_parse_simulation_reads_the_factors_from_the_path():
    assert rerun.parse_simulation(SIM) == {
        "run_key": "simplified__fp0.05__sp0__r01.5__coi2__med__sim0__s7", "simulation_model": "simplified",
        "total_nodes": 200.0, "founder_probability": 0.05, "superinfection_probability": 0.0,
        "secondary_case_rate": 1.5, "founder_coi_rate": 2.0, "asymptomatic_probability": 0.5, "panel_name": "med",
        "panel_alleles_per_locus": 10.0, "simulation_index": 0, "seed": 7}
    assert rerun.parse_simulation(OTHER)["simulation_model"] == "biological"
    with pytest.raises(rerun.RerunError):
        rerun.parse_simulation("somewhere/else/sim_0_s1")


def test_one_manifest_per_model_and_simulation(layout):
    root, inputs, base = layout
    written = rerun.manifests_from_rerun(root, inputs, base / "manifests")
    assert sorted(p.relative_to(base / "manifests").as_posix() for p in written) == sorted(
        [f"full/{SIM}.json", f"partial/{SIM}.json", f"null/{SIM}.json", f"full/{OTHER}.json"])


def test_manifests_open_with_their_chains_inputs_truth_and_factors(layout):
    root, inputs, base = layout
    rerun.manifests_from_rerun(root, inputs, base / "manifests")
    full = Analysis.open(base / "manifests" / "full" / f"{SIM}.json")
    assert list(full.runs) == ["chain1", "chain2"]
    assert full.runs["chain2"].path == (root / "full" / SIM / "chain2").resolve()
    assert full.input == (inputs / f"{SIM}_to_run.json.gz").resolve()
    assert full.truth == (inputs / f"{SIM}.json.gz").resolve()
    assert full.symptomatic_idp == (base / "symptomatic_idp.csv").resolve()
    assert full.factors["missingness_condition"] == "full" and full.factors["fit_type"] == "posterior"
    assert full.factors["run_key"] == "simplified__fp0.05__sp0__r01.5__coi2__med__sim0__s7"
    assert full.label == "simplified__fp0.05__sp0__r01.5__coi2__med__sim0__s7 full"


def test_the_null_fit_is_the_full_fits_null_model(layout):
    root, inputs, base = layout
    rerun.manifests_from_rerun(root, inputs, base / "manifests")
    full = Analysis.open(base / "manifests" / "full" / f"{SIM}.json")
    assert full.null_model == (base / "manifests" / "null" / f"{SIM}.json").resolve()
    null = Analysis.open(full.null_model)
    assert null.factors["fit_type"] == "null_prior" and null.factors["missingness_condition"] == "full"
    assert null.null_model is None
    partial = Analysis.open(base / "manifests" / "partial" / f"{SIM}.json")
    assert partial.input == (inputs / f"{SIM}_to_run_partially_obs.json.gz").resolve()
    assert partial.factors["missingness_condition"] == "missing_50" and partial.factors["fit_type"] == "posterior"
    assert partial.null_model is None


def test_a_chain_without_a_store_is_left_out_and_so_is_a_missing_null_fit(layout):
    root, inputs, base = layout
    rerun.manifests_from_rerun(root, inputs, base / "manifests")
    other = Analysis.open(base / "manifests" / "full" / f"{OTHER}.json")
    assert list(other.runs) == ["chain1"]
    assert other.null_model is None


def test_the_command_writes_the_manifests(layout, capsys):
    root, inputs, base = layout
    assert main(["manifests-from-rerun", str(root), "--inputs", str(inputs), "--out", str(base / "m")]) == 0
    assert "4 manifests" in capsys.readouterr().out
    assert (base / "m" / "null" / f"{SIM}.json").is_file()
