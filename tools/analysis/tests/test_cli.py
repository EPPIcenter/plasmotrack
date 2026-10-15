"""The `pt-analysis` commands: `info`, `derive`, and the summaries `edges`, `map` and `credible`."""

import shutil

from conftest import GOLDEN
from test_analysis import write_manifest

from pt_analysis.cli import main


def test_info_lists_runs_windows_per_phase_and_draws(tmp_path, capsys):
    manifest = write_manifest(tmp_path / "analysis.json")
    assert main(["info", str(manifest)]) == 0
    out = capsys.readouterr().out
    assert out == f"""Analysis: {manifest}
input: {GOLDEN.parents[2] / 'examples' / 'example_input.json'}
Runs: 2

single: {GOLDEN / 'single-chain'}
  burn_in: 1 Window (0), iterations 0-20
  sampling: 2 Windows (1-2), iterations 21-70
  Draws: 10

replica: {GOLDEN / 'replica-exchange'}
  burn_in: 1 Window (0), iterations 0-20
  sampling: 2 Windows (1-2), iterations 21-70
  Draws: 10
"""


def test_info_on_a_missing_run_directory_fails_naming_the_path(tmp_path, capsys):
    manifest = write_manifest(tmp_path / "analysis.json", runs={"single": "gone"})
    assert main(["info", str(manifest)]) == 1
    assert str(tmp_path / "gone") in capsys.readouterr().err


def test_derive_a_run_then_again(store, capsys):
    assert main(["derive", str(store.path)]) == 0
    assert capsys.readouterr().out == f"{store.path}: derived from burn_in Windows [0], sampling Windows [1, 2]\n"
    assert main(["derive", str(store.path)]) == 0
    assert capsys.readouterr().out == f"{store.path}: up to date\n"


def test_derive_every_run_of_an_analysis(store, tmp_path, capsys):
    shutil.copytree(GOLDEN / "replica-exchange", tmp_path / "chain2")
    manifest = write_manifest(tmp_path / "analysis.json", runs={"chain1": str(store.path), "chain2": "chain2"})
    assert main(["derive", str(manifest)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert [line.split(": ")[0] for line in out] == [str(store.path), str(tmp_path / "chain2"),
                                                     str(tmp_path / "analysis.derived")]
    assert out[2] == f"{tmp_path / 'analysis.derived'}: derived convergence from Chains chain1, chain2"
    assert (tmp_path / "chain2" / "derived" / "manifest.json").is_file()
    assert (tmp_path / "analysis.derived" / "convergence.parquet").is_file()
    assert main(["derive", str(manifest)]) == 0
    assert capsys.readouterr().out.splitlines()[2] == f"{tmp_path / 'analysis.derived'}: up to date"


def test_derive_on_a_directory_that_is_not_a_store_fails(tmp_path, capsys):
    assert main(["derive", str(tmp_path)]) == 1
    assert "not a Run store" in capsys.readouterr().err


def test_edges_prints_the_top_edges(store, capsys):
    assert main(["edges", str(store.path), "--top", "1"]) == 0
    weight = {"single-chain": "0.8914", "replica-exchange": "0.7757"}[store.path.name]
    assert capsys.readouterr().out == f"parent  child  weight\n1       3      {weight}\n"


def test_map_prints_the_map_graph_psi_map_and_the_size_range(store, capsys):
    assert main(["map", str(store.path)]) == 0
    out = capsys.readouterr().out.splitlines()
    assert out[0].split() == ["infection", "set", "weight", "members"]
    assert [line.split()[:2] for line in out[1:4]] == [["1", "0"], ["2", "0"], ["3", "4"]]
    assert out[1].split()[3:] == ["{1_copy}"] and out[3].split()[3:] == ["{1}"]
    assert out[4] == "psi_MAP: 0.3333"
    lo = {"single-chain": "0.3333", "replica-exchange": "0"}[store.path.name]
    assert out[5] == f"parent-set size range at 0.8: [{lo}, 0.3333]"


def test_credible_prints_an_infections_collection(store, capsys):
    assert main(["credible", str(store.path), "3", "--level", "0.8"]) == 0
    out = capsys.readouterr().out.splitlines()
    size = {"single-chain": 1, "replica-exchange": 2}[store.path.name]
    assert out[0] == f"Infection 3: {size} {'set' if size == 1 else 'sets'} with mass at least 0.8"
    assert out[1].split() == ["rank", "set", "weight", "cumulative", "class", "members"]
    assert out[2].split()[:2] == ["1", "4"] and out[2].split()[4:] == ["observed_1", "{1}"]
    assert len(out) == 2 + size
    if size == 2:
        assert out[3].split()[:2] == ["2", "0"] and out[3].split()[4:] == ["unobserved_only", "{3_copy}"]


def test_credible_for_an_unknown_infection_fails(store, capsys):
    assert main(["credible", str(store.path), "99"]) == 1
    assert "99" in capsys.readouterr().err
