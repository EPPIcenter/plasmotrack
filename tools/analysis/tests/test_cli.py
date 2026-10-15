"""`pt-analysis info analysis.json` lists each Run's closed Windows per Phase and its Draws."""

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
