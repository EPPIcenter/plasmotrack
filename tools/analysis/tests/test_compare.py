"""Comparison of Analyses (U5): several Analyses of one input opened together, the first primary, and their
convergence tables side by side. Fixture: the golden single-chain store alone as Analysis "single" and the
replica-exchange store alone as Analysis "replica", both fitted to the same input."""

import json
import shutil

import pytest
from conftest import GOLDEN
from test_analysis import write_manifest

from pt_analysis.compare import DifferentInputs, open_comparison


def one_run_analysis(tmp_path, store, label, **sections):
    shutil.copytree(GOLDEN / store, tmp_path / store)
    return write_manifest(tmp_path / f"{label}.json", label=label, runs={"chain1": store}, **sections)


@pytest.fixture
def manifests(tmp_path):
    return [one_run_analysis(tmp_path, "single-chain", "single"),
            one_run_analysis(tmp_path, "replica-exchange", "replica")]


def test_analyses_open_in_order_with_the_first_primary(manifests):
    comparison = open_comparison(manifests)
    assert list(comparison.analyses) == ["single", "replica"]
    assert comparison.primary.label == "single"
    assert comparison.input_sha256 == "cd84e24b6713ef19bc3e27d4a45b9d579634262c808ee3530375e1ac409774ca"


def tamper_input_hash(store_dir):
    manifest = json.loads((store_dir / "manifest.json").read_text())
    manifest["input_sha256"]["input"] = "0" * 64
    (store_dir / "manifest.json").write_text(json.dumps(manifest))


def test_analyses_of_different_inputs_are_refused_naming_each(manifests, tmp_path):
    tamper_input_hash(tmp_path / "replica-exchange")
    with pytest.raises(DifferentInputs) as e:
        open_comparison(manifests)
    message = str(e.value)
    assert "single" in message and "replica" in message
    assert "cd84e24b6713" in message and "000000000000" in message


def test_a_comparison_manifest_lists_the_analyses_relative_to_itself(manifests, tmp_path):
    listing = tmp_path / "lists" / "comparison.json"
    listing.parent.mkdir()
    listing.write_text(json.dumps({"analyses": ["../replica.json", "../single.json"]}))
    assert list(open_comparison([listing]).analyses) == ["replica", "single"]


def test_repeated_labels_get_a_suffix(manifests):
    assert list(open_comparison([manifests[0], manifests[0]]).analyses) == ["single", "single (2)"]


def by_name(rows):
    return {(r["group"], r["name"]): r["analyses"] for r in rows}


def test_convergence_matches_quantities_by_name_with_each_analysis_side_by_side(manifests):
    from pt_analysis.compare import compare_convergence

    rows = compare_convergence(open_comparison(manifests).analyses)
    # The 42 quantities of the two-Chain Analysis of these stores (test_server), each once.
    assert len(rows) == 42
    assert len(by_name(rows)) == 42
    assert all(set(r["analyses"]) == {"single", "replica"} for r in rows)
    log_posterior = by_name(rows)[("log posterior", "log posterior")]
    # One Chain of 10 Draws each: ESS as the Chain's own (test_server), split R-hat undefined.
    assert log_posterior == {"single": {"ess": pytest.approx(10.0), "rhat": None},
                             "replica": {"ess": pytest.approx(7.44641539328075), "rhat": None}}
    assert rows[0]["group"] == "log posterior" and rows[-1]["group"] == "functionals"


def test_without_timings_there_is_no_ess_per_second(manifests):
    from pt_analysis.compare import compare_convergence, sampling_seconds

    comparison = open_comparison(manifests)
    assert [sampling_seconds(a) for a in comparison.analyses.values()] == [None, None]
    rows = compare_convergence(comparison.analyses)
    assert all("ess_per_second" not in cell for r in rows for cell in r["analyses"].values())


def test_timings_give_ess_per_second_of_sampling_summed_over_chains(tmp_path):
    from pt_analysis.compare import compare_convergence, sampling_seconds

    timed = one_run_analysis(tmp_path, "single-chain", "single",
                             timings={"chain1": {"exit": 0, "seed": 1, "sampling_seconds": 4.0, "samples": 50}})
    untimed = one_run_analysis(tmp_path, "replica-exchange", "replica")
    comparison = open_comparison([timed, untimed])
    assert sampling_seconds(comparison.analyses["single"]) == 4.0
    log_posterior = by_name(compare_convergence(comparison.analyses))[("log posterior", "log posterior")]
    assert log_posterior["single"]["ess_per_second"] == pytest.approx(10.0 / 4.0)
    assert "ess_per_second" not in log_posterior["replica"]


def test_timings_missing_a_chain_count_as_absent(tmp_path):
    from pt_analysis.compare import sampling_seconds

    path = one_run_analysis(tmp_path, "single-chain", "single", timings={"other": {"sampling_seconds": 4.0}})
    assert sampling_seconds(open_comparison([path]).primary) is None


def test_serve_opens_every_analysis_named_the_first_primary(manifests, monkeypatch):
    from pt_analysis import cli, server

    served = {}
    monkeypatch.setattr(server, "serve", lambda analysis, **kw: served.update(primary=analysis.label, **kw))
    assert cli.main(["serve", str(manifests[0]), str(manifests[1]), "--no-browser"]) == 0
    assert served["primary"] == "single"
    assert list(served["analyses"]) == ["single", "replica"]


def test_serve_refuses_analyses_of_different_inputs(manifests, tmp_path, capsys):
    from pt_analysis import cli

    tamper_input_hash(tmp_path / "replica-exchange")
    assert cli.main(["serve", str(manifests[0]), str(manifests[1]), "--no-browser"]) == 1
    assert "different inputs" in capsys.readouterr().err
