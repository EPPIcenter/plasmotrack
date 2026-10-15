import json
import os
import shutil
from pathlib import Path

import pytest

from pt_analysis.analysis import Analysis
from pt_analysis.store import Store

GOLDEN = Path(__file__).resolve().parents[3] / "test" / "data" / "golden"

# The two golden stores (same input) as Chains chain1 and chain2 of one Analysis, by their copies' directory names.
TWO_CHAINS = {"chain1": "single-chain", "chain2": "replica-exchange"}


@pytest.fixture(params=["single-chain", "replica-exchange"])
def golden_store(request):
    return GOLDEN / request.param


@pytest.fixture
def store(golden_store, tmp_path):
    """A copy of a golden store, so that `derived/` is written outside the repository's test data."""
    path = tmp_path / golden_store.name
    shutil.copytree(golden_store, path)
    return Store(path)


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


def copy_golden(directory):
    """Copy both golden stores into `directory` (once), so that `derived/` is written outside the repository's test
    data. Returns the two-Chain `runs` of a manifest in `directory`."""
    for name in TWO_CHAINS.values():
        if not (directory / name).exists():
            shutil.copytree(GOLDEN / name, directory / name)
    return dict(TWO_CHAINS)


def open_golden(directory, name="analysis.json", runs=None, **sections):
    """The golden stores copied into `directory` and opened as one Analysis whose manifest is `directory/name`:
    `runs` (default both stores as chain1 and chain2) and any further manifest `sections`."""
    two = copy_golden(directory)
    return Analysis.open(write_manifest(directory / name, runs=two if runs is None else runs, **sections))


@pytest.fixture
def golden_analysis(tmp_path):
    """The two golden stores as Chains chain1 and chain2 of one Analysis labelled `golden`, copied into `tmp_path`."""
    return open_golden(tmp_path, label="golden")
