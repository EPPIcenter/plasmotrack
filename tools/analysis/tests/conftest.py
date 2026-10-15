import shutil
from pathlib import Path

import pytest

from pt_analysis.store import Store

GOLDEN = Path(__file__).resolve().parents[3] / "test" / "data" / "golden"


@pytest.fixture(params=["single-chain", "replica-exchange"])
def golden_store(request):
    return GOLDEN / request.param


@pytest.fixture
def store(golden_store, tmp_path):
    """A copy of a golden store, so that `derived/` is written outside the repository's test data."""
    path = tmp_path / golden_store.name
    shutil.copytree(golden_store, path)
    return Store(path)
