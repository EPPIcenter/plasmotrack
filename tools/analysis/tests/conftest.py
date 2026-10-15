from pathlib import Path

import pytest

GOLDEN = Path(__file__).resolve().parents[3] / "test" / "data" / "golden"


@pytest.fixture(params=["single-chain", "replica-exchange"])
def golden_store(request):
    return GOLDEN / request.param
