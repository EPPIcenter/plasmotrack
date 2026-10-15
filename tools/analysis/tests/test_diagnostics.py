"""Multi-chain ESS and split R-hat, the one implementation the package and the benchmark scripts share."""

import math

import pytest
from conftest import GOLDEN

from pt_analysis.diagnostics import multichain_ess, split_rhat
from pt_analysis.store import Store


def log_posterior(name):
    return [r["prior"] + r["likelihood"] for r in Store(GOLDEN / name).rows("run_scalars")]


# Printed by mixing_benchmark.multichain_ess and split_rhat before they moved into the package, on the golden stores'
# log posterior: single-chain then replica-exchange as two Chains, and each alone.
def test_ess_and_rhat_of_the_golden_log_posterior_match_the_benchmark_before_the_move():
    chains = [log_posterior("single-chain"), log_posterior("replica-exchange")]
    assert multichain_ess(chains) == pytest.approx(3.58385705875339, rel=1e-12)
    assert split_rhat(chains) == pytest.approx(1.2643097719105674, rel=1e-12)
    assert multichain_ess(chains[:1]) == pytest.approx(10.0, rel=1e-12)
    assert multichain_ess(chains[1:]) == pytest.approx(7.44641539328075, rel=1e-12)


def test_a_single_chain_has_no_split_rhat():
    assert split_rhat([log_posterior("single-chain")]) is None


def test_fewer_than_four_draws_have_no_ess():
    assert multichain_ess([[1.0, 2.0, 3.0], [1.0, 3.0, 2.0]]) is None


def test_a_constant_series_has_no_ess_or_rhat():
    assert multichain_ess([[1.0] * 10, [1.0] * 10]) is None
    assert split_rhat([[1.0] * 10, [1.0] * 10]) is None


def test_a_series_constant_up_to_rounding_has_no_ess_or_rhat():
    """Unobserved-source inclusion of an Infection whose every set contains its Latent parent: 1 up to rounding."""
    chains = [[1.0, 1.0 + 2.2e-16, 1.0 - 1.1e-16, 1.0] * 3, [1.0 - 4.4e-16, 1.0, 1.0 + 2.2e-16, 1.0] * 3]
    assert multichain_ess(chains) is None
    assert split_rhat(chains) is None


def test_chains_stuck_at_different_values_have_infinite_rhat():
    assert split_rhat([[0.0] * 10, [1.0] * 10]) == math.inf


def test_split_rhat_truncates_every_half_to_the_shortest_half():
    """Chains of 4 and 6 Draws split into halves of 2 and 3; every half keeps its first 2 Draws, so the second chain's
    100s are dropped. Every half is then [0, 1]: half means 0.5, so B = 0; W = var([0, 1]) = 0.5; with n = 2,
    R-hat = sqrt(((n - 1) / n * W + B / n) / W) = sqrt(0.5)."""
    assert split_rhat([[0, 1, 0, 1], [0, 1, 100, 0, 1, 100]]) == pytest.approx(math.sqrt(0.5))
