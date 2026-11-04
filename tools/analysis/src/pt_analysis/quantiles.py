"""Quantiles of pooled per-Draw values, as the manuscript's R pipeline takes them: type 7 (R's default), the p
quantile lying (n - 1) p of the way along the sorted values, interpolated linearly. SQL takes the same quantile with
DuckDB's `quantile_cont`."""

import math
from dataclasses import dataclass

PROBABILITIES = (0.025, 0.975)


@dataclass(frozen=True)
class Interval:
    """A quantity's mean over Draws and its 95% MCMC interval."""

    mean: float
    lower: float
    upper: float


def _sorted_quantile(ordered, p):
    h = (len(ordered) - 1) * p
    below = math.floor(h)
    above = min(below + 1, len(ordered) - 1)
    return ordered[below] + (h - below) * (ordered[above] - ordered[below])


def quantile(values, p):
    """The type-7 p quantile of `values` (in any order)."""
    return _sorted_quantile(sorted(values), p)


def quantiles(values, ps):
    """The type-7 quantiles of `values` at each of `ps`, sorting once."""
    ordered = sorted(values)
    return [_sorted_quantile(ordered, p) for p in ps]


def interval(values):
    """The mean of `values` and their 2.5% and 97.5% type-7 quantiles; None without values."""
    if not values:
        return None
    return Interval(sum(values) / len(values), *quantiles(values, PROBABILITIES))
