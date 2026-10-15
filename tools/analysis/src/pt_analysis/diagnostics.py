"""Convergence diagnostics over Chains: multi-chain effective sample size and split R-hat.

A Chain here is one Run's series of a quantity over its Draws, in Sampling-iteration order. Both statistics treat a
series that is constant up to rounding (every value within `ROUNDING` of the others, relative to its magnitude) as
constant and return None for it: a sum of Parent-set weights that is 1 at every Draw varies by about 1e-16, and an ESS
of that noise means nothing.
"""

import math
import statistics

ROUNDING = 1e-12


def constant(chains):
    """Whether every value of every chain lies within ROUNDING (relative, at least absolute) of every other."""
    values = [v for c in chains for v in c]
    return max(values) - min(values) <= ROUNDING * max(1.0, max(abs(v) for v in values))


def multichain_ess(chains):
    """Effective sample size of several chains, counting disagreement between chains as variance (Vehtari et al.
    2021, without rank normalisation). Chains are truncated to the shortest; None below 4 Draws or for a constant
    series. Autocorrelations are taken against the pooled variance and summed with Geyer's initial monotone sequence,
    computing autocovariances (divisor n) only up to the lag where the sum stops."""
    m = len(chains)
    n = min(len(c) for c in chains)
    if n < 4:
        return None
    chains = [list(c[:n]) for c in chains]
    if constant(chains):
        return None
    means = [sum(c) / n for c in chains]
    centred = [[v - mu for v in c] for c, mu in zip(chains, means)]

    def mean_autocov(lag):
        return statistics.mean(sum(c[i] * c[i + lag] for i in range(n - lag)) / n for c in centred)

    within = mean_autocov(0) * n / (n - 1)
    between = n * statistics.variance(means) if m > 1 else 0.0
    pooled = (n - 1) / n * within + between / n
    if pooled <= 0:
        return None

    def rho(lag):
        return 1.0 - (within - mean_autocov(lag)) / pooled

    total, previous = 0.0, math.inf
    for k in range(0, n - 1, 2):
        pair = rho(k) + rho(k + 1)
        if pair <= 0:
            break
        pair = min(pair, previous)  # initial monotone sequence
        total += pair
        previous = pair
    tau = max(-1.0 + 2.0 * total, 1.0 / math.log10(max(m * n, 10)))
    return m * n / tau


def split_rhat(chains):
    """Split R-hat: each chain of length L gives its halves [:L//2] and [L//2:2*(L//2)], kept if at least 2 Draws long;
    every half is truncated to the shortest half's length n. None with fewer than 4 halves (fewer than 2 chains) or
    for a constant series; infinite when every half is constant but the halves differ. Otherwise
    sqrt(((n - 1) / n * W + B / n) / W), W the mean within-half variance and B n times the variance of half means."""
    halves = []
    for c in chains:
        h = len(c) // 2
        if h >= 2:
            halves += [list(c[:h]), list(c[h:2 * h])]
    if len(halves) < 4:
        return None
    n = min(len(h) for h in halves)
    halves = [h[:n] for h in halves]
    if constant(halves):
        return None
    means = [sum(h) / n for h in halves]
    within = statistics.mean(statistics.variance(h) for h in halves)
    between = n * statistics.variance(means)
    if within <= 0:
        # Constant within every half but not between them: stuck at different values.
        return None if between <= 0 else math.inf
    return math.sqrt(((n - 1) / n * within + between / n) / within)
