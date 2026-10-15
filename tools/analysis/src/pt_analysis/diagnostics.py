"""Convergence diagnostics over Chains: multi-chain effective sample size and split R-hat.

A Chain here is one Run's series of a quantity over its Draws, in Sampling-iteration order. Both statistics treat a
series that is constant up to rounding (every value within `ROUNDING` of the others, relative to its magnitude) as
constant and return None for it: a sum of Parent-set weights that is 1 at every Draw varies by about 1e-16, and an ESS
of that noise means nothing.
"""

import math
import statistics

import numpy as np

ROUNDING = 1e-12


def constant(chains):
    """Whether every value of every chain lies within ROUNDING (relative, at least absolute) of every other."""
    values = [v for c in chains for v in c]
    return max(values) - min(values) <= ROUNDING * max(1.0, max(abs(v) for v in values))


def multichain_ess(chains):
    """Effective sample size of several chains, counting disagreement between chains as variance (Vehtari et al.
    2021, without rank normalisation). Chains are truncated to the shortest; None below 4 Draws or for a constant
    series. Autocorrelations are taken against the pooled variance and summed with Geyer's initial monotone sequence.
    Autocovariances (divisor n) come for every lag at once from an FFT of each chain, averaged over chains."""
    m = len(chains)
    n = min(len(c) for c in chains)
    if n < 4:
        return None
    chains = [list(c[:n]) for c in chains]
    if constant(chains):
        return None
    x = np.asarray(chains, dtype=float)
    means = x.mean(axis=1)
    centred = x - means[:, None]
    size = 1 << (2 * n - 1).bit_length()  # zero-padded past 2n - 1: no wrap-around between lags
    spectrum = np.fft.rfft(centred, size, axis=1)
    autocov = (np.fft.irfft(spectrum * spectrum.conj(), size, axis=1)[:, :n] / n).mean(axis=0)

    within = float(autocov[0]) * n / (n - 1)
    between = n * statistics.variance(means.tolist()) if m > 1 else 0.0
    pooled = (n - 1) / n * within + between / n
    if pooled <= 0:
        return None
    rho = 1.0 - (within - autocov) / pooled

    total, previous = 0.0, math.inf
    for k in range(0, n - 1, 2):
        pair = float(rho[k] + rho[k + 1])
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
