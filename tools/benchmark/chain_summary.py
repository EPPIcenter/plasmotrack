#!/usr/bin/env python3
"""Convergence summary of independent chains written to RUN_DIR/chain*/ by the model CLI.

    chain_summary.py RUN_DIR [--chains 1,2,3] [--pairs A:B ...]

Reports, per group of quantities, the multi-chain effective sample size and split R-hat (as mixing_benchmark.py
analyze does, without timings); for each pair A:B, how often each chain had A infected before B and how often that
order changed between logged draws; and, per infection, its most probable parent pooled over chains and the largest
difference between chains in any of its parent probabilities. Standard library only.
"""

import argparse
import gzip
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mixing_benchmark as mb  # noqa: E402


def multichain_ess(chains):
    """mixing_benchmark.multichain_ess, computing autocovariances only up to the lag where Geyer's sum stops."""
    m = len(chains)
    n = min(len(c) for c in chains)
    if n < 4:
        return None
    chains = [c[:n] for c in chains]
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

    total, previous = 0.0, float("inf")
    for k in range(0, n - 1, 2):
        pair = rho(k) + rho(k + 1)
        if pair <= 0:
            break
        pair = min(pair, previous)
        total += pair
        previous = pair
    tau = max(-1.0 + 2.0 * total, 1.0 / math.log10(max(m * n, 10)))
    return m * n / tau


def diagnostics(traces):
    print(f"{'group':<24}{'series':>8}{'ESS median':>12}{'ESS min':>10}{'R-hat>1.1':>11}{'R-hat max':>11}")
    best = {}
    for groups in traces:
        for key, series in groups.get("_edges", {}).items():
            best.setdefault(key.split("<-")[0], defaultdict(float))[key] += sum(series) / len(series)
    top_edges = {max(scores, key=scores.get) for scores in best.values()}
    for groups in traces:
        groups["top parent edge"] = {k: v for k, v in groups.pop("_edges", {}).items() if k in top_edges}
    for group in traces[0]:
        names = sorted(set().union(*(set(t.get(group, {})) for t in traces)))
        totals, rhats = [], []
        for name in names:
            chains = [t[group][name] for t in traces if name in t.get(group, {})]
            if (total := multichain_ess(chains)) is None:
                continue
            totals.append(total)
            if (r := mb.split_rhat(chains)) is not None:
                rhats.append(r)
        if totals:
            above = f"{sum(r > 1.1 for r in rhats) / len(rhats):.0%}" if rhats else "-"
            rmax = f"{max(rhats):.3g}" if rhats else "-"
            print(f"{group:<24}{len(totals):>8}{statistics.median(totals):>12.1f}{min(totals):>10.1f}{above:>11}{rmax:>11}")


def pair_orders(chains, pairs):
    for chain in chains:
        with gzip.open(chain / "parameters" / "infection_order.csv.gz", "rt") as f:
            rows = [line.strip() for line in f if line.strip()][2:]  # header, then the state before burn-in
        out = []
        for a, b in pairs:
            first = [r.split(",").index(a) < r.split(",").index(b) for r in rows]
            flips = sum(x != y for x, y in zip(first, first[1:]))
            out.append(f"{a} before {b}: {sum(first)}/{len(first)} (changes {flips})")
        print(f"{chain.name}: " + " | ".join(out))


def parent_probabilities(path, child):
    """Mean over draws of each parent's marginal probability, '(none)' for the latent parent alone."""
    latent = f"{{{child}_copy}}"
    sums, draws = defaultdict(float), set()
    for row in mb.read_rows(path)[1:]:
        members, prob, it = row.rsplit(",", 2)
        prob = float(prob)
        draws.add(it)
        if members == latent:
            sums["(none)"] += prob
        for member in members.strip("{}").split(";"):
            if member and not member.endswith("_copy") and member != "S":
                sums[member] += prob
    return {k: v / len(draws) for k, v in sums.items()}


def agreement(chains, show):
    children = sorted(f.name[: -len("_ps.csv.gz")] for f in (chains[0] / "parent_sets").glob("*_ps.csv.gz"))
    rows = []
    for child in children:
        per_chain = [parent_probabilities(c / "parent_sets" / f"{child}_ps.csv.gz", child) for c in chains]
        keys = set().union(*per_chain)
        pooled = {k: statistics.mean(m.get(k, 0.0) for m in per_chain) for k in keys}
        spread = {k: max(m.get(k, 0.0) for m in per_chain) - min(m.get(k, 0.0) for m in per_chain) for k in keys}
        top = max((k for k in keys if k != "(none)"), key=pooled.get, default="-")
        rows.append((max(spread.values()), child, top, pooled.get(top, 0.0), pooled.get("(none)", 0.0)))
    rows.sort(reverse=True)
    print(f"{'child':<10}{'top parent':<12}{'P(top)':>8}{'P(none)':>9}{'spread':>8}")
    for worst, child, top, p, none in rows[:show]:
        print(f"{child:<10}{top:<12}{p:>8.3f}{none:>9.3f}{worst:>8.3f}")
    spreads = sorted(r[0] for r in rows)
    print(f"largest between-chain difference in a parent probability, per infection: median {spreads[len(spreads) // 2]:.3f}, "
          f"90th percentile {spreads[int(0.9 * len(spreads))]:.3f}, max {spreads[-1]:.3f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("--chains", help="comma-separated chain numbers (default: all)")
    parser.add_argument("--pairs", nargs="*", default=[], help="infection pairs A:B")
    parser.add_argument("--show", type=int, default=15, help="infections listed, least agreement first")
    args = parser.parse_args()
    run = Path(args.run_dir)
    chains = sorted((p for p in run.glob("chain*") if p.is_dir()), key=lambda p: int(p.name[5:]))
    if args.chains:
        wanted = {f"chain{c}" for c in args.chains.split(",")}
        chains = [c for c in chains if c.name in wanted]
    print(f"{len(chains)} chains: {', '.join(c.name for c in chains)}\n")
    diagnostics([mb.run_traces(c) for c in chains])
    if args.pairs:
        print()
        pair_orders(chains, [p.split(":") for p in args.pairs])
    print()
    agreement(chains, args.show)


if __name__ == "__main__":
    main()
