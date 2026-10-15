"""The chain summary of an Analysis: convergence diagnostics per group, how often each Chain had one Infection before
another, and how far the Chains agree on each Infection's parents. Reads the derived tables (deriving them first when
they are not current) and, for pairs, the `ordering` Trace.
"""

# The SQL below names DuckDB relations held in local variables (`FROM edges`); DuckDB finds them by name in the
# calling frame (a replacement scan), so those variables are used even though nothing else reads them.

import statistics
from dataclasses import dataclass

from pt_analysis.derive import GROUPS, derive_analysis, derived, derived_analysis

# The key, beside the observed parents' ids, of the probability that an Infection's set holds no observed parent.
UNOBSERVED_ONLY = "(unobserved only)"


class SummaryError(Exception):
    """The summary cannot be computed from what it was given; the message says why."""


@dataclass(frozen=True)
class GroupDiagnostics:
    """A convergence group's pooled rows: how many quantities have an ESS, the median and smallest ESS, the share of
    quantities with an R-hat above 1.1 and the largest R-hat (None when no quantity has an R-hat)."""

    group: str
    series: int
    ess_median: float
    ess_min: float
    rhat_above: float | None
    rhat_max: float | None


@dataclass(frozen=True)
class Agreement:
    """An Infection's most probable observed parent pooled over Chains ('-' if it has none), that parent's and the
    unobserved-only class's pooled probability, and the largest between-Chain difference in any of those."""

    infection: str
    top: str
    p_top: float
    p_unobserved_only: float
    spread: float


def diagnostics(analysis):
    """Per convergence group with at least one ESS, in group order, its pooled diagnostics."""
    convergence = derived_analysis(analysis, "convergence")
    rows = convergence.filter("chain IS NULL").project('"group", ess, rhat').fetchall()
    out = []
    for group in GROUPS:
        ess = [e for g, e, _ in rows if g == group and e is not None]
        rhat = [r for g, e, r in rows if g == group and e is not None and r is not None]
        if ess:
            out.append(GroupDiagnostics(group, len(ess), statistics.median(ess), min(ess),
                                        sum(r > 1.1 for r in rhat) / len(rhat) if rhat else None,
                                        max(rhat) if rhat else None))
    return out


def pair_order(before):
    """(Draws with A before B, Draws, switches between adjacent Draws) from per-Draw booleans `before`."""
    return sum(before), len(before), sum(x != y for x, y in zip(before, before[1:]))


def ordered_before(store, a, b):
    """Per Draw, in Sampling-iteration order, whether Infection `a` precedes `b` in the Ordering."""
    ordering = store.table("ordering")
    if ordering is None:
        return []
    rows = store.connection.execute("""
        SELECT x.rank < y.rank FROM ordering x JOIN ordering y USING (iteration)
        WHERE x.infection = $a AND y.infection = $b ORDER BY iteration""", {"a": a, "b": b}).fetchall()
    return [before for (before,) in rows]


def parent_probabilities(store):
    """{infection: {parent: pooled Edge marginal, UNOBSERVED_ONLY: pooled unobserved-only weight}} for one Run, over
    its observed parents; {} without Draws."""
    draws = store.draws()
    if not draws:
        return {}
    edges = derived(store, "edges")
    infection_draws = derived(store, "infection_draws")
    latent = store.structure("infections")
    unobserved = "unobserved_only" if "unobserved_only" in infection_draws.columns else "0"
    out = {infection: {UNOBSERVED_ONLY: p} for infection, p in store.connection.sql(f"""
        SELECT infection, sum({unobserved}) / {draws} FROM infection_draws GROUP BY infection""").fetchall()}
    for child, parent, p in store.connection.sql(f"""
            SELECT child, parent, sum(prob) / {draws} FROM edges
            WHERE parent NOT IN (SELECT latent_parent FROM latent) GROUP BY ALL""").fetchall():
        out.setdefault(child, {UNOBSERVED_ONLY: 0.0})[parent] = p
    return out


def agreement(chains):
    """Per Infection of any Chain, from per-Chain `parent_probabilities`: the probabilities pooled as means over
    Chains (missing counts as 0), the top observed parent by pooled probability, and the largest between-Chain spread
    (max minus min) of any parent's or the unobserved-only probability. Sorted by spread then Infection, descending."""
    rows = []
    for infection in set().union(*chains):
        per_chain = [c.get(infection, {}) for c in chains]
        keys = set().union(*per_chain)
        pooled = {k: statistics.mean(c.get(k, 0.0) for c in per_chain) for k in keys}
        spread = max(max(c.get(k, 0.0) for c in per_chain) - min(c.get(k, 0.0) for c in per_chain) for k in keys)
        top = max(sorted(k for k in keys if k != UNOBSERVED_ONLY), key=pooled.get, default="-")
        rows.append(Agreement(infection, top, pooled.get(top, 0), pooled.get(UNOBSERVED_ONLY, 0), spread))
    return sorted(rows, key=lambda r: (r.spread, r.infection), reverse=True)


def spread_quantiles(rows):
    """(median, 90th percentile, max) of the spreads, taking the element at index n // 2 and int(0.9 n) of the sorted
    spreads, as chain_summary.py did."""
    spreads = sorted(r.spread for r in rows)
    return spreads[len(spreads) // 2], spreads[int(0.9 * len(spreads))], spreads[-1]


def print_summary(analysis, pairs=(), show=15):
    """Print the chain summary of `analysis`: the diagnostics table, a line per Chain for `pairs` of Infection ids
    (A, B), and the `show` Infections whose Chains agree least."""
    derive_analysis(analysis)
    known = {infection for store in analysis.runs.values()
             for (infection,) in store.structure("infections").project("infection").fetchall()}
    for a, b in pairs:
        for infection in (a, b):
            if infection not in known:
                raise SummaryError(f"no Infection {infection} in the Analysis")

    print(f"{len(analysis.runs)} Chains: {', '.join(analysis.runs)}\n")
    print(f"{'group':<30}{'series':>8}{'ESS median':>12}{'ESS min':>10}{'R-hat>1.1':>11}{'R-hat max':>11}")
    for d in diagnostics(analysis):
        above = f"{d.rhat_above:.0%}" if d.rhat_above is not None else "-"
        rmax = f"{d.rhat_max:.3g}" if d.rhat_max is not None else "-"
        print(f"{d.group:<30}{d.series:>8}{d.ess_median:>12.1f}{d.ess_min:>10.1f}{above:>11}{rmax:>11}")

    if pairs:
        print()
        for label, store in analysis.runs.items():
            orders = []
            for a, b in pairs:
                count, draws, switches = pair_order(ordered_before(store, a, b))
                orders.append(f"{a} before {b}: {count}/{draws} (changes {switches})")
            print(f"{label}: " + " | ".join(orders))

    rows = agreement([parent_probabilities(store) for store in analysis.runs.values()])
    print()
    print(f"{'infection':<12}{'top parent':<12}{'P(top)':>8}{'P(unobserved only)':>20}{'spread':>8}")
    for r in rows[:show]:
        print(f"{r.infection:<12}{r.top:<12}{r.p_top:>8.3f}{r.p_unobserved_only:>20.3f}{r.spread:>8.3f}")
    if rows:
        median, p90, largest = spread_quantiles(rows)
        print(f"largest between-Chain difference in a parent probability, per Infection: median {median:.3f}, "
              f"90th percentile {p90:.3f}, max {largest:.3f}")
