"""The pair inspector's data (U3): why do Infections A and B sit in one order?

For a pair (A, B) of an Analysis: per Chain, the fraction of Draws at which A precedes B in the Ordering, how often
that order changes between adjacent Draws, and the per-Draw indicator series; both Edge marginals per Chain (B <- A
and A <- B); the co-inclusion of A and B in each third Infection's Parent sets; genotypes per Locus (Observed from the
input, Latent posterior mode with its agreement, true Latent genotype when the Analysis has truth); and each
Infection's posterior infection time and the lag between their observation times.
"""

# The SQL below names DuckDB relations held in local variables (`FROM edges`); DuckDB finds them by name in the
# calling frame (a replacement scan), so those variables are used even though nothing else reads them.

from pt_analysis import inputs
from pt_analysis.derive import derived, table
from pt_analysis.summary import ordered_before, pair_order


class PairError(Exception):
    """The pair cannot be inspected: an id that is no Infection of the Analysis, or the same id twice."""


class UnknownInfection(PairError):
    """An id that is no Infection of the Analysis."""


def order_summary(iterations, before):
    """The order of a pair in one Chain from per-Draw booleans `before` (A precedes B) at Sampling `iterations`: the
    indicator series as 0/1, the number of Draws and of Draws with A first, their fraction (None without Draws) and
    the switches, the adjacent Draws at which the order changes."""
    count, draws, switches = pair_order(before)
    return {"iteration": list(iterations), "before": [int(b) for b in before], "draws": draws, "before_count": count,
            "probability": count / draws if draws else None, "switches": switches}


def pooled(chains, key, draws):
    """The mean over every Chain's Draws of a per-Chain mean `key`, weighting each Chain by its Draws; None without
    Draws."""
    total = sum(draws[c["chain"]] for c in chains)
    return sum(c[key] * draws[c["chain"]] for c in chains) / total if total else None


def order_iterations(store, a):
    """The Sampling iterations of the Draws at which Infection `a` has a place in the Ordering, in order."""
    ordering = store.table("ordering")
    if ordering is None:
        return []
    rows = store.connection.execute("SELECT iteration FROM ordering WHERE infection = $a ORDER BY iteration",
                                    {"a": a}).fetchall()
    return [iteration for (iteration,) in rows]


def pair_order_chains(analysis, a, b):
    """The order of A and B per Chain (`order_summary`, with `chain` the Chain's label) and P(A before B) pooled over
    every Chain's Draws (None without Draws). Reads the `ordering` Trace as `summary.ordered_before` does."""
    chains = [{"chain": label, **order_summary(order_iterations(store, a), ordered_before(store, a, b))}
              for label, store in analysis.runs.items()]
    draws = sum(c["draws"] for c in chains)
    return {"chains": chains, "probability": sum(c["before_count"] for c in chains) / draws if draws else None}


def edge_marginals(analysis, a, b):
    """Both Edge marginals of the pair per Chain, the mean over Draws of the `edges` weight: `a_to_b` with A the
    parent of B (B <- A) and `b_to_a`; 0 where the edge has no row at any Draw. Also pooled over Chains."""
    chains, draws = [], {}
    for label, store in analysis.runs.items():
        draws[label] = store.draws()
        edges = derived(store, "edges")
        found = dict(store.connection.execute("""
            SELECT parent, sum(prob) FROM edges
            WHERE (child = $b AND parent = $a) OR (child = $a AND parent = $b) GROUP BY parent""",
                                              {"a": a, "b": b}).fetchall())
        n = max(draws[label], 1)
        chains.append({"chain": label, "a_to_b": found.get(a, 0.0) / n, "b_to_a": found.get(b, 0.0) / n})
    return {"chains": chains, "a_to_b": pooled(chains, "a_to_b", draws), "b_to_a": pooled(chains, "b_to_a", draws)}


def co_inclusion(analysis, a, b):
    """Per third Infection with a Parent set holding both A and B in some Chain, the summed pooled weight of those
    sets per Chain (from `set_summary`, so a mean over the Chain's Draws) and over Chains; largest first."""
    per_infection, draws = {}, {}
    for label, store in analysis.runs.items():
        draws[label] = store.draws()
        sets = derived(store, "set_summary")
        members = table(store, "parent_set_members")
        rows = store.connection.execute("""
            WITH both_sets AS (SELECT infection, set_id FROM members WHERE member IN ($a, $b)
                               GROUP BY ALL HAVING count(DISTINCT member) = 2)
            SELECT s.infection, sum(s.weight) FROM sets s JOIN both_sets USING (infection, set_id)
            WHERE s.infection NOT IN ($a, $b) GROUP BY s.infection""", {"a": a, "b": b}).fetchall()
        for infection, weight in rows:
            per_infection.setdefault(infection, {})[label] = weight
    out = []
    for infection, weights in per_infection.items():
        chains = [{"chain": label, "weight": weights.get(label, 0.0)} for label in analysis.runs]
        out.append({"infection": infection, "weight": pooled(chains, "weight", draws), "chains": chains})
    return sorted(out, key=lambda r: (-r["weight"], r["infection"]))


def latent_modes(analysis, infection):
    """{locus: (posterior mode bitstring, agreement)} of an Infection's Latent genotype pooled over Chains: the
    bitstring most Draws of all Chains had (ties to the smaller bitstring) and the fraction of Draws that had it."""
    counts = {}
    for store in analysis.runs.values():
        modes = derived(store, "genotype_mode")
        for locus, alleles, draws in store.connection.execute(
                "SELECT locus, alleles, draws FROM modes WHERE infection = $i", {"i": infection}).fetchall():
            by_alleles = counts.setdefault(locus, {})
            by_alleles[alleles] = by_alleles.get(alleles, 0) + draws
    out = {}
    for locus, by_alleles in counts.items():
        mode = min(by_alleles, key=lambda alleles: (-by_alleles[alleles], alleles))
        out[locus] = (mode, by_alleles[mode] / sum(by_alleles.values()))
    return out


def genotype_rows(analysis, a, b):
    """One row per Locus, in the input's order: Observed A (None where unobserved), Latent A's pooled mode and
    agreement, Latent B's, Observed B, whether the two modes differ, and with truth the true Latent genotypes.
    Also the number of differing Loci and whether truth is shown."""
    node = inputs.nodes(analysis.input)
    modes = {i: latent_modes(analysis, i) for i in (a, b)}
    observed = {i: inputs.genotypes(node.get(i, {}).get("observed_genotype")) for i in (a, b)}
    truth = None
    if analysis.truth is not None:
        true_nodes = inputs.nodes(analysis.truth)
        truth = {i: inputs.genotypes(true_nodes.get(i, {}).get("flat_latent_genotype")) for i in (a, b)}
    order = inputs.loci(analysis.input)
    names = order + sorted((set(modes[a]) | set(modes[b])) - set(order))
    rows = []
    for locus in names:
        latent_a, agreement_a = modes[a].get(locus, (None, None))
        latent_b, agreement_b = modes[b].get(locus, (None, None))
        row = {"locus": locus, "observed_a": observed[a].get(locus), "latent_a": latent_a, "agreement_a": agreement_a,
               "latent_b": latent_b, "agreement_b": agreement_b, "observed_b": observed[b].get(locus),
               "differs": latent_a != latent_b}
        if truth is not None:
            row.update(true_a=truth[a].get(locus), true_b=truth[b].get(locus))
        rows.append(row)
    return {"rows": rows, "differing": sum(r["differs"] for r in rows), "truth": truth is not None}


def quantile(values, p):
    """The p quantile of `values` by linear interpolation between order statistics at (n - 1) p (R's type 7)."""
    values = sorted(values)
    h = (len(values) - 1) * p
    lo = int(h)
    return values[lo] + (h - lo) * (values[min(lo + 1, len(values) - 1)] - values[lo])


def infection_time(analysis, infection, observation_time):
    """An Infection's posterior infection time: its observation time minus the median, 97.5% and 2.5% quantiles of
    its `infection_duration` pooled over every Chain's Draws, as {observation_time, median, lower, upper} (None
    without Draws)."""
    durations = []
    for store in analysis.runs.values():
        trace = table(store, "infections")
        durations += [d for (d,) in store.connection.execute(
            "SELECT infection_duration FROM trace WHERE infection = $i", {"i": infection}).fetchall()]
    if not durations:
        return {"observation_time": observation_time, "median": None, "lower": None, "upper": None}
    return {"observation_time": observation_time, "median": observation_time - quantile(durations, 0.5),
            "lower": observation_time - quantile(durations, 0.975),
            "upper": observation_time - quantile(durations, 0.025)}


def infection_times(analysis, a, b):
    """Each Infection's `infection_time` and the lag, B's observation time minus A's, from the input (which gives
    every Infection an observation time)."""
    node = inputs.nodes(analysis.input)
    observed = {i: node[i]["observation_time"] for i in (a, b)}
    return {"a": infection_time(analysis, a, observed[a]), "b": infection_time(analysis, b, observed[b]),
            "lag": observed[b] - observed[a]}


def infections(analysis):
    """The Infection ids of every Chain of `analysis`."""
    return {infection for store in analysis.runs.values()
            for (infection,) in store.structure("infections").project("infection").fetchall()}


def inspect_pair(analysis, a, b):
    """Everything the pair inspector shows for Infections `a` and `b` of `analysis`, as a dict:

        order         {chains: [{chain, iteration, before, draws, before_count, probability, switches}], probability}
        edges         {chains: [{chain, a_to_b, b_to_a}], a_to_b, b_to_a}
        co_inclusion  [{infection, weight, chains: [{chain, weight}]}]
        genotypes     {rows: [{locus, observed_a, latent_a, agreement_a, latent_b, agreement_b, observed_b, differs
                      [, true_a, true_b]}], differing, truth}
        times         {a: {observation_time, median, lower, upper}, b: likewise, lag}
    """
    if a == b:
        raise PairError(f"a pair needs two Infections, not {a} twice")
    known = infections(analysis)
    for infection in (a, b):
        if infection not in known:
            raise UnknownInfection(f"no Infection {infection} in the Analysis")
    return {"a": a, "b": b, "order": pair_order_chains(analysis, a, b), "edges": edge_marginals(analysis, a, b),
            "co_inclusion": co_inclusion(analysis, a, b), "genotypes": genotype_rows(analysis, a, b),
            "times": infection_times(analysis, a, b)}
