"""The pair inspector's data (U3): why do Infections A and B sit in one order?

For a pair (A, B) of an Analysis: per Chain, the fraction of Draws at which A precedes B in the Ordering and how often
that order changes between adjacent Draws (from the derived `order_pairs`), and the per-Draw indicator series (read
from the `ordering` Trace on demand); both Edge marginals per Chain (B <- A
and A <- B); the co-inclusion of A and B in each third Infection's Parent sets; genotypes per Locus (Observed from the
input, Latent posterior mode with its agreement, true Latent genotype when the Analysis has truth); and each
Infection's posterior infection time and the lag between their observation times.
"""

from pt_analysis import inputs, times
from pt_analysis.chains import infection_ids, pooled_mean
from pt_analysis.derive import derived, table
from pt_analysis.sql import Scope, literal
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


def order_iterations(store, a):
    """The Sampling iterations of the Draws at which Infection `a` has a place in the Ordering, in order."""
    ordering = store.table("ordering")
    if ordering is None:
        return []
    return [iteration for (iteration,) in ordering.filter(f"infection = {literal(a)}").project("iteration")
            .order("iteration").fetchall()]


def order_counts(store, a, b):
    """(Draws with both A and B in the Ordering, Draws with A first, switches) from the Run's derived `order_pairs`,
    which holds each pair once, ids in text order; (0, 0, 0) when the pair has no row."""
    i, j = sorted((a, b))
    found = derived(store, "order_pairs").filter(f"i = {literal(i)} AND j = {literal(j)}").project(
        "draws, before, switches").fetchone()
    if found is None:
        return 0, 0, 0
    draws, before, switches = found
    return draws, before if a == i else draws - before, switches


def pair_order_chains(analysis, a, b):
    """The order of A and B per Chain, with `chain` the Chain's label: the per-Draw indicator series of
    `order_summary`, read from the `ordering` Trace as `summary.ordered_before` does, and the counts and probability
    from `order_pairs`. Also P(A before B) pooled over every Chain's Draws (None without Draws)."""
    chains = []
    for label, store in analysis.runs.items():
        draws, before, switches = order_counts(store, a, b)
        series = order_summary(order_iterations(store, a), ordered_before(store, a, b))
        chains.append({**series, "chain": label, "draws": draws, "before_count": before,
                       "probability": before / draws if draws else None, "switches": switches})
    draws = sum(c["draws"] for c in chains)
    return {"chains": chains, "probability": sum(c["before_count"] for c in chains) / draws if draws else None}


def edge_marginals(analysis, a, b):
    """Both Edge marginals of the pair per Chain, the mean over Draws of the `edges` weight: `a_to_b` with A the
    parent of B (B <- A) and `b_to_a`; 0 where the edge has no row at any Draw. Also pooled over Chains."""
    chains, draws = [], {}
    for label, store in analysis.runs.items():
        draws[label] = store.draws()
        with Scope(store.connection) as s:
            found = dict(s.rows(f"""
                SELECT parent, sum(prob) FROM {s.name(derived(store, "edges"))}
                WHERE (child = $b AND parent = $a) OR (child = $a AND parent = $b) GROUP BY parent""",
                                {"a": a, "b": b}))
        n = max(draws[label], 1)
        chains.append({"chain": label, "a_to_b": found.get(a, 0.0) / n, "b_to_a": found.get(b, 0.0) / n})
    return {"chains": chains, "a_to_b": pooled_mean(chains, "a_to_b", draws),
            "b_to_a": pooled_mean(chains, "b_to_a", draws)}


def co_inclusion(analysis, a, b):
    """Per third Infection with a Parent set holding both A and B in some Chain, the summed pooled weight of those
    sets per Chain (from `set_summary`, so a mean over the Chain's Draws) and over Chains; largest first."""
    per_infection, draws = {}, {}
    for label, store in analysis.runs.items():
        draws[label] = store.draws()
        with Scope(store.connection) as s:
            rows = s.rows(f"""
                WITH both_sets AS (SELECT infection, set_id FROM {s.name(table(store, "parent_set_members"))}
                                   WHERE member IN ($a, $b) GROUP BY ALL HAVING count(DISTINCT member) = 2)
                SELECT s.infection, sum(s.weight) FROM {s.name(derived(store, "set_summary"))} s
                     JOIN both_sets USING (infection, set_id)
                WHERE s.infection NOT IN ($a, $b) GROUP BY s.infection""", {"a": a, "b": b})
        for infection, weight in rows:
            per_infection.setdefault(infection, {})[label] = weight
    out = []
    for infection, weights in per_infection.items():
        chains = [{"chain": label, "weight": weights.get(label, 0.0)} for label in analysis.runs]
        out.append({"infection": infection, "weight": pooled_mean(chains, "weight", draws), "chains": chains})
    return sorted(out, key=lambda r: (-r["weight"], r["infection"]))


def latent_modes(analysis, infection):
    """{locus: (posterior mode bitstring, agreement)} of an Infection's Latent genotype pooled over Chains: the
    bitstring most Draws of all Chains had (ties to the smaller bitstring) and the fraction of Draws that had it."""
    return latent_modes_of(analysis, [infection]).get(infection, {})


def latent_modes_of(analysis, infections):
    """{infection: latent_modes(analysis, infection)} for every one of `infections` with a `genotype_mode` row, read
    with one query per Chain."""
    ids = sorted(set(infections))
    counts = {}
    for store in analysis.runs.values():
        with Scope(store.connection) as s:
            rows = s.rows(f"""SELECT infection, locus, alleles, draws FROM {s.name(derived(store, "genotype_mode"))}
                              WHERE list_contains($ids, infection)""", {"ids": ids})
        for infection, locus, alleles, draws in rows:
            by_alleles = counts.setdefault(infection, {}).setdefault(locus, {})
            by_alleles[alleles] = by_alleles.get(alleles, 0) + draws
    out = {}
    for infection, loci in counts.items():
        modes = out[infection] = {}
        for locus, by_alleles in loci.items():
            mode = min(by_alleles, key=lambda alleles: (-by_alleles[alleles], alleles))
            modes[locus] = (mode, by_alleles[mode] / sum(by_alleles.values()))
    return out


def differing_loci(modes_a, modes_b):
    """The number of Loci at which two Infections' Latent modes (`latent_modes`) differ, a Locus without a mode
    counting as None, as `genotype_rows` counts them."""
    return sum(modes_a.get(locus, (None,))[0] != modes_b.get(locus, (None,))[0] for locus in modes_a.keys() | modes_b)


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


def infection_times(analysis, a, b):
    """Each Infection's observation time and posterior infection time ({observation_time, median, lower, upper},
    from `times.infection_times`; None without Draws) and the lag, B's observation time minus A's, from the input
    (which gives every Infection an observation time)."""
    observed = times.observation_times(analysis)
    intervals = times.infection_times(analysis, (a, b))
    empty = {"median": None, "lower": None, "upper": None}
    return {"a": {"observation_time": observed[a], **intervals.get(a, empty)},
            "b": {"observation_time": observed[b], **intervals.get(b, empty)}, "lag": times.lag(observed, a, b)}


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
    known = infection_ids(analysis)
    for infection in (a, b):
        if infection not in known:
            raise UnknownInfection(f"no Infection {infection} in the Analysis")
    return {"a": a, "b": b, "order": pair_order_chains(analysis, a, b), "edges": edge_marginals(analysis, a, b),
            "co_inclusion": co_inclusion(analysis, a, b), "genotypes": genotype_rows(analysis, a, b),
            "times": infection_times(analysis, a, b)}
