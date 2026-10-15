"""The infection card's data (U10): what the posterior says about one Infection.

For an Infection of an Analysis: its expected false positive and false negative counts, duration and infection time
(observation time minus duration) per Chain with pooled quantiles; its COI per Draw (the derived `coi` table) and the
pooled distribution; its Latent genotype as posterior presence per Locus and Allele beside the Observed genotype; its
Parent sets per Chain and pooled, with members, composition class and the credible collection; its parent and child
Edge marginals; and its q, entropy and effective parent-set size. Pooled values weight each Chain by its Draws.
"""

from pt_analysis import functionals, inputs, times
from pt_analysis.chains import infection_ids, pooled_mean
from pt_analysis.derive import derived, derived_analysis, table
from pt_analysis.parent_sets import LEVEL, pooled_sets
from pt_analysis.quantiles import quantiles
from pt_analysis.sql import literal

# The card's quantities from the `infections` Trace: key, column, and the convergence table's group of its trace.
TRACE_QUANTITIES = [("false_positives", "expected_false_positives", "false positive rate"),
                    ("false_negatives", "expected_false_negatives", "false negative rate"),
                    ("duration", "infection_duration", "infection duration")]


class UnknownInfection(Exception):
    """An id that is no Infection of the Analysis."""


def series_summary(chains, group, name):
    """A quantity's per-Chain series and the 2.5%, 50% and 97.5% quantiles (R type 7) of its values pooled over every
    Chain's Draws (None without Draws)."""
    values = [v for c in chains for v in c["value"]]
    lower, median, upper = quantiles(values, (0.025, 0.5, 0.975)) if values else (None, None, None)
    return {"group": group, "name": name, "chains": chains, "lower": lower, "median": median, "upper": upper}


def quantities(analysis, infection, observation_time):
    """{false_positives, false_negatives, duration, infection_time}: each a `series_summary`. Infection time is the
    observation time minus the duration at each Draw, its quantiles the posterior infection time of
    `times.infection_times`; it has no trace in the convergence table (group None)."""
    per_chain = {key: [] for key, *_ in TRACE_QUANTITIES}
    series = []
    for label, store in analysis.runs.items():
        rows = table(store, "infections").filter(f"infection = {literal(infection)}").project(
            "iteration, expected_false_positives, expected_false_negatives, infection_duration").order(
            "iteration").fetchall()
        iterations = [r[0] for r in rows]
        for k, (key, *_) in enumerate(TRACE_QUANTITIES):
            per_chain[key].append({"chain": label, "iteration": iterations, "value": [r[k + 1] for r in rows]})
        series.append({"chain": label, "iteration": iterations, "value": [observation_time - r[3] for r in rows]})
    out = {key: series_summary(per_chain[key], group, infection) for key, _, group in TRACE_QUANTITIES}
    interval = times.infection_times(analysis, [infection]).get(infection, dict.fromkeys(("lower", "median", "upper")))
    out["infection_time"] = {"group": None, "name": None, "chains": series, **interval}
    return out


def coi(analysis, infection):
    """The COI per Draw per Chain from the derived `coi` table, and its distribution over every Chain's Draws:
    [{coi, draws, fraction}] by COI."""
    chains, counts = [], {}
    for label, store in analysis.runs.items():
        rows = derived(store, "coi").filter(f"infection = {literal(infection)}").project("iteration, coi").order(
            "iteration").fetchall()
        chains.append({"chain": label, "iteration": [r[0] for r in rows], "coi": [r[1] for r in rows]})
        for _, value in rows:
            counts[value] = counts.get(value, 0) + 1
    total = sum(counts.values())
    return {"chains": chains,
            "distribution": [{"coi": c, "draws": n, "fraction": n / total} for c, n in sorted(counts.items())]}


def genotype(analysis, infection):
    """Per Locus in the input's order, the number of Alleles (from `structure/loci`), the Observed genotype from the
    input (None where unobserved) and the posterior presence of each Allele pooled over Chains (from
    `genotype_summary`; 0 where no Draw has it)."""
    store = next(iter(analysis.runs.values()))
    n_alleles = dict(store.structure("loci").aggregate("locus, max(allele) + 1").fetchall())
    presence, draws = {}, {}
    for label, store in analysis.runs.items():
        draws[label] = store.draws()
        for locus, allele, p in derived(store, "genotype_summary").filter(f"infection = {literal(infection)}").project(
                "locus, allele, presence").fetchall():
            presence.setdefault((locus, allele), {})[label] = p
    total = sum(draws.values())
    observed = inputs.genotypes(inputs.nodes(analysis.input)[infection].get("observed_genotype"))
    order = inputs.loci(analysis.input)
    loci = []
    for locus in order + sorted(set(n_alleles) - set(order)):
        pooled_presence = [sum(p * draws[label] for label, p in presence.get((locus, a), {}).items()) / total
                           if total else 0.0 for a in range(n_alleles.get(locus, 0))]
        loci.append({"locus": locus, "alleles": n_alleles.get(locus, 0), "observed": observed.get(locus),
                     "presence": pooled_presence})
    return {"loci": loci}


def set_row(s, credible):
    return {"set_id": s.set_id, "members": list(s.members), "n_observed": s.n_observed, "class": s.composition,
            "weight": s.weight, "rank": s.rank, "cumulative_mass": s.cumulative_mass, "credible": credible}


def parent_sets(analysis, infection, level=LEVEL):
    """The Infection's Parent sets per Chain, ranked by the Chain's weights (`parent_sets.ranked_sets`), with members,
    composition class, cumulative mass and whether the set is in the Chain's credible collection at `level`. Also
    pooled over Chains (`parent_sets.pooled_sets`): each set's weight per Chain and pooled, ranked by the same MAP
    rule and with the credible collection by the same prefix rule; a pooled set's id is the first Chain's that has
    it."""
    chains = []
    for label, store in analysis.runs.items():
        sets = functionals.ranked_sets(store, infection=infection)
        credible = {s.set_id for s in functionals.credible_collection(store, infection, level)}
        chains.append({"chain": label, "sets": [set_row(s, s.set_id in credible) for s in sets]})
    pooled = []
    for s in pooled_sets(analysis, level, infection).get(infection, []):
        pooled.append({"set_id": next(iter(s.set_ids.values())), "members": list(s.members),
                       "n_observed": s.n_observed, "class": s.composition, "weight": s.weight, "rank": s.rank,
                       "cumulative_mass": s.cumulative_mass, "credible": s.credible,
                       "chains": [{"chain": label, "weight": w} for label, w in s.chains.items()]})
    return {"level": level, "pooled": pooled, "chains": chains}


def edge_marginals(analysis, infection):
    """The Infection's parent Edge marginals (it is the child; the unobserved source is its Latent parent id) and
    child Edge marginals (it is the parent), from the `edge_stats` table: pooled, spread, ESS and per Chain mean, with
    the trace quantity (group, name) of each; largest pooled first."""
    rows = derived_analysis(analysis, "edge_stats").filter(
        f"child = {literal(infection)} OR (parent = {literal(infection)} AND NOT unobserved)").project(
        "child, parent, unobserved, chain, mean, spread, ess").fetchall()
    edges = {}
    for child, parent, unobserved, chain, mean, spread, ess in rows:
        edge = edges.setdefault((child, parent), {
            "child": child, "parent": parent, "unobserved": unobserved,
            "group": "unobserved-source inclusion" if unobserved else "edge",
            "name": child if unobserved else f"{child}<-{parent}", "chains": []})
        if chain is None:
            edge.update(pooled=mean, spread=spread, ess=ess)
        else:
            edge["chains"].append({"chain": chain, "mean": mean})
    order = list(analysis.runs)
    for edge in edges.values():
        edge["chains"].sort(key=lambda c: order.index(c["chain"]))
    largest = lambda e: (-e["pooled"], e["child"], e["parent"])  # noqa: E731
    return (sorted((e for e in edges.values() if e["child"] == infection), key=largest),
            sorted((e for e in edges.values() if e["parent"] == infection), key=largest))


SUMMARY = ["q", "expected_observed_parents", "entropy", "effective_size"]


def summary(analysis, infection):
    """Per Chain the mean over its Draws of q, the expected observed-parent count, the entropy and the effective
    parent-set size from `infection_draws`, and each pooled over Chains (None without Draws)."""
    chains, draws = [], {}
    for label, store in analysis.runs.items():
        draws[label] = store.draws()
        means = derived(store, "infection_draws").filter(f"infection = {literal(infection)}").aggregate(
            ", ".join(f"avg({c})" for c in SUMMARY)).fetchone()
        chains.append({"chain": label, **dict(zip(SUMMARY, means))})
    with_draws = [c for c in chains if c["q"] is not None]
    return {"chains": chains, **{c: pooled_mean(with_draws, c, draws) for c in SUMMARY}}


def card(analysis, infection, level=LEVEL):
    """Everything the infection card shows for `infection`, as a dict:

        infection, observation_time
        quantities   {false_positives, false_negatives, duration, infection_time}: each {group, name, chains:
                     [{chain, iteration, value}], lower, median, upper}
        coi          {chains: [{chain, iteration, coi}], distribution: [{coi, draws, fraction}]}
        genotype     {loci: [{locus, alleles, observed, presence}]}
        parent_sets  {level, pooled: [set row + chains: [{chain, weight}]], chains: [{chain, sets: [set row]}]}, a
                     set row being {set_id, members, n_observed, class, weight, rank, cumulative_mass, credible}
        parents      [{child, parent, unobserved, group, name, pooled, spread, ess, chains: [{chain, mean}]}]
        children     likewise, with the Infection as parent
        summary      {chains: [{chain, q, expected_observed_parents, entropy, effective_size}], and each pooled}
    """
    if infection not in infection_ids(analysis):
        raise UnknownInfection(f"no Infection {infection} in the Analysis")
    observation_time = inputs.nodes(analysis.input)[infection]["observation_time"]
    parents, children = edge_marginals(analysis, infection)
    return {"infection": infection, "observation_time": observation_time,
            "quantities": quantities(analysis, infection, observation_time), "coi": coi(analysis, infection),
            "genotype": genotype(analysis, infection), "parent_sets": parent_sets(analysis, infection, level),
            "parents": parents, "children": children, "summary": summary(analysis, infection)}
