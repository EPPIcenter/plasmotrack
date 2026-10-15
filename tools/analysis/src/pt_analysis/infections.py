"""The Infections workspace's table: one row per Infection with what the posterior says about it, computed in DuckDB
from the derived tables of every Chain at once.

Per Infection, in input order:

    infection, observation_time, symptomatic       from the input
    infection_time_lower, _median, _upper          the posterior infection time (`times.infection_times`)
    q, q_ess, q_rhat                               unobserved-source inclusion pooled over every Draw of every Chain,
                                                   with its pooled ESS and split R-hat from the convergence table
    expected_observed_parents                      pooled mean of `infection_draws`
    multi_source, multi_source_observed            pooled means of `multi_source_augmented` and `_observed`
    posterior_entropy, posterior_effective_size    of the pooled Parent-set weights: -sum w ln w and its exponential
    entropy, effective_size                        pooled means of the per-Draw entropy and effective size
                                                   (`infection_draws`, as the infection card's summary has them)
    top_set, top_set_class, top_set_probability    the pooled MAP Parent set (`parent_sets` MAP rule): members
                                                   (observed parents in id order, then the Latent parent), class and
                                                   pooled weight
    top_parent, top_parent_marginal                the observed parent with the largest pooled Edge marginal (None
                                                   without one)
    out_degree                                     expected out-degree: pooled mean of `out_degree_draws`
    coi_median                                     median COI over every Draw (the lower median)
    worst_rhat, worst_rhat_group, worst_rhat_name  the largest pooled R-hat over the Infection's quantities (false
                                                   positive and negative rates, duration, q) and its parent edges
    epi:<column>                                   each epi column's value as text (`pt_analysis.epi`)

and per Chain, each Chain's mean q (`q_chains`). Pooling weights each Chain by its Draws, as everywhere else.

The two effective sizes answer different questions. The per-Draw one is how spread the Parent-set weights are given
one Draw's Ordering, genotypes and parameters; when each Draw allows a single set it is 1 even though the Draws
disagree on which. The pooled one is how spread the posterior over Parent sets is, the uncertainty a reader means.
"""

import math
from functools import reduce

import pyarrow as pa

from pt_analysis import inputs
from pt_analysis.derive import derive_analysis, derived, derived_analysis, table
from pt_analysis.epi import load_epi
from pt_analysis.sql import Scope
from pt_analysis.times import infection_times

# The convergence groups of an Infection's own quantities, each named by the Infection's id.
OWN_GROUPS = ("false positive rate", "false negative rate", "infection duration", "unobserved-source inclusion")


def _union(analysis, relation_of):
    """`relation_of(store)` of every Chain with the Chain's label as column `chain`, as one relation."""
    return reduce(lambda a, b: a.union(b), (relation_of(store).project(f"'{label}' AS chain, *")
                                            for label, store in analysis.runs.items()))


def _by_infection(rows):
    return {r[0]: r[1:] for r in rows}


def infections(analysis):
    """The table's payload: {chains, columns: {column: [value per Infection]}, q_chains: {Chain: [mean q]},
    epi_columns: [{name, type}]}. See the module docstring."""
    derive_analysis(analysis)
    chains = list(analysis.runs)
    draws = {label: store.draws() for label, store in analysis.runs.items()}
    total = sum(draws.values())
    nodes = inputs.nodes(analysis.input)
    ids = list(nodes)
    connection = next(iter(analysis.runs.values())).connection

    with Scope(connection) as s:
        per_draw = s.name(_union(analysis, lambda st: derived(st, "infection_draws").project(
            "infection, q, expected_observed_parents, entropy, effective_size, multi_source_augmented, "
            "multi_source_observed")))
        means = _by_infection(s.rows(f"""
            SELECT infection, avg(q), avg(expected_observed_parents), avg(multi_source_augmented),
                   avg(multi_source_observed), avg(entropy), avg(effective_size)
            FROM {per_draw} GROUP BY infection"""))
        q_chains = {(chain, infection): q for chain, infection, q in s.rows(
            f"SELECT chain, infection, avg(q) FROM {per_draw} GROUP BY ALL")}

        # Parent sets pooled over Chains, matched by members; each Chain's weight counts by its Draws.
        weights = {label: n / total if total else 0.0 for label, n in draws.items()}
        sets = s.name(_union(analysis, lambda st: derived(st, "set_summary").project("infection, set_id, weight")))
        classes = s.name(_union(analysis, lambda st: derived(st, "set_classes").project('infection, set_id, "class"')))
        members = s.name(_union(analysis, lambda st: table(st, "parent_set_members")))
        latent = s.name(_union(analysis, lambda st: st.structure("infections").project("infection, latent_parent")))
        share = s.name(pa.table({"chain": list(weights), "share": list(weights.values())}))
        top_sets = _by_infection(s.rows(f"""
            WITH named AS (
                SELECT m.chain, m.infection, m.set_id,
                       list(m.member ORDER BY m.member = l.latent_parent, m.member) AS members
                FROM {members} m JOIN {latent} l USING (chain, infection) GROUP BY ALL),
            pooled AS (
                SELECT x.infection, n.members, any_value(c."class") AS "class", sum(x.weight * w.share) AS weight
                FROM {sets} x JOIN named n USING (chain, infection, set_id)
                     JOIN {classes} c USING (chain, infection, set_id) JOIN {share} w USING (chain)
                GROUP BY x.infection, n.members)
            SELECT infection, arg_min(members, (-weight, members)), arg_min("class", (-weight, members)),
                   max(weight), coalesce(-sum(weight * ln(weight)) FILTER (WHERE weight > 0), 0)
            FROM pooled GROUP BY infection"""))

        out_degree = dict(s.rows(f"""
            SELECT infection, avg(out_degree)
            FROM {s.name(_union(analysis, lambda st: derived(st, "out_degree_draws")))} GROUP BY infection"""))
        coi = dict(s.rows(f"""
            SELECT infection, quantile_disc(coi, 0.5)
            FROM {s.name(_union(analysis, lambda st: derived(st, "coi")))} GROUP BY infection"""))

        edge_stats = s.name(derived_analysis(analysis, "edge_stats"))
        top_parent = _by_infection(s.rows(f"""
            SELECT child, arg_min(parent, (-mean, parent)), max(mean)
            FROM {edge_stats} WHERE chain IS NULL AND NOT unobserved GROUP BY child"""))

        convergence = s.name(derived_analysis(analysis, "convergence"))
        groups = ", ".join(f"'{g}'" for g in OWN_GROUPS)
        mixing = s.name(s.sql(f"""
            SELECT CASE WHEN "group" = 'edge' THEN split_part(name, '<-', 1) ELSE name END AS infection,
                   "group", name, ess, rhat
            FROM {convergence} WHERE chain IS NULL AND ("group" IN ({groups}) OR "group" = 'edge')"""))
        q_mixing = _by_infection(s.rows(f"""
            SELECT infection, ess, rhat FROM {mixing} WHERE "group" = 'unobserved-source inclusion'"""))
        worst = _by_infection(s.rows(f"""
            SELECT infection, max(rhat), arg_max("group", rhat), arg_max(name, rhat)
            FROM {mixing} WHERE rhat IS NOT NULL AND NOT isnan(rhat) GROUP BY infection"""))

    times = infection_times(analysis)
    blank = dict.fromkeys(("lower", "median", "upper"))
    epi = load_epi(analysis)
    epi_columns = list(epi.columns) if epi is not None else []

    def column(f):
        return [f(i) for i in ids]

    def at(source, k, default=None):
        return lambda i: source[i][k] if i in source else default

    columns = {
        "infection": ids,
        "observation_time": column(lambda i: nodes[i].get("observation_time")),
        "symptomatic": column(lambda i: nodes[i].get("symptomatic", True)),
        "infection_time_lower": column(lambda i: times.get(i, blank)["lower"]),
        "infection_time_median": column(lambda i: times.get(i, blank)["median"]),
        "infection_time_upper": column(lambda i: times.get(i, blank)["upper"]),
        "q": column(at(means, 0)),
        "q_ess": column(at(q_mixing, 0)),
        "q_rhat": column(at(q_mixing, 1)),
        "expected_observed_parents": column(at(means, 1)),
        "multi_source": column(at(means, 2)),
        "multi_source_observed": column(at(means, 3)),
        "posterior_entropy": column(at(top_sets, 3)),
        "posterior_effective_size": column(lambda i: None if i not in top_sets else math.exp(top_sets[i][3])),
        "entropy": column(at(means, 4)),
        "effective_size": column(at(means, 5)),
        "top_set": column(at(top_sets, 0)),
        "top_set_class": column(at(top_sets, 1)),
        "top_set_probability": column(at(top_sets, 2)),
        "top_parent": column(at(top_parent, 0)),
        "top_parent_marginal": column(at(top_parent, 1)),
        "out_degree": column(lambda i: out_degree.get(i)),
        "coi_median": column(lambda i: coi.get(i)),
        "worst_rhat": column(at(worst, 0)),
        "worst_rhat_group": column(at(worst, 1)),
        "worst_rhat_name": column(at(worst, 2)),
    }
    for name in epi_columns:
        values = epi.values(name)
        columns[f"epi:{name}"] = column(lambda i: values.get(i))
    return {"chains": chains, "columns": columns,
            "q_chains": {label: [q_chains.get((label, i)) for i in ids] for label in chains},
            "epi_columns": [{"name": name, "type": epi.columns[name]} for name in epi_columns]}

