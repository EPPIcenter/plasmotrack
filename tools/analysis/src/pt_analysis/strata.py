"""Functionals of an Analysis overall and per stratum of an epi column (U9), pooled over its Chains.

A stratum is the set of Infections with one value of the column; Infections with no value (no row in the epi data,
or an empty cell) form one more stratum, whose value is None. Without a column there is one stratum, every
Infection, also with value None. Strata partition the Infections, so each functional's stratum values, weighted by
the strata's Infection counts, average back to its overall value, Draw by Draw.

Each functional is a mean over the stratum's Infections, all of them, including those with no weight in it. Which
Infections count follows the manuscript's Zanzibar summaries (`scripts/summarize_zanzibar_functionals.R`):

- Out-degree, psi_edge and psi_MAP, counts an Infection as a source: the expected number of observed children of the
  stratum's Infections (from `out_degree_draws`), and the number of children whose modal set names them.
- Everything else counts an Infection as a recipient, from its own Parent sets: unobserved-source inclusion q, the
  composition classes, the multi-source probabilities (`infection_draws`) and the parent-set size range.

Per-Draw functionals are named as in `functional_draws`: `psi_edge`, `mean_q`, `mean_<class>` per composition class,
`mean_multi_source_augmented` and `mean_multi_source_observed`. A functional's estimate is the mean of its per-Draw
values over every Draw of every Chain; its interval is their 2.5% and 97.5% quantiles, type 7 (R's default), as in
the manuscript's pipeline. psi_MAP and the size range come from the pooled Parent-set weights and have no interval.
"""

import math
from dataclasses import dataclass, field

import pyarrow as pa

from pt_analysis.derive import class_labels, derived
from pt_analysis.epi import EpiError, load_epi
from pt_analysis.parent_sets import LEVEL, pooled_sets
from pt_analysis.quantiles import Interval, interval
from pt_analysis.sql import Scope

MULTI_SOURCE = ("multi_source_augmented", "multi_source_observed")


class StrataError(Exception):
    """A per-Draw functional or a stratum value the Analysis does not have."""


@dataclass(frozen=True)
class Stratum:
    """The functionals of one stratum. `value` is None for every Infection (no column) or for those with no value."""

    value: str | None
    n: int
    psi_map: float | None = None
    size_range: tuple[float, float] | None = None
    psi_edge: Interval | None = None
    mean_q: Interval | None = None
    classes: dict[str, Interval] = field(default_factory=dict)  # composition class label -> its mean weight
    multi_source_augmented: Interval | None = None
    multi_source_observed: Interval | None = None


def class_key(label):
    """Composition classes by observed-parent count, each count's sets without the unobserved source first,
    `unobserved_only` last, as `derive.class_labels` orders them."""
    if label == "unobserved_only":
        return (math.inf, False)
    count, _, rest = label.removeprefix("observed_").partition("_")
    return (int(count), rest == "unobserved")


def strata(analysis, column=None, epi=None):
    """[(value, [Infection, ...])] in the column's order, the Infections with no value last; one stratum of every
    Infection without a column. Raises EpiError when the Analysis has no epi data or the data has no `column`."""
    if column is None:
        store = next(iter(analysis.runs.values()))
        return [(None, sorted(i for (i,) in store.structure("infections").select("infection").fetchall()))]
    epi = epi or load_epi(analysis)
    if epi is None:
        raise EpiError(f"the Analysis {analysis.label} has no epi data; its manifest names no 'epi' CSV")
    values = epi.values(column)
    return [(value, [i for i, v in values.items() if v == value]) for value in epi.stratum_values(column)]


def labels(analysis):
    """The composition class labels of every Chain, in `class_key` order."""
    return sorted({label for store in analysis.runs.values() for label in class_labels(store)}, key=class_key)


def per_draw(analysis, groups):
    """Per Chain, per Draw and stratum, every per-Draw functional: {chain: [(iteration, stratum index, {name:
    value})]} in Sampling-iteration then stratum order. `groups` lists each stratum's Infections."""
    classes = labels(analysis)
    names = ["psi_edge", "mean_q", *(f"mean_{c}" for c in classes), *(f"mean_{m}" for m in MULTI_SOURCE)]
    n = [len(infections) for infections in groups]
    membership = pa.table({"infection": [i for infections in groups for i in infections],
                           "stratum": [k for k, infections in enumerate(groups) for _ in infections]})
    result = {}
    for label, store in analysis.runs.items():
        idraws = derived(store, "infection_draws")
        present = set(idraws.columns)
        sums = ", ".join(f'sum("{c}") AS "{c}"' if c in present else f'0.0 AS "{c}"'
                         for c in ["q", *classes, *MULTI_SOURCE])
        with Scope(store.connection) as sc:
            members = sc.name(membership)
            rows = sc.rows(f"""
                WITH r AS (SELECT d.iteration, s.stratum, {sums}
                           FROM {sc.name(idraws)} d JOIN {members} s USING (infection) GROUP BY ALL),
                     o AS (SELECT d.iteration, s.stratum, sum(d.out_degree) AS out_degree
                           FROM {sc.name(derived(store, "out_degree_draws"))} d JOIN {members} s USING (infection)
                           GROUP BY ALL)
                SELECT o.iteration, o.stratum, o.out_degree, r.* EXCLUDE (iteration, stratum)
                FROM o LEFT JOIN r USING (iteration, stratum) ORDER BY o.iteration, o.stratum""")
        result[label] = [(iteration, k, dict(zip(names, (0.0 if v is None else v / n[k] for v in values))))
                         for iteration, k, *values in rows]
    return result, names


def map_and_size_range(analysis, groups, level=LEVEL):
    """Per stratum, (psi_MAP, (psi_lo, psi_hi)): the mean over its Infections of the number of children whose modal
    set names them, and of their own k_lo and k_hi, the min and max observed-parent count in their credible
    collection at `level`. From the Parent sets pooled over every Draw of every Chain (`parent_sets.pooled_sets`);
    None where no Chain has Draws."""
    sets = pooled_sets(analysis, level)
    if not sets:
        return [(None, None) for _ in groups]
    children = {}
    for ranked in sets.values():
        for parent in ranked[0].parents:
            children[parent] = children.get(parent, 0) + 1
    counts = {i: [s.n_observed for s in ranked if s.credible] for i, ranked in sets.items()}
    result = []
    for infections in groups:
        n = len(infections)
        result.append((sum(children.get(i, 0) for i in infections) / n,
                       (sum(min(counts[i]) for i in infections) / n, sum(max(counts[i]) for i in infections) / n)))
    return result


def stratified_functionals(analysis, column=None, epi=None, level=LEVEL):
    """[Stratum] for each value of epi column `column` (see `strata`), or the one stratum of every Infection."""
    found = strata(analysis, column, epi)
    groups = [infections for _, infections in found]
    draws, names = per_draw(analysis, groups)
    points = map_and_size_range(analysis, groups, level)
    pooled = [{name: [] for name in names} for _ in found]
    for rows in draws.values():
        for _, k, values in rows:
            for name, value in values.items():
                pooled[k][name].append(value)
    result = []
    for k, (value, infections) in enumerate(found):
        p = {name: interval(values) for name, values in pooled[k].items()}
        result.append(Stratum(
            value=value, n=len(infections), psi_map=points[k][0], size_range=points[k][1], psi_edge=p["psi_edge"],
            mean_q=p["mean_q"],
            classes={name.removeprefix("mean_"): p[name] for name in names[2:-2]},
            multi_source_augmented=p["mean_multi_source_augmented"],
            multi_source_observed=p["mean_multi_source_observed"]))
    return result


def stratum_series(analysis, name, column=None, value=None, epi=None):
    """Per-Draw functional `name` of one stratum in every Chain: {chain: {"iteration": [...], "value": [...]}} in
    Sampling-iteration order. The stratum is the Infections with `value` of epi column `column` (None: those with no
    value), or every Infection without a column."""
    if column is None and value is not None:
        raise StrataError(f"stratum '{value}' needs a column")
    found = strata(analysis, column, epi)
    groups = [infections for v, infections in found if v == value]
    if not groups:
        raise StrataError(f"no Infection has {column} '{value}'")
    draws, names = per_draw(analysis, groups)
    if name not in names:
        raise StrataError(f"no per-Draw functional '{name}'; they are {', '.join(names)}")
    return {chain: {"iteration": [iteration for iteration, _, _ in rows], "value": [v[name] for _, _, v in rows]}
            for chain, rows in draws.items()}
