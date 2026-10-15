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

# The SQL below names DuckDB relations held in local variables (`FROM idraws`); DuckDB finds them by name in the
# calling frame (a replacement scan), so those variables are used even though nothing else reads them.

import math
from dataclasses import dataclass, field

import pyarrow as pa

from pt_analysis.derive import class_labels, derived
from pt_analysis.epi import EpiError, load_epi

LEVEL = 0.8
# A prefix closes when its mass reaches the level within this tolerance, as `functionals` and the R scripts do.
TOLERANCE = 1e-12
PROBABILITIES = (0.025, 0.975)
MULTI_SOURCE = ("multi_source_augmented", "multi_source_observed")


class StrataError(Exception):
    """A per-Draw functional or a stratum value the Analysis does not have."""


@dataclass(frozen=True)
class Interval:
    """A functional's mean over Draws and its 95% MCMC interval."""

    mean: float
    lower: float
    upper: float


def quantile(ordered, p):
    """Type-7 quantile of sorted values: (n - 1) p of the way along them, interpolated linearly."""
    h = (len(ordered) - 1) * p
    below = math.floor(h)
    above = min(below + 1, len(ordered) - 1)
    return ordered[below] + (h - below) * (ordered[above] - ordered[below])


def interval(values):
    """The mean of `values` and their 2.5% and 97.5% type-7 quantiles; None without values."""
    if not values:
        return None
    ordered = sorted(values)
    return Interval(sum(values) / len(values), *(quantile(ordered, p) for p in PROBABILITIES))


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
        odraws = derived(store, "out_degree_draws")
        members = store.connection.from_arrow(membership)
        present = set(idraws.columns)
        sums = ", ".join(f'sum("{c}") AS "{c}"' if c in present else f'0.0 AS "{c}"'
                         for c in ["q", *classes, *MULTI_SOURCE])
        rows = store.connection.sql(f"""
            WITH r AS (SELECT d.iteration, s.stratum, {sums}
                       FROM idraws d JOIN members s USING (infection) GROUP BY ALL),
                 o AS (SELECT d.iteration, s.stratum, sum(d.out_degree) AS out_degree
                       FROM odraws d JOIN members s USING (infection) GROUP BY ALL)
            SELECT o.iteration, o.stratum, o.out_degree, r.* EXCLUDE (iteration, stratum)
            FROM o LEFT JOIN r USING (iteration, stratum) ORDER BY o.iteration, o.stratum""").fetchall()
        result[label] = [(iteration, k, dict(zip(names, (0.0 if v is None else v / n[k] for v in values))))
                         for iteration, k, *values in rows]
    return result, names


def pooled_sets(analysis, level=LEVEL):
    """Per Infection, from Parent-set weights pooled over every Draw of every Chain: the observed parents of its
    modal set, and the min and max observed-parent count in its credible collection at `level` (the smallest prefix
    of its sets by pooled weight with mass at least `level`, or all of them). {infection: (parents, k_lo, k_hi)}.

    Sets are matched across Chains by their members, since set ids are each Run's own; ties in weight go to the
    members in text order."""
    pooled = None
    total = sum(store.draws() for store in analysis.runs.values())
    for store in analysis.runs.values():
        draws = store.draws()
        if not draws:
            continue
        summary = derived(store, "set_summary")
        classes = derived(store, "set_classes")
        members = store.table("parent_set_members")
        latent = store.structure("infections")
        weighted = store.connection.sql(f"""
            WITH m AS (SELECT infection, set_id, list_sort(list(member)) AS members,
                              list_sort(list(member) FILTER (WHERE member <> l.latent_parent)) AS parents
                       FROM members JOIN latent l USING (infection) GROUP BY ALL)
            SELECT s.infection, m.members::VARCHAR AS key, coalesce(m.parents, []) AS parents, c.n_observed,
                   s.weight * {draws} AS mass
            FROM summary s JOIN m USING (infection, set_id) JOIN classes c USING (infection, set_id)""")
        pooled = weighted if pooled is None else pooled.union(weighted)
    if pooled is None:
        return {}
    rows = pooled.query("pooled", f"""
        WITH sets AS (SELECT infection, key, any_value(parents) AS parents, any_value(n_observed) AS n_observed,
                             sum(mass) / {total} AS weight FROM pooled GROUP BY infection, key),
             ranked AS (SELECT *, row_number() OVER w AS rank,
                               coalesce(sum(weight) OVER (w ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING), 0)
                                   AS mass_before
                        FROM sets WINDOW w AS (PARTITION BY infection ORDER BY weight DESC, key))
        SELECT infection, any_value(parents) FILTER (WHERE rank = 1),
               min(n_observed) FILTER (WHERE mass_before < {level} - {TOLERANCE}),
               max(n_observed) FILTER (WHERE mass_before < {level} - {TOLERANCE})
        FROM ranked GROUP BY infection""").fetchall()
    return {infection: (parents, lo, hi) for infection, parents, lo, hi in rows}


def map_and_size_range(analysis, groups, level=LEVEL):
    """Per stratum, (psi_MAP, (psi_lo, psi_hi)): the mean over its Infections of the number of children whose modal
    set names them, and of their own k_lo and k_hi. None where no Chain has Draws."""
    sets = pooled_sets(analysis, level)
    if not sets:
        return [(None, None) for _ in groups]
    children = {}
    for parents, _, _ in sets.values():
        for parent in parents:
            children[parent] = children.get(parent, 0) + 1
    result = []
    for infections in groups:
        n = len(infections)
        result.append((sum(children.get(i, 0) for i in infections) / n,
                       (sum(sets[i][1] for i in infections) / n, sum(sets[i][2] for i in infections) / n)))
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
            value=value, n=len(infections), psi_map=points[k][0], size_range=points[k][1], psi_edge=p["psi_edge"], mean_q=p["mean_q"],
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
