"""Flows between strata of an epi column, lags and spatial scale (U9), pooled over an Analysis's Chains.

Every quantity here is a function of the pooled Edge marginals: Q[i, j], the weight of child j's Parent sets that hold
observed parent i, summed over every Draw of every Chain and divided by the number of Draws (zero at a Draw where the
store has no row), and q_j, the same for the unobserved source (j's Latent parent). Sums of Q are expected edge
counts. Strata are those of `pt_analysis.strata`: one per value of the column, plus one of the Infections with no
value. Definitions follow the manuscript's Zanzibar results (`scripts/summarize_zanzibar_functionals.R`):

- Attribution: per recipient stratum, the attributed parent mass of its Infections (expected observed parents plus
  q) split by source stratum, with the unobserved source as one more source, normalised within the recipient.
- Enrichment: per source x recipient cell, the cell's share of all edge mass over its share of the candidate pairs
  (the input's `allowed_parents`); no value where the cell has no candidate pair.
- Net flow: per stratum, expected out-edges minus in-edges, counting only edges between different strata.
- Lag: child minus parent observation time, from the input (`pt_analysis.times`). Reassigning beyond a lag moves the
  mass of every edge with lag above the cutoff to the unobserved source, Draw by Draw: psi_edge drops by the moved
  mass over n and mean q rises by it, the manuscript's accounting, except that each Infection's q is clamped at 1 at
  every Draw that gains moved mass. q is a probability; unclamped, an Infection whose sets hold both a moved parent
  and its Latent parent would count that weight twice and could pass 1. So mean q rises by at most the moved mass
  over n. Per stratum,
  psi_edge counts the stratum's Infections as sources and q as recipients, as in `pt_analysis.strata`.
- Spatial scale: edge mass by where the source lives relative to the recipient, for a unit column nested in a parent
  column: `same_unit` (same unit and same parent unit, so a unit name reused in two parent units is two units),
  `same_parent` (same parent unit, another unit), `other` (anything else, a missing value included).
"""

import math
import statistics
from dataclasses import dataclass

import pyarrow as pa

from pt_analysis.derive import derived, table
from pt_analysis.epi import EpiError, load_epi
from pt_analysis.inputs import nodes
from pt_analysis.quantiles import Interval, interval
from pt_analysis.sql import Scope
from pt_analysis.strata import strata
from pt_analysis import times

BINS = 30
RELATIONS = ("same_unit", "same_parent", "other")


class FlowsError(Exception):
    """A lag cutoff, bin width or pair of location columns that cannot be used."""


@dataclass(frozen=True)
class AttributionCell:
    """Of recipient stratum `recipient`'s attributed parent mass, the part from source stratum `source`, or from the
    unobserved source when `unobserved` (then `source` is None)."""

    recipient: str | None
    source: str | None
    unobserved: bool
    mass: float
    share: float | None


@dataclass(frozen=True)
class EnrichmentCell:
    source: str | None
    recipient: str | None
    expected_edges: float
    share: float | None  # of all edge mass
    candidates: int
    candidate_share: float | None  # of all candidate pairs
    enrichment: float | None  # share / candidate_share; None without candidate pairs


@dataclass(frozen=True)
class NetFlow:
    value: str | None
    n: int
    out_edges: float  # expected edges from the stratum to another stratum
    in_edges: float  # expected edges into the stratum from another stratum
    net: float


@dataclass(frozen=True)
class LagBin:
    x0: float
    x1: float
    mass: float  # expected edges with lag in [x0, x1)


@dataclass(frozen=True)
class Reassigned:
    """psi_edge and mean q of a stratum before and after reassigning mass beyond the lag, each a mean over Draws with
    its 95% MCMC interval."""

    value: str | None
    n: int
    psi_edge: Interval | None
    psi_edge_moved: Interval | None
    mean_q: Interval | None
    mean_q_moved: Interval | None


@dataclass(frozen=True)
class Reassignment:
    max_lag: float
    moved_mass: float  # expected edges with lag above max_lag
    moved_share: float | None  # of all edge mass
    overall: Reassigned
    strata: list[Reassigned]


@dataclass(frozen=True)
class SpatialShare:
    relation: str
    mass: float
    share: float | None


def _draws(analysis):
    return sum(store.draws() for store in analysis.runs.values())


def _marked_edges(store, s):
    """The Run's derived `edges` with `latent`: whether the parent is the child's Latent parent; named in scope `s`."""
    return s.name(s.sql(f"""
        SELECT e.iteration, e.child, e.parent, e.prob, e.parent = l.latent_parent AS latent
        FROM {s.name(derived(store, "edges"))} e JOIN {s.name(store.structure("infections"))} l
             ON l.infection = e.child"""))


def pooled(analysis):
    """({(parent, child): Q} for observed parents, {infection: q}), pooled over every Draw of every Chain."""
    total = _draws(analysis)
    edges, q = {}, {}
    if not total:
        return edges, q
    for store in analysis.runs.values():
        with Scope(store.connection) as s:
            rows = s.rows(f"SELECT child, parent, latent, sum(prob) FROM {_marked_edges(store, s)} GROUP BY ALL")
        for child, parent, latent, mass in rows:
            if latent:
                q[child] = q.get(child, 0.0) + mass / total
            else:
                edges[(parent, child)] = edges.get((parent, child), 0.0) + mass / total
    return edges, q


def lags(analysis):
    """lag(parent, child): the child's minus the parent's observation time (`times.lag`), None without both."""
    observed = times.observation_times(analysis)
    return lambda parent, child: times.lag(observed, parent, child)


def candidate_pairs(analysis):
    """[(parent, child)] from the input's `allowed_parents`, both ends Infections of the input."""
    found = nodes(analysis.input)
    return [(str(p), child) for child, node in found.items() for p in node.get("allowed_parents", [])
            if str(p) in found]


def _stratum_of(found):
    return {i: value for value, infections in found for i in infections}


def attribution(analysis, column=None, epi=None):
    """[AttributionCell] per recipient stratum, sources in stratum order then the unobserved source."""
    found = strata(analysis, column, epi)
    of = _stratum_of(found)
    edges, q = pooled(analysis)
    mass = {}
    for (parent, child), m in edges.items():
        if parent in of and child in of:
            mass[(of[child], of[parent])] = mass.get((of[child], of[parent]), 0.0) + m
    cells = []
    for recipient, infections in found:
        unobserved = sum(q.get(i, 0.0) for i in infections)
        row = [(source, False, mass.get((recipient, source), 0.0)) for source, _ in found] + [(None, True, unobserved)]
        total = sum(m for _, _, m in row)
        cells += [AttributionCell(recipient, source, u, m, m / total if total else None) for source, u, m in row]
    return cells


def enrichment(analysis, column=None, epi=None):
    """[EnrichmentCell] for every source x recipient pair of strata, source-major in stratum order."""
    found = strata(analysis, column, epi)
    of = _stratum_of(found)
    edges, _ = pooled(analysis)
    mass, count = {}, {}
    for (parent, child), m in edges.items():
        if parent in of and child in of:
            mass[(of[parent], of[child])] = mass.get((of[parent], of[child]), 0.0) + m
    for parent, child in candidate_pairs(analysis):
        if parent in of and child in of:
            count[(of[parent], of[child])] = count.get((of[parent], of[child]), 0) + 1
    total_mass, total_count = sum(mass.values()), sum(count.values())
    cells = []
    for source, _ in found:
        for recipient, _ in found:
            m, c = mass.get((source, recipient), 0.0), count.get((source, recipient), 0)
            share = m / total_mass if total_mass else None
            candidate_share = c / total_count if total_count else None
            cells.append(EnrichmentCell(source, recipient, m, share, c, candidate_share,
                                        share / candidate_share if share is not None and candidate_share else None))
    return cells


def net_flow(analysis, column=None, epi=None, min_infections=1):
    """([NetFlow] for strata with at least `min_infections` Infections, in stratum order; the number hidden)."""
    found = strata(analysis, column, epi)
    of = _stratum_of(found)
    edges, _ = pooled(analysis)
    out, into = {}, {}
    for (parent, child), m in edges.items():
        if parent in of and child in of and of[parent] != of[child]:
            out[of[parent]] = out.get(of[parent], 0.0) + m
            into[of[child]] = into.get(of[child], 0.0) + m
    shown = [NetFlow(v, len(i), out.get(v, 0.0), into.get(v, 0.0), out.get(v, 0.0) - into.get(v, 0.0))
             for v, i in found if len(i) >= min_infections]
    return shown, len(found) - len(shown)


def nice_width(span, bins=BINS):
    """A bin width of 1, 2 or 5 times a power of ten giving about `bins` bins over `span`."""
    if span <= 0:
        return 1.0
    raw = span / bins
    power = 10 ** math.floor(math.log10(raw))
    return next(power * f for f in (1, 2, 5, 10) if power * f >= raw)


def lag_histogram(analysis, width=None):
    """[LagBin] of edge mass by lag, contiguous from the bin holding the smallest lag to the one holding the
    largest; bins start at multiples of `width` (default `nice_width` of the lag range)."""
    edges, _ = pooled(analysis)
    lag = lags(analysis)
    pairs = [(lag(p, c), m) for (p, c), m in edges.items() if lag(p, c) is not None]
    if not pairs:
        return []
    lo, hi = min(x for x, _ in pairs), max(x for x, _ in pairs)
    if width is None:
        width = nice_width(hi - lo)
    if not width > 0:
        raise FlowsError(f"a bin width is positive, not {width}")
    first = math.floor(lo / width)
    mass = [0.0] * (math.floor(hi / width) - first + 1)
    for x, m in pairs:
        mass[math.floor(x / width) - first] += m
    return [LagBin((first + k) * width, (first + k + 1) * width, m) for k, m in enumerate(mass)]


def median_duration(analysis):
    """The median over Infections of each one's posterior mean infection duration (the `infections` Trace) over
    every Draw of every Chain; None without Draws."""
    total = _draws(analysis)
    if not total:
        return None
    sums = {}
    for store in analysis.runs.values():
        for infection, s in table(store, "infections").aggregate("infection, sum(infection_duration)",
                                                                 "infection").fetchall():
            sums[infection] = sums.get(infection, 0.0) + s
    return statistics.median(s / total for s in sums.values()) if sums else None


def _reassigned_draws(analysis, groups, max_lag):
    """Per stratum, per Draw of every Chain: (psi_edge, psi_edge moved, mean q, mean q moved), q moved clamped at 1
    per Infection."""
    lag = lags(analysis)
    edges, _ = pooled(analysis)
    beyond = [(p, c) for (p, c) in edges if (x := lag(p, c)) is not None and x > max_lag]
    moved = pa.table({"parent": [p for p, _ in beyond], "child": [c for _, c in beyond]},
                     schema=pa.schema([("parent", pa.string()), ("child", pa.string())]))
    membership = pa.table({"infection": [i for infections in groups for i in infections],
                           "stratum": [k for k, infections in enumerate(groups) for _ in infections]})
    n = [len(infections) for infections in groups]
    result = [[] for _ in groups]
    for store in analysis.runs.values():
        with Scope(store.connection) as s:
            members = s.name(membership)
            rows = s.rows(f"""
                WITH e AS (SELECT m.*, f.child IS NOT NULL AND NOT m.latent AS beyond
                           FROM {_marked_edges(store, s)} m LEFT JOIN {s.name(moved)} f USING (parent, child)),
                     draws AS (SELECT DISTINCT iteration FROM e),
                     grid AS (SELECT d.iteration, s.stratum
                              FROM draws d CROSS JOIN (SELECT DISTINCT stratum FROM {members}) s),
                     source AS (SELECT e.iteration, s.stratum, sum(e.prob) AS out_mass,
                                       coalesce(sum(e.prob) FILTER (WHERE e.beyond), 0) AS out_moved
                                FROM e JOIN {members} s ON s.infection = e.parent WHERE NOT e.latent GROUP BY ALL),
                     child AS (SELECT iteration, child, coalesce(sum(prob) FILTER (WHERE latent), 0) AS q,
                                      coalesce(sum(prob) FILTER (WHERE beyond), 0) AS moved
                               FROM e GROUP BY ALL),
                     recipient AS (SELECT c.iteration, s.stratum, sum(c.q) AS q_mass,
                                          sum(CASE WHEN c.moved > 0 THEN least(c.q + c.moved, 1.0) ELSE c.q END)
                                              AS q_moved
                                   FROM child c JOIN {members} s ON s.infection = c.child GROUP BY ALL)
                SELECT g.stratum, coalesce(o.out_mass, 0), coalesce(o.out_moved, 0), coalesce(r.q_mass, 0),
                       coalesce(r.q_moved, 0)
                FROM grid g LEFT JOIN source o USING (iteration, stratum)
                     LEFT JOIN recipient r USING (iteration, stratum)
                ORDER BY g.iteration, g.stratum""")
        for k, out_mass, out_moved, q_mass, q_moved in rows:
            result[k].append((out_mass / n[k], (out_mass - out_moved) / n[k], q_mass / n[k], q_moved / n[k]))
    return result


def _reassigned(value, infections, draws):
    psi, psi_moved, q, q_moved = (interval([d[k] for d in draws]) for k in range(4))
    return Reassigned(value, len(infections), psi, psi_moved, q, q_moved)


def reassign(analysis, max_lag, column=None, epi=None):
    """psi_edge and mean q overall and per stratum (none without a column) with the mass of every edge of lag above
    `max_lag` moved to the unobserved source."""
    if not isinstance(max_lag, (int, float)) or math.isnan(max_lag):
        raise FlowsError(f"a lag cutoff is a number, not {max_lag!r}")
    everyone = strata(analysis)
    found = strata(analysis, column, epi) if column is not None else []
    draws = _reassigned_draws(analysis, [i for _, i in everyone] + [i for _, i in found], max_lag)
    edges, _ = pooled(analysis)
    lag = lags(analysis)
    total = sum(edges.values())
    moved = sum(m for (p, c), m in edges.items() if (x := lag(p, c)) is not None and x > max_lag)
    return Reassignment(
        max_lag=max_lag, moved_mass=moved, moved_share=moved / total if total else None,
        overall=_reassigned(None, everyone[0][1], draws[0]),
        strata=[_reassigned(v, i, d) for (v, i), d in zip(found, draws[1:])])


def spatial_scale(analysis, unit, parent, max_lag=None, epi=None):
    """[SpatialShare] for each relation in RELATIONS: the edge mass, of edges with lag at most `max_lag` (every edge
    when None), whose source is in the recipient's unit, in its parent unit only, or elsewhere."""
    if unit == parent:
        raise FlowsError(f"the unit and parent unit columns are both '{unit}'; name two nested columns")
    epi = epi or load_epi(analysis)
    if epi is None:
        raise EpiError(f"the Analysis {analysis.label} has no epi data; its manifest names no 'epi' CSV")
    units, parents = epi.values(unit), epi.values(parent)
    edges, _ = pooled(analysis)
    lag = lags(analysis)
    mass = dict.fromkeys(RELATIONS, 0.0)
    for (source, child), m in edges.items():
        if max_lag is not None and ((x := lag(source, child)) is None or x > max_lag):
            continue
        same_parent = parents.get(source) is not None and parents.get(source) == parents.get(child)
        same_unit = same_parent and units.get(source) is not None and units.get(source) == units.get(child)
        mass["same_unit" if same_unit else "same_parent" if same_parent else "other"] += m
    total = sum(mass.values())
    return [SpatialShare(r, mass[r], mass[r] / total if total else None) for r in RELATIONS]
