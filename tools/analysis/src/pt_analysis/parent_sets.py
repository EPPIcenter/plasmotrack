"""Parent sets ranked by pooled weight, one Run's or pooled over an Analysis's Chains, with the two rules every
summary of them shares:

- The MAP rule: an Infection's sets rank by weight descending, ties going to the members tuple (observed parents in
  id order, then the Latent parent when it is a member) in ascending order. The modal set is rank 1.
- The credible collection at a level: the smallest prefix of the ranking whose mass reaches the level, within
  `TOLERANCE` (as the manuscript's R scripts do), or every set when the mass never reaches it.

A Run's weight of a set is its weight summed over the Run's Draws divided by their number (`set_summary`). Pooled over
Chains, it is summed over every Draw of every Chain divided by the number of Draws in all; sets are matched across
Chains by their members, since set ids are each Run's own.
"""

from dataclasses import dataclass, field

from pt_analysis.derive import derived, table
from pt_analysis.sql import Scope

LEVEL = 0.8
TOLERANCE = 1e-12


def map_key(weight, members):
    """The MAP rule's sort key of a set: weight descending, then members ascending."""
    return (-weight, tuple(members))


def credible_flags(weights, level=LEVEL):
    """For weights in rank order, whether each set is in the credible collection at `level`: those before the
    prefix's mass reaches it, within `TOLERANCE`."""
    flags, mass = [], 0.0
    for w in weights:
        flags.append(mass < level - TOLERANCE)
        mass += w
    return flags


@dataclass(frozen=True)
class RankedSet:
    """One of an Infection's Parent sets in one Run with its weight and its place in the Infection's ranking."""

    infection: str
    set_id: int
    members: tuple[str, ...]  # observed parents in id order, then the Latent parent if it is a member
    n_observed: int
    composition: str  # the composition class label of `set_classes`
    weight: float
    rank: int
    cumulative_mass: float


def _rank(rows, level):
    """Rows (infection, members, weight, ...) grouped by Infection and ranked by the MAP rule: [(row, rank,
    cumulative mass, in the credible collection at `level`)] in Infection then rank order."""
    by_infection = {}
    for row in rows:
        by_infection.setdefault(row[0], []).append(row)
    out = []
    for infection in sorted(by_infection):
        ranked = sorted(by_infection[infection], key=lambda r: map_key(r[2], r[1]))
        mass = 0.0
        for rank, (row, credible) in enumerate(zip(ranked, credible_flags([r[2] for r in ranked], level)), 1):
            mass += row[2]
            out.append((row, rank, mass, credible))
    return out


def _set_rows(store, infection=None, infections=None):
    """The Run's Parent sets unranked: [(infection, members, weight, set_id, n_observed, class)], members as in
    RankedSet. With `infection`, only that Infection's sets; with `infections`, only theirs."""
    with Scope(store.connection) as s:
        rows = s.rows(f"""
            SELECT s.infection, list(m.member ORDER BY m.member = i.latent_parent, m.member) AS members, s.weight,
                   s.set_id, c.n_observed, c."class"
            FROM {s.name(derived(store, "set_summary"))} s JOIN {s.name(derived(store, "set_classes"))} c
                     USING (infection, set_id)
                 JOIN {s.name(table(store, "parent_set_members"))} m USING (infection, set_id)
                 JOIN {s.name(store.structure("infections"))} i USING (infection)
            WHERE ($infection IS NULL OR s.infection = $infection)
              AND ($all OR list_contains($infections, s.infection))
            GROUP BY ALL""", {"infection": infection, "all": infections is None,
                              "infections": sorted(infections or [])})
    return [(inf, tuple(ms), w, set_id, n, cls) for inf, ms, w, set_id, n, cls in rows]


def ranked_sets(store, level=None, infection=None):
    """The Run's Parent sets with their weights (`set_summary`), classes and members, in Infection then rank order
    by the MAP rule. With `level`, only each Infection's credible collection at that level. With `infection`, only
    that Infection's sets."""
    ranked = _rank(_set_rows(store, infection), LEVEL if level is None else level)
    return [RankedSet(inf, set_id, members, n, cls, w, rank, mass)
            for (inf, members, w, set_id, n, cls), rank, mass, credible in ranked if level is None or credible]


@dataclass(frozen=True)
class PooledSet:
    """One of an Infection's Parent sets pooled over the Chains of an Analysis."""

    infection: str
    members: tuple[str, ...]  # as in RankedSet
    n_observed: int
    composition: str
    weight: float  # summed over every Draw of every Chain, divided by the number of Draws in all
    rank: int
    cumulative_mass: float
    credible: bool  # in the credible collection at the level `pooled_sets` was given
    chains: dict = field(default_factory=dict)  # Chain label -> the set's weight in that Chain (0 where it has none)
    set_ids: dict = field(default_factory=dict)  # Chain label -> the set's id there, for Chains that have it

    @property
    def parents(self):
        """The set's observed parents, in id order."""
        return self.members[:self.n_observed]


def collection_entries(analysis):
    """[(parent, child, entry)] for every observed parent that appears in any of a child's Parent sets, pooled over
    the Chains. `entry` is the pooled mass of the sets ranked above the first one containing the parent, so the parent
    is inside the child's credible collection at every level above `entry` and outside it at every level up to it.
    The unobserved source is not an edge. Ordered by entry, then child, then parent.

    The same numbers as `pooled_sets(analysis, 1)` gives (each set's pooled weight sums its Chains in Analysis order,
    the cumulative mass adds down the ranking), computed in one query over every Chain's `set_summary` rather than
    by building every PooledSet."""
    draws = {label: store.draws() for label, store in analysis.runs.items()}
    total = sum(draws.values())
    if not total:
        return []
    labels = list(analysis.runs)
    connection = next(iter(analysis.runs.values())).connection
    with Scope(connection) as s:
        def union(relation_of):
            return " UNION ALL ".join(f"SELECT {i} AS chain, * FROM {s.name(relation_of(store))}"
                                      for i, store in enumerate(analysis.runs.values()))
        # Each Chain's weight times its Draws, added in Analysis order as `pooled_sets` adds them.
        pooled = " + ".join(f"coalesce(max(weight) FILTER (WHERE chain = {i}), 0) * {draws[label]}"
                            for i, label in enumerate(labels))
        members = union(lambda st: table(st, "parent_set_members"))
        latent = union(lambda st: st.structure("infections").project("infection, latent_parent"))
        summary = union(lambda st: derived(st, "set_summary").project("infection, set_id, weight"))
        # A set's members joined by a separator below every character of an id, so the keys sort as the members
        # tuples do (the MAP rule's tie-break) and match a set across Chains.
        rows = s.rows(f"""
            WITH m AS MATERIALIZED (SELECT m.*, m.member = l.latent_parent AS latent
                       FROM ({members}) m JOIN ({latent}) l USING (chain, infection)),
            x AS ({summary}),
            named AS MATERIALIZED (  -- observed parents in id order, then the Latent parent when it is a member
                SELECT chain, infection, set_id,
                       concat_ws(chr(1), array_to_string(list_sort(list(member) FILTER (WHERE NOT latent)), chr(1)),
                                 any_value(member) FILTER (WHERE latent)) AS key
                FROM m SEMI JOIN x USING (chain, infection, set_id) GROUP BY chain, infection, set_id),
            sets AS (
                SELECT infection, key, ({pooled}) / {total} AS weight
                FROM x JOIN named USING (chain, infection, set_id) GROUP BY infection, key),
            ranked AS MATERIALIZED (
                SELECT infection, key, weight,
                       row_number() OVER (PARTITION BY infection ORDER BY weight DESC, key) AS rank
                FROM sets),
            firsts AS (  -- the rank of the first set containing each observed parent
                SELECT infection, m.member AS parent, min(rank) AS rank
                FROM named JOIN ranked USING (infection, key) JOIN m USING (chain, infection, set_id)
                WHERE NOT m.latent GROUP BY ALL)
            SELECT w.infection, w.weights, f.firsts
            FROM (SELECT infection, list(weight ORDER BY rank) AS weights FROM ranked GROUP BY infection) w
                 LEFT JOIN (SELECT infection, list((parent, rank)) AS firsts FROM firsts GROUP BY infection) f
                 USING (infection)""")
    out = []
    for child, weights, firsts in rows:
        mass, entries = 0.0, []
        for w in weights:  # the entry of the set at each rank: the mass down to it, less its weight
            mass += w
            entries.append(mass - w)
        out.extend((parent, child, entries[rank - 1]) for parent, rank in firsts or [])
    out.sort(key=lambda r: (r[2], r[1], r[0]))
    return out


def pooled_sets(analysis, level=LEVEL, infection=None, infections=None):
    """{Infection: [PooledSet] in rank order by the MAP rule}, pooled over every Draw of every Chain, with each set's
    membership of the credible collection at `level`. With `infection`, only that Infection; with `infections`, only
    those. {} when no Chain has Draws."""
    draws = {label: store.draws() for label, store in analysis.runs.items()}
    total = sum(draws.values())
    if not total:
        return {}
    found = {}  # Chains in Analysis order, so each set's pooled weight sums its Chains in that order
    for label, store in analysis.runs.items():
        for inf, members, weight, set_id, n_observed, composition in _set_rows(store, infection, infections):
            entry = found.setdefault((inf, members), {"set": (n_observed, composition), "chains": {}, "set_ids": {}})
            entry["chains"][label] = weight
            entry["set_ids"][label] = set_id
    rows = [(inf, members, sum(e["chains"][c] * draws[c] for c in e["chains"]) / total, e)
            for (inf, members), e in found.items()]
    out = {}
    for (inf, members, weight, e), rank, mass, credible in _rank(rows, level):
        out.setdefault(inf, []).append(PooledSet(
            inf, members, *e["set"], weight, rank, mass, credible,
            chains={label: e["chains"].get(label, 0.0) for label in analysis.runs}, set_ids=e["set_ids"]))
    return out
