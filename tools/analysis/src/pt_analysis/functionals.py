"""Transmission-network summaries of one Run, read from its derived tables (`pt_analysis.derive`), which are derived
first if they are not current. Pooled weights are sums over Draws divided by the number of Draws.
"""

from dataclasses import dataclass

from pt_analysis.derive import derived
from pt_analysis.parent_sets import LEVEL, TOLERANCE, RankedSet, ranked_sets  # noqa: F401  (this module's interface)
from pt_analysis.sql import Scope


@dataclass(frozen=True)
class Edge:
    """An entry of the posterior transmission matrix: the pooled probability that `parent` is in `child`'s set."""

    parent: str
    child: str
    weight: float


def credible_collection(store, infection, level=LEVEL):
    """The credible collection of `infection`'s sets at `level` (`pt_analysis.parent_sets`)."""
    return ranked_sets(store, level=level, infection=infection)


def map_graph(store):
    """The MAP graph: every Infection's modal set under the pooled weights, by the MAP rule of
    `pt_analysis.parent_sets` (ties go to the members in order)."""
    return {s.infection: s for s in ranked_sets(store) if s.rank == 1}


def psi_map(store):
    """ψ_MAP, the MAP graph's mean out-degree: the mean over Infections of the modal set's observed-parent count."""
    graph = map_graph(store)
    return sum(s.n_observed for s in graph.values()) / len(graph) if graph else None


def size_range(store, level=LEVEL):
    """The parent-set size range [ψ_lo, ψ_hi]: the min and max observed-parent count in each Infection's credible
    collection at `level`, each averaged over Infections. None when the Run has no Draws."""
    counts = {}
    for s in ranked_sets(store, level=level):
        counts.setdefault(s.infection, []).append(s.n_observed)
    if not counts:
        return None
    return (sum(min(c) for c in counts.values()) / len(counts), sum(max(c) for c in counts.values()) / len(counts))


def top_edges(store, n=None, unobserved=False):
    """Entries of the posterior transmission matrix by weight, descending, then by child and parent: the pooled Edge
    marginals. Edges from the unobserved source (a Latent parent) only with `unobserved`. All of them unless `n`."""
    with Scope(store.connection) as s:
        edges, infections = s.name(derived(store, "edges")), s.name(store.structure("infections"))
        rows = s.rows(f"""
            SELECT parent, child, sum(prob) / (SELECT count(DISTINCT iteration) FROM {edges}) AS weight
            FROM {edges}
            WHERE $unobserved OR parent NOT IN (SELECT latent_parent FROM {infections})
            GROUP BY ALL ORDER BY weight DESC, child, parent
            LIMIT $n""", {"unobserved": unobserved, "n": n})
    return [Edge(*row) for row in rows]
