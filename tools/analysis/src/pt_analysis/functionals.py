"""Transmission-network summaries of one Run, read from its derived tables (`pt_analysis.derive`), which are derived
first if they are not current. Pooled weights are sums over Draws divided by the number of Draws.
"""

# The SQL below names DuckDB relations held in local variables (`FROM parent_sets`); DuckDB finds them by name in
# the calling frame (a replacement scan), so those variables are used even though nothing else reads them.

from dataclasses import dataclass

from pt_analysis.derive import derived, table

LEVEL = 0.8
# A prefix closes when its mass reaches the level within this tolerance, as the manuscript's R scripts do.
TOLERANCE = 1e-12


@dataclass(frozen=True)
class RankedSet:
    """One of an Infection's Parent sets with its pooled weight and its place in the Infection's ranking."""

    infection: str
    set_id: int
    members: tuple[str, ...]  # observed parents in id order, then the Latent parent if it is a member
    n_observed: int
    composition: str  # the composition class label of `set_classes`
    weight: float
    rank: int
    cumulative_mass: float


@dataclass(frozen=True)
class Edge:
    """An entry of the posterior transmission matrix: the pooled probability that `parent` is in `child`'s set."""

    parent: str
    child: str
    weight: float


def ranked_sets(store, level=None, infection=None):
    """Rows of `set_summary` with their classes and members, in Infection then rank order. With `level`, only each
    Infection's credible collection at that level: the smallest prefix by rank whose mass is at least `level`, or
    every set when the mass never reaches it. With `infection`, only that Infection's sets."""
    summary = derived(store, "set_summary")
    classes = derived(store, "set_classes")
    members = table(store, "parent_set_members")
    infections = store.structure("infections")
    rows = store.connection.execute("""
        WITH prefixed AS (
            SELECT *, coalesce(lag(cumulative_mass) OVER (PARTITION BY infection ORDER BY rank), 0) AS mass_before
            FROM summary)
        SELECT s.infection, s.set_id, list(m.member ORDER BY m.member = i.latent_parent, m.member) AS members,
               c.n_observed, c."class", s.weight, s.rank, s.cumulative_mass
        FROM prefixed s JOIN classes c USING (infection, set_id)
             JOIN members m USING (infection, set_id) JOIN infections i USING (infection)
        WHERE ($level IS NULL OR s.mass_before < $level - $tolerance)
          AND ($infection IS NULL OR s.infection = $infection)
        GROUP BY ALL ORDER BY s.infection, s.rank""",
        {"level": level, "tolerance": TOLERANCE, "infection": infection}).fetchall()
    return [RankedSet(inf, set_id, tuple(ms), n, cls, w, rank, cum) for inf, set_id, ms, n, cls, w, rank, cum in rows]


def credible_collection(store, infection, level=LEVEL):
    """The smallest prefix of `infection`'s sets by pooled weight whose mass is at least `level`."""
    return ranked_sets(store, level=level, infection=infection)


def map_graph(store):
    """The MAP graph: every Infection's modal set under the pooled weights (ties go to the lower set_id)."""
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
    edges = derived(store, "edges")
    infections = store.structure("infections")
    rows = store.connection.execute("""
        SELECT parent, child, sum(prob) / (SELECT count(DISTINCT iteration) FROM edges) AS weight
        FROM edges
        WHERE $unobserved OR parent NOT IN (SELECT latent_parent FROM infections)
        GROUP BY ALL ORDER BY weight DESC, child, parent
        LIMIT $n""", {"unobserved": unobserved, "n": n}).fetchall()
    return [Edge(*row) for row in rows]
