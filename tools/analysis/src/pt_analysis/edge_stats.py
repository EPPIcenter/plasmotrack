"""Edge statistics (U2): for every (child, parent) of an Analysis, its pooled Edge marginal, each Chain's mean, the
between-Chain spread, each Chain's switch count and the multi-chain ESS. The unobserved source is a parent like any
other, under the child's Latent parent id, and is flagged as such.

A Chain's mean of an edge is its `edges` weight summed over the Chain's Draws and divided by their number: a Draw
without a row has weight 0 there. The pooled marginal is the mean over every Draw of every Chain. The spread is the
largest Chain mean minus the smallest, a Chain with no row of the edge counting as 0. The switch count of a Chain is
the number of adjacent Draws at which the weight crossed one half: one Draw above 0.5, the other at or below it.

The statistics are the Analysis-level derived table `edge_stats` (`derive_analysis` writes it after `convergence`,
whose ESS it reuses), one row per edge and Chain with that Chain's mean and switch count, then a pooled row (chain
null) with the pooled marginal, the summed switch count, the spread and the ESS:

    edge_stats    child, parent, unobserved, chain, mean, switches, spread, ess
"""

from dataclasses import dataclass, replace

import pyarrow as pa

from pt_analysis.derive import analysis_dir, derived, derived_analysis, table
from pt_analysis.sql import Scope


@dataclass(frozen=True)
class EdgeStats:
    """An edge over the Chains of an Analysis. `unobserved` says whether the parent is the unobserved source (the
    child's Latent parent). `switches` sums the Chains' switch counts; `chains` maps each Chain with Draws to its
    ChainEdge, ChainEdge(0.0, 0) where the Chain has no row of the edge. `ess` is the pooled multi-chain ESS from the
    convergence table, None where it has none."""

    child: str
    parent: str
    unobserved: bool
    pooled: float
    spread: float
    switches: int
    ess: float | None
    chains: dict

    @property
    def group(self):
        """The convergence table's group of this edge's trace quantity: an edge from the unobserved source is the
        child's unobserved-source inclusion q."""
        return "unobserved-source inclusion" if self.unobserved else "edge"

    @property
    def name(self):
        """The convergence table's name of this edge's trace quantity: `child<-parent`, or the child for q."""
        return self.child if self.unobserved else f"{self.child}<-{self.parent}"


@dataclass(frozen=True)
class ChainEdge:
    """One Chain's view of an edge: its mean weight over the Chain's Draws and its switch count."""

    mean: float
    switches: int


def chain_edge_stats(connection, edges, draws):
    """{(child, parent): ChainEdge} for one Chain, from its `edges` relation (iteration, child, parent, prob) on
    `connection` and the Sampling iterations of its Draws. Computed in DuckDB over every Draw of every edge with a
    row in the Chain."""
    if not draws:
        return {}
    with Scope(connection) as s:
        edges = s.name(edges)
        rows = s.rows(f"""
            WITH pairs AS (SELECT DISTINCT child, parent FROM {edges}),
                 dense AS (SELECT child, parent, iteration, coalesce(e.prob, 0) AS prob
                           FROM pairs CROSS JOIN {s.name(pa.table({"iteration": pa.array(draws, pa.int64())}))}
                                LEFT JOIN {edges} e USING (child, parent, iteration)),
                 flagged AS (SELECT child, parent, prob, prob > 0.5 AS above,
                                    lag(prob > 0.5) OVER (PARTITION BY child, parent ORDER BY iteration) AS before
                             FROM dense)
            SELECT child, parent, avg(prob), count(*) FILTER (WHERE above <> before)
            FROM flagged GROUP BY ALL""")
    return {(child, parent): ChainEdge(mean, switches) for child, parent, mean, switches in rows}


def edge_stats(per_chain, draws, latent_parents, ess):
    """Every edge of any Chain, in (child, parent) order, as EdgeStats.

    `per_chain` maps each Chain's label to its `chain_edge_stats`, `draws` each label to its number of Draws,
    `latent_parents` is the set of Latent parent ids, and `ess` maps a trace quantity (group, name) to its pooled ESS.
    Chains without Draws are left out."""
    chains = [label for label in per_chain if draws.get(label)]
    total = sum(draws[label] for label in chains)
    keys = sorted({key for label in chains for key in per_chain[label]})
    out = []
    for child, parent in keys:
        views = {label: per_chain[label].get((child, parent), ChainEdge(0.0, 0)) for label in chains}
        means = [v.mean for v in views.values()]
        stats = EdgeStats(child, parent, parent in latent_parents,
                          pooled=sum(views[label].mean * draws[label] for label in chains) / total,
                          spread=max(means) - min(means), switches=sum(v.switches for v in views.values()),
                          ess=None, chains=views)
        out.append(replace(stats, ess=ess.get((stats.group, stats.name))))
    return out


SCHEMA = pa.schema([("child", pa.string()), ("parent", pa.string()), ("unobserved", pa.bool_()),
                    ("chain", pa.string()), ("mean", pa.float64()), ("switches", pa.int64()),
                    ("spread", pa.float64()), ("ess", pa.float64())])


def compute(analysis):
    """The `edge_stats` table of an Analysis, as a pyarrow Table. Reads each Run's `edges` and the Analysis's
    `convergence.parquet`, which `derive_analysis` writes first."""
    per_chain, draws, latent_parents = {}, {}, set()
    for label, store in analysis.runs.items():
        iterations = [i for (i,) in table(store, "run_scalars").aggregate("iteration").fetchall()]
        per_chain[label] = chain_edge_stats(store.connection, derived(store, "edges"), iterations)
        draws[label] = len(iterations)
        latent_parents |= {p for (p,) in store.structure("infections").project("latent_parent").fetchall()}
    store = next(iter(analysis.runs.values()))
    convergence = store.connection.read_parquet(str(analysis_dir(analysis) / "convergence.parquet"))
    ess = {(group, name): e for group, name, e in convergence.filter("chain IS NULL").project('"group", name, ess')
           .fetchall()}
    rows = []
    for e in edge_stats(per_chain, draws, latent_parents, ess):
        edge = {"child": e.child, "parent": e.parent, "unobserved": e.unobserved}
        rows += [{**edge, "chain": label, "mean": c.mean, "switches": c.switches, "spread": None, "ess": None}
                 for label, c in e.chains.items()]
        rows.append({**edge, "chain": None, "mean": e.pooled, "switches": e.switches, "spread": e.spread,
                     "ess": e.ess})
    return pa.Table.from_pylist(rows, schema=SCHEMA)


def edge_table(analysis):
    """Every edge of `analysis` as EdgeStats, in (child, parent) order, from the derived `edge_stats` table (deriving
    first when it is not current)."""
    rows = derived_analysis(analysis, "edge_stats").order("child, parent").fetchall()
    out = {}
    for child, parent, unobserved, chain, mean, switches, spread, ess in rows:
        key = (child, parent)
        if key not in out:
            out[key] = {"chains": {}}
        if chain is None:
            out[key].update(child=child, parent=parent, unobserved=unobserved, pooled=mean, spread=spread,
                            switches=switches, ess=ess)
        else:
            out[key]["chains"][chain] = ChainEdge(mean, switches)
    order = list(analysis.runs)
    return [EdgeStats(**{**e, "chains": {label: e["chains"][label] for label in order if label in e["chains"]}})
            for e in out.values()]
