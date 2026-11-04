"""Derived tables: Parquet tables computed once from a Run's store and cached under `<run>/derived/`.

`derived/manifest.json` records the closed Windows per Phase the tables came from and the package version that wrote
them. `derive` recomputes every table when either differs from the store's closed Windows or this package's version,
and does nothing otherwise. Tables are computed from the sampling Phase, so every row comes from Draws;
`parent_set_members` is read over both Phases.

    set_classes   infection, set_id, n_observed, has_unobserved, class
    edges         iteration, child, parent, prob
    set_summary   infection, set_id, weight, rank, cumulative_mass
    infection_draws   iteration, infection, q, expected_observed_parents, entropy, effective_size, one column per
                      composition class present in the Run, multi_source_augmented, multi_source_observed
    functional_draws  iteration, psi_edge, mean_q, mean_entropy, mean_<class> per composition class
    genotype_summary  infection, locus, allele, presence
    genotype_mode     infection, locus, alleles, draws, fraction, rank
    out_degree_draws  iteration, infection, out_degree
    coi               iteration, infection, coi
    order_pairs       i, j, draws, before, switches

An Analysis's own derived tables go in `<manifest stem>.derived/` beside its manifest, with a manifest recording each
Run's directory and closed Windows; `derive_analysis` derives every Run first, then rewrites them when that changes:

    convergence   group, name, chain, ess, rhat
    edge_stats    child, parent, unobserved, chain, mean, switches, spread, ess   (see `pt_analysis.edge_stats`)
"""

# Each table is computed by a function of the store and a `Scope` (`pt_analysis.sql`) that names the relations its SQL
# reads; `derive` writes the table before the scope closes.

import json
import os
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from pt_analysis import __version__
from pt_analysis.diagnostics import multichain_ess, split_rhat
from pt_analysis.sql import Scope
from pt_analysis.store import BOTH_PHASES

# Columns the derivations read, for a table with no rows in the closed Windows (a Run still in burn-in has no
# sampling `parent_sets`).
COLUMNS = {
    "parent_sets": {"iteration": "BIGINT", "infection": "VARCHAR", "set_id": "BIGINT", "prob": "DOUBLE"},
    "parent_set_members": {"infection": "VARCHAR", "set_id": "BIGINT", "member": "VARCHAR"},
    "run_scalars": {"iteration": "BIGINT", "prior": "DOUBLE", "likelihood": "DOUBLE", "mean_coi": "DOUBLE",
                    "mean_strains_tx": "DOUBLE", "parent_set_size_prob": "DOUBLE"},
    "infections": {"iteration": "BIGINT", "infection": "VARCHAR", "expected_false_positives": "DOUBLE",
                   "expected_false_negatives": "DOUBLE", "infection_duration": "DOUBLE"},
    "allele_frequencies": {"iteration": "BIGINT", "locus": "VARCHAR", "allele": "BIGINT", "frequency": "DOUBLE"},
    "latent_genotypes": {"iteration": "BIGINT", "infection": "VARCHAR", "locus": "VARCHAR", "alleles": "VARCHAR"},
    "ordering": {"iteration": "BIGINT", "infection": "VARCHAR", "rank": "BIGINT"},
}


def table(store, name):
    """`store.table(name)`, or an empty relation with the columns the derivations read."""
    relation = store.table(name)
    if relation is None:
        columns = ", ".join(f'NULL::{t} AS "{c}"' for c, t in COLUMNS[name].items())
        relation = store.connection.sql(f"SELECT {columns} WHERE false")
    return relation


# The composition class of a Parent set: `unobserved_only` when it has no observed parent (its only member is the
# Latent parent; there is never an empty set), else `observed_<n>`, suffixed `_unobserved` when the Latent parent is a
# member too.
CLASS_LABEL = """CASE WHEN n_observed = 0 THEN 'unobserved_only'
                     ELSE 'observed_' || n_observed || CASE WHEN has_unobserved THEN '_unobserved' ELSE '' END END"""


def set_classes(store, s):
    """Every Parent set of the Run with its observed-parent count, whether the unobserved source (the Infection's
    Latent parent) is a member, and its composition class."""
    members = s.name(table(store, "parent_set_members"))
    infections = s.name(store.structure("infections"))
    return s.sql(f"""
        WITH sets AS (
            SELECT m.infection, m.set_id,
                   count(*) FILTER (WHERE m.member <> i.latent_parent) AS n_observed,
                   bool_or(m.member = i.latent_parent) AS has_unobserved
            FROM {members} m JOIN {infections} i USING (infection)
            GROUP BY ALL)
        SELECT *, {CLASS_LABEL} AS "class" FROM sets ORDER BY infection, set_id""")


def edges(store, s):
    """Per Draw, each entry of the conditional transmission matrix: the summed weight of the child's Parent sets that
    contain the parent. The unobserved source is the child's Latent parent id. A parent with no row at a Draw has
    weight 0 there: the store drops sets below the Run's floor."""
    return s.sql(f"""
        SELECT p.iteration, p.infection AS child, m.member AS parent, sum(p.prob) AS prob
        FROM {s.name(table(store, "parent_sets"))} p JOIN {s.name(table(store, "parent_set_members"))} m
             USING (infection, set_id)
        GROUP BY ALL ORDER BY iteration, child, parent""")


def set_summary(store, s):
    """Each Infection's Parent sets by pooled weight: the sum of the set's weight over Draws divided by the number of
    Draws (a set with no row at a Draw has weight 0 there), ranked from 1 by weight descending then set_id, with the
    cumulative weight down the ranking. Sets with no weight at any Draw are left out."""
    parent_sets = s.name(table(store, "parent_sets"))
    draws = store.draws()
    return s.sql(f"""
        WITH pooled AS (SELECT infection, set_id, sum(prob) / {draws} AS weight FROM {parent_sets} GROUP BY ALL),
        ranked AS (SELECT *, row_number() OVER (PARTITION BY infection ORDER BY weight DESC, set_id) AS rank
                   FROM pooled)
        SELECT *, sum(weight) OVER (PARTITION BY infection ORDER BY rank) AS cumulative_mass
        FROM ranked ORDER BY infection, rank""")


def class_labels(store):
    """The composition class labels of the Run's Parent sets, by observed-parent count, each count's sets without the
    unobserved source first, `unobserved_only` last."""
    with Scope(store.connection) as s:
        rows = s.rows(f"""
            SELECT DISTINCT "class", n_observed, has_unobserved FROM {s.name(set_classes(store, s))}
            ORDER BY n_observed = 0, n_observed, has_unobserved""")
    return [label for label, *_ in rows]


def infection_draws(store, s):
    """Per Draw and Infection, functions of its Parent-set weights: unobserved-source inclusion q (the weight of sets
    containing the Latent parent), the expected observed-parent count, the entropy of the weights and its exponential
    (the effective parent-set size), the weight of each composition class present in the Run (a column per class
    label, 0 where the Infection has no set of that class), and the multi-source probabilities: the weight of sets
    with more than one member (augmented) and with more than one observed parent (observed only)."""
    parent_sets = s.name(table(store, "parent_sets"))
    classes = s.name(set_classes(store, s))
    columns = "".join(f',\n               coalesce(sum(p.prob) FILTER (WHERE c."class" = \'{label}\'), 0) AS "{label}"'
                      for label in class_labels(store))
    return s.sql(f"""
        SELECT p.iteration, p.infection,
               coalesce(sum(p.prob) FILTER (WHERE c.has_unobserved), 0) AS q,
               sum(p.prob * c.n_observed) AS expected_observed_parents,
               coalesce(-sum(p.prob * ln(p.prob)) FILTER (WHERE p.prob > 0), 0) AS entropy,
               exp(entropy) AS effective_size{columns},
               coalesce(sum(p.prob) FILTER (WHERE c.n_observed + c.has_unobserved::INTEGER > 1), 0)
                   AS multi_source_augmented,
               coalesce(sum(p.prob) FILTER (WHERE c.n_observed > 1), 0) AS multi_source_observed
        FROM {parent_sets} p JOIN {classes} c USING (infection, set_id)
        GROUP BY ALL ORDER BY iteration, infection""")


def functional_draws(store, s):
    """Per Draw, means over Infections of `infection_draws`: the edge-sum psi_edge (the mean expected observed-parent
    count, the conditional transmission matrix's observed entries summed and divided by the number of Infections),
    mean q, mean entropy and the mean weight of each composition class, as `mean_<class label>`."""
    draws = s.name(infection_draws(store, s))
    columns = "".join(f', avg("{label}") AS "mean_{label}"' for label in class_labels(store))
    return s.sql(f"""
        SELECT iteration, avg(expected_observed_parents) AS psi_edge, avg(q) AS mean_q, avg(entropy) AS mean_entropy
               {columns}
        FROM {draws} GROUP BY iteration ORDER BY iteration""")


def out_degree_draws(store, s):
    """Per Draw and Infection, its expected number of observed children: the summed Edge marginals from it to every
    other Infection, 0 where it is no Infection's parent at that Draw. The Infection as a source, where
    `infection_draws` has it as a recipient; both average to psi_edge over Infections."""
    parent_sets = s.name(table(store, "parent_sets"))
    observed = s.name(edges(store, s))
    infections = s.name(store.structure("infections"))
    return s.sql(f"""
        WITH draws AS (SELECT DISTINCT iteration FROM {parent_sets}),
             grid AS (SELECT d.iteration, i.infection FROM draws d CROSS JOIN {infections} i)
        SELECT g.iteration, g.infection, coalesce(sum(o.prob), 0) AS out_degree
        FROM grid g LEFT JOIN {observed} o ON o.iteration = g.iteration AND o.parent = g.infection
        GROUP BY ALL ORDER BY iteration, infection""")


def genotype_summary(store, s):
    """Per Infection, Locus and Allele (every Allele of `structure/loci`), the posterior presence: the fraction of
    Draws whose Latent genotype has the Allele. Character k of a `latent_genotypes` bitstring is Allele k."""
    latent = s.name(table(store, "latent_genotypes"))
    loci = s.name(store.structure("loci"))
    draws = store.draws()
    return s.sql(f"""
        SELECT g.infection, g.locus, l.allele,
               count(*) FILTER (WHERE substr(g.alleles, l.allele::INTEGER + 1, 1) = '1') / {max(draws, 1)} AS presence
        FROM {latent} g JOIN {loci} l USING (locus)
        GROUP BY ALL ORDER BY infection, locus, allele""")


def genotype_mode(store, s):
    """Per Infection and Locus, every Latent genotype bitstring a Draw had, with the number and fraction of Draws that
    had it, ranked from 1 by Draws descending then bitstring. Rank 1 is the posterior mode and its fraction the
    agreement; the other rows let Chains be pooled exactly by summing `draws`."""
    latent = s.name(table(store, "latent_genotypes"))
    draws = store.draws()
    return s.sql(f"""
        WITH counts AS (SELECT infection, locus, alleles, count(*) AS draws FROM {latent} GROUP BY ALL)
        SELECT *, draws / {max(draws, 1)} AS fraction,
               row_number() OVER (PARTITION BY infection, locus ORDER BY draws DESC, alleles) AS rank
        FROM counts ORDER BY infection, locus, rank""")


def coi(store, s):
    """Per Draw and Infection, the COI: the largest Allele count (the number of `1` characters of the bitstring) over
    the Loci of the Infection's Latent genotype at that Draw."""
    return s.sql(f"""
        SELECT iteration, infection, max(length(replace(alleles, '0', ''))) AS coi
        FROM {s.name(table(store, "latent_genotypes"))} GROUP BY ALL ORDER BY iteration, infection""")


def order_pairs(store, s):
    """Per pair of Infections i < j (ids in text order) placed in the Ordering at the same Draws, the number of those
    Draws, of the Draws with i before j, and the switches: adjacent Draws at which the order changes. P(i before j)
    is before / draws and P(j before i) its complement; the per-Draw indicator of one pair is read from `ordering` on
    demand (`pt_analysis.pair`)."""
    ordering = s.name(table(store, "ordering"))
    return s.sql(f"""
        WITH pairs AS (SELECT x.iteration, x.infection AS i, y.infection AS j, x.rank < y.rank AS before
                       FROM {ordering} x JOIN {ordering} y USING (iteration) WHERE x.infection < y.infection),
             flagged AS (SELECT i, j, before, lag(before) OVER (PARTITION BY i, j ORDER BY iteration) AS previous
                         FROM pairs)
        SELECT i, j, count(*) AS draws, count(*) FILTER (WHERE before) AS before,
               count(*) FILTER (WHERE before <> previous) AS switches
        FROM flagged GROUP BY ALL ORDER BY i, j""")


TABLES = {"set_classes": set_classes, "edges": edges, "set_summary": set_summary, "infection_draws": infection_draws,
          "functional_draws": functional_draws, "genotype_summary": genotype_summary, "genotype_mode": genotype_mode,
          "out_degree_draws": out_degree_draws}
TABLES["coi"] = coi  # U10 infection card
TABLES["order_pairs"] = order_pairs  # U3 pair inspector


def source(store):
    """What the derived tables of `store` are computed from, as `derived/manifest.json` records it."""
    windows = {phase: [w.number for w in store.windows((phase,))] for phase in BOTH_PHASES}
    return {"windows": windows, "package_version": __version__}


def derived_dir(store):
    return Path(store.path) / "derived"


def is_current(store):
    """Whether `derived/` holds every table, computed from the store's closed Windows by this package version."""
    out = derived_dir(store)
    try:
        recorded = json.loads((out / "manifest.json").read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return False
    return recorded == source(store) and all((out / f"{name}.parquet").is_file() for name in TABLES)


def derive(store):
    """Write the derived tables of `store` under `<run>/derived/` unless they are current. True if they were written.
    Each table is written beside its final name and renamed into place, and the manifest is written last, so an
    interrupted derive leaves no manifest and is redone."""
    if is_current(store):
        return False
    out = derived_dir(store)
    out.mkdir(exist_ok=True)
    manifest = out / "manifest.json"
    manifest.unlink(missing_ok=True)
    for name, compute in TABLES.items():
        partial = out / f".{name}.parquet"
        with Scope(store.connection) as s:
            compute(store, s).write_parquet(str(partial))
        os.replace(partial, out / f"{name}.parquet")
    partial = out / ".manifest.json"
    partial.write_text(json.dumps(source(store), indent=2) + "\n")
    os.replace(partial, manifest)
    return True


def derived(store, name):
    """Derived table `name` of `store` as a DuckDB relation, deriving the tables first if they are not current."""
    derive(store)
    return store.connection.read_parquet(str(derived_dir(store) / f"{name}.parquet"))


# --- Analysis level -------------------------------------------------------------------------------------------------

# Groups of the convergence table, in display order.
GROUPS = ["log posterior", "global scalars", "false positive rate", "false negative rate", "infection duration",
          "allele frequency", "edge", "unobserved-source inclusion", "functionals"]

CONVERGENCE_SCHEMA = pa.schema([("group", pa.string()), ("name", pa.string()), ("chain", pa.string()),
                                ("ess", pa.float64()), ("rhat", pa.float64())])


def chain_series(store):
    """Every quantity the convergence table covers, for one Run: {(group, name): [value per Draw]} in
    Sampling-iteration order. Edges from observed parents are named `child<-parent` and are 0 at Draws where the
    child's sets containing the parent are all below the floor; unobserved-source inclusion is named by Infection."""
    derive(store)
    with Scope(store.connection) as s:
        return _chain_series(store, s)


def _chain_series(store, s):
    scalars = s.name(table(store, "run_scalars"))
    infections = s.name(table(store, "infections"))
    frequencies = s.name(table(store, "allele_frequencies"))
    edges = s.name(derived(store, "edges"))
    infection_draws = s.name(derived(store, "infection_draws"))
    functional_draws = s.name(derived(store, "functional_draws"))
    latent = s.name(store.structure("infections"))
    queries = {
        "log posterior": f"SELECT 'log posterior' AS name, iteration, prior + likelihood AS value FROM {scalars}",
        "global scalars": f"""
            UNPIVOT (SELECT iteration, mean_coi, mean_strains_tx, parent_set_size_prob FROM {scalars})
            ON mean_coi, mean_strains_tx, parent_set_size_prob INTO NAME name VALUE value""",
        "false positive rate": f"""
            SELECT infection AS name, iteration, expected_false_positives AS value FROM {infections}""",
        "false negative rate": f"""
            SELECT infection AS name, iteration, expected_false_negatives AS value FROM {infections}""",
        "infection duration": f"SELECT infection AS name, iteration, infection_duration AS value FROM {infections}",
        "allele frequency": f"""
            SELECT locus || '[' || allele || ']' AS name, iteration, frequency AS value FROM {frequencies}""",
        "edge": f"""
            WITH observed AS (SELECT * FROM {edges} WHERE parent NOT IN (SELECT latent_parent FROM {latent})),
                 draws AS (SELECT DISTINCT iteration FROM {scalars}),
                 pairs AS (SELECT DISTINCT child, parent FROM observed)
            SELECT child || '<-' || parent AS name, iteration, coalesce(o.prob, 0) AS value
            FROM pairs CROSS JOIN draws LEFT JOIN observed o USING (child, parent, iteration)""",
        "unobserved-source inclusion": f"SELECT infection AS name, iteration, q AS value FROM {infection_draws}",
        "functionals": f"""
            UNPIVOT (SELECT iteration, psi_edge, mean_q, mean_entropy FROM {functional_draws})
            ON psi_edge, mean_q, mean_entropy INTO NAME name VALUE value""",
    }
    series = {}
    for group, query in queries.items():
        rows = s.rows(f"SELECT name, list(value::DOUBLE ORDER BY iteration) FROM ({query}) GROUP BY name")
        series.update({(group, name): values for name, values in rows})
    return series


def convergence(analysis):
    """The convergence table of an Analysis, as a pyarrow Table: per group and quantity name (sorted), one row per
    Chain with that Chain's ESS, then a pooled row (chain null) with the multi-chain ESS and split R-hat over every
    Chain with Draws. R-hat is null on Chain rows. ESS and R-hat are null where `diagnostics` gives None, and for a
    Chain without Draws. An edge a Chain never gives weight to is 0 at its every Draw."""
    per_chain = {label: chain_series(store) for label, store in analysis.runs.items()}
    draws = {label: store.draws() for label, store in analysis.runs.items()}
    rows = []
    for group in GROUPS:
        names = sorted({name for series in per_chain.values() for g, name in series if g == group})
        for name in names:
            missing = lambda label: [0.0] * draws[label] if group == "edge" else []  # noqa: E731
            chains = {label: series.get((group, name)) or missing(label) for label, series in per_chain.items()}
            for label, values in chains.items():
                rows.append({"group": group, "name": name, "chain": label, "ess": multichain_ess([values]),
                             "rhat": None})
            pooled = [values for values in chains.values() if values]
            rows.append({"group": group, "name": name, "chain": None,
                         "ess": multichain_ess(pooled) if pooled else None,
                         "rhat": split_rhat(pooled) if pooled else None})
    return pa.Table.from_pylist(rows, schema=CONVERGENCE_SCHEMA)


ANALYSIS_TABLES = {"convergence": convergence}


def edge_stats(analysis):
    """U2's per-edge statistics over the Chains, in `pt_analysis.edge_stats`; written after `convergence`, whose ESS
    it reads."""
    from pt_analysis import edge_stats

    return edge_stats.compute(analysis)


ANALYSIS_TABLES["edge_stats"] = edge_stats


def analysis_dir(analysis):
    """Where an Analysis's derived tables go: `<manifest stem>.derived/` beside its manifest."""
    return analysis.path.parent / f"{analysis.path.stem}.derived"


def analysis_source(analysis):
    """What the Analysis-level tables are computed from: each Run's directory and closed Windows per Phase."""
    runs = {label: {"path": str(store.path), "windows": source(store)["windows"]}
            for label, store in analysis.runs.items()}
    return {"runs": runs, "package_version": __version__}


def derive_analysis(analysis):
    """Derive every Run of `analysis`, then write the Analysis-level tables under `analysis_dir` unless they are
    current: recorded from the same Runs, closed Windows and package version. True if they were written. Written the
    way `derive` writes, manifest last."""
    for store in analysis.runs.values():
        derive(store)
    out = analysis_dir(analysis)
    manifest = out / "manifest.json"
    try:
        current = json.loads(manifest.read_text()) == analysis_source(analysis)
    except (FileNotFoundError, json.JSONDecodeError):
        current = False
    if current and all((out / f"{name}.parquet").is_file() for name in ANALYSIS_TABLES):
        return False
    out.mkdir(exist_ok=True)
    manifest.unlink(missing_ok=True)
    for name, compute in ANALYSIS_TABLES.items():
        partial = out / f".{name}.parquet"
        pq.write_table(compute(analysis), partial)
        os.replace(partial, out / f"{name}.parquet")
    partial = out / ".manifest.json"
    partial.write_text(json.dumps(analysis_source(analysis), indent=2) + "\n")
    os.replace(partial, manifest)
    return True


def derived_analysis(analysis, name):
    """Analysis-level derived table `name` as a DuckDB relation, deriving first if it is not current."""
    derive_analysis(analysis)
    store = next(iter(analysis.runs.values()))
    return store.connection.read_parquet(str(analysis_dir(analysis) / f"{name}.parquet"))
