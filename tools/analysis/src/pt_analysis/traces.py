"""Per-Chain series of one quantity of the convergence table, read on demand for the trace panel, and ESS and R-hat
over a range of Draws.

A quantity is a (group, name) pair as in `derive.GROUPS` and the convergence table. Quantities of the first six
groups are read from the store's Traces (`run_scalars`, `infections`, `allele_frequencies`), so their burn-in rows
can be read too; edges, unobserved-source inclusion and functionals come from the derived tables, which hold Draws
only. An edge with no row at a Draw has weight 0 there, as in the convergence table.
"""

import re

from pt_analysis.derive import GROUPS, derived, table
from pt_analysis.diagnostics import multichain_ess, split_rhat
from pt_analysis.store import BOTH_PHASES, SAMPLING

# Groups read from the store's Traces, which have burn-in rows.
TRACE_GROUPS = GROUPS[:6]

GLOBAL_SCALARS = ("mean_coi", "mean_strains_tx", "parent_set_size_prob")
INFECTION_COLUMNS = {"false positive rate": "expected_false_positives",
                     "false negative rate": "expected_false_negatives",
                     "infection duration": "infection_duration"}
FUNCTIONALS = ("psi_edge", "mean_q", "mean_entropy")


class TraceError(Exception):
    """A quantity that is not in the convergence table's groups or is malformed for its group."""


class UnknownQuantity(TraceError):
    """A well-formed quantity with no rows in any Chain."""


def _rows(store, query, params, relations):
    """Rows (iteration, phase, value) of `query` in Sampling-iteration order, as a dict of lists, or None if empty.
    `relations` maps the names the query reads from to DuckDB relations."""
    connection = store.connection
    for alias, relation in relations.items():
        relation.create_view(alias, replace=True)
    rows = connection.execute(f"SELECT iteration, phase, value::DOUBLE FROM ({query}) ORDER BY iteration",
                              params).fetchall()
    if not rows:
        return None
    iteration, phase, value = (list(column) for column in zip(*rows))
    return {"iteration": iteration, "phase": phase, "value": value}


def series(store, group, name, phases=SAMPLING):
    """One Run's series of quantity (group, name): {"iteration", "phase", "value"} lists in Sampling-iteration order,
    over `phases` for the Trace groups and over the Draws otherwise; None if the Run has no rows of it."""
    if group not in GROUPS:
        raise TraceError(f"unknown group '{group}'; the groups are {', '.join(GROUPS)}")
    if group in TRACE_GROUPS:
        if group in ("log posterior", "global scalars"):
            relations = {"pt_trace_source": table_over(store, "run_scalars", phases)}
            if group == "log posterior":
                if name != "log posterior":
                    raise TraceError("the log posterior group has one quantity, 'log posterior'")
                query, params = "SELECT iteration, phase, prior + likelihood AS value FROM pt_trace_source", []
            else:
                if name not in GLOBAL_SCALARS:
                    raise TraceError(f"the global scalars are {', '.join(GLOBAL_SCALARS)}")
                query, params = f'SELECT iteration, phase, "{name}" AS value FROM pt_trace_source', []
        elif group in INFECTION_COLUMNS:
            relations = {"pt_trace_source": table_over(store, "infections", phases)}
            column = INFECTION_COLUMNS[group]
            query = f"SELECT iteration, phase, {column} AS value FROM pt_trace_source WHERE infection = ?"
            params = [name]
        else:
            match = re.fullmatch(r"(.+)\[(\d+)\]", name)
            if not match:
                raise TraceError(f"an Allele frequency is named <locus>[<allele>], not '{name}'")
            relations = {"pt_trace_source": table_over(store, "allele_frequencies", phases)}
            query = "SELECT iteration, phase, frequency AS value FROM pt_trace_source WHERE locus = ? AND allele = ?"
            params = [match[1], int(match[2])]
        return _rows(store, query, params, relations)

    if group == "edge":
        child, sep, parent = name.partition("<-")
        if not sep or not child or not parent:
            raise TraceError(f"an edge is named <child><-<parent>, not '{name}'")
        relations = {"pt_trace_source": derived(store, "edges"), "pt_trace_draws": table(store, "run_scalars")}
        found = _rows(store, "SELECT iteration, 'sampling' AS phase, prob AS value FROM pt_trace_source "
                             "WHERE child = ? AND parent = ?", [child, parent], relations)
        if found is None:
            return None
        query = """
            SELECT d.iteration, 'sampling' AS phase, coalesce(e.prob, 0) AS value
            FROM (SELECT DISTINCT iteration FROM pt_trace_draws) d
            LEFT JOIN (SELECT * FROM pt_trace_source WHERE child = ? AND parent = ?) e USING (iteration)"""
        return _rows(store, query, [child, parent], relations)
    if group == "unobserved-source inclusion":
        relations = {"pt_trace_source": derived(store, "infection_draws")}
        query = "SELECT iteration, 'sampling' AS phase, q AS value FROM pt_trace_source WHERE infection = ?"
        return _rows(store, query, [name], relations)
    if name not in FUNCTIONALS:
        raise TraceError(f"the functionals are {', '.join(FUNCTIONALS)}")
    relations = {"pt_trace_source": derived(store, "functional_draws")}
    return _rows(store, f"SELECT iteration, 'sampling' AS phase, {name} AS value FROM pt_trace_source", [], relations)


def table_over(store, name, phases):
    """Trace `name` over `phases`, or the empty relation `derive.table` gives when it has no rows there."""
    relation = store.table(name, phases)
    return relation if relation is not None else table(store, name).select("*, NULL::VARCHAR AS phase")


def stratum_traces(analysis, name, stratum):
    """A functional of one stratum (`pt_analysis.strata`) in every Chain, as `chain_traces` gives series. `stratum`
    is (epi column or None for every Infection, value or None for the Infections with no value)."""
    from pt_analysis.epi import EpiError
    from pt_analysis.strata import StrataError, stratum_series

    try:
        found = stratum_series(analysis, name, *stratum)
    except (EpiError, StrataError) as e:
        raise TraceError(str(e)) from None
    return {label: {**s, "phase": ["sampling"] * len(s["iteration"])} for label, s in found.items()}, False


def chain_traces(analysis, group, name, burn_in=False, stratum=None):
    """Quantity (group, name) in every Chain of `analysis`: ({label: series}, whether burn-in rows are included).
    Burn-in rows are included only when asked for and the group is read from the store's Traces. A Chain with no rows
    of an edge another Chain has is 0 at each of its Draws; with no rows of another quantity, its series is empty.
    Raises UnknownQuantity when no Chain has rows. With `stratum`, a functional of one stratum (`stratum_traces`)."""
    if stratum is not None:
        if group != "functionals":
            raise TraceError(f"only functionals have strata, not {group}")
        return stratum_traces(analysis, name, stratum)
    with_burn_in = burn_in and group in TRACE_GROUPS
    phases = BOTH_PHASES if with_burn_in else SAMPLING
    found = {label: series(store, group, name, phases) for label, store in analysis.runs.items()}
    if all(s is None for s in found.values()):
        raise UnknownQuantity(f"no Chain has rows of {group} '{name}'")
    empty = {"iteration": [], "phase": [], "value": []}
    traces = {}
    for label, s in found.items():
        if s is None and group == "edge":
            draws = analysis.runs[label].table("run_scalars")
            iterations = sorted(r[0] for r in draws.aggregate("iteration").fetchall()) if draws else []
            s = {"iteration": iterations, "phase": ["sampling"] * len(iterations), "value": [0.0] * len(iterations)}
        traces[label] = s or empty
    return traces, with_burn_in


def range_diagnostics(analysis, group, name, first=None, last=None, stratum=None):
    """ESS and R-hat of quantity (group, name) over the Draws with Sampling iteration in [first, last] (either end
    open when None): {"chains": [{"chain", "draws", "ess"}], "ess", "rhat"}, the pooled statistics over the Chains
    with Draws in the range. Statistics are None where `diagnostics` gives None. `stratum` as for `chain_traces`."""
    traces, _ = chain_traces(analysis, group, name, stratum=stratum)
    chains, pooled = [], []
    for label, s in traces.items():
        values = [v for i, v in zip(s["iteration"], s["value"])
                  if (first is None or i >= first) and (last is None or i <= last)]
        chains.append({"chain": label, "draws": len(values), "ess": multichain_ess([values]) if values else None})
        if values:
            pooled.append(values)
    return {"chains": chains, "ess": multichain_ess(pooled) if pooled else None,
            "rhat": split_rhat(pooled) if pooled else None}
