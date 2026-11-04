"""The fit summary (U7): one row per Analysis with exactly the 70 columns of the thesis pipeline's per-fit summary,
`simulation_fit_summary.rds` (`remote/R/final_analysis/R`), so the manuscript's R table and figure scripts read the
rerun's rows unchanged. Character columns are strings, counts int64, the rest double.

- Simulation factors (`run_key` to `fit_type`) come from the manifest's `factors` object, which the rerun loader
  (`pt_analysis.rerun`) fills from the simulation's path; null where the manifest has none.
- Realised truth, edge ranking and parent-set columns are `metrics.fit_summary_metrics`; null without truth.
- Functionals are per-Draw means over the Analysis's Infections (`strata.per_draw` over the one stratum of every
  Infection, and `infection_draws` for entropy and effective size), pooled over every Draw of every Chain: the
  estimate is their mean, `_lower` and `_upper` their 2.5% and 97.5% type-7 quantiles (`quantiles.interval`).
  `average_outdegree` is psi_edge; `unobserved_source` is mean q; `observed_only`, `mixed_observed_unobserved` and
  `unobserved_only` sum the composition classes without the unobserved source, with it, and with no observed parent.
- `mean_no_parent_probability` and its interval are 0: it is the weight of the empty set, which a store never has.
- Errors are estimate minus truth, with its absolute value and square.

`edge_rows` gives the per-edge rows behind the cohort view's calibration curve: each scored or true ordered pair of
Infections with its Edge marginal and whether it is a true edge.
"""

import math

import pyarrow as pa
import pyarrow.parquet as pq

from pt_analysis import metrics, strata
from pt_analysis.quantiles import interval
from pt_analysis.derive import derived
from pt_analysis.parent_sets import LEVEL

FACTORS = [
    ("run_key", pa.string()), ("simulation_model", pa.string()), ("total_nodes", pa.float64()),
    ("founder_probability", pa.float64()), ("superinfection_probability", pa.float64()),
    ("secondary_case_rate", pa.float64()), ("founder_coi_rate", pa.float64()),
    ("asymptomatic_probability", pa.float64()), ("panel_name", pa.string()), ("panel_alleles_per_locus", pa.float64()),
    ("simulation_index", pa.int64()), ("seed", pa.int64()), ("missingness_condition", pa.string()),
    ("fit_type", pa.string()),
]
TRUTH = [("n_true_cases", pa.int64()), ("n_observed_cases", pa.int64()), ("true_average_outdegree", pa.float64()),
         ("true_unobserved_source_proportion", pa.float64())]
# (estimate column, interval prefix) of each functional, in column order.
FUNCTIONALS = [
    ("estimated_average_outdegree", "average_outdegree"),
    ("estimated_unobserved_source_proportion", "unobserved_source"),
    *((name, name) for name in (
        "mean_observed_only_probability", "mean_mixed_observed_unobserved_probability",
        "mean_unobserved_only_probability", "mean_no_parent_probability", "mean_multi_source_probability_augmented",
        "mean_multi_source_probability_observed_only", "mean_parent_set_entropy", "mean_effective_parent_sets")),
]
RANKING = [(f"{kind}_{name}", t) for kind in ("directed", "undirected")
           for name, t in (("n_edge_universe", pa.int64()), ("n_true_edges", pa.int64()),
                           ("n_distinct_scores", pa.int64()), ("average_precision", pa.float64()),
                           ("pr_auc", pa.float64()), ("top_k_precision", pa.float64()))]
PARENT_SETS = ["mean_credible_parent_set_size_80", "median_credible_parent_set_size_80",
               "credible_parent_set_coverage_80", "mean_top_parent_set_probability"]
# (error column prefix, estimate, truth)
ERRORS = [("average_outdegree", "estimated_average_outdegree", "true_average_outdegree"),
          ("unobserved_source", "estimated_unobserved_source_proportion", "true_unobserved_source_proportion")]

SCHEMA = pa.schema([
    *FACTORS, *TRUTH,
    *((c, pa.float64()) for estimate, prefix in FUNCTIONALS
      for c in (estimate, f"{prefix}_lower", f"{prefix}_upper")),
    *RANKING, *((c, pa.float64()) for c in PARENT_SETS),
    *((f"{prefix}_{kind}", pa.float64()) for prefix, _, _ in ERRORS for kind in ("error", "abs_error", "sq_error")),
])
# Columns that identify a fit within a cohort: a simulation's full, partial and null fits share a run_key.
FIT_KEY = ("run_key", "missingness_condition", "fit_type")
EDGE_SCHEMA = pa.schema([*((k, pa.string()) for k in FIT_KEY), ("parent", pa.string()), ("child", pa.string()),
                         ("score", pa.float64()), ("truth", pa.bool_())])
# Parquet key-value metadata of a row file: the Analysis manifest it summarises, for the cohort's `manifest` column.
MANIFEST_KEY = b"pt_analysis.manifest"


def per_draw_functionals(analysis):
    """{functional estimate column: [per-Draw value over every Draw of every Chain]}, Chains in manifest order."""
    found = strata.strata(analysis)
    infections = found[0][1]
    n = len(infections)
    draws, names = strata.per_draw(analysis, [infections])
    classes = [name for name in names if name.startswith("mean_observed_") or name == "mean_unobserved_only"]
    values = {estimate: [] for estimate, _ in FUNCTIONALS}
    for label, store in analysis.runs.items():
        entropy = {it: (h / n, e / n) for it, h, e in derived(store, "infection_draws").aggregate(
            "iteration, sum(entropy), sum(effective_size)", "iteration").fetchall()}
        for iteration, _, v in draws[label]:
            observed = [c for c in classes if c != "mean_unobserved_only"]
            h, e = entropy.get(iteration, (0.0, 0.0))
            for column, value in (
                    ("estimated_average_outdegree", v["psi_edge"]),
                    ("estimated_unobserved_source_proportion", v["mean_q"]),
                    ("mean_observed_only_probability", sum(v[c] for c in observed if not c.endswith("_unobserved"))),
                    ("mean_mixed_observed_unobserved_probability",
                     sum(v[c] for c in observed if c.endswith("_unobserved"))),
                    ("mean_unobserved_only_probability", v.get("mean_unobserved_only", 0.0)),
                    ("mean_no_parent_probability", 0.0),
                    ("mean_multi_source_probability_augmented", v["mean_multi_source_augmented"]),
                    ("mean_multi_source_probability_observed_only", v["mean_multi_source_observed"]),
                    ("mean_parent_set_entropy", h), ("mean_effective_parent_sets", e)):
                values[column].append(value)
    return values


def factor_values(analysis):
    """The manifest's simulation factors as the schema's types; null for a factor it does not give."""
    given = analysis.factors or {}
    result = {}
    for name, t in FACTORS:
        value = given.get(name)
        if value is not None:
            value = str(value) if t == pa.string() else int(value) if t == pa.int64() else float(value)
        result[name] = value
    return result


def summarize(analysis, level=LEVEL):
    """The Analysis's fit-summary row: {column: value} in schema order."""
    row = dict.fromkeys(SCHEMA.names)
    row.update(factor_values(analysis))
    for (estimate, prefix), (_, values) in zip(FUNCTIONALS, per_draw_functionals(analysis).items(), strict=True):
        found = interval(values)
        if found is not None:
            row[estimate], row[f"{prefix}_lower"], row[f"{prefix}_upper"] = found.mean, found.lower, found.upper
    truth = metrics.analysis_truth(analysis)
    if truth is not None:
        row.update({k: v for k, v in metrics.fit_summary_metrics(analysis, truth, level).items() if k in row})
        for prefix, estimate, true in ERRORS:
            if row[estimate] is not None and row[true] is not None:
                error = row[estimate] - row[true]
                row[f"{prefix}_error"], row[f"{prefix}_abs_error"], row[f"{prefix}_sq_error"] = \
                    error, abs(error), error * error
    else:
        sets = metrics.parent_set_summary(metrics.parent_set_rows(analysis, None, level), level)
        row.update({k: v for k, v in sets.items() if k in row})
    return row


def table(row, manifest=None):
    """One fit-summary row as an Arrow table with `SCHEMA`, recording the manifest's path in its metadata."""
    values = {name: [None if isinstance(v, float) and math.isnan(v) else v]
              for name, v in zip(SCHEMA.names, (row[n] for n in SCHEMA.names))}
    result = pa.table(values, schema=SCHEMA)
    return result.replace_schema_metadata({MANIFEST_KEY: str(manifest)}) if manifest else result


def write_row(row, path, manifest=None):
    pq.write_table(table(row, manifest), path)


def edge_rows(analysis, truth=None):
    """Per-edge rows of the Analysis for the calibration curve, as a table with `EDGE_SCHEMA`: every ordered pair of
    distinct Infections with a positive Edge marginal or a true edge between them, with its score (the pooled Edge
    marginal clipped to [0, 1] and rounded as the ranking scores them, 0 when unscored) and whether it is a true
    edge. Pairs neither scored nor true are left out: their count is the row's `directed_n_edge_universe` less these
    rows. Without truth, every `truth` is false."""
    truth = truth if truth is not None else metrics.analysis_truth(analysis)
    ids = set(metrics.infections(analysis))
    scores = {k: metrics.normalise(s) for k, s in metrics.edge_scores(analysis).items()
              if k[0] in ids and k[1] in ids and k[0] != k[1]}
    true = set(truth.edges) if truth is not None else set()
    keys = sorted({k for k, s in scores.items() if s > 0} | true)
    factors = factor_values(analysis)
    return pa.table({**{k: [factors[k]] * len(keys) for k in FIT_KEY},
                     "parent": [p for p, _ in keys], "child": [c for _, c in keys],
                     "score": [scores.get(k, 0.0) for k in keys], "truth": [k in true for k in keys]},
                    schema=EDGE_SCHEMA)


def write_edges(analysis, path):
    """Write the Analysis's per-edge rows (`edge_rows`) to `path`; return their number."""
    rows = edge_rows(analysis)
    pq.write_table(rows, path)
    return rows.num_rows
