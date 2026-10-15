"""Ranking and parent-set metrics against truth (U6, part 2), computed by `pt_analysis.metrics`.

    GET /api/ranking        {truth: {n_true_cases, n_observed_cases, true_average_outdegree,
                                     true_unobserved_source_proportion},
                             fits: [{fit, label, directed, undirected}]}
                            404 when the manifest names no truth file
    GET /api/parent-sets[?level=0.8]
                            {level, truth, fits: [{fit, label, summary}],
                             rows: [{infection, effective_size, credible_size, covered, top_probability,
                                     null_effective_size, null_credible_size, null_covered, null_top_probability}]}

`fits` holds the posterior, then the null-model Analysis and the pairwise scores (IBD; undirected only, `directed`
null) when the manifest names them; a baseline it does not name is not in the list. A ranking is its metrics by
fit-summary name (`average_precision`, `pr_auc`, `top_k_precision`, `n_edge_universe`, `n_true_edges`,
`n_distinct_scores`) and `curve`, a list of {score, recall, precision, tp, fp} by score descending, thinned to the
points a line through the whole curve needs. Parent-set rows are per Infection in id order; the `null_` columns are
there only with a null-model Analysis, and `covered` is null without truth. Results are computed once per level.
"""

from dataclasses import asdict

from starlette.routing import Route

from pt_analysis import metrics
from pt_analysis.analysis import AnalysisError
from pt_analysis.endpoints import error, json


def cached(request, key, compute):
    cache = getattr(request.app.state, "metrics_cache", None)
    if cache is None:
        cache = request.app.state.metrics_cache = {}
    if key not in cache:
        cache[key] = compute()
    return cache[key]


def ranking_json(ranking):
    if ranking is None:
        return None
    return {**ranking.metrics(), "curve": [asdict(p) for p in metrics.drawn_curve(ranking)]}


def ranking_payload(analysis):
    truth = metrics.analysis_truth(analysis)
    fits = metrics.ranking_fits(analysis, truth)
    return {"truth": {name: getattr(truth, name) for name in ("n_true_cases", "n_observed_cases",
                                                               "true_average_outdegree",
                                                               "true_unobserved_source_proportion")},
            "fits": [{"fit": f.fit, "label": f.label, "directed": ranking_json(f.directed),
                      "undirected": ranking_json(f.undirected)} for f in fits]}


def parent_sets_payload(analysis, level):
    truth = metrics.analysis_truth(analysis)
    fits = metrics.parent_set_fits(analysis, truth, level)
    rows = {}
    for fit in fits:
        prefix = "" if fit.fit == "posterior" else f"{fit.fit}_"
        for r in fit.rows:
            row = rows.setdefault(r.infection, {"infection": r.infection})
            row.update({f"{prefix}{k}": v for k, v in asdict(r).items() if k != "infection"})
    return {"level": level, "truth": truth is not None,
            "fits": [{"fit": f.fit, "label": f.label, "summary": f.summary} for f in fits],
            "rows": [rows[i] for i in sorted(rows)]}


async def ranking(request):
    analysis = request.app.state.analysis
    if analysis.truth is None:
        return error("the Analysis manifest names no truth file, so there is nothing to rank against", 404)
    try:
        return json(cached(request, "ranking", lambda: ranking_payload(analysis)))
    except (OSError, KeyError, ValueError, AnalysisError, metrics.MetricsError) as e:
        return error(str(e))


async def parent_sets(request):
    try:
        level = float(request.query_params.get("level", metrics.LEVEL))
    except ValueError:
        return error(f"level is a probability, not '{request.query_params['level']}'")
    if not 0 < level <= 1:
        return error(f"level must be in (0, 1], not {level}")
    analysis = request.app.state.analysis
    try:
        return json(cached(request, ("parent-sets", level), lambda: parent_sets_payload(analysis, level)))
    except (OSError, KeyError, ValueError, AnalysisError, metrics.MetricsError) as e:
        return error(str(e))


ROUTES = [Route("/api/ranking", ranking), Route("/api/parent-sets", parent_sets)]
