"""The network view's truth overlay and the error table (U6, part 1): where a fit's edges disagree with the truth of
a simulated input.

An edge set is one fit's scores of (parent, child) pairs: the posterior's pooled Edge marginals, the null-model
Analysis's when the manifest names one, and the pairwise scores (IBD) when it names a scores CSV. The first two are
directed; the pairwise scores are undirected, a pair's score being the larger of its rows in either order, and a pair
is true when a true edge joins it either way. Scores, truth and pooling are `pt_analysis.metrics`'.

The overlay draws every true edge with its score and the false edges scoring at least a floor. The error table calls
an edge at a threshold: a true edge scoring below it is missed, a false edge scoring at or above it is a false edge,
or a direction error when its reverse is true. A true edge missed because its false reverse is called is one error,
the direction error, so it has no missed row. Rows rank by how confidently wrong they are: 1 - score for a missed
edge, the score for a false one.
"""

from dataclasses import dataclass

from pt_analysis import inputs, metrics, times

SOURCES = ("posterior", "null", "ibd")
DEFAULT_FLOOR = 0.01
DEFAULT_THRESHOLD = 0.5


class UnknownSource(LookupError):
    """An edge set the Analysis does not have."""


@dataclass(frozen=True)
class EdgeSet:
    source: str  # "posterior", "null" or "ibd"
    label: str
    directed: bool
    scores: dict  # {(parent, child): score}; undirected sets hold each pair once, ids in order


def edge_sets(analysis, only=None):
    """The posterior's edge set, then the null-model Analysis's and the pairwise scores' when the manifest names
    them; with `only`, just that source's."""
    sets = []
    if only in (None, "posterior"):
        sets.append(EdgeSet("posterior", analysis.label, True, metrics.edge_scores(analysis)))
    if only in (None, "null") and analysis.null_model:
        null = metrics.null_analysis(analysis)
        sets.append(EdgeSet("null", null.label, True, metrics.edge_scores(null)))
    if only in (None, "ibd") and analysis.pairwise_scores:
        scores = {}
        for (a, b), s in metrics.read_pairwise_scores(analysis.pairwise_scores).items():
            if a != b:
                k = metrics.pair(a, b)
                scores[k] = max(scores.get(k, s), s)
        sets.append(EdgeSet("ibd", "IBD", False, scores))
    return sets


def match(analysis):
    """How the truth file's Infection ids meet the Analysis's: the number of true edges with both ends in it, the
    Analysis's Infections without a truth node, and the truth's Infections the Analysis does not observe."""
    ids = [str(i) for i in metrics.infections(analysis)]
    truth = metrics.load_truth(analysis.truth, ids)
    simulated = list(inputs.nodes(analysis.truth))
    inside = set(ids)
    return {"n_true_edges": len(truth.edges),
            "unmatched_infections": [i for i in ids if i not in set(simulated)],
            "unobserved_truth_infections": [i for i in simulated if i not in inside]}


def _edges(edge_set, ids, true_edges, observed):
    """Every true edge and every scored pair of the set, both ends among `ids`, as rows {parent, child, marginal, true,
    reverse_true, reverse_marginal, lag}. `true_edges` None marks no edge true or false."""
    inside = set(ids)
    scores = {k: metrics.normalise(s) for k, s in edge_set.scores.items()
              if k[0] in inside and k[1] in inside and k[0] != k[1]}
    known = true_edges is not None
    truth = true_edges or frozenset()
    if edge_set.directed:
        keys = scores.keys() | truth
        rows = [(p, c, scores.get((p, c), 0.0), (p, c) in truth, (c, p) in truth, scores.get((c, p), 0.0))
                for p, c in keys]
    else:
        true_pairs = {}
        for p, c in sorted(truth):
            true_pairs.setdefault(metrics.pair(p, c), (p, c))
        keys = scores.keys() | true_pairs.keys()
        rows = [(*true_pairs.get(k, k), scores.get(k, 0.0), k in true_pairs, None, None) for k in keys]
    return [{"parent": p, "child": c, "marginal": m, "true": t if known else None,
             "reverse_true": r if known else None, "reverse_marginal": rm, "lag": times.lag(observed, p, c)}
            for p, c, m, t, r, rm in rows]


def _context(analysis):
    ids = [str(i) for i in metrics.infections(analysis)]
    truth = metrics.load_truth(analysis.truth, ids) if analysis.truth else None
    return ids, times.observation_times(analysis), truth


def overlay(analysis, floor=DEFAULT_FLOOR):
    """{floor, truth, sources}: per edge set its true edges whatever their score and its other edges scoring at least
    `floor`, by score descending. `truth` is `match` without a truth file None, and then no edge is marked."""
    ids, observed, truth = _context(analysis)
    sources = []
    for s in edge_sets(analysis):
        rows = [r for r in _edges(s, ids, truth and truth.edges, observed) if r["true"] or r["marginal"] >= floor]
        rows.sort(key=lambda r: (-r["marginal"], r["parent"], r["child"]))
        sources.append({"source": s.source, "label": s.label, "directed": s.directed, "edges": rows})
    return {"floor": floor, "truth": match(analysis) if truth else None, "sources": sources}


def errors(analysis, source="posterior", threshold=DEFAULT_THRESHOLD):
    """{source, label, directed, threshold, truth, sources, rows}: the error rows of one edge set at `threshold`, each
    {kind, parent, child, marginal, reverse_true, reverse_marginal, lag, confidence, parent_true_infection_time,
    child_true_infection_time}, most confidently wrong first. A true p->c below `threshold` whose reverse c->p is
    called and false is counted once, as c->p's direction row; p->c has no missed row.
    Raises UnknownSource for an edge set the Analysis does not have. Needs a truth file."""
    if source not in SOURCES:
        raise UnknownSource(f"source is one of {', '.join(SOURCES)}, not '{source}'")
    [edge_set] = edge_sets(analysis, only=source) or [None]
    if edge_set is None:
        missing = "null-model Analysis" if source == "null" else "pairwise scores"
        raise UnknownSource(f"the Analysis manifest names no {missing}")
    ids, observed, truth = _context(analysis)
    rows = []
    for r in _edges(edge_set, ids, truth.edges, observed):
        if r["true"] and r["marginal"] < threshold:
            if r["reverse_true"] is False and r["reverse_marginal"] >= threshold:
                continue  # the called reverse is a direction error; one error, its row
            rows.append({"kind": "missed", **r, "confidence": 1 - r["marginal"]})
        elif not r["true"] and r["marginal"] >= threshold:
            rows.append({"kind": "direction" if r["reverse_true"] else "false", **r, "confidence": r["marginal"]})
    true_time = {i: n.get("infection_time") for i, n in inputs.nodes(analysis.truth).items()}
    for r in rows:
        del r["true"]
        r["parent_true_infection_time"] = true_time.get(r["parent"])
        r["child_true_infection_time"] = true_time.get(r["child"])
    rows.sort(key=lambda r: (-r["confidence"], r["parent"], r["child"]))
    return {"source": source, "label": edge_set.label, "directed": edge_set.directed, "threshold": threshold,
            "truth": match(analysis), "sources": available_sources(analysis), "rows": rows}


def available_sources(analysis):
    """The edge sets the Analysis has, in `SOURCES` order, without computing their scores."""
    return [s for s, named in zip(SOURCES, (True, analysis.null_model, analysis.pairwise_scores)) if named]
