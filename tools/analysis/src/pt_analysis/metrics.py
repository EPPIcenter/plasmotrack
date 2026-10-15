"""Accuracy of an Analysis against the truth of a simulated input (U6): edge ranking and parent-set metrics.

Definitions follow the manuscript supplement and the thesis pipeline's per-fit summary
(`remote/R/final_analysis/R/edge_metrics.R`, `credible_parent_sets.R`, `truth_from_json.R`), with one change: edge
scores and set weights are pooled as sums over all Draws of all Chains divided by the number of Draws (the
manuscript's estimator), not means over the Draws where a row is present. Metric names are the fit summary's
column names, so the per-fit summary maps them directly.
"""

import csv
import math
import statistics
from dataclasses import dataclass, field

from pt_analysis import functionals, inputs
from pt_analysis.analysis import Analysis
from pt_analysis.derive import derived
from pt_analysis.functionals import LEVEL, TOLERANCE

PAIRWISE_COLUMNS = ("a", "b", "score")


class MetricsError(Exception):
    """A file the metrics read cannot be used."""


@dataclass(frozen=True)
class Truth:
    """The true network of a simulated input, restricted to the Analysis's Infections."""

    infections: tuple[str, ...]  # the Analysis's Infections, the observed cases
    edges: frozenset[tuple[str, str]]  # (parent, child) with both ends among the Infections
    # Per Infection, its true Parent set as tokens comparable with a fitted set whose Latent parent's `_copy` is
    # stripped: the true parents among the Infections, plus the Infection's own id (the unobserved source) when it
    # is a root or has a true parent outside the Infections.
    parent_sets: dict[str, frozenset[str]]
    n_true_cases: int  # Infections of the simulation, observed or not
    n_observed_cases: int
    true_average_outdegree: float  # true edges with both ends observed per observed Infection
    true_unobserved_source_proportion: float  # share of Infections whose true set includes the unobserved source


def load_truth(path, infections):
    """The truth file of a simulation (`network: [{src, dest}]`, `nodes: [{id, ...}]`; gzipped when named `.gz`),
    restricted to `infections`, the Analysis's Infection ids."""
    raw = inputs.read_json(path)
    observed = tuple(str(i) for i in infections)
    inside = set(observed)
    network = [(str(e["src"]), str(e["dest"])) for e in raw.get("network") or []]
    edges = frozenset((s, d) for s, d in network if s in inside and d in inside and s != d)
    parent_sets = {}
    for child in observed:
        parents = [s for s, d in network if d == child]
        tokens = {s if s in inside else child for s in parents}
        parent_sets[child] = frozenset(tokens or {child})
    n = len(observed)
    return Truth(
        infections=observed, edges=edges, parent_sets=parent_sets, n_true_cases=len(raw.get("nodes") or []),
        n_observed_cases=n, true_average_outdegree=len(edges) / n if n else None,
        true_unobserved_source_proportion=sum(c in s for c, s in parent_sets.items()) / n if n else None)


def infections(analysis):
    """The Analysis's Infection ids, from its first Chain's structure."""
    store = next(iter(analysis.runs.values()))
    return [i for (i,) in store.structure("infections").order("infection").project("infection").fetchall()]


def analysis_truth(analysis):
    """The truth file the Analysis's manifest names, restricted to its Infections; None when it names none."""
    return load_truth(analysis.truth, infections(analysis)) if analysis.truth else None


# --- edge ranking ----------------------------------------------------------------------------------------------------

def normalise(score):
    """A score clipped to [0, 1] and rounded to 10 decimals, so that sums equal up to rounding tie."""
    return round(min(1.0, max(0.0, float(score))), 10)


@dataclass(frozen=True)
class CurvePoint:
    """The precision-recall curve at one distinct score: every pair scoring at least `score` is called an edge."""

    score: float
    tp: int
    fp: int
    recall: float
    precision: float


@dataclass(frozen=True)
class Ranking:
    """How well scores rank the true edges of a universe of pairs."""

    n_edge_universe: int
    n_true_edges: int
    n_distinct_scores: int
    average_precision: float | None  # sum over breakpoints of the recall gained times the precision there
    pr_auc: float | None  # trapezoid area under the curve with (recall 0, precision 1) prepended
    top_k_precision: float | None  # precision of the pairs scoring at least the k-th highest, k = n_true_edges
    curve: list[CurvePoint] = field(default_factory=list)  # by score, descending

    def metrics(self):
        """The metrics without the curve, by their fit-summary names less the `directed_`/`undirected_` prefix."""
        return {"n_edge_universe": self.n_edge_universe, "n_true_edges": self.n_true_edges,
                "n_distinct_scores": self.n_distinct_scores, "average_precision": self.average_precision,
                "pr_auc": self.pr_auc, "top_k_precision": self.top_k_precision}


def rank_groups(groups):
    """The Ranking of a universe given as {normalised score: (true pairs, false pairs)} at that score."""
    positives = sum(t for t, _ in groups.values())
    universe = sum(t + f for t, f in groups.values())
    if not positives:
        return Ranking(universe, 0, len(groups), None, None, None, [])
    curve, tp, fp = [], 0, 0
    for score in sorted(groups, reverse=True):
        t, f = groups[score]
        tp, fp = tp + t, fp + f
        curve.append(CurvePoint(score, tp, fp, tp / positives, tp / (tp + fp)))
    average_precision, auc, recall, precision = 0.0, 0.0, 0.0, 1.0
    for p in curve:
        average_precision += (p.recall - recall) * p.precision
        auc += (p.recall - recall) * (p.precision + precision) / 2
        recall, precision = p.recall, p.precision
    # Every pair scoring at least the k-th highest score: the curve up to the breakpoint where k pairs are reached.
    k = min(positives, universe)
    top_k = next(p for p in curve if p.tp + p.fp >= k)
    return Ranking(universe, positives, len(groups), average_precision, auc, top_k.precision, curve)


def rank(scores, truth):
    """The Ranking of pairs with these scores, of which those with a true `truth` are true edges."""
    groups = {}
    for score, true in zip(scores, truth, strict=True):
        t, f = groups.get(normalise(score), (0, 0))
        groups[normalise(score)] = (t + 1, f) if true else (t, f + 1)
    return rank_groups(groups)


def drawn_curve(ranking):
    """The curve's points a line through all of them needs: the first, the last, and those where the true-positive
    count changes on either side. Between two kept points at one count, the dropped points lie on the vertical
    segment joining them (same recall, precision falling), so the line is the same."""
    curve = ranking.curve
    return [p for i, p in enumerate(curve)
            if i in (0, len(curve) - 1) or curve[i - 1].tp != p.tp or curve[i + 1].tp != p.tp]


def add(groups, score, true, count=1):
    t, f = groups.get(score, (0, 0))
    groups[score] = (t + count, f) if true else (t, f + count)


def universe_groups(scores, true_keys, size):
    """Groups for `rank_groups` of a universe of `size` pairs: those in `scores` at their score, the rest at 0."""
    groups = {}
    for key, score in scores.items():
        add(groups, score, key in true_keys)
    unscored_true = len(true_keys - scores.keys())
    unscored = size - len(scores)
    if unscored:
        add(groups, 0.0, True, unscored_true)
        add(groups, 0.0, False, unscored - unscored_true)
    return groups


def directed_ranking(scores, infections, true_edges):
    """Ranking of every ordered pair (parent, child) of distinct `infections` by `scores`, a missing score being 0,
    against `true_edges`."""
    inside = set(infections)
    scored = {(p, c): normalise(s) for (p, c), s in scores.items() if p in inside and c in inside and p != c}
    true = {(p, c) for p, c in true_edges if p in inside and c in inside and p != c}
    n = len(inside)
    return rank_groups(universe_groups(scored, true, n * (n - 1)))


def pair(a, b):
    return (a, b) if a <= b else (b, a)


def undirected_ranking(scores, infections, true_edges, symmetric=False):
    """Ranking of every unordered pair of distinct `infections`, true when a true edge joins them either way. With
    `symmetric`, `scores` holds one score per pair (in either order; the largest counts); otherwise it holds directed
    scores and a pair's score is the sum of its two directions clipped to 1."""
    inside = set(infections)
    scored = {}
    for (a, b), s in scores.items():
        if a in inside and b in inside and a != b:
            key, s = pair(a, b), normalise(s)
            scored[key] = max(scored.get(key, 0.0), s) if symmetric else normalise(scored.get(key, 0.0) + s)
    true_pairs = {pair(a, b) for a, b in true_edges if a in inside and b in inside and a != b}
    n = len(inside)
    return rank_groups(universe_groups(scored, true_pairs, n * (n - 1) // 2))


def edge_scores(analysis):
    """The posterior transmission matrix's observed entries: {(parent, child): Edge marginal}, each the sum of its
    `edges` weight over every Draw of every Chain divided by the number of Draws. Edges with no weight at any Draw,
    and edges from the unobserved source, are left out."""
    ids = set(infections(analysis))
    totals, draws = {}, 0
    for store in analysis.runs.values():
        draws += store.draws()
        rows = derived(store, "edges").aggregate("parent, child, sum(prob)").fetchall()
        for parent, child, weight in rows:
            if parent in ids:
                totals[(parent, child)] = totals.get((parent, child), 0.0) + weight
    return {key: total / draws for key, total in totals.items()} if draws else {}


def edge_rankings(analysis, truth):
    """The directed and undirected Rankings of the Analysis's Edge marginals against `truth`."""
    scores = edge_scores(analysis)
    return {"directed": directed_ranking(scores, truth.infections, truth.edges),
            "undirected": undirected_ranking(scores, truth.infections, truth.edges)}


def read_pairwise_scores(path):
    """A pairwise scores CSV (IBD) with columns a, b, score, one row per pair of Infection ids: {(a, b): score}.
    Rows for the same pair in the same order keep the largest score."""
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in PAIRWISE_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise MetricsError(f"{path}: pairwise scores need columns a, b, score; missing {', '.join(missing)}")
        scores = {}
        for row in reader:
            key = (row["a"].strip(), row["b"].strip())
            try:
                score = float(row["score"])
            except ValueError:
                raise MetricsError(f"{path}: score of {key[0]}, {key[1]} is not a number: {row['score']!r}") from None
            scores[key] = max(scores.get(key, score), score)
    return scores


def pairwise_ranking(analysis, truth):
    """The undirected Ranking of the Analysis's pairwise scores CSV against `truth`: a pair's score is the largest
    of its rows in either order, 0 without one. Pairs with an id outside the Analysis are left out."""
    return undirected_ranking(read_pairwise_scores(analysis.pairwise_scores), truth.infections, truth.edges,
                              symmetric=True)


@dataclass(frozen=True)
class RankingFit:
    """One set of edge scores ranked against truth: the posterior, the null-model posterior or the pairwise
    scores (undirected only)."""

    fit: str  # "posterior", "null" or "ibd"
    label: str
    directed: Ranking | None
    undirected: Ranking


def null_analysis(analysis):
    """The null-model Analysis the manifest names, or None."""
    return Analysis.open(analysis.null_model) if analysis.null_model else None


def ranking_fits(analysis, truth):
    """The posterior's Rankings, then the null-model Analysis's and the pairwise scores' when the manifest names
    them, each against `truth`."""
    fits = [RankingFit("posterior", analysis.label, **edge_rankings(analysis, truth))]
    null = null_analysis(analysis)
    if null is not None:
        fits.append(RankingFit("null", null.label, **edge_rankings(null, truth)))
    if analysis.pairwise_scores:
        fits.append(RankingFit("ibd", "IBD", None, pairwise_ranking(analysis, truth)))
    return fits


# --- parent-set metrics ----------------------------------------------------------------------------------------------

def pooled_sets(analysis):
    """Each Infection's Parent sets with their weights pooled over Chains, {infection: [(members, weight)]} by
    weight descending then members: a set's weight summed over every Draw of every Chain divided by the number of
    Draws. Sets are matched across Chains by their members, as set ids are a Run's own."""
    totals, draws = {}, 0
    for store in analysis.runs.values():
        n = store.draws()
        draws += n
        for s in functionals.ranked_sets(store):
            key = (s.infection, s.members)
            totals[key] = totals.get(key, 0.0) + s.weight * n
    if not draws:
        return {}
    sets = {}
    for (infection, members), total in totals.items():
        sets.setdefault(infection, []).append((members, total / draws))
    for ranked in sets.values():
        ranked.sort(key=lambda s: (-s[1], s[0]))
    return sets


def credible_size(weights, level=LEVEL):
    """The size of the credible collection of sets with these weights, in rank order: the smallest prefix whose
    mass reaches `level` (within the tolerance of `functionals.credible_collection`), or all of them."""
    mass = 0.0
    for k, w in enumerate(weights, start=1):
        mass += w
        if mass >= level - TOLERANCE:
            return k
    return len(weights)


def truth_tokens(members):
    """A fitted set's members as tokens comparable with `Truth.parent_sets`: the Latent parent `<id>_copy` as `<id>`."""
    return frozenset(m.removesuffix("_copy") for m in members)


@dataclass(frozen=True)
class InfectionParentSets:
    """One Infection's pooled Parent-set distribution summarised."""

    infection: str
    effective_size: float  # exp of the entropy of the pooled weights
    credible_size: int  # sets in the credible collection at the level
    covered: bool | None  # whether the true set is in the credible collection; None without truth
    top_probability: float  # the modal set's pooled weight


def parent_set_rows(analysis, truth=None, level=LEVEL):
    """Per Infection, in id order, its effective parent-set size, credible collection size at `level`, whether
    `truth`'s set is in that collection, and its modal set's pooled weight."""
    rows = []
    for infection, sets in sorted(pooled_sets(analysis).items()):
        weights = [w for _, w in sets]
        k = credible_size(weights, level)
        covered = None
        if truth is not None:
            true_set = truth.parent_sets.get(infection)
            covered = any(truth_tokens(members) == true_set for members, _ in sets[:k])
        entropy = -sum(w * math.log(w) for w in weights if w > 0)
        rows.append(InfectionParentSets(infection, math.exp(entropy), k, covered, max(weights)))
    return rows


def level_suffix(level):
    """The fit summary's suffix for a credible level: `_80` at 0.8, else `_a<percent>` (`_a095` at 0.95)."""
    return "_80" if abs(level - 0.8) < 1e-9 else f"_a{round(100 * level):03d}"


def parent_set_summary(rows, level=LEVEL):
    """The fit summary's parent-set columns over `rows`: mean and median credible collection size and coverage at
    `level`, and the mean modal-set weight. Coverage is None without truth."""
    suffix = level_suffix(level)
    sizes = [r.credible_size for r in rows]
    covered = [r.covered for r in rows if r.covered is not None]
    mean = lambda xs: sum(xs) / len(xs) if xs else None  # noqa: E731
    return {f"mean_credible_parent_set_size{suffix}": mean(sizes),
            f"median_credible_parent_set_size{suffix}": statistics.median(sizes) if sizes else None,
            f"credible_parent_set_coverage{suffix}": mean(covered),
            "mean_top_parent_set_probability": mean([r.top_probability for r in rows])}


@dataclass(frozen=True)
class ParentSetFit:
    """Parent-set metrics of one fit: the posterior or the null-model posterior."""

    fit: str  # "posterior" or "null"
    label: str
    rows: list[InfectionParentSets]
    summary: dict


def parent_set_fits(analysis, truth=None, level=LEVEL):
    """The posterior's parent-set metrics, then the null-model Analysis's when the manifest names one."""
    fits = []
    for fit, a in (("posterior", analysis), ("null", null_analysis(analysis))):
        if a is not None:
            rows = parent_set_rows(a, truth, level)
            fits.append(ParentSetFit(fit, a.label, rows, parent_set_summary(rows, level)))
    return fits


# --- the fit summary's columns ---------------------------------------------------------------------------------------

def fit_summary_metrics(analysis, truth, level=LEVEL):
    """The fit summary's truth, ranking and parent-set columns for the Analysis's posterior, by column name."""
    rankings = edge_rankings(analysis, truth)
    rows = parent_set_rows(analysis, truth, level)
    columns = {"n_true_cases": truth.n_true_cases, "n_observed_cases": truth.n_observed_cases,
               "true_average_outdegree": truth.true_average_outdegree,
               "true_unobserved_source_proportion": truth.true_unobserved_source_proportion}
    for kind in ("directed", "undirected"):
        columns.update({f"{kind}_{name}": value for name, value in rankings[kind].metrics().items()})
    return columns | parent_set_summary(rows, level)
