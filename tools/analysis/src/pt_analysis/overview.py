"""The Overview workspace's health check and headline posterior, and the status bar's summary, each in one payload.

`status` is cheap: closed Windows, what the manifest has, mixing counts from the convergence table and Run events of
note. `overview` reads the derived tables once: mixing by group, the attention lists, the headline functionals (ψ_edge,
mean q and composition, from `strata.per_draw`; ψ_MAP needs the pooled Parent sets and stays with
`/api/functionals`), per-Infection q and effective parent-set size, the top edges and a network thumbnail.
"""

import json

from pt_analysis import sampler, summary
from pt_analysis.derive import GROUPS, analysis_dir, analysis_source, derived, derived_analysis, source
from pt_analysis.edge_stats import edge_table
from pt_analysis.quantiles import interval
from pt_analysis.strata import strata as all_strata
from pt_analysis.strata import per_draw
from pt_analysis.times import observation_times

# R-hat and ESS thresholds of the design (`web/DESIGN.md`, Visual system): above 1.1 and under 100 warn, above 1.2
# and under 20 are critical.
RHAT_WARN, RHAT_CRITICAL = 1.1, 1.2
ESS_WARN, ESS_CRITICAL = 100, 20

# The sampler panel's default target band of acceptance rates (`web/src/lib/sampler.ts`, DEFAULT_BAND).
BAND = (0.15, 0.5)

# Run events worth a reader's attention; burn-in ending and a Run ending normally are routine.
NOTABLE = ("non_finite_value", "interrupted", "run_resumed")

# How many rows each attention list and the top-edge list hold, and the thumbnail's edge threshold.
LIST_LENGTH = 10
THUMBNAIL_THRESHOLD = 0.5


def _pooled(analysis):
    """The convergence table's pooled rows: (group, name, ess, rhat)."""
    return derived_analysis(analysis, "convergence").filter("chain IS NULL").project('"group", name, ess, rhat') \
        .fetchall()


def _events(analysis):
    return [{"chain": label, "iteration": e["iteration"], "phase": e["phase"], "type": e["type"]}
            for label, store in analysis.runs.items() for e in sampler.events(store, NOTABLE)]


def _is_current(analysis):
    try:
        recorded = json.loads((analysis_dir(analysis) / "manifest.json").read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return False
    return recorded == analysis_source(analysis)


def status(analysis, analyses=1):
    """What the status bar shows: each Chain's closed Windows per Phase and whether the derived tables are current,
    the manifest's optional inputs, the mixing counts, the Run events of note, and the attention badges of the
    workspaces: quantities above R-hat 1.1 for Convergence, sampler types outside the default band for
    Samplers."""
    current = _is_current(analysis)
    rows = _pooled(analysis)
    rhat = [r for _, _, _, r in rows if r is not None]
    ess = [e for _, _, e, _ in rows if e is not None]
    mixing = {"quantities": len(rows),
              "rhat_above_1_1": sum(r > RHAT_WARN for r in rhat), "rhat_above_1_2": sum(r > RHAT_CRITICAL for r in rhat),
              "ess_below_100": sum(e < ESS_WARN for e in ess), "ess_below_20": sum(e < ESS_CRITICAL for e in ess)}
    return {
        "derive": {"current": current, "chains": {label: source(store)["windows"]
                                                  for label, store in analysis.runs.items()}},
        "capabilities": {"epi": analysis.epi is not None, "truth": analysis.truth is not None,
                         "cohort": analysis.cohort is not None, "pairwise_scores": analysis.pairwise_scores is not None,
                         "analyses": analyses},
        "mixing": mixing,
        "events": _events(analysis),
        "attention": {"convergence": mixing["rhat_above_1_1"],
                      "samplers": len({s["type"] for s in _samplers(analysis, BAND)})},
    }


def _groups(rows):
    out = []
    for group in GROUPS:
        mine = [(e, r) for g, _, e, r in rows if g == group]
        if not mine:
            continue
        ess = [e for e, _ in mine if e is not None]
        rhat = [r for _, r in mine if r is not None]
        out.append({"group": group, "n": len(mine), "ess": ess, "rhat": rhat,
                    "ess_min": min(ess, default=None), "ess_median": sorted(ess)[len(ess) // 2] if ess else None,
                    "rhat_max": max(rhat, default=None),
                    "rhat_above": sum(r > RHAT_WARN for r in rhat)})
    return out


def _worst(rows):
    """The worst-mixing quantities: largest R-hat first, then smallest ESS."""
    ranked = sorted((r for r in rows if r[3] is not None),
                    key=lambda r: (-r[3], r[2] if r[2] is not None else float("inf")))
    return [{"group": g, "name": n, "ess": e, "rhat": r} for g, n, e, r in ranked[:LIST_LENGTH]]


def _samplers(analysis, band):
    """Per Chain, the sampler types whose acceptance over the sampling Windows lies outside `band`."""
    out = []
    for label, store in analysis.runs.items():
        sampling = {w.number for w in store.closed_windows() if w.phase == "sampling"}
        by_type = sampler.acceptance_by_type(store)
        totals = {}
        for kind, window, accepted, proposed, low, high in zip(
                by_type["type"], by_type["window"], by_type["accepted"], by_type["proposed"], by_type["rate_min"],
                by_type["rate_max"]):
            if window not in sampling:
                continue
            t = totals.setdefault(kind, [0, 0, None, None])
            t[0] += accepted
            t[1] += proposed
            t[2] = low if t[2] is None or (low is not None and low < t[2]) else t[2]
            t[3] = high if t[3] is None or (high is not None and high > t[3]) else t[3]
        for kind, (accepted, proposed, low, high) in sorted(totals.items()):
            if proposed and not band[0] <= accepted / proposed <= band[1]:
                out.append({"chain": label, "type": kind, "rate": accepted / proposed, "rate_min": low,
                            "rate_max": high})
    return out


def _infection_means(analysis, columns):
    """{column: {Infection: mean over every Draw of every Chain}} from `infection_draws`."""
    sums = {c: {} for c in columns}
    counts = {}
    for store in analysis.runs.values():
        rows = derived(store, "infection_draws").aggregate(
            "infection, count(*), " + ", ".join(f'sum("{c}")' for c in columns), "infection").fetchall()
        for infection, n, *values in rows:
            counts[infection] = counts.get(infection, 0) + n
            for c, v in zip(columns, values):
                sums[c][infection] = sums[c].get(infection, 0.0) + (v or 0.0)
    return {c: {i: s / counts[i] for i, s in sums[c].items()} for c in columns}


def _headline(analysis, edges):
    groups = [infections for _, infections in all_strata(analysis)]
    draws, names = per_draw(analysis, groups)
    pooled = {name: [] for name in names}
    for rows in draws.values():
        for _, _, values in rows:
            for name, value in values.items():
                pooled[name].append(value)
    intervals = {name: interval(values) for name, values in pooled.items()}
    asdict = lambda i: None if i is None else {"mean": i.mean, "lower": i.lower, "upper": i.upper}  # noqa: E731
    classes = {name.removeprefix("mean_"): asdict(intervals[name]) for name in names[2:-2]}
    means = _infection_means(analysis, ["q", "effective_size"])
    observed = [e for e in edges if not e.unobserved]
    top = sorted(observed, key=lambda e: (-e.pooled, e.child, e.parent))[:LIST_LENGTH]
    return {"psi_edge": asdict(intervals["psi_edge"]), "mean_q": asdict(intervals["mean_q"]), "classes": classes,
            "q": list(means["q"].values()), "effective_size": list(means["effective_size"].values()),
            "top_edges": [{"parent": e.parent, "child": e.child, "marginal": e.pooled} for e in top]}


def _thumbnail(analysis, edges):
    q = _infection_means(analysis, ["q"])["q"]
    times = observation_times(analysis)
    nodes = [{"id": i, "observation_time": t, "q": q.get(i)} for i, t in times.items()]
    shown = sorted((e for e in edges if not e.unobserved and e.pooled >= THUMBNAIL_THRESHOLD),
                   key=lambda e: -e.pooled)
    return {"threshold": THUMBNAIL_THRESHOLD, "nodes": nodes,
            "edges": [{"parent": e.parent, "child": e.child, "marginal": e.pooled} for e in shown]}


def _disagreement(analysis):
    chains = [summary.parent_probabilities(store) for store in analysis.runs.values()]
    return [{"infection": a.infection, "top": a.top, "p_top": a.p_top, "p_unobserved_only": a.p_unobserved_only,
             "spread": a.spread} for a in summary.agreement(chains)[:LIST_LENGTH]]


def overview(analysis, band=BAND):
    """Everything the Overview workspace draws but ψ_MAP. See the module docstring."""
    rows = _pooled(analysis)
    edges = edge_table(analysis)
    observed = [e for e in edges if not e.unobserved]
    spread = sorted(observed, key=lambda e: (-e.spread, e.child, e.parent))[:LIST_LENGTH]
    return {
        "groups": _groups(rows),
        "worst": _worst(rows),
        "edge_spread": [{"parent": e.parent, "child": e.child, "pooled": e.pooled, "spread": e.spread,
                         "switches": e.switches, "chains": {c: s.mean for c, s in e.chains.items()}} for e in spread],
        "disagreement": _disagreement(analysis),
        "band": list(band),
        "samplers": _samplers(analysis, band),
        "events": _events(analysis),
        "headline": _headline(analysis, edges),
        "thumbnail": _thumbnail(analysis, edges),
    }
